import daemon
import os
import time
import logging
import configparser
from impl.dual_logging import get_last_log_ID
from utils.condor_utils import *


class log_manager:
    
    def __init__(self,cfg):
        #  dT time after which a job is considered historic (Folder is zipped and moved to "old")
        self.dT = float(cfg['Time variables']['Historic after'])
        # dt time after which we considered a job UNSUCCESSFULL after being removed from the condor_q (job_control_daemon gets killed)
        self.dt = float(cfg['Time variables']['Unsuccessfull after'])
        # dH time that a job can be on hold before being UNSUCCESSFULL
        self.dH = float(cfg['Time variables']['On hold untill'])
        # how often we run the daemon
        self.running_freq = float(cfg['Time variables']['Running_freq'])
        # global log file
        self.logfile = cfg['Paths']['Global logfile']
        # global control file for cleaning up and debug purposes
        self.cfg_log_file = cfg['Paths']['Cfg_log_file']
        # log daemon logs
        self.log_daemon_logs = cfg['Paths']['Log_daemon_logs']
        # Directory for historical jobs (completed and failed)
        self.completed_dir = cfg['Paths']['Completed_dir']
        self.failed_dir = cfg['Paths']['Failed_dir']
        # Archive for log files
        self.logs_folder = cfg['Paths']['Logs_folder']
        # dictionary of possible condor job status
        self.job_status_dict = cfg['Job status']
        # List possible epilogue of condor jobs
        self.end_status = [self.job_status_dict['Submission_err'],self.job_status_dict['Unsuccessfull'], self.job_status_dict['Done'], self.job_status_dict['Killed_by_log_daem']]
        
        self.parser = None
        self.log = self.get_log_file(self.log_daemon_logs,'%(asctime)s - %(levelname)s - %(message)s')
        
        
    def get_log_file(self,log_daemon_logs,formatt):
        formatter = logging.Formatter(formatt)
        handler = logging.FileHandler(log_daemon_logs)        
        handler.setFormatter(formatter)
    
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        logger.addHandler(handler)
    
        return logger
        
    # read files is called outside the class to allow easier unittesting    
    def read_files(self): 
    	# Read config file
    	self.parser = configparser.ConfigParser()
    	try:
    		self.parser.read(self.cfg_log_file)
    	except:
    		raise ImportError("Could not read global cfg file")
    		
    	# Find last ID in config file
    	if self.parser.sections():
    		last_ID_cfg = int(self.parser.sections()[-1])
    	else: last_ID_cfg = 0
    	# Find last ID in log file
    	last_ID_log = int(get_last_log_ID())

        # Get daemons. Read daemons before updating config!
    	self.log.info("Get job daemons")
    	self.daemons = get_job_daemons("job_control_daemon.py")
            
    	# Create new sections for the newly added IDs
    	if last_ID_log > last_ID_cfg:
    		id_range = range(last_ID_cfg+1,last_ID_log+1)
    		self.add_id_sections(id_range)
    		# Check for free running daemons
            # NOTE: done here to avoid checking on not up to date cfg file
    		self.log.info("Find and kill daemons for failed or untracked jobs")
    		self.kill_untracked_daemons(self.daemons)
            
    
    def update_log_file(self):
        	# Wake up to work 
        # Get job status
        	self.log.info("Reading condor queue")
        	self.all_jobs = get_jobs_status()
            
        	parser = self.parser
        	for i in parser.sections():
        		if parser[i]['Submission date']!= '-':
        			job_age = get_job_age(parser[i]['Submission date'],"%Y-%m-%d %H:%M:%S")
        			# Update the status for recent sections (date - submission date < dT)
        			if job_age < self.dT:
        				self.log.info("Checking section {}".format(i))
        				# ideal status
        				self.log.info("Update ideal status")
        				self.update_ideal_status(parser[i])
        				# condor status
        				if parser[i]['Condor status'] not in self.end_status:
        					self.log.info("Update condor status")
        					self.update_job_status(parser[i],self.all_jobs)
        				# job control daemon
        				self.log.info("Update daemon status")
        				self.update_job_daemon_status(parser[i],self.daemons)
        				# kill daemons for unsuccessful jobs
        				self.kill_running_daemons(parser[i])
        			# Clean up data for historic jobs
        			else:
        				self.cleanup_workdir(parser[i],self.completed_dir,self.failed_dir)
        		else:
        			parser[i]['Condor status'] = self.job_status_dict['submission_err']
        			self.log.info("Submission error. Going to archive the job as failed")
        			self.cleanup_workdir(parser[i],self.completed_dir,self.failed_dir)
            
        # Clean up historical logs
        	self.log.info("Zip historical log files")
        	self.archive_logs()
        	
        	# Dump the new configuration file
        	with open(self.cfg_log_file, 'w') as configfile:
        			self.parser.write(configfile)
                    
        # Sleep
        	self.log.info("Going to sleep for {} s\n\n".format(self.running_freq))
        	time.sleep(self.running_freq)
        	self.log.info("Waking up to work")
        
    def update_ideal_status(self,pars_sec):
        cfg = configparser.ConfigParser()
        cfg.read(pars_sec['Simulation settings'])
        pars_sec['Status'] = cfg['DEFAULT']['status']
        self.log.debug("status: {}".format(pars_sec['Status']))
    

    def update_job_status(self,pars_sec, all_jobs):
        if pars_sec['Condor id'] in all_jobs:  # job is in condor_q (running, idle, done, hold)
            job_id = pars_sec['Condor id']
            pars_sec['Condor status'] = str(all_jobs[job_id])
            self.log.info("Job in condor queue")
            if job_on_hold(all_jobs,job_id):
                if 'On hold since' not in pars_sec:
                    self.log.debug("Job was put on hold")
                    pars_sec['On hold since'] = time.strftime("%Y-%m-%d %H:%M:%S")
                if condor_check_run() == 0: # no problems with condor
                    # try to release
                    self.log.debug("Try to release job")
                    ret = release_condor_job(job_id)
                    if ret == 0:
                        pars_sec['On hold since'] = 'Released successfully'
                        self.log.debug("Job released with exit code 0")
                    else: # if it was too long on hold, remove and mark as UNSUCCESSFULL
                        if get_job_age(pars_sec['On hold since'],"%Y-%m-%d %H:%M:%S") > self.dH:
                            self.log.debug("Tried to release for {} sand failed. Job will be removed from queue".format(self.dH))
                            remove_condor_job(job_id) 
                            pars_sec['Condor status'] = self.job_status_dict['killed_by_log_daem']
        else:            
            if pars_sec['Status'] == 'FINISHED':
                self.log.info("Job not in the queue because of successful termination")
                pars_sec['Condor status'] = self.job_status_dict['done']
            else:
                 if pars_sec['Condor status'] != self.job_status_dict['checking']:
                     pars_sec['Condor status'] = self.job_status_dict['checking']
                     pars_sec['Last checked'] = time.strftime("%Y-%m-%d %H:%M:%S")
                     self.log.info("Job not in the queue but not terminated. Wait.")
                 else:
                     if get_job_age(pars_sec['Last checked'],"%Y-%m-%d %H:%M:%S") > self.dt: # if the job has been checked for longer then an hour
                         self.log.info("Job not in the queue and still not terminated after {} s. Marked as failed.".format(self.dt))
                         pars_sec['Condor status'] = self.job_status_dict['unsuccessfull']
        self.log.debug("Condor status: {}".format(pars_sec['Condor status']))
                         
    def update_job_daemon_status(self,pars_sec,daemons):
        if pars_sec['Work_dir'] in daemons:
            self.log.info("Daemon is running")
            pid = daemons[pars_sec['Work_dir']]
            pars_sec['Job control daemon'] = 'Running with pid {}'.format(pid)
        elif pars_sec['Job control daemon'] != 'Daemon killed': 
            self.log.info("Daemon not running. Assumed successful termination")
            pars_sec['Job control daemon'] = 'Daemon successfully finished'
    
    def is_daemon_to_kill(self,pars_sec):
        is_daemon = 'Running' in pars_sec['Job control daemon']
        if pars_sec['Condor status'] in self.end_status and is_daemon:
            return True
        else:
            return False
    
    def kill_running_daemons(self,pars_sec,test=False):
        if self.is_daemon_to_kill(pars_sec):
            pid = pars_sec['Job control daemon'].split(" ")[-1]
            if not test:
                ret = kill_process(pid)
            pars_sec['Job control daemon'] = 'Daemon killed'
            self.log.info("Daemon for job {} killed by program".format(pars_sec['Condor id']))
            
            if test: return pid
            
            
    def kill_untracked_daemons(self,daemons,test=False):
        # kill daemons not associated to any job
        a = list() # for unittest
        workdirs = [self.parser[p]['Work_dir'] for p in self.parser.sections()]
        for d in daemons.items():
            if d[0] not in workdirs:
                pid = d[1]
                if not test:
                    ret = kill_process(pid)
                    self.log.info("Daemon killed because not connected to any tracked job.\nWorkdir: {}".format(d[0]))
                else: a.append(pid)
                    
        if test: return a
                    
    def archive_logs(self):
        old_dir = self.logs_folder+"/old"
        if not os.path.isdir(old_dir):
            os.makedirs(old_dir)
        os.chdir(self.logs_folder)
        files = [f for f in os.listdir() if ".py" in f]
        for f in files:
            date_str = f.split(".py_")[1].split(".")[0]
            age = get_job_age(date_str, "%Y_%m_%d_%H_%M_%S")
            if age > self.dT:
                name = "old/"+f.split(".log")[0]+".zip"
                zip_files(name,[f])
                os.remove(f)
    
    
    def cleanup_workdir(self,pars_sec,completed_dir,failed_dir):
        if pars_sec['Status'] == 'ARCHIVED':
            return
        
        if not os.path.isdir(failed_dir):
            os.makedirs(failed_dir)
        if not os.path.isdir(completed_dir):
            os.makedirs(completed_dir)
            
        self.log.info("Job considered historical. Going to zip output and work dir")
        
        workdir = pars_sec['Work_dir']
        data = workdir.split("work")[0] 
        workdir = workdir.split("/rungate")[0]
        self.log.debug("Working dir: {}".format(workdir))
        # get output dir
        fname =  workdir.split("work")[1]
        output = data+"output"+fname
        self.log.debug("Output dir: {}".format(output))
        
        base_out = data + "--output" # destination file for output dir
        base_work = data + "--work" # destination file for work dir
        
        self.log.info("Zipping {0} in {1}. Removing {0}".format(output,base_out+".zip"))
        zip_and_clean_folder(base_out,output)
        self.log.info("Zipping {0} in {1}. Removing {0}".format(workdir,base_work+".zip"))
        zip_and_clean_folder(base_work,workdir)
        
        # make zip file for the two previouslzy zipped files   
        if pars_sec['Status'] == 'FINISHED':
            os.chdir(data)
            zip_files(completed_dir+fname+'.zip',['--output.zip','--work.zip'])
            self.log.info("Archiving BOTH in {}".format(completed_dir+fname+'.zip'))
        else:
            os.chdir(data)
            zip_files(failed_dir+fname+'.zip',['--output.zip','--work.zip'])
            self.log.info("Archiving in {}".format(failed_dir+fname+'.zip'))
            
        os.remove(base_out+'.zip')    
        os.remove(base_work+'.zip')
        pars_sec['Status'] = 'ARCHIVED' 
            
        
            
    def add_id_sections(self,id_range):
    	# TODO: we are reading the whole file. Better to read obly a buffer with the last N lines?
        nw = 1 # how many lines after IdealID is the working_dir in the log file
        nd = 3 # same for submission date
        ns = 5 # same for user settings
        nc = 6 # condor id
        
        with open(self.logfile,'r') as f:
            lines = f.readlines()
            
        ID_lines = [l for l in lines if ("IdealID:" in l.split(" ") and int(l.split(" ")[1]) in id_range)]
        
        #new_ids = [line.split(" ")[1] for line in ID_lines]
		
		# Get workdir for missing IDs
        work_dir_lines = [lines[lines.index(l)+nw][:-1] for l in ID_lines]
        work_dirs = [line.split(" ")[2] for line in work_dir_lines]

        dates = list()
        settings = list()
        condor_ids = list()
        for line in ID_lines:
            # Get submission date for missing IDs
            date_line = lines[lines.index(line)+nd][:-1]
            if date_line.split(" ")[2] != 'error':
                dates.append(date_line.split(" ")[3]+" "+date_line.split(" ")[4])
                # Get settings for missing IDs
                settings.append(lines[lines.index(line)+ns][:-1])
                # Get condor cluster ID
                condor_ids.append(lines[lines.index(line)+nc][:-1].split(" ")[-1])
            else:
                dates.append('-')
                settings.append('-')
    
     
        for i in range(len(id_range)):
            self.config_template(self.parser,id_range[i],work_dirs[i],dates[i],condor_ids[i],settings[i])
            self.log.info("Added section with ID: {}".format(id_range[i]))
        
    		
    def config_template(self,config,idealID,workdir,date,condor_id,settings):
    
        config[idealID] = {'Submission date': date,
    						 'Work_dir': workdir,
                             'Simulation settings': settings,
                             'Status': '',
    						 'Condor id': condor_id,
                             'Condor status': '',
                             'Job control daemon': ''}



if __name__ == '__main__':
    
    cfg_parser = configparser.ConfigParser()
    cfg_parser.read("/opt/IDEAL-1.1test/cfg/log_daemon.cfg")
    
    #with daemon.DaemonContext():
	# Create log manager class:
    manager = log_manager(cfg_parser)
	
    while True:  # To stop run bin/stop_log_daemon
		# Read main log file and update config file with new entries
        manager.read_files()
        manager.update_log_file()

    
