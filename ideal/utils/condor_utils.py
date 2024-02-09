import subprocess
import os
import time
import shutil
from zipfile import ZipFile

def shell_output_ret(shell_command):
    output = subprocess.getstatusoutput(shell_command) # Byte object
    ret = output[0]
    out = output[1]
    return ret, out
    
def shell_output(shell_command):
	out = subprocess.check_output(shell_command, shell = True)
	return str(out)
	
def condor_check_run():
    # TODO: check if all machines are connected (nr of expected machine should be known)
    ret, out = shell_output_ret("ps -ef | grep condor")
    if "condor_master" in out and "condor_schedd" in out:
        return 0
    else: 
        return -1
    
def remove_condor_job(job_id):
    os.system("condor_rm {}".format(job_id))
    
def release_condor_job(job_id):
    ret = os.system("condor_release {}".format(job_id))
    return ret
    
def get_pids(daemon):
    out = shell_output("ps -ef | grep " + daemon)
    lines = out.split("\\n")
    daemons = list()
    for l in lines[:-3]:
        cont = [p for p in l.split(" ") if p!='']
        #print(cont)
        if 'python' in cont or 'python3' in cont:
            pid = cont[1] 
            if not pid.isnumeric():
                raise AssertionError("pid is not a number")
            daemons.append(pid)
			   
    return daemons

def kill_process(pid):
    ret = os.system("kill -9 {}".format(pid))  
#    if ret != 0:
#        # TODO: don't kill, update logs
#        raise Exception('Could not kill daemon')
    return ret
       
def get_job_daemons(job_daemon):
    out = shell_output("ps -ef | grep " + job_daemon)
    lines = out.split("\\n")
    daemons = dict()
    for l in lines[:-3]:
        cont = [p for p in l.split(" ") if p!='']
        if 'python' in cont or 'python3' in cont:
            pid = cont[1]
            if not pid.isnumeric():
                raise AssertionError("pid is not a number") 
            wdir = cont[-1]
            if '/rungate.' not in wdir:
                raise AssertionError("got wrong or non-existing working directory")
            daemons[wdir] = pid
    
    return daemons
    
def get_jobs_status():
    jobs_status = dict()
    out = shell_output("condor_q -all -wide")
    lines = out.split("\\n")
    header = [i for i in lines if i.startswith('OWNER')][0]
    start_job_lines = lines.index(header)+1
    #job = lines[start_job_lines]
    for job in lines[start_job_lines:]:
        if job == '': break
        job_info = [j for j in job.split(" ") if j!='']
        #job_id = job_info[-1].split(".")[0]
        job_id = job_info[1].split('+')[-1]
        jobs_status[job_id] = dict()
        jobs_status[job_id]["IDs"] = job_info[-1].split(".")[1]
        jobs_status[job_id]["RUN"] = job_info[5]
        jobs_status[job_id]["IDLE"] = job_info[6]
        jobs_status[job_id]["TOTAL"] = job_info[7]
        jobs_status[job_id]["DONE"] = job_info[4]
        if 'HOLD' in header:
            jobs_status[job_id]["HOLD"] = job_info[7]
        else:
            jobs_status[job_id]["HOLD"] = "_"
        
    return jobs_status

def job_on_hold(all_jobs,job_id):
    if all_jobs[job_id]['HOLD']!='_':
        return True
    else:
        return False
    
def get_job_age(date_str,tformat):
    now = time.time()
    date = time.strptime(date_str, tformat)
    age = now - time.mktime(date)
    
    return age
    
def condor_id(run_dagmann):
    ret, out = shell_output_ret(run_dagmann)
    dagmann_id = out.splitlines()[1].split(" ")[-1][:-1]
    #dagmann_id = int(dagmann_id) # + 1
    
    return ret, dagmann_id

def zip_dir_tree(base_name,form,root_dir):
    shutil.make_archive(base_name,form,root_dir)
    
def clean_dir_tree(d):
    shutil.rmtree(d, ignore_errors=True)

def zip_files(destin_fname,original_file_vec):
    myzipfile = ZipFile(destin_fname, mode='a')
    for d in original_file_vec:
        myzipfile.write(d)
    myzipfile.close()       
    
def zip_and_clean_folder(destin_fname,original_dir,form='zip'):
    zip_dir_tree(destin_fname,form,original_dir)
    clean_dir_tree(original_dir)
    
def change_folder_access_rights(path,gid,recursive=True, new_mode = 0o775):
    # Note: new_mode value is specified in octal notation, where the first digit specifies
    # the access rights for the owner, the second digit specifies the access rights for the group,
    # and the third digit specifies the access rights for everyone else. Each digit is a combination
    # of the following values: 4 for read access, 2 for write access, and 1 for execute access.
    # The prefix 0o indicates an octal integer
    try:
        if not recursive or os.path.isfile(path):
            os.chown(path, -1, gid)
            os.chmod(path,new_mode)
        else:
            for root, dirs, files in os.walk(path):
                os.chown(root, -1, gid)
                os.chmod(root,new_mode)
                for item in dirs:
                    os.chown(os.path.join(root, item), -1, gid)
                    os.chmod(os.path.join(root, item),new_mode)
                for item in files:
                    os.chown(os.path.join(root, item), -1, gid)
                    os.chmod(os.path.join(root, item),new_mode)
    except OSError as e:
        raise e 
    
    
