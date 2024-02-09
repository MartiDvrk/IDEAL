#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#   Copyright (C): MedAustron GmbH, ACMIT Gmbh and Medical University Vienna
#   This software is distributed under the terms
#   of the GNU Lesser General  Public Licence (LGPL)
#   See LICENSE for further details
# -----------------------------------------------------------------------------


import sys,os
import logging
from glob import glob

import impl.dual_logging as duall
from impl.system_configuration import get_sysconfig, system_configuration
from impl.idc_details import IDC_details
from impl.idc_enum_types import MCStatType
from impl.job_executor import job_executor
from impl.hlut_conf import hlut_conf
from impl.version import version_info
import impl.dicom_functions as dcm
from job_control_daemon import check_accuracy_for_beam, dose_monitoring_config, update_user_logs, periodically_check_statistical_accuracy

#global logger

class ideal_simulation():  
    def __init__(self,username,RP_path,n_particles=0,uncertainty=0,time_limit=0,debug=False,score_on_full_CT=False,
                    beams_to_simulate=None,ct_protocol=None,phantom=None,nvoxels=None,beamline_override=None,
                    padding_material="",material_overrides=None,sysconfig="",n_cores=0, condor_memory = 0):
        # username
        self.username = username
        # dicom plan
        self.dicom_planfile = RP_path
        # stopping criteria (input at least one)
        self.number_of_primaries_per_beam = n_particles
        self.percent_uncertainty_goal = uncertainty
        self.time_limit_in_minutes = time_limit
        # optional
        self.beams = beams_to_simulate
        self.ctprotocol = ct_protocol
        self.phantom = phantom
        self.nvoxels = nvoxels
        self.padding_material = padding_material
        self.material_overrides = material_overrides
        self.score_dose_on_full_CT = score_on_full_CT
        self.beamline_override = beamline_override
        self.number_of_cores = n_cores
        self.condor_memory = condor_memory
        # Initialize simulation object with the given inputs
        self.current_details = self.create_sim_object()
        # otput dir
        self.outputdir = self.current_details.output_job
        self.jobId = self.outputdir.split("/")[-1]
        # Job configuration 
        self.cfg = None
        # Simulation statistics
        self.stats = list()
  
    def get_plan_roi_names(self):
         return self.current_details.roinames
                
    def get_plan_beam_names(self):
         return self.current_details.beam_names
     
    def get_plan_nvoxels(self):
         nx,ny,nz = self.current_details.GetNVoxels()
         return nx,ny,nz

    def get_plan_resolution(self):         
         sx,sy,sz = self.current_details.GetDoseResolution()
         return sx,sy,sz
    
    def create_sim_object(self):
        #want_logfile = "default"
        prefix="\n * "
        sysconfig = system_configuration.getInstance()
        self.sysconfig = sysconfig
        logfilename = self.username + '_' + __name__ + '_' + duall.timestamp()
        sysconfig.set_logger(logfilename) 
        sysconfig.override('username',self.username)
        sysconfig.override("log file path",logfilename)
        if self.condor_memory:
            sysconfig.override('condor memory request default [MB]',self.condor_memory)

        logger = sysconfig.logger
        all_phantom_specs      = sysconfig["phantom_defs"]
        all_override_materials = sysconfig['ct override list']
        njobs = sysconfig['number of cores'] if self.number_of_cores else self.number_of_cores
        if self.number_of_cores:
            sysconfig.override('number of cores',self.number_of_cores)
        #username = sysconfig["username"]
        material_overrides = dict()
        if self.material_overrides is None:
            logger.debug("did not get any material override")
        else:
            logger.debug("got material override(s):{}{}".format(prefix,prefix.join(self.material_overrides)))
            for override in self.material_overrides:
                if len(override.split(":")) != 2:
                    raise RuntimeError("got bad material override '{}': should be of the form 'ROI:MATERIAL'".format(override))
                    sys.exit(1)
                roi,mat = override.split(":")
                if roi in material_overrides.keys():
                    raise RuntimeError("got multiple material overrides for the same ROI '{}'".format(roi))
                    sys.exit(2)
                if mat.upper() not in [m.upper() for m in all_override_materials.keys()]:
                    raise RuntimeError("the material '{}' is unknown/unrecognized/unsupported/forbidden".format(mat))
                # it's OK, maybe. :)
                material_overrides[roi] = mat
        matdb = os.path.join(sysconfig['commissioning'], sysconfig['materials database'])
        logger.debug("material database is {}".format(matdb))
        ##############################################################################################################
        current_details = IDC_details()
        rp = str(self.dicom_planfile)
        current_details.SetPlanFilePath(str(self.dicom_planfile))
        if bool(self.padding_material):
            mat=self.padding_material
            if mat not in [m.upper() for m in all_override_materials.keys()]:
                raise RuntimeError("the padding material '{}' is unknown/unrecognized/unsupported/forbidden".format(mat))
            current_details.dosepad_material = mat
        if bool(self.score_dose_on_full_CT):
            current_details.score_dose_on_full_CT = True
        if bool(self.beamline_override):
            logger.info(f"beamline override {self.beamline_override}")
            current_details.beamline_override = self.beamline_override
            for beam in current_details.beams:
                beam.TreatmentMachineName = self.beamline_override
        else:
            logger.debug("NO beamline override")
        if self.phantom is not None:
            phantoms = [ spec for label,spec in all_phantom_specs.items() if self.phantom.lower() in label.lower() ]
            exact_phantoms = [ spec for label,spec in all_phantom_specs.items() if self.phantom.lower() == label.lower() ]
            if len(phantoms)==1:
                logger.info("Running plan with phantom {}".format(phantoms[0]))
                current_details.UpdatePhantomGEO(phantoms[0])
            elif len(exact_phantoms)==1:
                logger.info("Running plan with phantom {}".format(exact_phantoms[0]))
                current_details.UpdatePhantomGEO(exact_phantoms[0])
            else:
                logger.error("unknown or ambiguous phantom name '{0}', available phantoms: '{1}'".format(self.phantom, "\n".join([label for label,spec in all_phantom_specs.items()])))
                raise RuntimeError("unknown or ambiguous phantom name '{0}', available phantoms: '{1}'".format(self.phantom, "\n".join([label for label,spec in all_phantom_specs.items()])))
                sys.exit(4)
            if self.material_overrides is not None:
                raise RuntimeError("for a geometrical phantom (no CT) you cannot specify override materials")
        elif current_details.run_with_CT_geometry:
            logger.debug("Running plan with CT")
        else:
            logger.error("Plan requires a phantom geometry, please specify phantom with the -p option, see -P to get a list of available phantom options")
            raise RuntimeError("Plan requires a phantom geometry, please specify phantom.")
            sys.exit(3)
        if self.nvoxels:
            logger.info("got nvoxels override: {}".format(self.nvoxels))
            for idim,nvoxel in enumerate(self.nvoxels):
                if nvoxel<1:
                    raise RuntimeError("number of voxels should be positive, got nvoxels[{}]={}".format(idim,nvoxel))
                if nvoxel>1000:
                    raise RuntimeError("number of voxels should be less than or equal to 1000, got nvoxels[{}]={}".format(idim,nvoxel))
                current_details.UpdateDoseGridResolution(idim,nvoxel)
        if self.phantom is None:
            logger.debug("no phantom => we will try to use the CT")
            # if the user did not specify a protocol, then 'hlut_conf' will try to detect it automagically
            # TODO: maybe wrap the 'SetHLUT' call in a try-except block?
            current_details.SetHLUT(self.ctprotocol)
            for roi,mat in material_overrides.items():
                if roi not in current_details.roinames:
                    logger.error("unknown roi name '{}', the structure set used by the current plan consists of:{}{}".format(
                        roi,prefix,prefix.join(current_details.roinames)))
                    raise RuntimeError("unknown roi name '{}'".format(roi))
                current_details.SetHUOverride(roi,mat)
                
           
                #rsids = beam.RangeShifterIDs
                #rmids = beam.RangeModulatorIDs
#        else:
#                 11th Oct 2022: removed option to override RangeShifter and RM as it is only accesible by gui
#            #TODOmfa: check overrides; only possible from gui; 
#            beam.RangeShifterIDs  = current_details.RSOverrides.get(beam.Name,beam.RangeShifterIDs)
#            rmids = self.details.RMOverrides.get(beam.Name,beam.RangeModulatorIDs)
#            rsflag="(as PLANNED)" if rsids == beam.RangeShifterIDs else "(OVERRIDE)"
#            rmflag="(as PLANNED)" if rmids == beam.RangeModulatorIDs else "(OVERRIDE)"
        if self.beams is not None:
            selection = dict([(name,False) for name in current_details.beam_names])
            for name in self.beams:
                if name not in selection.keys():
                    logger.error("{} is not a beam name in plan {}; known beams:{}{}".format(name,rp,prefix,prefix.join(current_details.beam_names)))
                    raise RuntimeError("{} is not a beam name in plan {}; known beams:{}{}".format(name,rp,prefix,prefix.join(current_details.beam_names)))
                    sys.exit(4)
                selection[name]=True
            current_details.SetBeamSelection(selection)
        statset=False
        if self.percent_uncertainty_goal > 0:
            logger.debug("simulation goal is {} % average uncertainty".format(self.percent_uncertainty_goal))
            current_details.SetStatistics(MCStatType.Xpct_unc_in_target,self.percent_uncertainty_goal)
            statset=True
        if self.time_limit_in_minutes > 0:
            logger.debug("simulation goal is {} minutes per job, for {} jobs.".format(self.time_limit_in_minutes,njobs))
            current_details.SetNJobs(njobs)
            current_details.SetStatistics(MCStatType.Nminutes_per_job,self.time_limit_in_minutes)
            statset=True
        if self.number_of_primaries_per_beam > 0:
            logger.debug("simulation goal is {} ions per beam".format(self.number_of_primaries_per_beam))
            current_details.SetStatistics(MCStatType.Nions_per_beam,self.number_of_primaries_per_beam)
            statset=True
        if not statset:
            logger.error("at least positive simulation goal should be set")
            raise RuntimeError("at least one positive simulation goal should be set")
            sys.exit(1)
            
        return current_details
        
        
    def start_simulation(self):
        logger = self.sysconfig.logger
        jobexec = job_executor.create_condor_job_executor(self.current_details)
        self.workdir = jobexec.template_gate_work_directory
        ret, condor_id =jobexec.launch_subjobs()
        self.condor_id = condor_id
        self.submission_date = jobexec.submission_date
        self.settings = jobexec.settings
        if ret!=0:
            logger.error("Something went wrong when submitting the job, got return code {}".format(ret))
            sys.exit(ret)
        # set job configuration
        self.cfg = dose_monitoring_config(self.workdir,self.username,daemonize=False,uncertainty_goal_percent=self.percent_uncertainty_goal,
                                          minimum_number_of_primaries=self.number_of_primaries_per_beam,time_out_minutes=self.time_limit_in_minutes)
    
    def check_accuracy(self,sim_time_minutes,input_stop=False):
        cfg = self.cfg
        current_dict = dict()
        stop = False
        for beamname,dosemhd in zip(cfg.beamname_list,cfg.dose_mhd_list):
            status = f"RUNNING GATE FOR BEAM={beamname}"
            current_dict[beamname]=dict()
            print(f"checking {dosemhd} for beam={beamname}")
            dose_files = glob(os.path.join(cfg.workdir,"tmp","output.*.*",dosemhd))
            if len(dose_files) == 0:
                print(f"looks like simulation for {dosemhd} did not start yet (zero dose files)")
                continue

            dc = check_accuracy_for_beam(cfg,beamname,dosemhd,dose_files)
            current_dict[beamname]['n_particles']=dc.tot_n_primaries
            current_dict[beamname]['average uncertainty']=dc.mean_unc_pct
            current_dict['simulation time in minutes'] = sim_time_minutes
            self.stats.append(current_dict)
            
            tmsg = f"Tsim = {sim_time_minutes} minutes (timeout = {cfg.time_out_minutes} minutes)"
            nmsg = f"Nsim = {dc.tot_n_primaries} primaries (minimum = {cfg.min_num_primaries})"
            umsg = f"Average Uncertainty = {dc.mean_unc_pct} pct (goal = {dc.cfg.unc_goal_pct} pct)"
            msg = ""
            # Maybe the following logic tree can be compactified, but for now I prefer to spell it out very explicitly
            if sim_time_minutes > cfg.time_out_minutes > 0:
                stop = True
                msg = "STOP: time is up: " + tmsg
            elif cfg.min_num_primaries > 0:
                if dc.tot_n_primaries < cfg.min_num_primaries:
                    stop = False
                    msg = "CONTINUE: not yet enough primaries: " + nmsg
                elif dc.cfg.unc_goal_pct > 0:
                    if dc.mean_unc_pct < dc.cfg.unc_goal_pct:
                        stop = True
                        msg = "STOP: uncertainty goal reached: " + umsg
                    else:
                        stop = False
                        msg = "CONTINUE: uncertainty goal NOT YET reached: " + umsg
                else:
                    stop = True
                    msg = "STOP: desired number of primaries reached: " + nmsg
            elif dc.cfg.unc_goal_pct > 0:
                if dc.mean_unc_pct < dc.cfg.unc_goal_pct:
                    stop = True
                    msg = "STOP: uncertainty goal reached: " + umsg
                else:
                    stop = False
                    msg = "CONTINUE: uncertainty goal NOT YET reached: " + umsg
            else:
                stop = False
                msg = "CONTINUE: time out not yet reached: " + tmsg
            print(f"{dosemhd} {tmsg} {nmsg} {umsg}")
            print(msg)
            update_user_logs(cfg.user_cfg,status,section=beamname,changes={"job control daemon status":msg})
            if stop:
                self.soft_stop_simulation(cfg)
                
        if input_stop:
            self.soft_stop_simulation(cfg)
            
        return stop, self.stats
        
    def soft_stop_simulation(self,cfg):
        # create stop file
        print("Going to stop GATE simulation")
        for beamname,dosemhd in zip(cfg.beamname_list,cfg.dose_mhd_list):
            with open(os.path.join(cfg.workdir,"STOP_"+dosemhd),"w") as stopfd:
                stopfd.write("{msg}\n")
            cfg.dose_mhd_list.remove(dosemhd)
            cfg.beamname_list.remove(beamname)
    
    def start_job_control_daemon(self):
        syscfg = system_configuration.getInstance()
        ret=os.system( "{bindir}/job_control_daemon.py -l {username} -t {timeout} -n {minprim} -u {uncgoal} -p {poll} -d -w '{workdir}'".format(
        bindir=syscfg['bindir'],
        username=syscfg['username'],
        # DONE: change this into Nprim, Unc, TimeOut settings
        #goal=self.details.mc_stat_thr,
        timeout=self.current_details.mc_stat_thr[MCStatType.Nminutes_per_job],
        minprim=self.current_details.mc_stat_thr[MCStatType.Nions_per_beam],
        uncgoal=self.current_details.mc_stat_thr[MCStatType.Xpct_unc_in_target],
        poll=syscfg['stop on script actor time interval [s]'],  # 300s
        workdir=self.workdir))
        
    
    def periodically_check_accuracy(self,frequency):
        cfg = self.cfg
        cfg.polling_interval_seconds = frequency
        periodically_check_statistical_accuracy(cfg)


# Initialize sysconfig
def initialize_sysconfig(filepath = '', username = '', debug=False):
    sysconfig = None
    try:
        sysconfig = get_sysconfig(filepath = filepath, username = username, debug=False)
    except Exception as e:
        print(f"OOPS, sorry! Problems getting system configuration, error message = '{e}'")  
        
    return sysconfig  
            
# Functions to enable queries
def verify_dicom_input_files(dicom_planfile):
    return dcm.dicom_files(dicom_planfile).check_all_dcm()

def get_version():
    return version_info
    
def list_available_phantoms():
    sysconfig = system_configuration.getInstance()
    all_phantom_specs = sysconfig["phantom_defs"]
    return all_phantom_specs.keys()
    
def get_ct_protocols():
    all_hluts=hlut_conf.getInstance()
    return all_hluts.keys()
    
def list_available_override_materials():
    sysconfig = system_configuration.getInstance()
    all_override_materials = sysconfig['ct override list']
    return all_override_materials.keys()
    
def list_available_beamline_names():
    sysconfig = system_configuration.getInstance()
    blmap=dict()
    for src_prop in glob(os.path.join(sysconfig['beamlines'],"*","*_*_source_properties.txt")):
        b=os.path.basename(src_prop)
        d=os.path.basename(os.path.dirname(src_prop))
        p=b[len(d)+1:-len("_source_properties.txt")]
        if d in blmap:
            blmap[d].append(p)
        else:
            blmap[d] = [p]
            
    return blmap
    

        
