# -----------------------------------------------------------------------------
#   Copyright (C): MedAustron GmbH, ACMIT Gmbh and Medical University Vienna
#   This software is distributed under the terms
#   of the GNU Lesser General  Public Licence (LGPL)
#   See LICENSE for further details
# -----------------------------------------------------------------------------

import pydicom
import itk
import configparser
import os, sys
import re
import numpy as np
from datetime import datetime
from impl.dual_logging import timestamp
from impl.idc_enum_types import MCStatType, MCPriorityType
from impl.gate_hlut_cache import generate_hlut_cache, hlut_cache_dir
from impl.system_configuration import system_configuration
from impl.hlut_conf import hlut_conf
from utils.bounding_box import bounding_box
from utils.roi_utils import region_of_interest, list_roinames
from utils.ct_dicom_to_img import ct_image_from_dicom
from utils.beamset_info import beamset_info
from utils.crop import crop_image
from impl.dicom_dose_template import write_dicom_dose_template
import logging
logger=logging.getLogger(__name__)

class IDC_details:
    def __init__(self,guimain=None):
        self._gui_main = guimain
        self.NA = "NA"
        self.subscribers = []
        self.ctphantom_subscribers = []
        self.ctprotocol_subscribers = []
        self.Reset()
    def Reset(self):
        syscfg = system_configuration.getInstance()
        self._HaveCT = False
        self._CT = True
        self._NeedDosePadding = False
        self._PHANTOM = False
        self._PhantomISOmm = dict()
        self._phantom_specs = None
        self.gamma_key = "gamma index parameters dta_mm dd_percent thr_percent def"
        self.gamma_parameters = syscfg[self.gamma_key]
        self.do_gamma = bool(self.gamma_parameters)
        self._beamline_override = None
        self.HUoverride = dict()
        self.humin = np.int16(-1024)
        self.humax = np.int16(3000)
        self.beam_selection = dict()
        self.rp_filepath = None
        self.rp_dataset = None
        self.ct_info = None
        self._ctprotocol_name = None
        self._mass_mhd = ""
        self.njobs = syscfg['number of cores']
        #self.mc_stat_q = MCStatType.Nions_per_beam
        #self.mc_stat_thr = 1000000
        self.mc_stat_thr = list(MCStatType.default_values)
        self.dosegrid_spacing = np.ones(3,dtype=float)
        self.dosegrid_nvoxels = np.ones(3,dtype=int)
        self.dosegrid_size = np.ones(3,dtype=float)
        self.dosegrid_center = np.zeros(3,dtype=float)
        self.dosegrid_origins = dict()
        self.dosemask = None
        self._score_dose_on_full_CT = False
        self.dosepad_material = "G4_AIR" # if dose grid sticks out of CT, pad the CT with this material
        self.bs_info = None
        self.rd_plan_info = None
        self.structure_set = None
        self.structure_set_filename = None
        self.roinames = list()
        self.roinumbers = list()
        self.roicontoursets = list()
        self.roitypes = list()
        self.external_roiname = ""
        self._tmpdir_job = None
        self.output_job = None
        self.output_job_2nd = None
        self._priority = MCPriorityType.condor_priorities[MCPriorityType.labels.index("Normal")]
        self._rs_overrides = dict()
        self._rm_overrides = dict()
    def override_gamma_parameters(self,*args):
        if 0 == len(args):
            self.do_gamma = False
        else:
            # TODO: check that these are actually four floats?
            self.gamma_parameters = args
    def set_job_dirs(self,jobname):
        syscfg = system_configuration.getInstance()
        self._tmpdir_job = os.path.join(syscfg['tmpdir jobs'],jobname)
        os.mkdir(self._tmpdir_job)
        self.output_job = os.path.join(syscfg['first output dicom'],jobname)
        os.mkdir(self.output_job)
        self.output_job_2nd = "" if not syscfg['second output dicom'] else os.path.join(syscfg['second output dicom'],jobname)
    #@property
    #def mc_stat_type(self):
    #    return self.mc_stat_q
    @property
    def priority(self):
        return self._priority
    @priority.setter
    def priority(self,newvalue):
        assert(type(newvalue)==int)
        assert(newvalue in MCPriorityType.condor_priorities)
        self._priority = newvalue
    @property
    def ctprotocol_name(self):
        return self._ctprotocol_name
    def SetPlanFilePath(self,rpfilepath):
        if not bool(rpfilepath):
            return
        if not os.path.exists(rpfilepath):
            logger.error("got non-existing RP DICOM file path {}".format(rpfilepath))
        self.Reset()
        #self.rp_filepath = os.path.realpath(rpfilepath)
        #self.rp_filepath = rpfilepath
        self.rp_filepath = os.path.join(os.path.realpath(os.path.dirname(rpfilepath)),os.path.basename(rpfilepath))
        logger.debug("going to get new file {}".format(rpfilepath))
        try:
            self.rp_dataset = pydicom.dcmread(self.rp_filepath)
            self.bs_info = beamset_info(self.rp_filepath)
            self._PhantomISOmm = dict([(beamname,np.array(self.bs_info[beamname].IsoCenter)) for beamname in self.bs_info.beam_names])
            self.beam_selection = dict([(name,True) for name in self.bs_info.beam_names])
        except IOError as fnfe:
            logger.error("OOPSIE: {}".format(fnfe))
            self.Reset()
            raise
        # try to get CT image
        rpdir = os.path.dirname(self.rp_filepath)
        try:
            patient_name = self.bs_info.patient_info["Patient Name"]
            badchars=re.compile("[^a-zA-Z0-9_]")
            sanitized_patient_name = re.sub(badchars,"_",patient_name)
            jobname="_".join([self.username,sanitized_patient_name,self.bs_info.sanitized_plan_label,timestamp()])
            logger.debug("jobname is {}".format(jobname))
            logger.debug("RP plan directory is {}".format(rpdir))
            self.set_job_dirs(jobname)
            logger.debug("getting plan dose info")
            self.rd_plan_info = self.bs_info.plan_dose
            logger.debug("got plan dose info")
            ss_ref_uid = self.rp_dataset.ReferencedStructureSetSequence[0].ReferencedSOPInstanceUID
            logger.debug("going to try to find the file with structure set with UID '{}'".format(ss_ref_uid))
            nskip=0
            ndcmfail=0
            nwrongtype=0
            nfiles=len([s for s in os.listdir(rpdir)])
            for s in os.listdir(rpdir):
                if s[-4:].lower() != '.dcm':
                    nskip+=1
                    logger.debug("no .dcm suffix: {}".format(s))
                    continue
                try:
                    logger.debug(s)
                    ds = pydicom.dcmread(os.path.join(rpdir,s))
                    dcmtype = ds.SOPClassUID.name
                except:
                    ndcmfail+=1
                    continue
                if dcmtype == "RT Structure Set Storage" and ss_ref_uid == ds.SOPInstanceUID:
                    logger.debug("found structure set for CT: {}".format(s))
                    ct_series_uid = ds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].SeriesInstanceUID
                    self.structure_set = ds
                    self.structure_set_filename = s
                    break
                else:
                    nwrongtype+=1
                    logger.debug("rejected structure set for CT: {}".format(s))
                    logger.debug("because it as a wrong SOP class ID: {}".format(dcmtype))
                    logger.debug("AND/OR because it has the wrong SOP Instance UID: {} != {}".format(ds.SOPInstanceUID,ss_ref_uid))
            if self.structure_set is None:
                raise RuntimeError("could not find structure set with UID={}; skipped {} with wrong suffix, got {} with 'dcm' suffix but pydicom could not read it, got {} with wrong class UID and/or instance UID. It could well be that this is a commissioning plan without CT and structure set data.".format(ss_ref_uid,nskip,ndcmfail,nwrongtype))
            self.ct_info = ct_image_from_dicom(rpdir,uid=ct_series_uid)
            logger.debug("image spacing is {}".format(self.ct_info.img.GetSpacing()))
            logger.debug("image size is {}".format(self.ct_info.img.GetLargestPossibleRegion().GetSize()))
            logger.debug("image origin is {}".format(self.ct_info.img.GetOrigin()))
            # the following code should be encapsulated in a "roi_info" class/object
            logger.debug("checking out structure set with {} ROIs".format(len(self.structure_set.StructureSetROISequence)))
            for i,roi in enumerate(self.structure_set.StructureSetROISequence):
                try:
                    logger.debug("{}. ROI number {}".format(i,roi.ROINumber))
                    logger.debug("{}. ROI name   {}".format(i,roi.ROIName))
                    roinumber = str(roi.ROINumber) # NOTE: roi numbers are *strings*
                    roiname = str(roi.ROIName)
                    contourset = None
                    roitype = None
                    if i<len( self.structure_set.ROIContourSequence ):
                        # normally this works
                        ci = self.structure_set.ROIContourSequence[i]
                        if str(ci.ReferencedROINumber) == roinumber:
                            contourset = ci
                    if not contourset:
                        logger.debug("(nr={},name={}) looks like this is a messed up structure set...".format(roinumber,roiname))
                        for ci in self.structure_set.ROIContourSequence:
                            if str(ci.ReferencedROINumber) == roinumber:
                                logger.debug("(nr={},name={}) contour found, phew!".format(roinumber,roiname))
                                contourset = ci
                                break
                    if not contourset:
                        logger.warn("ROI nr={} name={} does not have a contour, skipping it".format(roinumber,roiname))
                    if i<len( self.structure_set.RTROIObservationsSequence ):
                        # normally this works
                        obsi = self.structure_set.RTROIObservationsSequence[i]
                        if str(obsi.ReferencedROINumber) == roinumber:
                            roitype = str(obsi.RTROIInterpretedType)
                    if not roitype:
                        logger.debug("(nr={},name={}) looks like this is a messed up structure set...".format(roinumber,roiname))
                        for obsi in self.structure_set.RTROIObservationsSequence:
                            if str(obsi.ReferencedROINumber) == roinumber:
                                roitype = str(obsi.RTROIInterpretedType)
                                logger.debug("(nr={},name={}) type={} found, phew!".format(roinumber,roiname,roitype))
                                break
                    if not roitype:
                        logger.warn("ROI nr={} name={} does not have a type, skipping it".format(roinumber,roiname))
                    if bool(roitype) and bool(contourset):
                        self.roinumbers.append(roinumber)
                        self.roinames.append(roiname)
                        self.roicontoursets.append(contourset)
                        self.roitypes.append(roitype)
                except Exception as e:
                    logger.error("something went wrong with {}th ROI in the structure set: {}".format(i,e))
                    logger.error("skipping that for now, keep fingers crossed")
            dose_roinr = str(self.bs_info.target_ROI_number)
            dose_roiname = "NOT FOUND" if dose_roinr not in self.roinumbers else self.roinames[self.roinumbers.index(dose_roinr)]
            self.bs_info.target_ROI_name = dose_roiname
            iext = [i for i,typ in enumerate(self.roitypes) if typ.upper() == "EXTERNAL"]
            if len(iext) == 0:
                raise RuntimeError("Missing external ROI for CT image.")
            if len(iext) > 1:
                logger.warn("More than 1 external ROI for CT image: {}".format(
                            ", ".join(["'{}'".format(self.roinames[i]) for i in iext]) ) )
                logger.warn("Picking the first one, cross thumbs...")
            self.external_roiname = self.roinames[iext[0]]
            logger.debug(f"Successfully read plan file {rpfilepath}, found image and structure info.")
            self._HaveCT = True
            self._CT = True
            self._PHANTOM = False
            self._phantom_specs = None
            self.SetDefaultDoseGridSettings()
        except Exception as e:
            logger.info("Failed to find structure set and CT information: {}".format(e))
            self._HaveCT = False
            self._CT = False
            self._PHANTOM = True
            self._phantom_specs = None
            #self.Reset()
        if self.bs_info is not None:
            for s in self.subscribers:
                logger.debug("updating widget of type {}".format(type(s)))
                s.PlanUpdate()
        self.SetGeometry(0 if self.have_CT else 1)
        if self._gui_main:
            self._gui_main.update()
    def calculate_ram_request_mb(self,beamname):
        syscfg = system_configuration.getInstance()
        mb_min = syscfg['condor memory request minimum [MB]']
        mb_max = syscfg['condor memory request maximum [MB]']
        ## MFA 12/21/2022
        mb_default = syscfg['condor memory request default [MB]']
        geo="ct" if self._CT else "phantom"
        if self._CT:
            # this is silly, needs cleanup/revisiting
            ct_bb_nvoxels = self.ct_bb.volume/np.prod(self.ct_info.voxel_size)
            dose_nvoxels = ct_bb_nvoxels
        else:
            dose_nvoxels = np.prod(self.dosegrid_nvoxels)
        radtype="proton" if "PROTON" == self.bs_info[beamname].RadiationType.upper() else "carbon"
        logger.debug("going to use memory fit for radtype={} and geo={}".format(radtype,geo))
        mb_fit = syscfg['condor memory fit {} {}'.format(radtype,geo)]
        mb_guess = 0.
        for k,v in mb_fit.items():
            if k=='offset':
                mb_guess+=v
            elif k=='dosegrid':
                mb_guess+=v*dose_nvoxels
            elif k=='ct':
                mb_guess+=v*ct_bb_nvoxels
            elif k=='nspots':
                mb_guess+=v*self.bs_info[beamname].nspots
            else:
                logger.error("don't know with memory fit item key='{}' value='{}'".format(k,v))
        #print(f"{mb_guess=}")
        mb_guess = mb_default
        mb = min(mb_max,max(mb_min,mb_guess))
        logger.debug("RAM fit gives guess {} MB for this beam, minmb={} and maxmb={}, so {} is used".format(mb_guess,mb_min,mb_max,mb))
        return mb
    def DoseGridSticksPartlyOutsideOfCTVolume(self):
        return self._NeedDosePadding
    def set_gui_main(self,guimain):
        self._gui_main = guimain
    @property
    def username(self):
        syscfg = system_configuration.getInstance()
        return syscfg['username']
    @property
    def tmpdir_job(self):
        return self._tmpdir_job
    @property
    def input_dicom(self):
        syscfg = system_configuration.getInstance()
        return syscfg['input dicom']
    @property
    def logging_directory(self):
        syscfg = system_configuration.getInstance()
        return syscfg['logging']
    @property
    def min_dose_res(self):
        if self.run_with_CT_geometry:
            return self.ct_info.voxel_size
        syscfg = system_configuration.getInstance()
        return syscfg['minimum dose grid resolution [mm]']
    @property
    def max_dose_res(self):
        syscfg = system_configuration.getInstance()
        return syscfg['maximum dose grid resolution [mm]']
    @property
    def dose_step(self):
        syscfg = system_configuration.getInstance()
        return syscfg['dose grid resolution stepsize [mm]']
    @property
    def beam_names(self):
        return self.bs_info.beam_names if self.bs_info else list()
    @property
    def beam_numbers(self):
        return self.bs_info.beam_numbers if self.bs_info else list()
    @property
    def score_dose_on_full_CT(self):
        return self._score_dose_on_full_CT
    @score_dose_on_full_CT.setter
    def score_dose_on_full_CT(self,override):
        assert((override is True) or (override is False))
        self._score_dose_on_full_CT = override
    @property
    def have_CT(self):
        return self._HaveCT
    @property
    def run_with_CT_geometry(self):
        assert(self._CT != self._PHANTOM)
        return self._CT
    def PhantomISOinMM(self,beamname):
        return self._PhantomISOmm[beamname].copy()
    @property
    def have_dose_grid(self):
        logger.debug("CT is {}".format(self._CT))
        logger.debug("phantom specs is {}".format(self._phantom_specs is not None))
        return self._CT or self._phantom_specs is not None
    @property
    def PhantomSpecs(self):
        return self._phantom_specs
    def SetHLUT(self,kw=None):
        all_hluts = hlut_conf.getInstance()
        need_update = False
        # TODO: the "hlut match" methods may throw a KeyError, should we catch that here or leave that ot the caller?
        if kw:
            ctprotocol = all_hluts.hlut_match_keyword(kw)
        elif self.have_CT:
            ctprotocol = all_hluts.hlut_match_dicom(self.ct_info.slices[0])
        else:
            logger.warn("attempt to set HLUT without keyword nor CT?")
            return
        if self._ctprotocol_name == ctprotocol:
            logger.debug(f"HLUT was and remains '{ctprotocol}'")
        else:
            self._ctprotocol_name = ctprotocol
            logger.debug(f"HLUT protocol set to '{ctprotocol}'")
            assert( self.ctprotocol_name)
            self.UpdateHURange()
            for widget in self.ctprotocol_subscribers:
                widget.update_ctprotocol()
    @property
    def beamline_override(self):
        return self._beamline_override
    @beamline_override.setter
    def beamline_override(self,override):
        """
        For all beams in the beamset, use this beam model instead of the one specified in the treatment plan.
        TODO: allow beam-specific overrides?
        """
        syscfg = system_configuration.getInstance()
        if syscfg['role'] == 'clinical':
            raise RuntimeError("As user {} you have a 'clinical' role, therefore you cannot override the beamline model!".format(syscfg['username']))
        elif os.path.isdir(os.path.join(syscfg['beamlines'],override)):
            logger.info(f"storing beamline override to '{override}'")
            self._beamline_override = override
        else:
            raise ValueError(f"{override} not a supported beamline model name")
    def HaveHUOverride(self,roiname):
        answer = roiname in self.HUoverride.keys()
        return answer
    def RemoveHUOverride(self,roiname):
        if roiname in self.HUoverride.keys():
            oldvalue=self.HUoverride.pop(roiname)
            logger.debug('REMOVED material={} override for ROI={}'.format(oldvalue,roiname))
        else:
            logger.warn('got request to remove HU override for ROI, but there is no such override.'.format(roiname))
    def SetHUOverride(self,roiname,matname):
        assert(roiname in self.roinames)
        syscfg = system_configuration.getInstance()
        assert(matname in syscfg['ct override list'].keys())
        self.HUoverride[roiname] = matname
    def UpdateHURange(self):
        all_hluts = hlut_conf.getInstance()
        assert(self.ctprotocol_name)
        cache_hu2mat,dummy = all_hluts[self.ctprotocol_name].get_hu2mat_files()
        if not os.path.exists(cache_hu2mat):
            raise RuntimeError(f"BUG: Could not find cached HU-to-material file {cache_hu2mat}. Please complain to developer.")
        logger.debug("going to read HU range from cache file {}".format(cache_hu2mat))
        with open(cache_hu2mat,"r") as hu2mat:
            lines = hu2mat.readlines()
            # assume table is sorted properly
            humin=float(lines[0].strip().split()[0])
            humax=float(lines[-1].strip().split()[1])
        assert(humin>-2000)
        assert(humin<0)
        assert(humax>0)
        assert(humax<=32000) # max value for signed short is 32767, but we need some room for material overrides.
        self.humin=humin
        self.humax=humax
        logger.debug("HU min (~air) is {}, HU max is {}".format(self.humin,self.humax))
    def GetROIBoundingBox(self):
        logger.debug("bounding box based on {} ROIs in HU override list".format(len(self.HUoverride.keys())))
        syscfg = system_configuration.getInstance()
        dosegrid_air_margin = float(syscfg["air box margin [mm]"])
        self.roi_bb = bounding_box()
        for roiname in set(list(self.HUoverride.keys())+[self.external_roiname]):
            if roiname == "!HUMAX":
                logger.debug("ignoring !HUMAX")
                continue
            roi_id,margin = (roiname[1:],dose_air_margin) if roiname[0]=="!" else (roiname,0.)
            roi = region_of_interest(ds=self.structure_set,roi_id=roi_id)
            logger.debug("merging roi_bb with id={} BB={}, using margin={}".format(roi_id,roi.bb,margin))
            self.roi_bb.should_contain(roi.bb.mincorner-margin)
            self.roi_bb.should_contain(roi.bb.maxcorner+margin)
        logger.debug("ROI bounding box is now {}".format(self.roi_bb))
    def GetCTBoundingBox(self):
        if self.score_dose_on_full_CT:
            self.ct_bb = bounding_box(img=self.ct_info.img)
        else:
            logger.debug("Get combined bounding box containing the dose matrix and all ROIs for material overrides (including external)")
            logger.debug("(Note: if the dose matrix sticks outside of the CT, then so will the bounding box!)")
            self.ct_bb = bounding_box()
            self.ct_bb.should_contain(self.dosegrid_center-0.5*self.dosegrid_size)
            self.ct_bb.should_contain(self.dosegrid_center+0.5*self.dosegrid_size)
            logger.debug("BB around dose matrix, with air margin: {}".format(self.ct_bb))
            self.GetROIBoundingBox()
            self.ct_bb.should_contain(self.roi_bb.mincorner)
            self.ct_bb.should_contain(self.roi_bb.maxcorner)
        return bounding_box(bb=self.ct_bb)
    def WritePreProcessingConfigFile(self,submitdir,mhd,hu2mat,hudensity):
        syscfg = system_configuration.getInstance()
        logger.debug("going to write preprocessing config file")
        parser=configparser.RawConfigParser()
        parser.optionxform = lambda option : option
        logger.debug("dicom section")
        parser.add_section('dicom')
        parser['dicom'].update({"directory":os.path.dirname(self.rp_filepath)})
        parser['dicom'].update({"RSfile":self.structure_set_filename})
        parser['dicom'].update({"CTuid":self.ct_info.uid})
        ct_bb = self.GetCTBoundingBox()
        logger.debug("ct bounding box section")
        parser.add_section('ct bounding box')
        parser['ct bounding box'].update({"min corner":" ".join([str(v) for v in ct_bb.mincorner])})
        parser['ct bounding box'].update({"max corner":" ".join([str(v) for v in ct_bb.maxcorner])})
        logger.debug("output section")
        parser.add_section('mhd files')
        #parser['mhd files'].update({'mhd resized':"ct_resized.mhd"})
        #ct_orig_mhd = os.path.join(os.path.dirname(mhd),"ct_orig.mhd")
        ct_orig_mhd = os.path.join(submitdir,"ct_orig.mhd")
        parser['mhd files'].update({'mhd with original ct':ct_orig_mhd})
        parser['mhd files'].update({'mhd with overrides':mhd})
        self._mass_mhd = mhd.replace(".mhd","_mass.mhd")
        parser['mhd files'].update({'mhd with mass':self._mass_mhd})
        self.dosemask = os.path.join(os.path.dirname(mhd),"dose_grid_mask.mhd")
        parser['mhd files'].update({'mhd dose grid mask':self.dosemask})
        #parser['mhd files'].update({'tar file':'ct.tar.gz'})
        parser.add_section('HUoverride')
        override_materials = set(["G4_AIR","G4_WATER"]+[self.dosepad_material]+list(self.HUoverride.values()))
        humap=dict()
        parser.add_section('density')
        # the density curve is needed for the mass image (for mass-weighted resampling)
        parser['density']["hlut_path"] = hudensity
        assert(os.path.exists(hu2mat))
        with open(hu2mat,"a") as fh:
            # Create an addition to the existing copy of the HU2mat table, continue where that table ends.
            # Choose HU intervals of length 2, and use the central values of those intervals for the overrides.
            # Hopefully this avoids any material confusion due to rounding & boundary issues...
            h1,hu,h2 = self.humax, int(np.ceil(self.humax+1)), int(np.ceil(self.humax+2))
            for material in override_materials:
                parser['density'][str(hu)] = str(syscfg['ct override list'][material])
                humap[material] = str(hu)
                fh.write("{} {} {}\n".format(h1,h2,material))
                h1=h2
                hu=h1+1
                h2=h1+2
        parser['HUoverride'].update({"!HUMAX":str(int(self.humax))})
        parser['HUoverride'].update({"!DOSEPAD":str(int(humap[self.dosepad_material]))})
        parser['HUoverride'].update({"!"+self.external_roiname:humap["G4_AIR"]})
        parser['HUoverride'].update(dict([(k,humap[v]) for k,v in self.HUoverride.items()]))
        ########################
        #resized_ct_img = self.CropAndPadCTImage(humap["G4_AIR"])
        #itk.imwrite(resized_ct_img,"ct_resized.mhd")
        logger.debug("writing original CT to {}".format(ct_orig_mhd))
        itk.imwrite(self.ct_info.img,ct_orig_mhd)
        ########################
        parser.add_section('dose grid')
        assert(len(self.dosegrid_size)==3)
        parser['dose grid'].update({'dose grid center':" ".join([str(v) for v in self.dosegrid_center])})
        parser['dose grid'].update({'dose grid size':" ".join([str(v) for v in self.dosegrid_size])})
        parser['dose grid'].update({'dose grid nvoxels':" ".join([str(v) for v in self.dosegrid_nvoxels])})
        parser['dose grid'].update({'dose grid air margin': str(syscfg["air box margin [mm]"])})
        with open(os.path.join(submitdir,"preprocessor.cfg"),"w") as fp:
            parser.write(fp)
        if self.score_dose_on_full_CT:
            self.ct_nvoxels = self.ct_info.size
        else:
            ibbmin,ibbmax = ct_bb.indices_in_image(self.ct_info.img)
            self.ct_nvoxels=ibbmax-ibbmin
        return ct_bb,self.ct_nvoxels
    def WritePostProcessingConfigFile(self,submitdir,qspecs,plan_dose_file="idc-PLAN"):
        assert(len(qspecs)>0)
        syscfg = system_configuration.getInstance()
        parser=configparser.RawConfigParser()
        parser.optionxform = lambda option : option
        parser['DEFAULT']["run gamma analysis"]       = str(syscfg["run gamma analysis"])
        parser['DEFAULT']["debug"]       = str(syscfg["debug"])
        parser['DEFAULT']["first output dicom"]       = self.output_job
        parser['DEFAULT']["second output dicom"]      = self.output_job_2nd
        parser['DEFAULT']["nFractions"]               = str(self.bs_info.Nfractions)
        parser['DEFAULT']["write mhd unscaled dose"]  = str(syscfg["write mhd unscaled dose"])
        parser['DEFAULT']["write mhd scaled dose"]    = str(syscfg["write mhd scaled dose"])
        parser['DEFAULT']["write mhd physical dose"]  = str(syscfg["write mhd physical dose"])
        parser['DEFAULT']["write mhd rbe dose"]       = str(syscfg["write mhd rbe dose"])
        parser['DEFAULT']["write dicom physical dose"]= str(syscfg["write dicom physical dose"])
        parser['DEFAULT']["write dicom rbe dose"]     = str(syscfg["write dicom rbe dose"])
        parser['DEFAULT']["write dicom re alanine"]     = str(syscfg["write dicom re alanine"])
        parser['DEFAULT']["write unresampled dose"]   = "yes" if self.score_dose_on_full_CT else "no"
        parser['DEFAULT']["dose grid size"]           = " ".join([str(val) for val in self.dosegrid_size])
        parser['DEFAULT']["dose grid resolution"]     = " ".join([str(val) for val in self.dosegrid_nvoxels])
        if self._CT:
            parser['DEFAULT']["sim dose resolution"]  = " ".join([str(val) for val in self.ct_nvoxels])
        else:
            parser['DEFAULT']["sim dose resolution"]  = parser['DEFAULT']["dose grid resolution"]
        if syscfg["write mhd plan dose"]:
            parser['DEFAULT']["mhd plan dose"]        = plan_dose_file + ".mhd"
        else:
            parser['DEFAULT']["mhd plan dose"]        = ""
        if syscfg["write dicom plan dose"]:
            template_path="plan_dose_template.dcm"
            write_dicom_dose_template(self.rp_dataset,"PLAN",template_path,phantom=self._PHANTOM)
            parser['DEFAULT']["dicom plan dose"]      = plan_dose_file + ".dcm"
            parser['DEFAULT']["plan dcm template"]    = template_path
        else:
            parser['DEFAULT']["dicom plan dose"]      = ""
            parser['DEFAULT']["plan dcm template"]    = ""
        parser['DEFAULT']["mass mhd"]                 = self._mass_mhd
        if self.do_gamma:
            parser['DEFAULT'][self.gamma_key] = self.gamma_parameters
        if self.run_with_CT_geometry:
            parser['DEFAULT']["apply external dose mask"] = "yes" if syscfg["remove dose outside external"] else "no"
            parser['DEFAULT']["external dose mask"] = self.dosemask
        for beamname,qspec in qspecs.items():
            origname=qspec["origname"]
            beam = self.bs_info[origname]
            rbe = str(syscfg["rbe factor protons"]) if "PROTON" == self.bs_info[origname].RadiationType.upper() else str(1.0)
            parser.add_section(beamname)
            parser[beamname].update(qspec)
            assert(int(qspec["nJobs"])>0)
            # calculate corrected msw for beam (MFA, 8/17/23)
            def_msw_scaling=syscfg['msw scaling']["default"]
            dose_corr_key=(beam.TreatmentMachineName+"_"+beam.RadiationType).lower()
            params_msw_scaling = syscfg['msw scaling'].get(dose_corr_key,def_msw_scaling)
            conversion = lambda msw, energy : msw*np.polyval(params_msw_scaling,energy)
            beam.msw_conv_func = conversion
            nTPS = self.calc_msw_tot_beam(beam, conversion)
            parser[beamname].update({"nTPS":str(nTPS)})
            beamnr=self.bs_info[origname].number
            template_path="dose_template_beam_{}_{}.dcm".format(beamnr,beamname)
            write_dicom_dose_template(self.rp_dataset,beamnr,template_path,phantom=self._PHANTOM)
            parser[beamname].update({"dcm template":template_path})
            orig = self.dosegrid_origin if self.run_with_CT_geometry else np.zeros(3,dtype=float) - 0.5 * self.dosegrid_size + 0.5*self.dosegrid_spacing
            parser[beamname].update({"dose grid origin":" ".join([str(val) for val in orig])})
            dosetype="PHYSICAL" if rbe=="1.0" else "EFFECTIVE"
            if self.do_gamma and self.bs_info.have_tps_dose(origname,sumtype="BEAM",dosetype=dosetype):
                fpath = self.bs_info.tps_dose(origname,sumtype="BEAM",dosetype=dosetype).filepath
                parser[beamname]["path to reference dose image for gamma index calculation"] = fpath
        parser['DEFAULT']["RBE"]=rbe
        dosetype="PHYSICAL" if rbe=="1.0" else "EFFECTIVE"
        if self.do_gamma and self.bs_info.have_tps_dose(None,sumtype="PLAN",dosetype=dosetype):
            fpath = self.bs_info.tps_dose(None,sumtype="PLAN",dosetype=dosetype).filepath
            parser[beamname][f"path to reference {dosetype} plan dose image for gamma index calculation"] = fpath
        with open(os.path.join(submitdir,"postprocessor.cfg"),"w") as fp:
            parser.write(fp)
            
    '''
    function to scale the msw of a single beam according to the scaling factors
    defined in the config file. Same concept is applied when writing the plan txt file
    '''
    def calc_msw_tot_beam(self, beam, conversion = lambda x: x):
        new_msw_tot = 0
        for i,l in enumerate(beam.layers):
            #k_e = np.polyval(params_msw_scaling,l.energy)
            for k, spot in enumerate(l.spots):
                new_msw_tot += conversion(spot.msw,l.energy) 
        return new_msw_tot
    def WriteUserSettings(self,qspecs,ymd_hms,condordir):
        ####################
        logger.debug("Experimental feature: writing cfg file with user specifications and semi-minimal logging info.")
        parser=configparser.RawConfigParser()
        parser.optionxform = lambda option : option
        ####################
        syscfg = system_configuration.getInstance()
        ####################
        parser['DEFAULT']["username"] = self.username
        parser['DEFAULT']["status"] = "submitted"
        parser['DEFAULT']["date and time of last update"] = datetime.now().ctime()
        parser['DEFAULT']["submission program"] = sys.argv[0]
        parser['DEFAULT']["submission date and time"] = datetime.now().ctime()
        #parser['DEFAULT']["statistics goal quantity"] = MCStatType.guilabels[self.mc_stat_q]
        parser['DEFAULT']["statistics goal threshold values (timeout_minutes, min_n_primaries, unc_goal_pct)"] = " ".join([str(v) for v in self.mc_stat_thr])
        parser['DEFAULT']["TPS dicom plan file path"] = self.rp_filepath
        parser['DEFAULT']["condor submit directory"] = condordir
        parser['DEFAULT']["dose grid resolution"] = " ".join([str(v) for v in self.dosegrid_nvoxels])
        parser['DEFAULT']["selected beams"] = " ".join(["'{}'".format(beamname) for beamname,yes in self.beam_selection.items() if yes])
        parser['DEFAULT']["deselected beams"] = " ".join(["'{}'".format(beamname) for beamname,yes in self.beam_selection.items() if not yes])
        parser['DEFAULT']["log file path"] = syscfg["log file path"]
        parser['DEFAULT']["output directory on cluster (1st dicom output)"] = self.output_job
        parser['DEFAULT']["output directory on file share (2nd dicom output)"] =  str(self.output_job_2nd) if bool(self.output_job_2nd) else "(None)"
        ####################
        parser.add_section("Patient")
        parser["Patient"].update(self.bs_info.patient_info)
        ####################
        parser.add_section("Plan")
        parser["Plan"].update(self.bs_info.plan_info)
        ####################
        parser.add_section("BS")
        parser["BS"].update(self.bs_info.bs_info)
        ####################
        if self.run_with_CT_geometry:
            parser.add_section("CT")
            parser["CT"]["CT protocol"]=self.ctprotocol_name
            parser["CT"]["material to use for padding CT in case it does not include the dosegrid"]=self.dosepad_material
            parser["CT"].update(dict([(str(k),str(v)) for k,v in self.ct_info.meta_data.items()]))
            parser["CT"].update(dict([("HU override for '{}'".format(roiname),matname) for roiname,matname in self.HUoverride.items() if roiname!="!HUMAX"]))
            parser.add_section("Passive Elements")
            for beam in self.bs_info.beams:
                if not self.beam_selection[beam.Name]:
                    continue
                rsids = self.RSOverrides.get(beam.Name,beam.RangeShifterIDs)
                rmids = self.RMOverrides.get(beam.Name,beam.RangeModulatorIDs)
                rstxt = "(none)" if len(rsids)==0 else " ".join(rsids)
                rmtxt = "(none)" if len(rmids)==0 else " ".join(rmids)
                parser["Passive Elements"]["'{}' range shifters".format(beam.Name)] = rstxt
                parser["Passive Elements"]["'{}' range modulators".format(beam.Name)] = rmtxt
        else:
            parser.add_section("PHANTOM")
            parser["PHANTOM"].update(dict({"CT":"None" if not self.have_CT else str(self.ct_info._uid)}))
            parser["PHANTOM"].update(dict([(str(k),str(v)) for k,v in self._phantom_specs.meta_data.items()]))
            parser.add_section("Passive Elements")
            for beam in self.bs_info.beams:
                if not self.beam_selection[beam.Name]:
                    continue
                rsids = self.RSOverrides.get(beam.Name,beam.RangeShifterIDs)
                rmids = self.RMOverrides.get(beam.Name,beam.RangeModulatorIDs)
                rsflag="(as PLANNED)" if rsids == beam.RangeShifterIDs else "(OVERRIDE)"
                rmflag="(as PLANNED)" if rmids == beam.RangeModulatorIDs else "(OVERRIDE)"
                rstxt = "(none)" if len(rsids)==0 else " ".join(rsids)
                rmtxt = "(none)" if len(rmids)==0 else " ".join(rmids)
                parser["Passive Elements"]["'{}' range shifters".format(beam.Name)] = "{} {}".format(rsflag,rstxt)
                parser["Passive Elements"]["'{}' range modulators".format(beam.Name)] = "{} {}".format(rmflag,rmtxt)
        ####################
        for beamname,qspec in qspecs.items():
            origname=qspec["origname"]
            nTPS=self.bs_info[origname].mswtot
            parser.add_section(origname)
            parser[origname].update(qspec)
            parser[origname]["nTPS"]=str(nTPS)
            parser[origname]["sanitized beam name (e.g. used in name of main mac file)"]=str(beamname)
        fpath = os.path.join(self.output_job,"user_logs_{}.cfg".format(ymd_hms.translate(str.maketrans(": -","___"))))
        ####################
        parser.add_section("Logs")
        parser["Logs"]['UI logs'] = syscfg["log file path"]
        parser["Logs"]['preprocessor logs'] = os.path.join(condordir,"preprocessor.log")
        parser["Logs"]['GateRTion log directory'] = os.path.join(condordir,"logs")
        parser["Logs"]['postprocessor logs'] = os.path.join(condordir,"postprocessor.log")
        logger.debug("user settings file: {}".format(fpath))
        with open(fpath,"w") as fp:
            parser.write(fp)
        for proc_cfg in ["postprocessor.cfg","preprocessor.cfg"]:
            if os.path.exists(proc_cfg):
                procparser=configparser.RawConfigParser()
                procparser.optionxform = lambda option : option
                with open(proc_cfg,"r") as fp:
                    procparser.read_file(fp)
                procparser.add_section("user logs file")
                procparser['user logs file']['path'] = fpath
                with open(proc_cfg,"w") as fp:
                    procparser.write(fp)
        return fpath
    def TargetROIName(self):
        return self.bs_info.target_ROI_name
    def ExternalROIName(self):
        return self.external_roiname
    def GetDoseResolution(self):
        return (self.dosegrid_spacing).tolist() # mm
    @property
    def dosegrid_changed(self):
        return (self.dosegrid_nvoxels != self.def_dosegrid_nvoxels).any()
    #@property
    #def dosegrid_center(self):
    #    return self.dosegrid_origin - 0.5 * self.dosegrid_spacing + 0.5 * self.dosegrid_size
    @property
    def dosegrid_origin(self):
        return self.dosegrid_center + 0.5 * self.dosegrid_spacing - 0.5 * self.dosegrid_size
    def phantom_dosegrid_center(self,beamname):
        return -1*self._PhantomISOmm[beamname] - 0.5 * self.dosegrid_size[0] + 0.5*self.dosegrid_spacing[0]
    def GetDoseCenter(self):
        return self.dosegrid_center.tolist() # mm
    def GetDoseSize(self):
        return (self.dosegrid_size).tolist() # mm
    def GetNVoxels(self):
        return self.dosegrid_nvoxels.tolist()
    def Subscribe(self,widget,special=None):
        if not special:
            self.subscribers.append(widget)
        elif special=="CTPHANTOM":
            self.ctphantom_subscribers.append(widget)
        elif special=="CTPROTOCOL":
            self.ctprotocol_subscribers.append(widget)
    def UpdateDoseGridResolution(self,idim,newval):
        assert(idim in [0,1,2])
        assert(newval>0)
        new_nvoxel = int(newval)
        self.dosegrid_nvoxels[idim] = new_nvoxel
        old_spacing = self.dosegrid_spacing[idim]
        new_spacing = self.dosegrid_size[idim]/new_nvoxel
        #self.dosegrid_origin[idim] += 0.5*(new_spacing-old_spacing)
        self.dosegrid_spacing[idim] = new_spacing
        #self.dosegrid_spacing[idim] = newval
        #new_nvoxels = int(np.ceil(self.def_dosegrid_size[idim]/newval))
        #self.dosegrid_nvoxels[idim] = new_nvoxels
        #self.dosegrid_size[idim] = new_nvoxels*newval
        ## if new_nvoxels==1, then dosegrid_center==dosegrid_origin
        #self.dosegrid_origin[idim] = self.dosegrid_center[idim] - 0.5*(new_nvoxels-1)*newval
    def SetDefaultDoseGridSettings(self):
        logger.debug("going to set default dose grid settings")
        if self._CT:
            if self.rd_plan_info:
                self.def_dosegrid_nvoxels = self.rd_plan_info.nvoxels
                self.dosegrid_spacing = self.rd_plan_info.spacing
                self.def_dosegrid_origin =  self.rd_plan_info.origin
            else:
                self._NeedDosePadding = False
                logger.warn("OOPS no TPS dose available??? Falling back to CT geometry")
                #self.def_dosegrid_nvoxels = self.ct_info.size
                #self.dosegrid_spacing = self.ct_info.voxel_size
                #self.def_dosegrid_origin = self.ct_info.origin 
                logger.debug("going to get ROI bb")
                self.GetROIBoundingBox()
                logger.debug("got ROI bb: {}".format(self.roi_bb))
                cropped_ct = crop_image(self.ct_info.img,*self.roi_bb.indices_in_image(self.ct_info.img))
                logger.debug("got cropped image")
                self.def_dosegrid_nvoxels = np.array(itk.size(cropped_ct))
                logger.debug("got size")
                self.dosegrid_spacing = np.array(itk.spacing(cropped_ct))
                logger.debug("got spacing")
                self.def_dosegrid_origin = np.array(itk.origin(cropped_ct))
                logger.debug("got origin")
                # TODO: raise an error, so that the GUI can show an error/warning dialog?
            self.dosegrid_nvoxels = self.def_dosegrid_nvoxels.copy()
            self.def_dosegrid_size = self.dosegrid_nvoxels * self.dosegrid_spacing
            self.dosegrid_size = self.def_dosegrid_size.copy()
            self.dosegrid_center = self.def_dosegrid_origin-0.5*self.dosegrid_spacing+0.5*self.dosegrid_size
            if not self.score_dose_on_full_CT:
                bb_dose = bounding_box(xyz=np.stack(( self.dosegrid_center-0.5*self.dosegrid_size, self.dosegrid_center+0.5*self.dosegrid_size )))
                bb_ct = bounding_box(img=self.ct_info.img)
                self._NeedDosePadding = not bb_ct.encloses(bb_dose)
                logger.debug("CT {} {} DOSEGRID {}".format(bb_ct,"DOES NOT ENCLOSE" if self._NeedDosePadding else "ENCLOSES", bb_dose))
        elif self._PHANTOM:
            self._NeedDosePadding = False
            if self._phantom_specs is None:
                logger.error("this should not happen: dose grid settings for phantoms cannot be set if phantom is not yet defined")
                return
            self.def_dosegrid_size = self._phantom_specs.dose_grid_size
            self.dosegrid_size = self.def_dosegrid_size.copy()
            self.dosegrid_spacing = self._phantom_specs.dose_voxel_size
            self.def_dosegrid_nvoxels = self._phantom_specs.dose_nvoxels
            self.dosegrid_nvoxels = self.def_dosegrid_nvoxels.copy()
            self.dosegrid_origins = dict([(beamname,-1*self._PhantomISOmm[beamname] - 0.5 * self.dosegrid_size + 0.5*self.dosegrid_spacing) for beamname in self.beam_names])
            self.dosegrid_centers = dict([(beamname,-1*self._PhantomISOmm[beamname]) for beamname in self.beam_names])
        logger.debug("done setting default dose grid settings")
    @property
    def uid(self):
        if self.bs_info:
            return self.bs_info.uid
        return self.NA
    def GetPatientInfo(self):
        logger.debug("get PATIENT info")
        labels = list()
        values = list()
        if self.rp_dataset:
            info = self.bs_info.patient_info
            for a in beamset_info.patient_attrs:
                labels.append(a)
                values.append(info[a])
        return values,labels
    def GetPlanInfo(self):
        logger.debug("get PLAN info")
        labels = list()
        values = list()
        if self.bs_info:
            # for a in ["PlanName","SOPInstanceUID","OperatorsName","ReviewerName","ReviewDate","ReviewTime","ReferringPhysicianName","PlanIntent"]:
            #    value = getattr(self.rp_dataset,a) if hasattr(self.rp_dataset,a) else self.NA
            info = self.bs_info.plan_info
            for a in beamset_info.plan_attrs:
                labels.append(a)
                values.append(info[a])
        return values,labels
    def GetCTInfo(self):
        logger.debug("get CT info")
        labels = list()
        values = list()
        if self.ct_info:
            for k,v in self.ct_info.meta_data.items():
                labels.append(k)
                values.append(v)
        else:
            labels.append("CT")
            values.append("NOT AVAILABLE (COMMISSIONING PLAN?)")
        return values,labels
    def GetBeamSetInfo(self):
        logger.debug("get BEAMSET info")
        labels = list()
        values = list()
        if self.bs_info.bs_info:
            info = self.bs_info.bs_info
            for a in beamset_info.bs_attrs:
                labels.append(a)
                values.append(info[a])
        return values,labels
    def GetAndClearWarnings(self):
        # TODO: apply same warning collection system to other objects, like ct image?
        return self.bs_info.GetAndClearWarnings()
    def NoDataText(self):
        return self.NA
    def BeamIsSelected(self,beamname):
        if beamname not in self.beam_selection.keys():
            logger.error("UNKNOWN BEAM NAME '{}'".format(beamname))
            # raise error?
            return False
        if not self.beam_selection[beamname] :
            logger.info("beam name {} was DESELECTED => no simulation".format(beamname))
            return False
        return True
    def SetBeamSelection(self,selection):
        assert(type(selection)==dict)
        if set(selection.keys()) != set(self.beam_selection.keys()):
            missing = [ k for k in self.beam_selection.keys() if k not in selection.keys() ]
            unknown = [ k for k in selection.keys() if k not in self.beam_selection.keys() ]
            logger.error("wrong set of beam names in beam selection call")
            logger.error("missing beam names: " + (",".join(missing) if missing else "(none)"))
            logger.error("unknown beam names: " + (",".join(unknown) if unknown else "(none)"))
            logger.error("expected beam names: " + ",".join(self.beam_selection.keys()))
            logger.error("got beam names: " + ",".join(selection.keys()))
            raise KeyError("wrong set of beam names for beam selection")
        self.beam_selection = selection
    #def GetNumberOfPrimariesPerCore(self,beamname,ncores):
    #    assert(ncores>0)
    #    # this will actually return THREE numbers:
    #    # 1. number of primaries per job
    #    # 2. number of jobs
    #    # 3. total number of primaries (product of the first two numbers)
    #    # For the uncertainty based MCStatType, the first and third number are constant upper limits
    #    # (maximum number of primaries to simulate in order to achieve accuracy goal)
    #    if beamname not in self.beam_selection.keys():
    #        logger.error("UNKNOWN BEAM NAME '{}'".format(beamname))
    #        return (0,0,0)
    #    if not self.beam_selection[beamname] :
    #        logger.info("beam name {} was DESELECTED => 0 primaries".format(beamname))
    #        return (0,0,0)
    #    logger.debug("type ncores is {}, value is {}".format(type(ncores),ncores))
    #    logger.debug("type mc_stat_thr is {}, value is {}".format(type(self.mc_stat_thr),self.mc_stat_thr))
    #    logger.debug("type mc_stat_q is {}".format(self.mc_stat_q))
    #    #if self.mc_stat_q == MCStatType.Nions_per_spot:
    #    #    nspots = self.bs_info[beamname].nspots
    #    #    n_per_spot_per_job = 1 + (self.mc_stat_thr - 1) // ncores
    #    #    nprim_per_job = nspots * n_per_spot_per_job
    #    #    njobs = 1 + (self.mc_stat_thr - 1) // n_per_spot_per_job
    #    #    nprim_total = njobs * nprim_per_job
    #    if self.mc_stat_q == MCStatType.Nions_per_beam:
    #        nprim_per_job = 1 + (self.mc_stat_thr-1) // ncores
    #        njobs = ncores
    #        nprim_total = njobs * nprim_per_job
    #    elif self.mc_stat_q == MCStatType.Xpct_unc_in_target or self.mc_stat_q == MCStatType.Nminutes_per_job:
    #        nprim_max = syscfg[MCStatType.cfglabels[MCStatType.Nions_per_beam]][2]
    #        nprim_per_job = 1 + (nprim_max-1) // ncores
    #        njobs = self.njobs if self.mc_stat_q == MCStatType.Nminutes_per_job else ncores
    #        nprim_total = njobs * nprim_per_job
    #    # TODO: enforce a maximum number of primaries per job?
    #    return (nprim_per_job,njobs,nprim_total)
    def UpdateRSOverride(self,name,rslist):
        logger.debug("RS override for beam={}: {}".format(name,"NONE" if len(rslist)==0 else "+".join(rslist)))
        self._rs_overrides[name] = list(rslist)
    def UpdateRMOverride(self,name,rmlist):
        logger.debug("RM override for beam={}: {}".format(name,"NONE" if len(rmlist)==0 else "+".join(rmlist)))
        self._rm_overrides[name] = list(rmlist)
    def ResetOverrides(self):
        self._rs_overrides = dict()
        self._rm_overrides = dict()
    @property
    def RSOverrides(self):
        return dict(self._rs_overrides)
    @property
    def RMOverrides(self):
        return dict(self._rm_overrides)
    def SetNJobs(self,value):
        self.njobs = int(value)
    def SetStatistics(self,i,value):
        if i>=0 and i<MCStatType.NTypes:
            #self.mc_stat_q = i
            #if self.mc_stat_q == MCStatType.Nions_per_spot or self.mc_stat_q == MCStatType.Nions_per_beam:
            if MCStatType.is_int[i]:
                self.mc_stat_thr[i] = int(value)
            else:
                self.mc_stat_thr[i] = float(value)
        else:
            raise RuntimeError(f"PROGRAMMING ERROR: wrong statistics enum value {i}")
        logger.debug("q={} ({}) value={}".format(i,MCStatType.guilabels[i],value))
    def SetGeometry(self,i):
        if i==0:
            assert(self.have_CT)
            logger.debug("CT geometry")
            self._CT = True
            self._PHANTOM = False
            self.SetDefaultDoseGridSettings()
        elif i==1:
            logger.debug("PHANTOM geometry")
            self._CT = False
            self._PHANTOM = True
            if self._phantom_specs is not None:
                self.SetDefaultDoseGridSettings()
        else:
            logger.error("illegal geometry setting: {}".format(igeo))
        for widget in self.ctphantom_subscribers:
            widget.update_ctphantom()
    def UpdatePhantomIsoOverride(self,name,new_xyz):
        # Are we paranoid enough?
        assert(name in self.beam_names)
        assert(type(new_xyz)==tuple)
        assert(len(new_xyz)==3)
        assert(type(new_xyz[0])==float)
        assert(type(new_xyz[1])==float)
        assert(type(new_xyz[2])==float)
        logger.debug("plan details: updating phantom ISO = {} mm (old value: {})".format(new_xyz,self._PhantomISOmm[name]))
        self._PhantomISOmm[name] = np.array(new_xyz)
        for widget in self.ctphantom_subscribers:
            widget.update_ctphantom()
    def UpdatePhantomGEO(self,newspecs):
        oldlabel = "None" if self._phantom_specs is None else self._phantom_specs.label
        newlabel = "None" if newspecs is None else newspecs.label
        logger.debug("plan details: updating phantom specs from {} to {}.".format(oldlabel,newlabel))
        self._phantom_specs = newspecs
        self.SetGeometry(1)
        #self.SetDefaultDoseGridSettings()
        #for widget in self.ctphantom_subscribers:
        #    widget.update_ctphantom()

# vim: set et softtabstop=4 sw=4 smartindent:
