#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#   Copyright (C): MedAustron GmbH, ACMIT Gmbh and Medical University Vienna
#   This software is distributed under the terms
#   of the GNU Lesser General  Public Licence (LGPL)
#   See LICENSE for further details
# -----------------------------------------------------------------------------

import os, sys, re
import pydicom
import configparser
import itk
import numpy as np
import logging
import tarfile
import shutil
from glob import glob
from datetime import datetime
from utils.gamma_index import get_gamma_index
from functools import reduce

if False:
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger()
else:
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(pathname)s - %(lineno)d - %(levelname)s - %(message)s')
    logfilename="postprocessor.log"
    fh = logging.FileHandler(logfilename)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

from utils.resample_dose import mass_weighted_resampling

def update_user_logs(user_cfg,status,section="DEFAULT",changes=dict()):
    if bool(user_cfg):
        parser=configparser.ConfigParser()
        logger.debug("going to read user logs/settings in {}".format(user_cfg))
        with open(user_cfg,"r") as fp:
            parser.read_file(fp)
        logger.debug("going to update user logs/settings in section {}".format(section))
        if bool(changes):
            parser[section].update(changes)
        parser['DEFAULT']["status"] = status
        parser['DEFAULT']["date and time of last update"] = datetime.now().ctime()
        with open(user_cfg,"w") as fp:
            parser.write(fp)
        logger.debug("finished updating user logs/settings in {}".format(user_cfg))

######################################################################################
# Write MHD image to DICOM (this should probably go to "utils")
######################################################################################
def image_2_dicom_dose(img_dose,dose_dcm_template,my_dose_dcm,physical=True):
    try:
        dose_dcm = pydicom.dcmread(dose_dcm_template)
        aimg_dose = itk.GetArrayFromImage(img_dose)
        img_size = img_dose.GetLargestPossibleRegion().GetSize()
        img_spacing = img_dose.GetSpacing()
        img_origin = img_dose.GetOrigin()
        logger.debug("template dose shape is {}".format(dose_dcm.pixel_array.shape))
        logger.debug("Gate output dose shape is {}".format(aimg_dose.shape))
        maxdose = np.max(aimg_dose)
        logger.debug("max dose in IDC image is {}, max dose in template DICOM is {}".format(maxdose,float(dose_dcm.DoseGridScaling)*np.max(dose_dcm.pixel_array)))
        nbitA = int(dose_dcm.BitsAllocated)
        nbitS = int(dose_dcm.BitsStored)
        if not nbitA == nbitS:
            raise RuntimeError("inconsistent template dicom dose file {}: nbitA={}, nbitS={}".format(dose_dcm_template,nbitA,nbitS))
        if maxdose > 0.:
            logger.debug("max dose in image is positive, good!")
        else:
            logger.error("max dose in image is NOT positive, BAD!")
            raise RuntimeError("max dose in image is NOT positive, BAD!")
        logger.debug("check that these are 16: nbitA={}, nbitS={}".format(nbitA,nbitS))
        dose_grid_scaling = maxdose/(2**nbitS-1)
        aimg_dose_int = np.round(aimg_dose / dose_grid_scaling).astype('uint16')
        #logger.debug("computed dose grid scale factor is {}, input DICOM DGSF is {}".format(dose_grid_scaling,float(dose_dcm.DoseGridScaling)))
        #logger.debug("{} nonzero voxels in input image, max is {}".format(np.sum(aimg_dose_int>0),np.max(aimg_dose_int)))
        #logger.debug("{} nonzero voxels in DICOM template, max is {}".format(np.sum(dose_dcm.pixel_array>0),np.max(dose_dcm.pixel_array)))
        dose_dcm.PixelData = aimg_dose_int.tobytes()
        dose_dcm.DoseGridScaling = dose_grid_scaling
        dose_dcm.DoseType = 'PHYSICAL' if physical else 'EFFECTIVE'
        dose_dcm.Rows = img_size[1]
        dose_dcm.Columns = img_size[0]
        dose_dcm.NumberOfFrames = img_size[2]
        dose_dcm.PixelSpacing[0] = img_spacing[1]
        dose_dcm.PixelSpacing[1] = img_spacing[0]
        dose_dcm.SliceThickness = img_spacing[2]
        logger.debug("going to resize grid frame offset vector")
        logger.debug("offset vector size is {} in the DICOM template".format(len(dose_dcm.GridFrameOffsetVector)))
        #dose_dcm.GridFrameOffsetVector.clear()
        #logger.debug("after the clearning, size is {}".format(len(dose_dcm.GridFrameOffsetVector)))
        logger.debug("going to create a new vector with grid offsets with spacing {} and size {}".format(float(img_spacing[2]),int(img_size[2])))
        new_gfov=pydicom.multival.MultiValue(pydicom.valuerep.DSfloat,np.arange(int(img_size[2]))*float(img_spacing[2]))
        logger.debug("going to override the old offset vector with the new vector")
        dose_dcm.GridFrameOffsetVector = new_gfov
        logger.debug("after the override, size is {}".format(len(dose_dcm.GridFrameOffsetVector)))
        dose_dcm.ImagePositionPatient = list(img_origin)
        #dose_dcm.InstanceNumber = ""
        #logger.debug("generating DCM names")
        #logger.debug("writing DICOM file with same UID like template, with minimal changes: {}".format(dcmA))
        #pydicom.write_file(my_dose_dcm.replace(".dcm","A.dcm"),dose_dcm,False)
        #logger.debug("writing DICOM file with same UID like template, with some DICOM standard compliance changes: {}".format(dcmB))
        #pydicom.write_file(my_dose_dcm.replace(".dcm","B.dcm"),dose_dcm,True)
        logger.debug("assigning new UID for this dose file. Old UID={}".format(dose_dcm.SOPInstanceUID))
        dose_dcm.SOPInstanceUID = pydicom.uid.generate_uid() # create a new unique UID
        logger.debug("new UID is: {}".format(dose_dcm.SOPInstanceUID))
        #logger.debug("writing DICOM file with new UID, with minimal changes: {}".format(dcmC))
        #pydicom.write_file(my_dose_dcm.replace(".dcm","C.dcm"),dose_dcm,True)
        #logger.debug("writing DICOM file with new UID, and with some DICOM standard compliance changes: {}".format(dcmD))
        #pydicom.write_file(my_dose_dcm.replace(".dcm","D.dcm"),dose_dcm,False)
        #logger.info("wrote A,B,C,D DICOM file: {}".format(my_dose_dcm))
        logger.debug("going to write to file: {}".format(my_dose_dcm))
        pydicom.write_file(my_dose_dcm,dose_dcm,False)
        logger.info("wrote DICOM file: {}".format(my_dose_dcm))
    except Exception as e:
        logger.info("something went wrong: {}".format(e))

def run_gamma_analysis(ref_dose_path,gamma_parameters,dose_sum_final,mhd_dose_final):
    ushort_imgref=itk.imread(ref_dose_path)
    aimgref=itk.array_from_image(ushort_imgref)*float(pydicom.dcmread(ref_dose_path).DoseGridScaling)
    imgref=itk.image_from_array(np.float32(aimgref))
    imgref.CopyInformation(ushort_imgref)
    npar = len(gamma_parameters)
    if not npar==4:
        raise ValueError(f"wrong number gamma index parameters ({npar}, should be 4)")
    dta_mm,dd_percent,dosethr,defgamma = gamma_parameters.tolist()
    g=get_gamma_index(ref=imgref,
                      target=dose_sum_final,
                      dta=dta_mm,
                      dd=dd_percent,
                      ddpercent=True,
                      threshold=dosethr,
                      defvalue=defgamma,
                      verbose=False,
                      threshold_percent=True)
    itk.imwrite(g,mhd_dose_final.replace(".mhd","_gamma.mhd"))
    itk.imwrite(imgref,mhd_dose_final.replace(".mhd","_tpsdose.mhd"))

def update_plan_dose(pdd,label,beam_dose_image):
    # 'pdd' is plan dose dictionary
    # label will be "unresampled", "Physical" or "RBE"
    if label in pdd:
        # add dose image to dose image in the dict
        a_plandose = itk.array_from_image(pdd[label])
        a_beamdose = itk.array_view_from_image(beam_dose_image)
        assert(a_plandose.shape==a_beamdose.shape)
        a_plandose += a_beamdose
        img_plandose = itk.image_from_array(a_plandose)
        img_plandose.CopyInformation(pdd[label])
        pdd[label] = img_plandose
    else:
        # copy dose image into dict
        a_plandose = itk.array_from_image(beam_dose_image)
        img_plandose = itk.image_from_array(a_plandose)
        img_plandose.CopyInformation(beam_dose_image)
        pdd[label] = img_plandose

def get_job_stats(mhd):
    """
    This function assumes that each output directory contains the output of
    exactly one GATE run, with one dose actor output (mhd file) and one stat
    actor output (txt file).
    """
    stats=glob(os.path.join(os.path.dirname(mhd),"stat*.txt"))
    gate_exit_value_txt=glob(os.path.join(os.path.dirname(mhd),"gate_exit_value.txt"))
    logger.debug("stat files: {}".format("\n".join(stats)))
    if not 1 == len(stats):
        raise RuntimeError("cannot retrieve number of primaries for {}".format(mhd) + \
                           "got {} associated stat actor files: {}".format(len(stats),"\n".join(stats)))
    if not 1 == len(gate_exit_value_txt):
        raise RuntimeError("cannot retrieve gate exit value for {}".format(mhd))
    gate_exit_value=0
    with open(gate_exit_value_txt[0],"r") as gev:
        line=gev.readline()
        gate_exit_value=int(line.strip())
    statfilepath=stats[0]
    stats_dict={"StatsFile":statfilepath}
    hash_key_val=re.compile(r'^#\s*(\b.*\S)\s*=\s*(\b.*\b)\s*$')
    with open(statfilepath,"r") as sf:
        stats_dict.update(dict([hash_key_val.search(line).groups() for line in sf if hash_key_val.search(line) is not None ]))
        n = int(stats_dict.get('NumberOfEvents',-1))
        if n>=0:
            logger.info("{} got {} primaries".format(mhd,n))
            return stats_dict,gate_exit_value
    raise RuntimeError("N primaries not found for {}".format(mhd))

def compress_jobdata(cfg,outputdirs,statfiles):
        try:
            logger.debug("start logging of tarball compression of {} output directories".format(len(outputdirs)))
            dir_tot=0
            tgz_tot=0
            basenames=set([os.path.basename(txt) for txt in statfiles])
            if len(basenames)!=1:
                logger.error("the statistics actor files all should have the same basename, but I got {} different ones: {}".format(len(basenames),basenames))
            else:
                # should we be paranoid and check that they have a .txt suffix?
                tgz_stats_name=basenames.pop().replace("txt","tar.gz")
                with tarfile.open(tgz_stats_name,"w:gz") as tgz_stats:
                    for s in statfiles:
                        tgz_stats.add(s)
            for d in outputdirs:
                cwd=os.path.realpath(os.curdir)
                dbase=os.path.basename(os.path.realpath(d))
                logger.debug("cwd={} d={}".format(cwd,d))
                if not os.path.isdir(d):
                    logger.debug("compression request for {} which is not a directory!".format(d))
                    continue
                if cfg.debug:
                    dtgz=dbase+".tar.gz"
                    logger.debug("going to change directory")
                    os.chdir(os.path.dirname(os.path.realpath(d)))
                    dir_size=sum([os.stat(os.path.join(dbase,f)).st_size for f in os.listdir(dbase)])
                    dir_tot+=dir_size
                    logger.debug("going to compress {}".format(dbase))
                    with tarfile.open(dtgz,"w:gz") as tgz:
                        logger.debug("created {}".format(dtgz))
                        tgz.add(dbase)
                        logger.debug("contents:\n"+"\n".join(tgz.getnames()))
                    tgz_size=os.stat(dtgz).st_size
                    tgz_tot+=tgz_size
                    logger.debug("directory size was approx {}, tar ball size is {}, saved {} bytes".format(dir_size,tgz_size,dir_size-tgz_size))
                    logger.debug("finished compressing {}".format(dbase))
                    MiB=1024.**2
                    GiB=1024.**3
                    logger.info("directories: {0} bytes = {1:.2f} GiB; tar balls: {2} bytes = {3:.1f} MiB; saved: {4} bytes = {5:.2f} GiB".format(
                                              dir_tot,    dir_tot/GiB,            tgz_tot,    tgz_tot/MiB, dir_tot-tgz_tot,(dir_tot-tgz_tot)/GiB))
                shutil.rmtree(dbase)
                logger.debug("removed directory {}".format(dbase))
                os.chdir(cwd)
        except Exception as e:
            logger.error("oopsie: '{}'".format(e))

def get_img_array(img_path):
    dose=itk.imread(img_path)
    return itk.GetArrayFromImage(dose)

######################################################################################
# Implementation details: accumulate the doses, apply rescaling and correction factors
######################################################################################
def post_processing(cfg,pdd,cul):
    # cfg=config
    # pdd=plan dose dictionary
    # cul=cleanup list
    if bool(cfg.user_cfg):
        update_user_logs(cfg.user_cfg,status=f"POSTPROCESSING beam '{cfg.origname}'")

    # filename definitions
    mhd_dose_sum = str(os.path.join(str(cfg.output_dicom1),cfg.dosemhd))
    mhd_dose_rescaled = mhd_dose_sum.replace(".mhd","-Rescaled.mhd")
    if bool(cfg.apply_external_dose_mask):
        mhd_dose_masked = mhd_dose_rescaled.replace(".mhd","-Masked.mhd")
    else:
        mhd_dose_masked = mhd_dose_rescaled.replace(".mhd","-Unmasked.mhd")
    mhd_dose_physical = mhd_dose_masked.replace(".mhd","-Physical.mhd")
    mhd_dose_rbe = mhd_dose_masked.replace(".mhd","-RBE.mhd")
    dcm_dose_physical = mhd_dose_physical.replace(".mhd",".dcm")
    dcm_dose_full_ct = dcm_dose_physical.replace(".dcm","_DEBUG_FULL_CT_GRID.dcm")
    dcm_dose_rbe = mhd_dose_rbe.replace(".mhd",".dcm")
    mhd_re_alanine = mhd_dose_sum.replace(".mhd","-RE_Alanine.mhd")
    dcm_re_alanine = mhd_re_alanine.replace(".mhd",".dcm")
    mhdlist = [str(os.path.join(d,cfg.dosemhd)) for d in os.listdir(os.curdir) if d[:7]=="output." and os.path.isdir(d) and os.path.exists(os.path.join(d,cfg.dosemhd))]
    re_num_list = [str(os.path.join(d,cfg.re_num_mhd)) for d in os.listdir(os.curdir) if d[:7]=="output." and os.path.isdir(d) and os.path.exists(os.path.join(d,cfg.re_num_mhd))]
    re_den_list = [str(os.path.join(d,cfg.re_den_mhd)) for d in os.listdir(os.curdir) if d[:7]=="output." and os.path.isdir(d) and os.path.exists(os.path.join(d,cfg.re_den_mhd))]
    if len(mhdlist)==cfg.nJobs:
        logger.info("got all {} dose files".format(cfg.nJobs))
    else:
        logger.warn("got {} dose files, actually {} were expected!".format(len(mhdlist),cfg.nJobs))
    if len(mhdlist)==0:
        logger.error("did not find any dose files named '{}' for beam '{}'".format(cfg.dosemhd,cfg.origname))
        return False
    # sum
    logger.debug("going to sum all {} doses and do some rescaling".format(len(mhdlist)))
    logger.debug("first dose file is {}".format(mhdlist[0]))
    logger.debug("type first dose file is {}".format(type(mhdlist[0])))
    dose0=itk.imread(mhdlist[0])
    logger.debug("dose distribution has orig={} spacing={} size={}".format(dose0.GetOrigin(),dose0.GetSpacing(),dose0.GetLargestPossibleRegion().GetSize()))
    adose=itk.GetArrayFromImage(dose0)
    statdict,retval=get_job_stats(mhdlist[0])
    nMC=int(statdict['NumberOfEvents'])
    nBADretval=0
    nBADzeronmc=0
    tCPUbrutto=0.
    tCPUnetto=0.
    statfiles=[statdict['StatsFile']]
    logger.debug("adding dose from {} primaries, Gate return value was {}".format(nMC,retval))
    if retval != 0:
        logger.error("return value {} means that something went WRONG, Gate did not terminate normally".format(retval))
        nBADretval += 1
    elif nMC <= 0:
        logger.error("ZERO ({}) primaries from mhd={}".format(nMCjob,mhd))
        nBADzeronmc += 1
    else:
        tCPUbrutto += float(statdict['ElapsedTime'])
        tCPUnetto += float(statdict['ElapsedTimeWoInit'])
    for mhd in mhdlist[1:]:
        logger.debug("next dose file is {}".format(mhd))
        try:
            dose=itk.imread(mhd)
            statdict,retval = get_job_stats(mhd)
            nMCjob = int(statdict['NumberOfEvents'])
            statfiles.append(statdict['StatsFile'])
            #nMCjob,statfile,retval = nprimaries(mhd)
            #statfiles.append(statfile)
            logger.debug("adding dose from {} primaries, Gate return value was {}".format(nMCjob,retval))
            # FIXME: such errors should be reported in the final result
            if retval != 0:
                nBADretval += 1
                raise RuntimeError("return value {} means that something went WRONG, Gate did not terminate normally".format(retval))
            # FIXME: such errors should be reported in the final result
            elif nMCjob <= 0:
                nBADzeronmc += 1
                raise RuntimeError("ZERO ({}) primaries from mhd={}".format(nMCjob,mhd))
            nMC += nMCjob
            tCPUbrutto += float(statdict['ElapsedTime'])
            tCPUnetto += float(statdict['ElapsedTimeWoInit'])
            assert bool(tuple(dose0.GetLargestPossibleRegion().GetSize()) == tuple(dose.GetLargestPossibleRegion().GetSize())), str("sizes {} and {} don't match".format(dose0.GetLargestPossibleRegion().GetSize(),dose.GetLargestPossibleRegion().GetSize()))
            assert bool(np.allclose(tuple(dose0.GetOrigin()), tuple(dose.GetOrigin()))), str("origins don't match")  # TODO: check that this sufficiently allows rounding differences
            assert bool(np.allclose(tuple(dose0.GetSpacing()),tuple(dose.GetSpacing()))), str("spacings don't match") # TODO: check that this sufficiently allows rounding differences
            adose += itk.GetArrayViewFromImage(dose)
            logger.debug("max dose (unscaled) is now {}".format(np.max(adose)))
        except Exception as e:
            # FIXME: such errors should be reported in the final result
            logger.error("something went wrong while processing {}: {}".format(mhd,e))
    if nMC <= 0:
        logger.error("failed to find any primaries for beam '{}', cannot scale any dose.".format(cfg.origname))
        return False
    logger.info("total simulated number of primaries is {}".format(nMC))
    # if relative efficiency analysis: create rel eff image
    logger.info("write_dicom_re_alanine = {}".format(cfg.write_dicom_re_alanine))
    if cfg.write_dicom_re_alanine:
        logger.info("Going to sum up separately numerator and denominator for re alanine")
        logger.info("found {} images".format(len(re_num_list)))
        num = get_img_array(re_num_list[0])
        denom = get_img_array(re_den_list[0])
        for mhd_n, mhd_d in zip(re_num_list[1:],re_den_list[1:]):
            print(mhd)
            num += get_img_array(mhd_n)
            denom += get_img_array(mhd_d)
        #sum_imgs = lambda mhd1, mhd2: get_img_array(mhd1) + get_img_array(mhd2)
        #num = reduce(sum_imgs,re_num_list)
        #denom = reduce(sum_imgs,re_den_list)
        logger.info("Sums are ready, going to perform division")
        re_alanine_array = np.divide(num, denom, out=np.zeros_like(num), where=denom!=0)
        logger.info("going to make an image")
        re_alanine_img = itk.GetImageFromArray(np.float32(re_alanine_array))
        re_alanine_img.CopyInformation(dose0) # assuming all image info (origine, resolution, spacing) are the same as the dose img
        #itk.imwrite(re_alanine_img,mhd_re_alanine)
        logger.info("Write RE image to dicom")
        image_2_dicom_dose(re_alanine_img,str(cfg.dcm_beam_in),str(dcm_re_alanine),physical=True)
        
    # now make an image
    dose_sum = itk.GetImageFromArray(np.float32(adose))
    dose_sum.CopyInformation(dose0)
    if cfg.write_mhd_unscaled_dose:
        itk.imwrite(dose_sum,mhd_dose_sum)
    # rescaling: get physical dose
    logger.info("scaling with number of fractions = {}".format(cfg.nFractions))
    adose *= cfg.nFractions
    scale_factor = cfg.dosecorrfactor*float(cfg.nTPS)/float(nMC)
    logger.info("scaling with number dose_correction_factor*nTPS/nMC = {}*{}/{} = {}".format(cfg.dosecorrfactor,cfg.nTPS,nMC,scale_factor))
    adose*=scale_factor
    dose_sum_rescaled = itk.GetImageFromArray(np.float32(adose))
    dose_sum_rescaled.CopyInformation(dose0)
    if cfg.write_mhd_scaled_dose:
        itk.imwrite(dose_sum_rescaled,mhd_dose_rescaled)
    if cfg.write_unresampled_dose:
        image_2_dicom_dose(dose_sum_rescaled,str(cfg.dcm_beam_in),str(dcm_dose_full_ct),physical=True)
        if cfg.dicom_plan_dose or cfg.mhd_plan_dose:
            update_plan_dose(pdd,"unresampled",dose_sum_rescaled)
    # override
    dose_spacing = cfg.dose_size/cfg.dose_nvoxels
    if cfg.mass_mhd:
        try:
            dose_resampled_ref = itk.GetImageFromArray(np.zeros(cfg.dose_nvoxels[::-1],dtype=np.float32))
            dose_resampled_ref.SetOrigin(cfg.dose_origin)
            dose_resampled_ref.SetSpacing(dose_spacing)
            mass_img=itk.imread(cfg.mass_mhd)
            logger.debug("dose_sum has dimsize={} mass has dimsize={}".format(np.array(itk.size(dose_sum_rescaled)),np.array(itk.size(mass_img))))
            logger.debug("going to resample from voxels with spacing {} to voxels with spacing {}".format(dose_sum_rescaled.GetSpacing(),dose_resampled_ref.GetSpacing()))
            t0=datetime.now()
            dose_physical = mass_weighted_resampling(dose_sum_rescaled,mass_img,dose_resampled_ref)
            t1=datetime.now()
            logger.debug("resampling took {} seconds".format((t1-t0).total_seconds()))
        except Exception as e:
            # whatever goes wrong, it should be reported in the log file
            logger.error(f"something when wrong during resampling: {e}")
            raise
    else:
        logger.debug("check: size=size {}".format("TRUE" if (cfg.dose_nvoxels==np.array(dose_sum_rescaled.GetLargestPossibleRegion().GetSize())).all() else "FALSE"))
        logger.debug("check: spacing=spacing {}".format("TRUE" if np.allclose(dose_spacing,dose_sum_rescaled.GetSpacing()) else "FALSE"))
        logger.debug("check: origin=origin {}".format("TRUE" if np.allclose(cfg.dose_origin,dose_sum_rescaled.GetOrigin()) else "FALSE"))
        dose_physical = dose_sum_rescaled
        dose_physical.SetOrigin(cfg.dose_origin)
    adose = itk.GetArrayFromImage(dose_physical)
    if cfg.apply_external_dose_mask:
        logger.debug("going to apply ROI mask from file {}".format(cfg.external_dose_mask))
        mask=itk.imread(str(cfg.external_dose_mask))
        logger.debug("succeeded reading mask image from file {}".format(cfg.external_dose_mask))
        amask=itk.GetArrayViewFromImage(mask)>0
        logger.debug("got mask as array: shape dose={} shape ROI mask = {}".format(adose.shape,amask.shape))
        amask_not = np.logical_not(amask)
        for iz in range(amask.shape[0]):
            logger.debug("iz={} #enables(iz)={} #disables(iz)={}".format(iz,np.sum(amask[iz,:,:]),np.sum(amask_not[iz,:,:])))
        one_percent=np.prod(amask.shape)/100.0
        n_in=np.sum(amask)
        n_out=np.sum(amask_not)
        logger.debug("total: mask enables/disables {0}/{1} voxels ({2:.2f}/{3:.2f} percent of the image)".format(n_in,n_out,n_in/one_percent,n_out/one_percent))
        adose*=amask
    new_dose_physical=itk.GetImageFromArray(np.float32(adose))
    new_dose_physical.CopyInformation(dose_physical)
    dose_physical = new_dose_physical
    if cfg.write_mhd_physical_dose:
        itk.imwrite(dose_physical,mhd_dose_physical)
    if cfg.write_dicom_physical_dose:
        image_2_dicom_dose(dose_physical,str(cfg.dcm_beam_in),str(dcm_dose_physical),physical=True)
        if cfg.dicom_plan_dose or cfg.mhd_plan_dose:
            update_plan_dose(pdd,"Physical",dose_physical)
    # rescaling: get RBE-corrected dose
    if np.isclose(cfg.RBE_factor,1.0):
        logger.debug("NOT multiplying dose with RBE factor because it is close to unity: {}".format(cfg.RBE_factor))
        dose_sum_final = dose_physical
        mhd_dose_final = mhd_dose_physical
    else:
        logger.debug("multiplying dose with RBE factor = {}".format(cfg.RBE_factor))
        adose*=cfg.RBE_factor
        dose_rbe = itk.GetImageFromArray(np.float32(adose))
        dose_rbe.CopyInformation(dose_physical)
        dose_sum_final = dose_rbe
        mhd_dose_final = mhd_dose_rbe
        if cfg.write_mhd_rbe_dose:
            itk.imwrite(dose_rbe,mhd_dose_rbe)
        if cfg.write_dicom_rbe_dose:
            image_2_dicom_dose(dose_rbe,str(cfg.dcm_beam_in),str(dcm_dose_rbe),physical=False)
            if cfg.dicom_plan_dose or cfg.mhd_plan_dose:
                update_plan_dose(pdd,"RBE",dose_rbe)
    if cfg.ref_dose_path:
        if cfg.gamma_analysis:
            try:
                logger.debug("going to run gamma analysis, using ref dose = {}".format(cfg.ref_dose_path))
                t0=datetime.now()
                run_gamma_analysis(cfg.ref_dose_path,cfg.gamma_parameters,dose_sum_final,mhd_dose_final)
                #ushort_imgref=itk.imread(cfg.ref_dose_path)
                #aimgref=itk.GetArrayFromImage(ushort_imgref)*float(pydicom.dcmread(cfg.ref_dose_path).DoseGridScaling)
                #imgref=itk.GetImageFromArray(np.float32(aimgref))
                #imgref.CopyInformation(ushort_imgref)
                #npar = len(cfg.gamma_parameters)
                #assert bool(npar==4),"wrong number gamma index parameters ({}, should be 4)".format(npar)
                #dta_mm,dd_percent,dosethr,defgamma = cfg.gamma_parameters.tolist()
                #g=get_gamma_index(ref=imgref,target=dose_sum_final,dta=dta_mm,dd=dd_percent, ddpercent=True,threshold=dosethr,defvalue=defgamma,verbose=False)
                #itk.imwrite(g,mhd_dose_final.replace(".mhd","_gamma.mhd"))
                #itk.imwrite(imgref,mhd_dose_final.replace(".mhd","_tpsdose.mhd"))
                t1=datetime.now()
                logger.debug("gamma index calculation took {} seconds".format((t1-t0).total_seconds()))
            except Exception as e:
                logger.error("something went wrong when attempting to compute the gamma index distribution: {}".format(e))
    else:
        logger.debug(f"NO gamma index calculation for '{os.path.basename(mhd_dose_final)}'")
    # update user settings/logs
    if bool(cfg.user_cfg):
        update_user_logs(cfg.user_cfg,status=f"FINISHED POSTPROCESSING beam '{cfg.origname}'",
                section=cfg.origname,
                changes={"total number of primaries":str(nMC),
                         "number of failed jobs":str(nBADretval),
                         "number of jobs with zero primaries":str(nBADzeronmc),
                         "CPU time [seconds] including init":str(tCPUbrutto),
                         "CPU time [seconds] excluding init":str(tCPUnetto),
                         "CPU time [hours] including init":str(tCPUbrutto/3600.),
                         "CPU time [hours] excluding init":str(tCPUnetto/3600.),
                         "number of primaries per second per core":str(nMC/tCPUnetto) })
    # update the clean up list
    outputdirs = [ os.path.realpath(os.path.dirname(mhd)) for mhd in mhdlist ]
    #logger.debug("going to compress {} output directories".format(len(outputdirs)))
    #compress_jobdata(outputdirs,statfiles)
    cul.append( (outputdirs,statfiles) )
    return True

class post_proc_config:
    def __init__(self,prsr,beamname):
        sec=prsr[beamname]
        self.user_cfg = prsr['user logs file']['path']
        self.beamname = beamname
        self.origname = sec['origname']
        #nMC=sec.getint("nmc")
        self.postproc_time = datetime.now()
        self.nJobs=sec.getint("njobs")
        #nMCtot=sec.getint("nmctot")
        self.nTPS=sec.getfloat("ntps")
        self.dosecorrfactor=sec.getfloat("dosecorrfactor")
        dosemhd=sec.get("dosemhd")
        self.dose2water=sec.getboolean("dose2water")
        self.dosemhd=dosemhd.replace(".mhd","-DoseToWater.mhd" if self.dose2water else "-Dose.mhd")
        # Alanine
        self.re_num_mhd = self.dosemhd.replace(".mhd","_RE_alanine-relEfficiency-letToG4_ALANINE-numerator.mhd")
        self.re_den_mhd = self.dosemhd.replace(".mhd","_RE_alanine-relEfficiency-letToG4_ALANINE-denominator.mhd")
        self.dose_origin=np.array([float(v) for v in sec.get("dose grid origin").split()])
        self.dose_size=np.array([float(v) for v in sec.get("dose grid size").split()])
        self.dose_nvoxels=np.array([int(v) for v in sec.get("dose grid resolution").split()])
        self.output_dicom1=sec.get("first output dicom")
        self.output_dicom2=sec.get("second output dicom")
        self.RBE_factor=sec.getfloat("rbe")
        self.nFractions=sec.getint("nfractions")
        self.dcm_beam_in=sec.get("dcm template")
        self.mass_mhd = sec.get("mass mhd","")
        # MFA 11/16/22
        self.gamma_analysis = sec.getboolean("run gamma analysis")
        self.debug = sec.getboolean("debug")
        
        self.write_mhd_unscaled_dose = sec.getboolean("write mhd unscaled dose")
        self.write_mhd_scaled_dose = sec.getboolean("write mhd scaled dose")
        self.write_mhd_physical_dose = sec.getboolean("write mhd physical dose")
        self.write_mhd_rbe_dose = sec.getboolean("write mhd rbe dose")
        self.write_dicom_physical_dose = sec.getboolean("write dicom physical dose")
        self.write_dicom_rbe_dose = sec.getboolean("write dicom rbe dose")
        # Alanine
        self.write_dicom_re_alanine = sec.getboolean("write dicom re alanine")
        self.dicom_plan_dose = sec.get("dicom plan dose")
        self.mhd_plan_dose = sec.get("mhd plan dose")
        self.dcm_plan_in = sec.get("plan dcm template")
        self.write_unresampled_dose = sec.getboolean("write unresampled dose")
        if self.write_dicom_physical_dose or self.write_dicom_rbe_dose:
            if not os.path.exists(self.dcm_beam_in):
                msg = "missing DICOM beam dose template file {self.dcm_beam_in}"
                logger.error(msg)
                raise RuntimeError(msg)
        if self.dicom_plan_dose:
            if not os.path.exists(self.dcm_plan_in):
                msg = "missing DICOM plan dose template file {self.dcm_plan_in}"
                logger.error(msg)
                raise RuntimeError(msg)
        self.apply_external_dose_mask = sec.getboolean("apply external dose mask",fallback=False)
        self.external_dose_mask = sec.get("external dose mask","")
        if self.apply_external_dose_mask and not bool(self.external_dose_mask):
            msg = "masking of dose with external requested, but no mask provided"
            logger.error(msg)
            raise RuntimeError(msg)
        self.ref_dose_path = sec.get("path to reference dose image for gamma index calculation","")
        self.gamma_parameters = np.array([float(v) for v in sec.get("gamma index parameters dta_mm dd_percent thr_percent def","").split()])
        self.ref_physical_plan_dose_path = sec.get("path to reference PHYSICAL plan dose image for gamma index calculation","")
        self.ref_effective_plan_dose_path = sec.get("path to reference EFFECTIVE plan dose image for gamma index calculation","")
        ## TODO own config for server
        self.send_result_to_url = sec.getboolean("send result")
        self.url_to_send_result_to = sec.get("url to send result", "")

######################################################################################
# MAIN
######################################################################################
if __name__ == '__main__':
    parser=configparser.ConfigParser()
    with open("postprocessor.cfg","r") as fp:
        parser.read_file(fp)
    ok = True
    plan_dose_dict = dict()
    cleanup_list = list()
    # MFA 11/21/22
#    api_cfg = configparser.ConfigParser()
#    with open("/opt/IDEAL-1.1test/cfg/api.cfg","r") as fp:
#        api_cfg.read_file(fp)
        
    for beamname in parser.sections():
        if beamname=='default' or beamname=='user logs file':
            continue
        t0 = datetime.now()
        cfg = post_proc_config(parser,beamname)
        # TODO: the post_processing now also includes the archiving (making a tarball of) the output directories of all subjobs. Maybe this needs to be separated.
        success = post_processing(cfg,plan_dose_dict,cleanup_list)
        t1 = datetime.now()
        dt = (t1-t0).total_seconds()
        if success:
            logger.info('SUCCESSFUL post processing (including the archiving of job data) of beam "{}" took {} seconds'.format(cfg.origname,dt))
        else:
            logger.error('post processing of beam "{}" FAILED after {} seconds'.format(cfg.origname,dt))
            ok = False
    if ok:
        update_user_logs(cfg.user_cfg,status=f"BEAM DOSES OK, COMPUTING PLAN DOSES")
        for label,img_dose in plan_dose_dict.items():
            physical = label.upper()!="RBE"
            if cfg.dicom_plan_dose != "":
                logger.debug(f"going to write {label} PLAN dose to DICOM")
                plan_dose_dcm = str(os.path.join(str(cfg.output_dicom1), cfg.dicom_plan_dose.replace("PLAN.dcm",f"PLAN-{label}.dcm")))
                image_2_dicom_dose(img_dose,cfg.dcm_plan_in,plan_dose_dcm,physical)
                mhd_gamma = plan_dose_dcm[:-4]+".mhd"
                logger.debug(f"finished writing {label} PLAN dose to DICOM")
            if cfg.mhd_plan_dose != "":
                logger.debug(f"going to write {label} PLAN dose to MHD")
                plan_dose_mhd = str(os.path.join(str(cfg.output_dicom1), cfg.mhd_plan_dose.replace("PLAN.mhd",f"PLAN-{label}.mhd")))
                itk.imwrite(img_dose,plan_dose_mhd)
                mhd_gamma = plan_dose_mhd
                logger.debug(f"finished writing {label} PLAN dose to MHD")
            if cfg.ref_physical_plan_dose_path != "" and label.upper() == "PHYSICAL" and cfg.gamma_analysis:
                logger.debug("start gamma index calculation PHYSICAL PLAN DOSE")
                t0=datetime.now()
                run_gamma_analysis(cfg.ref_physical_plan_dose_path,cfg.gamma_parameters,img_dose,mhd_gamma)
                t1=datetime.now()
                logger.debug("gamma index calculation PHYSICAL PLAN DOSE took {} seconds".format((t1-t0).total_seconds()))
            elif cfg.ref_effective_plan_dose_path != "" and label.upper() == "RBE" and cfg.gamma_analysis:
                logger.debug("start gamma index calculation EFFECTIVE PLAN DOSE")
                t0=datetime.now()
                run_gamma_analysis(cfg.ref_effective_plan_dose_path,cfg.gamma_parameters,img_dose,mhd_gamma)
                t1=datetime.now()
                logger.debug("gamma index calculation EFFECTIVE PLAN DOSE took {} seconds".format((t1-t0).total_seconds()))
        if cfg.output_dicom2:
            update_user_logs(cfg.user_cfg,status=f"DOSE POSTPROCESSING OK, COPYING DATA")
            try:
                shutil.copytree(cfg.output_dicom1,cfg.output_dicom2)
                logger.info("succeeded to copy {} as {}".format(cfg.output_dicom1,cfg.output_dicom2))
            except Exception as e:
                logger.warn("failed copy the first dicom output directory {} as {}".format(cfg.output_dicom1,cfg.output_dicom2))
                logger.warn("failure exception: '{}'".format(e))
        else:
            logger.info("no second copy of dicom output")
        update_user_logs(cfg.user_cfg,status=f"POSTPROCESSING OK, CLEANING UP")
        logger.info("going to clean up the 'tmp' directory")
        t0=datetime.now()
        tmp=os.path.join(os.curdir,'tmp')
        shutil.rmtree(tmp)
        t1=datetime.now()
        logger.info("cleaning up the 'tmp' directory {} successful and took {} seconds".format("was NOT" if os.path.exists(tmp) else "was", (t1-t0).total_seconds()))
        for outputdirs,statfiles in cleanup_list:
            compress_jobdata(cfg,outputdirs,statfiles)
        t2=datetime.now()
        logger.info("compressing all output directories took {} seconds".format((t2-t1).total_seconds()))
        # TODO i'm trying to send the files (i.e. put them on a specific folder for now)
#        if api_cfg['receiver'].getboolean('send result'):
#            try:
#                jobId = str(cfg.output_dicom1).split("/")[-1]
#                logger.info(f"JoId: {jobId}")
#                server_url = "/user/fava"
#                import shutil
#                shutil.make_archive(os.path.join(server_url, "thisfileissentoaserver"), 'zip', cfg.output_dicom1)
#                logger.info("Archive created")
#                import requests
#                with open(os.path.join(server_url, "thisfileissentoaserver.zip"), 'rb') as f:
#                    logger.info(f"Opened archive to send to {cfg.url_to_send_result_to}")
#                    r = requests.post(api_cfg['receiver']['url to send result']+"/"+jobId, files={"result": f})
#                    logger.info(f"{r.status_code} || {r.text}")
#                t3=datetime.now()
#                logger.info("sending the files to server took {} seconds".format((t3-t2).total_seconds()))
#            except Exception as e:
#                logger.warn("failed to transfer zipped output to server: '{}'".format(e))
#                
        # TODO end
        update_user_logs(cfg.user_cfg,status=f"FINISHED")
    else:
        update_user_logs(cfg.user_cfg,status=f"BEAM DOSE POST PROCESSING FAILED")
        logger.warn("NOT going to clean up the 'tmp' directory, to allow debugging of the reported errors")

# vim: set et softtabstop=4 sw=4 smartindent:
