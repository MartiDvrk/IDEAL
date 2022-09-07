import pydicom
import os
import itk
from impl.IDEAL_dictionary import *
from utils.dose_info import dose_info

class dicom_files:
    def __init__(self,rp_path):
        self.dcm_dir =  os.path.dirname(rp_path) # directory with all dicom files
        # RP
        print("Get RP file")
        self.rp_path = rp_path
        self.rp_data = pydicom.read_file(rp_path)
        self.uid = self.rp_data.SOPInstanceUID # same for all files
        # RD
        print("Get RD files")
        self.rds = dose_info.get_dose_files(self.dcm_dir,self.uid) #dictionary with (dcm_data,path)
        # RS
        print("Get RS file")
        self.rs_data = None
        self.rs_path = None
        self.get_RS_file()
        # CT
        print("Get CT files")
        self.ct_paths = self.get_CT_files() # list with the CT files paths
        
    
    def check_all_dcm(self):
        print("Checking RP file")
        check_RP(self.rp_path)
        
        print("Checking RS file")
        check_RS(self.rs_path)
        
        print("Checking RD files")
        for dp in self.rds.values():
            check_RD(dp.filepath)
            
        print("Checking CT files")
        for ct in self.ct_paths[1]:
            check_CT(ct)
        
        
    def get_RS_file(self):
        ss_ref_uid = self.rp_data.ReferencedStructureSetSequence[0].ReferencedSOPInstanceUID
        print("going to try to find the file with structure set with UID '{}'".format(ss_ref_uid))
        nskip=0
        ndcmfail=0
        nwrongtype=0
        nfiles=len([s for s in os.listdir(self.dcm_dir)])
        for s in os.listdir(self.dcm_dir):
            if s[-4:].lower() != '.dcm':
                nskip+=1
                print("no .dcm suffix: {}".format(s))
                continue
            try:
                #print(s)
                ds = pydicom.dcmread(os.path.join(self.dcm_dir,s))
                dcmtype = ds.SOPClassUID.name
            except:
                ndcmfail+=1
                continue
            if dcmtype == "RT Structure Set Storage" and ss_ref_uid == ds.SOPInstanceUID:
                print("found structure set for CT: {}".format(s))
                self.rs_data = ds
                self.rs_path = os.path.join(self.dcm_dir,s)
                break
            else:
                nwrongtype+=1
                #print("rejected structure set for CT: {}".format(s))
                #print("because it as a wrong SOP class ID: {}".format(dcmtype))
                #print("AND/OR because it has the wrong SOP Instance UID: {} != {}".format(ds.SOPInstanceUID,ss_ref_uid))
        if self.rs_data is None:
            raise RuntimeError("could not find structure set with UID={}; skipped {} with wrong suffix, got {} with 'dcm' suffix but pydicom could not read it, got {} with wrong class UID and/or instance UID. It could well be that this is a commissioning plan without CT and structure set data.".format(ss_ref_uid,nskip,ndcmfail,nwrongtype))

    def get_CT_files(self):
        #ids = sitk.ImageSeriesReader_GetGDCMSeriesIDs(ddir)
        dcmseries_reader = itk.GDCMSeriesFileNames.New(Directory=self.dcm_dir)
        ids = dcmseries_reader.GetSeriesUIDs()
        #print("got DICOM {} series IDs".format(len(ids)))
        flist=list()
        uid = self.rs_data.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].SeriesInstanceUID
        if uid:
            if uid in ids:
                try:
                    #flist = sitk.ImageSeriesReader_GetGDCMSeriesFileNames(ddir,uid)
                    flist = dcmseries_reader.GetFileNames(uid)
                    return uid,flist
                except:
                    logger.error('something wrong with series uid={} in directory {}'.format(uid,self.dcm_dir))
                    raise
        else:
            ctid = list()
            for suid in ids:
                #flist = sitk.ImageSeriesReader_GetGDCMSeriesFileNames(ddir,suid)
                flist = dcmseries_reader.GetFileNames(suid)
                f0 = pydicom.dcmread(flist[0])
                if not hasattr(f0,'SOPClassUID'):
                    logger.warn("weird, file {} has no SOPClassUID".format(os.path.basename(flist[0])))
                    continue
                descr = pydicom.uid.UID_dictionary[f0.SOPClassUID][0]
                if descr == 'CT Image Storage':
                    print('found CT series id {}'.format(suid))
                    ctid.append(suid)
                else:
                    print('not CT: series id {} is a "{}"'.format(suid,descr))
            if len(ctid)>1:
                raise ValueError('no series UID was given, and I found {} different CT image series: {}'.format(len(ctid), ",".join(ctid)))
            elif len(ctid)==1:
                uid = ctid[0]
                #flist = sitk.ImageSeriesReader_GetGDCMSeriesFileNames(ddir,uid)
                flist = dcmseries_reader.GetFileNames(uid)
                return flist
        
        
def check_RP(filepath):

	data = pydicom.read_file(filepath)
	dp = IDEAL_RP_dictionary()
	
	# keys used by IDEAL from RP file (maybe keys are enought?)
	genericTags = dp.RPGeneral
	ionBeamTags = dp.IonBeamSequence
	doseSeqTags = dp.DoseReferenceSequence
	refStructTags = dp.ReferencedStructureSetSequence
	fractionTags = dp.FractionGroupSequence
	icpTags = dp.IonControlPointSequence
	snoutTag = dp.SnoutID
	raShiTag = dp.RangeShifterID
	rangeModTag = dp.RangeModulatorID
	
	## --- Verify that all the tags are present and return an error if some are missing --- ##
		
	missing_keys = []
	
	# check first layer of the hierarchy
	loop_over_tags_level(genericTags, data, missing_keys)

	if "IonBeamSequence" in data:
	
		# check ion beam sequence
		loop_over_tags_level(ionBeamTags, data.IonBeamSequence[0], missing_keys)
		
		# check icp sequence
		if "IonControlPointSequence" not in missing_keys:			
			loop_over_tags_level(icpTags, data.IonBeamSequence[0].IonControlPointSequence[0], missing_keys)
			
		# check snout, rashi and rangMod
		if "NumberOfRangeModulators" not in missing_keys:
			if data.IonBeamSequence[0].NumberOfRangeModulators != 0:
				if "RangeModulatorSequence" not in data.IonBeamSequence[0]:
					missing_keys.append("RangeModulatorSequence")
				elif rangeModTag not in  data.IonBeamSequence[0].RangeModulatorSequence[0]:
					missing_keys.append(rangeModTag)			
			
		if "NumberOfRangeShifters" not in missing_keys:
			if data.IonBeamSequence[0].NumberOfRangeShifters != 0:
				if "RangeShifterSequence" not in data.IonBeamSequence[0]:
					missing_keys.append("RangeShifterSequence")
				elif rangeModTag not in  data.IonBeamSequence[0].RangeModulatorSequence[0]:
					missing_keys.append(raShiTag)
		
		if "SnoutSequence" not in missing_keys:
			if snoutTag not in  data.IonBeamSequence[0].SnoutSequence[0]:
				missing_keys.append("SnoutID")
			
				
	if "DoseReferenceSequence" in data:
			
		# check dose sequence
		loop_over_tags_level(doseSeqTags, data.DoseReferenceSequence[0], missing_keys)
	
	if "ReferencedStructureSetSequence" in data:
		
		# check reference structure sequence
		loop_over_tags_level(refStructTags, data.ReferencedStructureSetSequence[0], missing_keys)
		
	if "FractionGroupSequence" in data:
		
		# check fractions sequence
		loop_over_tags_level(fractionTags, data.FractionGroupSequence[0], missing_keys)

	if missing_keys:
		raise ImportError("DICOM RP file not conform. Missing keys: ",missing_keys)
	else: print("\033[92mRP file ok \033[0m")

		
def check_RS(filepath):
	
	data = pydicom.read_file(filepath) 
	ds = IDEAL_RS_dictionary()
	
	# keys and tags used by IDEAL from RS file
	genericTags = ds.RS
	structTags = ds.StructureSetROISequence 
	contourTags = ds.ROIContourSequence 
	observTags = ds.RTROIObservationsSequence
	
	## --- Verify that all the tags are present and return an error if some are missing --- ##
		
	missing_keys = []
	
	# check first layer of the hierarchy
	loop_over_tags_level(genericTags, data, missing_keys)
	
	if "StructureSetROISequence" in data:
	
		# check structure set ROI sequence
		loop_over_tags_level(structTags, data.StructureSetROISequence[0], missing_keys)
		
	if "ROIContourSequence" in data:
	
		# check ROI contour sequence
		loop_over_tags_level(contourTags, data.ROIContourSequence[0], missing_keys)
		
	if "RTROIObservationsSequence" in data:
	
		# check ROI contour sequence
		loop_over_tags_level(observTags, data.RTROIObservationsSequence[0], missing_keys)
		
	if missing_keys:
		raise ImportError("DICOM RS file not conform. Missing keys: ",missing_keys) 
	else: print("\033[92mRS file ok \033[0m")
	
def check_RD(filepath):
	
	data = pydicom.read_file(filepath) 
	dd = IDEAL_RD_dictionary()
	
	# keys and tags used by IDEAL from RD file
	genericTags = dd.RD 
	planSeqTag = dd.ReferencedRTPlanSequence
	refBeamTag = dd.ReferencedBeamNumber
	
	## --- Verify that all the tags are present and return an error if some are missing --- ##
		
	missing_keys = []
	
	# check first layer of the hierarchy
	loop_over_tags_level(genericTags, data, missing_keys)
	
	# check referenced RT Plan seq 
	if "ReferencedRTPlanSequence" in data:
	
		# check ROI contour sequence
		loop_over_tags_level(planSeqTag, data.ReferencedRTPlanSequence[0], missing_keys)
		
		if "DoseSummationType" in data:
			if data.DoseSummationType != "PLAN":
				# check also ReferencedFractionGroupSequence
				if "ReferencedFractionGroupSequence" not in data.ReferencedRTPlanSequence[0]:
					missing_keys.append("ReferencedFractionGroupSequence")
				elif refBeamTag not in data.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence[0].ReferencedBeamSequence[0]:
					missing_keys.append("ReferencedBeamNumber under ReferencedRTPlanSequence/ReferencedFractionGroupSequence/ReferencedBeamSequence")
		
	if missing_keys:
		raise ImportError("DICOM RD file not conform. Missing keys: ",missing_keys) 
	else: print("\033[92mRD file ok \033[0m")
	
def check_CT(filepath):
	
	data = pydicom.read_file(filepath) 
	dct = IDEAL_CT_dictionary()
	
	# keys and tags used by IDEAL from CT file
	genericTags = dct.CT
	
	## --- Verify that all the tags are present and return an error if some are missing --- ##
	missing_keys = []
	
	# check first layer of the hierarchy
	loop_over_tags_level(genericTags, data, missing_keys)
	
	if missing_keys:
		raise ImportError("DICOM CT file not conform. Missing keys: ",missing_keys) 
	else: print("\033[92mCT file ok \033[0m")
		
def loop_over_tags_level(tags, data, missing_keys):
	
	for key in tags:
		
		if key not in data:
			
			missing_keys.append(key)
		
	#return missing_keys

# function used in IDEAL code to check tags. Alternative to my approach.

def sequence_check(obj,attr,nmin=1,nmax=0,name="object"):
    print("checking that {} has attribute {}".format(name,attr))
    assert(hasattr(obj,attr))
    seq=getattr(obj,attr)
    print("{} has length {}, will check if it >={} and <={}".format(name,len(seq),nmin,nmax))
    assert(len(seq)>=nmin)
    assert(nmax==0 or len(seq)<=nmax)			
		

# ~ if __name__ == '__main__':
		
	# ~ filepath = input()
	# ~ RP_info(filepath)
