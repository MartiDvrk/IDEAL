import numpy as np
import pydicom

class IDEAL_RP_dictionary:
	
	def __init__(self):
		  
		self.RPGeneral = ["PatientID","PatientName","PatientBirthDate","PatientSex","RTPlanLabel",
		"SOPInstanceUID","ReferringPhysicianName","PlanIntent","RTPlanLabel",
		"SOPClassUID","SOPInstanceUID","IonBeamSequence",
		"FractionGroupSequence","ReferencedStructureSetSequence"]    
		# Optional:  "OperatorsName","ReviewerName","ReviewDate","ReviewTime","DoseReferenceSequence"
		
		self.IonBeamSequence = ["BeamNumber","IonControlPointSequence","FinalCumulativeMetersetWeight",
		"BeamName","RadiationType", #"RadiationAtomicNumber","RadiationMassNumber","RadiationChargeState", (only if RadiationType = 'ION')
		"TreatmentMachineName","NumberOfRangeModulators","NumberOfRangeShifters","PrimaryDosimeterUnit",
		"SnoutSequence"]  # "RangeModulatorSequence","RangeShifterSequence" optiona, depend on "NumberOfRangeModulators","NumberOfRangeShifters"
		
		self.DoseReferenceSequence = ["ReferencedROINumber"]   # Optional: "TargetPrescriptionDose"
		self.ReferencedStructureSetSequence = ["ReferencedSOPInstanceUID"]
		self.FractionGroupSequence = ["ReferencedBeamSequence","NumberOfFractionsPlanned"]
		self.IonControlPointSequence = ["PatientSupportAngle","IsocenterPosition","GantryAngle","SnoutPosition","NominalBeamEnergy",
		"NumberOfScanSpotPositions","ScanSpotMetersetWeights","ScanSpotPositionMap","CumulativeMetersetWeight",
		"ScanSpotTuneID","NumberOfPaintings"]
		self.SnoutID = "SnoutID"
		self.RangeShifterID = "RangeShifterID"
		self.RangeModulatorID =  "RangeModulatorID"
			
		
		
class IDEAL_RD_dictionary:
		
	def __init__(self):
		
		self.RD = ["NumberOfFrames","ReferencedRTPlanSequence","Rows","Columns","DoseGridScaling","PixelSpacing",
		"SliceThickness","ImagePositionPatient","DoseType","SOPClassUID","DoseSummationType","DoseUnits"]
		
		self.ReferencedRTPlanSequence = ["ReferencedSOPInstanceUID"] # "ReferencedFractionGroupSequence" only if "DoseSummationType" == PLAN
		self.ReferencedBeamNumber = "ReferencedBeamNumber"
		


class IDEAL_CT_dictionary:
	
	def __init__(self):
		
		self.CT = ["InstanceCreationDate","SeriesInstanceUID","ImagePositionPatient", #"SeriesDescription",
		"RescaleIntercept","RescaleSlope","InstanceCreationTime","ImagePositionPatient","PixelSpacing"] # Optional: "InstitutionName"
			
			

class IDEAL_RS_dictionary:
	
	def __init__(self):
		
		self.RS = ["SOPClassUID","SeriesInstanceUID","StructureSetROISequence", "ROIContourSequence",
		 "RTROIObservationsSequence","ReferencedFrameOfReferenceSequence"]	
		self.StructureSetROISequence = ["ROIName","ROINumber"]
		self.ROIContourSequence = ["ReferencedROINumber"]
		self.RTROIObservationsSequence = ["ReferencedROINumber","RTROIInterpretedType"]
		
