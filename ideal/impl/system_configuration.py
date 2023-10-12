# -----------------------------------------------------------------------------
#   Copyright (C): MedAustron GmbH, ACMIT Gmbh and Medical University Vienna
#   This software is distributed under the terms
#   of the GNU Lesser General  Public Licence (LGPL)
#   See LICENSE for further details
# -----------------------------------------------------------------------------

import os,sys
import logging
import getpass
import copy
from impl.idc_enum_types import MCStatType
from impl.phantom_specs import phantom_specs
from impl.dual_logging import get_dual_logging, create_logger, timestamp, get_logging_n
import configparser
from glob import glob
#logger=None

class system_configuration:
    """
    This is a 'singleton' class, only one system configuration object is
    supposed to exist.  The singleton system configuration object behaves like
    a dictionary, except that you get a (shallow) copy of an object from it,
    instead of a reference, when you "get" an "item".

    The system configuration is initialized once at the start of whatever
    program is using it, through the function call `get_sysconfig` (see below).
    After initialization it should then not change anymore.

    The implementation aims to avoid changing the system configuration
    *inadvertently*.  This is not banking software.
    """
    __instance = None
    @staticmethod
    def getInstance():
        """
        Get a reference to the one and only "system configuration" instance, if it exists.
        """
        if __class__.__instance is None:
            raise RuntimeError("system configuration used before initialized")
        return __class__.__instance
    def __init__(self,s=dict()):
        """
        Initialize the singleton object, just once.
        """
        if __class__.__instance is not None:
            raise RuntimeError("system configuration initialized more than once?")
        self.__settings = s
        self.logger = None
        __class__.__instance = self
    def __getitem__(self,name):
        """
        Get a copy of any (available) system configuration item.
        """
        if name in self.__settings:
            # return a copy, not a reference
            return copy.copy(self.__settings[name])
        raise KeyError(f"system config setting not found: '{name}'")
    def override(self,name,newvalue):
        if name in self.__settings:
            self.logger.warn("SYSCFG OVERRIDE {} = {}".format(name,newvalue))
            # new value
            self.__settings[name] = newvalue
        else:
            raise KeyError(f"system config setting not found: '{name}'")
    def set_logger(self,filepath):
        self.logger = get_logging_n(self,want_logfile="default", jobId = filepath)
    
def get_basedirs(syscfg,sysprsr,logger):
    dircfg = sysprsr['directories']
    problems = []
    dlist=["input dicom", "tmpdir jobs", "first output dicom", "second output dicom", "logging", "commissioning" ]
    for k,v in dircfg.items():
        if k not in dlist:
            problems.append("Unknown directory specification '{}' with value '{}'".format(k,v))
        else:
            syscfg[k] = v
    for d in dlist:
        if d not in syscfg.keys():
            problems.append("directory specification missing for: {}".format(d))
        elif d == "second output dicom" and not syscfg[d]:
            logger.debug("NOTE: {}={} in system.cfg means: NO SECOND COPY".format(d,syscfg[d]))
        elif not os.path.isdir(syscfg[d]):
            problems.append("ERROR: not an existing directory: {}={}".format(d,syscfg[d]))
        else:
            logger.debug("good directory {}={}".format(d,syscfg[d]))
    if problems:
        logger.error("ERROR in {}:\n{}".format(syscfg['sysconfig'],'\n'.join(problems)))
        raise IOError("ERRORs in {}, please fix:\n{}".format(syscfg['sysconfig'],'\n'.join(problems)))

def get_commissioning_dirs(syscfg,logger):
    problems = []
    for commd in ["CT", "CT/density","CT/composition","CT/cache","beamlines","phantoms"]:
        chkpath = os.path.join(syscfg['commissioning'],*commd.split("/"))
        if os.path.isdir(chkpath):
            syscfg[commd] = chkpath
            logger.debug("DEBUG: successfully found commissioning directory: {}".format(chkpath))
            if commd == "CT/composition":
                # TODO: do we still need this? Should be defined in CT/hlut.conf
                schneider = os.path.join(syscfg[commd],"Schneider2000MaterialsTable.txt")
                if os.path.exists(schneider):
                    syscfg["SchneiderMaterials"] = schneider
                else:
                    problems.append("ERROR: missing file {} in directory: {}".format(os.path.basename(schneider),
                                                                                     os.path.dirname(schneider) ) )
        else:
            problems.append("ERROR: missing commissioning directory: {}".format(chkpath))
    if problems:
        logger.error("ERROR in {}:\n{}".format(syscfg['sysconfig'],'\n'.join(problems)))
        raise IOError("ERRORs in {}, please fix:\n{}".format(syscfg['sysconfig'],'\n'.join(problems)))


def get_phantoms(syscfg,logger):
    syscfg["phantom_defs"] = dict()
    phdir = syscfg["phantoms"]
    for f in os.listdir(phdir):
        if f == "phantom-parameters.mac":
            continue
        logger.debug("checking out {} ...".format(f))
        if f[-4:] == ".mac":
            label = f[:-4]
            syscfg["phantom_defs"][label] = phantom_specs(phdir,label)
            logger.debug("got a full setup for phantom with label {}, yay!".format(label))
        elif f[-4:] == ".cfg":
            mac = f[:-4]+".mac"
            if not os.path.exists(os.path.join(phdir,mac)):
                raise IOError("found orphan cfg file {} in phantoms directory {}: mac file {} not found".format(f,phdir,mac))
            else:
                logger.debug("OK")
        else:
            raise IOError("found file {} of unknown type in phantoms directory {}, please remove or rename".format(f,phdir))

def get_condor_memory_req_fits(syscfg,sysprsr,logger):
    if sysprsr.has_section('condor memory'):
        parser = sysprsr['condor memory']
        for txt,defval in zip(["minimum","default","maximum"],[1500,2000,16000]):
            key='condor memory request {} [MB]'.format(txt)
            syscfg[key] = parser.getfloat(key,defval)
            logger.debug("{} = {}".format(key,syscfg[key]))
        for geo in ["ct","phantom"]:
            for radtype in ["proton","carbon"]:
                key='condor memory fit {} {}'.format(radtype,geo)
                # e.g.: condor memory fit proton ct = offset 500 ct 1e-4 dosegrid 1e-3 nspots 0.
                syscfg[key]=dict(offset=syscfg['condor memory request default [MB]'])
                values=[v.lower() for v in parser.get(key,"").split()]
                for var in ['offset','dosegrid','ct','nspots']:
                    if var in values:
                        i=values.index(var)
                        values.pop(i)
                        syscfg[key][var]=float(values.pop(i))
                if len(values)>0:
                    logger.error("ERROR: unrecognized entries in condor memory fit:".format(" ".join(values)))
                    raise RuntimeError("ERROR: unrecognized entries in condor memory fit:".format(" ".join(values)))
                logger.debug("{} = {}".format(key,syscfg[key]))


def get_mc_stats_settings(syscfg,sysprsr,logger):
    mc_stats_config = {
        #MCStatType.cfglabels[MCStatType.Nions_per_spot]         : [ 100  ,   1000 ,   1000000 , 100  ],
        MCStatType.cfglabels[MCStatType.Nminutes_per_job]        : [   5  ,     20 ,    60*24*7,   5  ], # maximum: one week! :-)
        MCStatType.cfglabels[MCStatType.Nions_per_beam]          : [1000  ,1000000 ,1000000000 ,1000  ],
        MCStatType.cfglabels[MCStatType.Xpct_unc_in_target]      : [   0.1,      1.,        99.,   0.1],
        "n top voxels for mean dose max"                         : 100,
        "dose threshold as fraction in percent of mean dose max" : 50.,
    }
    problems = []
    kdef = "nions per beam"
    if sysprsr.has_section('mc stats'):
        logger.debug("found mc stats section in config file")
        for k,v in sysprsr['mc stats'].items():
            if k.strip().lower()=='n top voxels for mean dose max':
                mc_stats_config[k]=int(v)
                continue
            elif k.strip().lower()=='dose threshold as fraction in percent of mean dose max':
                mc_stats_config[k]=float(v)
                continue
            if k in mc_stats_config.keys():
                logger.debug("found setting for {}".format(k))
                words = v.split()
                if len(words) < 4 or len(words)>5:
                    problems.append('MC stats settings should have 4 or 5 values: '
                                     + 'minval defval maxval stepsize [default], '
                                     + 'but I found {} values in the line for "{}"'.format(len(words),k))
                    continue
                try:
                    # parse the values
                    if k[0]=='n':
                        values=[int(val) for val in words[:4]]
                    else:
                        values=[float(val) for val in words[:4]]
                    # be paranoid
                    if values[0]>values[2]:
                        problems.append("min {} is larger than max {} for '{}'".format(values[0],values[2],k))
                    if values[0]>values[1] or values[1]>values[2]:
                        problems.append("default value {} is not between min {} and max {} for '{}'".format(values[1],values[0],values[2],k))
                    if values[3]<=0:
                        problems.append("step sizes should be positive, got {} for '{}'".format(values[3],k))
                    if values[0]<=0:
                        problems.append("got min={}<=0 for '{}', should be positive".format(values[0],k))
                    # store if OK
                    if not bool(problems):
                        mc_stats_config[k]=values[:]
                        if len(words)==5 and words[4]=='default':
                            kdef=k
                except ValueError as ve:
                    problems.append(ve)
            else:
                problems.append('Unknown mc stats setting "{}" (set to "{}")'.format(k,v))
        if problems:
            logger.error("ERROR in {}:\n{}".format(syscfg['sysconfig'],'\n'.join(problems)))
            raise IOError("ERRORs in {}, please fix:\n{}".format(syscfg['sysconfig'],'\n'.join(problems)))
    else:
        logger.debug('No "mc stats" section in system config file, will use defaults.')
    # store, freeze, never change again
    for k,v in mc_stats_config.items():
        if k in MCStatType.cfglabels:
            syscfg[k]=tuple(v+[k==kdef])
        else:
            syscfg[k]=v

def get_simulation_install(syscfg,sysprsr,logger):
    if not sysprsr.has_section('simulation'):
        raise RuntimeError("missing gate environment shell setting in system configuration")
    simulation = sysprsr['simulation']
    simulation_options = ['gate shell environment',
                          'proton physics list',
                          'ion physics list',
                          'air box margin [mm]',
                          'number of cores',
                          'minimum dose grid resolution [mm]',
                          'rbe factor protons',
                          'remove dose outside external',
                          'gamma index parameters dta_mm dd_percent thr_percent def',
                          'stop on script actor time interval [s]',
                          'htcondor next job start delay [s]',
                          'run gamma analysis',
                          'write mhd unscaled dose',
                          'write mhd scaled dose',
                          'write mhd physical dose',
                          'write mhd rbe dose',
                          'write mhd plan dose',
                          'write dicom physical dose',
                          'write dicom rbe dose',
                          'write dicom plan dose',
                          'write dicom re alanine']
    for k,v in simulation.items():
        if k not in simulation_options:
            msg="unknown option in simulation section of {}: '{}'; recognized options are:\n * {}".format(syscfg['sysconfig'],k,'\n * '.join(simulation_options))
            logger.error(msg)
            raise RuntimeError(msg)
    syscfg['gate_env.sh'] = simulation['gate shell environment']
    syscfg['proton physics list'] = simulation.get('proton physics list','QBBC_EMZ')
    syscfg['ion physics list'] = simulation.get('ion physics list','QBBC_EMZ')
    syscfg["air box margin [mm]"] = simulation.getfloat('air box margin [mm]',10.0)
    syscfg['number of cores'] = simulation.getint('number of cores',10)
    syscfg['rbe factor protons'] = simulation.getfloat('rbe factor protons',1.1)
    syscfg["minimum dose grid resolution [mm]"] = simulation.getfloat("minimum dose grid resolution [mm]")
    # TODO: introduce a new section "output options"?
    syscfg['remove dose outside external'] = simulation.getboolean('remove dose outside external',False)
    syscfg["gamma index parameters dta_mm dd_percent thr_percent def"] = simulation.get("gamma index parameters dta_mm dd_percent thr_percent def","")
    syscfg['stop on script actor time interval [s]'] = simulation.getint('stop on script actor time interval [s]',300)
    syscfg['htcondor next job start delay [s]'] = simulation.getfloat('htcondor next job start delay [s]',1.)
    # TODO: check that SoS actor time interval and next job start delay are not crazy
    syscfg['run gamma analysis']=simulation.getboolean('run gamma analysis',False)
    syscfg['write mhd unscaled dose']=simulation.getboolean('write mhd unscaled dose',False)
    syscfg['write mhd scaled dose']=simulation.getboolean('write mhd scaled dose',False)
    syscfg['write mhd physical dose']=simulation.getboolean('write mhd physical dose',False)
    syscfg['write mhd rbe dose']=simulation.getboolean('write mhd rbe dose',False)
    syscfg['write mhd plan dose']=simulation.getboolean('write mhd plan dose',False)
    syscfg['write dicom physical dose']=simulation.getboolean('write dicom physical dose',True)
    syscfg['write dicom rbe dose']=simulation.getboolean('write dicom rbe dose',False)
    syscfg['write dicom plan dose']=simulation.getboolean('write dicom plan dose',False)
    syscfg['write dicom re alanine']=simulation.getboolean('write dicom re alanine',False)
    # TODO: paranoid GATE version test (should be GateRTion 1.0)
    # TODO: silly density tolerance check (positive, less than 1.0)
    # TODO: check that the physics list is actually an existing one

# taken from G4NistMaterialBuilder, geant4 version 10.03.p03
def get_all_nist_materials():
    return [ "G4_WATER", "G4_H", "G4_He", "G4_Li", "G4_Be", "G4_B", "G4_C", "G4_N", "G4_O", "G4_F", "G4_Ne", "G4_Na", "G4_Mg", "G4_Al", "G4_Si", "G4_P", "G4_S",
"G4_Cl", "G4_Ar", "G4_K", "G4_Ca", "G4_Sc", "G4_Ti", "G4_V", "G4_Cr", "G4_Mn", "G4_Fe", "G4_Co", "G4_Ni", "G4_Cu", "G4_Zn", "G4_Ga", "G4_Ge", "G4_As",
"G4_Se", "G4_Br", "G4_Kr", "G4_Rb", "G4_Sr", "G4_Y", "G4_Zr", "G4_Nb", "G4_Mo", "G4_Tc", "G4_Ru", "G4_Rh", "G4_Pd", "G4_Ag", "G4_Cd", "G4_In", "G4_Sn",
"G4_Sb", "G4_Te", "G4_I", "G4_Xe", "G4_Cs", "G4_Ba", "G4_La", "G4_Ce", "G4_Pr", "G4_Nd", "G4_Pm", "G4_Sm", "G4_Eu", "G4_Gd", "G4_Tb", "G4_Dy", "G4_Ho",
"G4_Er", "G4_Tm", "G4_Yb", "G4_Lu", "G4_Hf", "G4_Ta", "G4_W", "G4_Re", "G4_Os", "G4_Ir", "G4_Pt", "G4_Au", "G4_Hg", "G4_Tl", "G4_Pb", "G4_Bi", "G4_Po",
"G4_At", "G4_Rn", "G4_Fr", "G4_Ra", "G4_Ac", "G4_Th", "G4_Pa", "G4_U", "G4_Np", "G4_Pu", "G4_Am", "G4_Cm", "G4_Bk", "G4_Cf", "G4_A-150_TISSUE", "G4_ACETONE",
"G4_ACETYLENE", "G4_ADENINE", "G4_ADIPOSE_TISSUE_ICRP", "G4_AIR", "G4_ALANINE", "G4_ALUMINUM_OXIDE", "G4_AMBER", "G4_AMMONIA", "G4_ANILINE", "G4_ANTHRACENE",
"G4_B-100_BONE", "G4_BAKELITE", "G4_BARIUM_FLUORIDE", "G4_BARIUM_SULFATE", "G4_BENZENE", "G4_BERYLLIUM_OXIDE", "G4_BGO", "G4_BLOOD_ICRP", "G4_BONE_COMPACT_ICRU",
"G4_BONE_CORTICAL_ICRP", "G4_BORON_CARBIDE", "G4_BORON_OXIDE", "G4_BRAIN_ICRP", "G4_BUTANE", "G4_N-BUTYL_ALCOHOL", "G4_C-552", "G4_CADMIUM_TELLURIDE", "G4_CADMIUM_TUNGSTATE",
"G4_CALCIUM_CARBONATE", "G4_CALCIUM_FLUORIDE", "G4_CALCIUM_OXIDE", "G4_CALCIUM_SULFATE", "G4_CALCIUM_TUNGSTATE", "G4_CARBON_DIOXIDE", "G4_CARBON_TETRACHLORIDE",
"G4_CELLULOSE_CELLOPHANE", "G4_CELLULOSE_BUTYRATE", "G4_CELLULOSE_NITRATE", "G4_CERIC_SULFATE", "G4_CESIUM_FLUORIDE", "G4_CESIUM_IODIDE", "G4_CHLOROBENZENE",
"G4_CHLOROFORM", "G4_CONCRETE", "G4_CYCLOHEXANE", "G4_1,2-DICHLOROBENZENE", "G4_DICHLORODIETHYL_ETHER", "G4_1,2-DICHLOROETHANE", "G4_DIETHYL_ETHER", "G4_N,N-DIMETHYL_FORMAMIDE",
"G4_DIMETHYL_SULFOXIDE", "G4_ETHANE", "G4_ETHYL_ALCOHOL", "G4_ETHYL_CELLULOSE", "G4_ETHYLENE", "G4_EYE_LENS_ICRP", "G4_FERRIC_OXIDE", "G4_FERROBORIDE",
"G4_FERROUS_OXIDE", "G4_FERROUS_SULFATE", "G4_FREON-12", "G4_FREON-12B2", "G4_FREON-13", "G4_FREON-13B1", "G4_FREON-13I1", "G4_GADOLINIUM_OXYSULFIDE", "G4_GALLIUM_ARSENIDE",
"G4_GEL_PHOTO_EMULSION", "G4_Pyrex_Glass", "G4_GLASS_LEAD", "G4_GLASS_PLATE", "G4_GLUCOSE", "G4_GLUTAMINE", "G4_GLYCEROL", "G4_GUANINE", "G4_GYPSUM", "G4_N-HEPTANE",
"G4_N-HEXANE", "G4_KAPTON", "G4_LANTHANUM_OXYBROMIDE", "G4_LANTHANUM_OXYSULFIDE", "G4_LEAD_OXIDE", "G4_LITHIUM_AMIDE", "G4_LITHIUM_CARBONATE", "G4_LITHIUM_FLUORIDE",
"G4_LITHIUM_HYDRIDE", "G4_LITHIUM_IODIDE", "G4_LITHIUM_OXIDE", "G4_LITHIUM_TETRABORATE", "G4_LUNG_ICRP", "G4_M3_WAX", "G4_MAGNESIUM_CARBONATE", "G4_MAGNESIUM_FLUORIDE",
"G4_MAGNESIUM_OXIDE", "G4_MAGNESIUM_TETRABORATE", "G4_MERCURIC_IODIDE", "G4_METHANE", "G4_METHANOL", "G4_MIX_D_WAX", "G4_MS20_TISSUE", "G4_MUSCLE_SKELETAL_ICRP",
"G4_MUSCLE_STRIATED_ICRU", "G4_MUSCLE_WITH_SUCROSE", "G4_MUSCLE_WITHOUT_SUCROSE", "G4_NAPHTHALENE", "G4_NITROBENZENE", "G4_NITROUS_OXIDE", "G4_NYLON-8062", "G4_NYLON-6-6",
"G4_NYLON-6-10", "G4_NYLON-11_RILSAN", "G4_OCTANE", "G4_PARAFFIN", "G4_N-PENTANE", "G4_PHOTO_EMULSION", "G4_PLASTIC_SC_VINYLTOLUENE", "G4_PLUTONIUM_DIOXIDE", "G4_POLYACRYLONITRILE",
"G4_POLYCARBONATE", "G4_POLYCHLOROSTYRENE", "G4_POLYETHYLENE", "G4_MYLAR", "G4_PLEXIGLASS", "G4_POLYOXYMETHYLENE", "G4_POLYPROPYLENE", "G4_POLYSTYRENE", "G4_TEFLON",
"G4_POLYTRIFLUOROCHLOROETHYLENE", "G4_POLYVINYL_ACETATE", "G4_POLYVINYL_ALCOHOL", "G4_POLYVINYL_BUTYRAL", "G4_POLYVINYL_CHLORIDE", "G4_POLYVINYLIDENE_CHLORIDE",
"G4_POLYVINYLIDENE_FLUORIDE", "G4_POLYVINYL_PYRROLIDONE", "G4_POTASSIUM_IODIDE", "G4_POTASSIUM_OXIDE", "G4_PROPANE", "G4_lPROPANE", "G4_N-PROPYL_ALCOHOL", "G4_PYRIDINE",
"G4_RUBBER_BUTYL", "G4_RUBBER_NATURAL", "G4_RUBBER_NEOPRENE", "G4_SILICON_DIOXIDE", "G4_SILVER_BROMIDE", "G4_SILVER_CHLORIDE", "G4_SILVER_HALIDES", "G4_SILVER_IODIDE",
"G4_SKIN_ICRP", "G4_SODIUM_CARBONATE", "G4_SODIUM_IODIDE", "G4_SODIUM_MONOXIDE", "G4_SODIUM_NITRATE", "G4_STILBENE", "G4_SUCROSE", "G4_TERPHENYL", "G4_TESTIS_ICRP",
"G4_TETRACHLOROETHYLENE", "G4_THALLIUM_CHLORIDE", "G4_TISSUE_SOFT_ICRP", "G4_TISSUE_SOFT_ICRU-4", "G4_TISSUE-METHANE", "G4_TISSUE-PROPANE", "G4_TITANIUM_DIOXIDE", "G4_TOLUENE",
"G4_TRICHLOROETHYLENE", "G4_TRIETHYL_PHOSPHATE", "G4_TUNGSTEN_HEXAFLUORIDE", "G4_URANIUM_DICARBIDE", "G4_URANIUM_MONOCARBIDE", "G4_URANIUM_OXIDE", "G4_UREA", "G4_VALINE",
"G4_VITON", "G4_WATER_VAPOR", "G4_XYLENE", "G4_GRAPHITE", "G4_lH2", "G4_lN2", "G4_lO2", "G4_lAr", "G4_lBr", "G4_lKr", "G4_lXe", "G4_PbWO4", "G4_Galactic", "G4_GRAPHITE_POROUS",
"G4_LUCITE", "G4_BRASS", "G4_BRONZE", "G4_STAINLESS-STEEL", "G4_CR39", "G4_OCTADECANOL", "G4_KEVLAR", "G4_DACRON", "G4_NEOPRENE", "G4_CYTOSINE", "G4_THYMINE", "G4_URACIL",
"G4_DEOXYRIBOSE", "G4_DNA_DEOXYRIBOSE", "G4_DNA_PHOSPHATE", "G4_DNA_ADENINE", "G4_DNA_GUANINE", "G4_DNA_CYTOSINE", "G4_DNA_THYMINE", "G4_DNA_URACIL", "G4_BODY" ]

def get_all_gate_materials(syscfg,logger):
    materials_parser = configparser.ConfigParser()
    materials_parser.optionxform = str # keys (material names) should be case sensitive
    matdb = os.path.join(syscfg['commissioning'], syscfg['materials database'])
    if not os.path.exists(matdb):
        raise RuntimeError("ERROR: missing file {} in directory: {}".format(syscfg['materials database'],
                                                                            syscfg['commissioning']) )
    with open(matdb,"r") as fp:
        materials_parser.read_file(fp)
    logger.debug("gate materials: '{}'".format("', '".join(list(materials_parser['Materials'].keys()))))
    return list(materials_parser['Materials'].keys())

def get_materials(syscfg,sysprsr,logger):
    hutol = 0.001
    minimal_nistlist = ["G4_WATER","G4_AIR"]
    if not sysprsr.has_section('materials'):
        raise RuntimeError("ERROR: 'materials' section missing in system config file {}".format(syscfg['sysconfig']))
    materials_section = sysprsr['materials']
    gatedb = materials_section.get('materials database',"GateMaterials.db")
    syscfg['materials database'] = gatedb
    syscfg['ct override list'] = dict()
    syscfg['hu density tolerance [g/cm3]'] = materials_section.getfloat('hu density tolerance [g/cm3]',hutol)
    hutol = syscfg['hu density tolerance [g/cm3]']
    all_gate_materials = get_all_gate_materials(syscfg,logger)
    all_nist_materials = get_all_nist_materials()
    ALL_GATE_MATERIALS = [ k.upper() for k in all_gate_materials ]
    ALL_NIST_MATERIALS = [ k.upper() for k in all_nist_materials ]
    for key in materials_section.keys():
        KEY=key.upper()
        if 'HU DENSITY' in KEY:
            continue
        if 'MATERIALS DATABASE' == KEY:
            continue
        if KEY in ALL_NIST_MATERIALS:
            material = all_nist_materials[ALL_NIST_MATERIALS.index(KEY)]
            dbname='NIST'
        elif KEY in ALL_GATE_MATERIALS:
            material = all_gate_materials[ALL_GATE_MATERIALS.index(KEY)]
            dbname=os.path.basename(gatedb)
        else:
            raise KeyError("Unknown entry {} in ct override list section, not in NIST nor in {}".format(key,os.path.basename(gatedb)))
        density = materials_section.getfloat(key)
        logger.debug("allowing {} (defined by {}) in the CT override list".format(material,dbname))
        syscfg['ct override list'][material] = density
    logger.debug('I found a "materials" section in {}, got HUtol={} and override list = {}'.format(syscfg['sysconfig'],hutol,syscfg['ct override list'].keys()))

def get_tmp_correction_factors(syscfg,sysprsr,logger):
    syscfg['(tmp) correction factors']=dict(default=1.0)
    if '(tmp) correction factors' not in sysprsr.sections():
        logger.debug("no '(tmp) correction factors' section in sysconfig file, correction factor is always 1.0")
        return
    tmp_section = sysprsr['(tmp) correction factors']
    src_props = glob(os.path.join(syscfg['beamlines'],"*","*_*_source_properties.txt"))
    for key in tmp_section.keys():
        value = tmp_section.getfloat(key)
        logger.debug("got correction factor {} for: '{}'".format(value,key))
        k_e_y = "_".join(key.split()).lower()
        if k_e_y == "default":
            def_value = tmp_section.getfloat(key)
            if def_value > 0.:
                syscfg['(tmp) correction factors'][k_e_y] == def_value
            else:
                raise ValueError("ERROR in {}: default correction factor '{}' is not positve".format(syscfg['sysconfig'],def_value))
            logger.debug("overriding default correction factor to: {}".format(syscfg['(tmp) correction factors']["default"]))
        else:
            for src_prop in src_props:
                if os.path.basename(src_prop).lower() == (k_e_y+"_source_properties.txt"):
                    # got a match
                    corr_value = tmp_section.getfloat(key)
                    if corr_value > 0.:
                        logger.debug("BINGO: {} matches {}".format(key,src_prop))
                        syscfg['(tmp) correction factors'][k_e_y] = corr_value
                    else:
                        raise ValueError("ERROR in {}: correction factor '{}' for beamline/radtype '{}' is not positve".format(syscfg['sysconfig'],corr_value,key))
                    break
                else:
                    logger.debug("{} does not match {}".format(key,src_prop))
            if k_e_y not in  syscfg['(tmp) correction factors'].keys():
                raise KeyError("ERROR in {}: correction factor for beamline/radtype '{}' that does not match any source properties file in the beamlines directory {}".format(syscfg['sysconfig'],key,syscfg['beamlines']))

def get_msw_scaling(syscfg,sysprsr,logger):
    syscfg['msw scaling']=dict(default=[1.0,0.0])
    if 'msw scaling' not in sysprsr.sections():
        logger.debug("no 'msw scaling' section in sysconfig file, no msw scaling applied")
        return
    tmp_section = sysprsr['msw scaling']
    src_props = glob(os.path.join(syscfg['beamlines'],"*","*_*_source_properties.txt"))
    for key in tmp_section.keys():
        value = tmp_section.get(key)
        logger.debug("got msw sclaing slpe and offset {} for: '{}'".format(value,key))
        k_e_y = "_".join(key.split()).lower()
        if k_e_y == "default":
            s = tmp_section.get(key)
            def_value = [float(i) for i in s.split()]
            syscfg['msw scaling'][k_e_y] == def_value
            logger.debug("overriding default  msw sclaing slpe and offset to: {}".format(syscfg['msw scaling']["default"]))
        else:
            for src_prop in src_props:
                if os.path.basename(src_prop).lower() == (k_e_y+"_source_properties.txt"):
                    # got a match
                    s = tmp_section.get(key)
                    corr_value = [float(i) for i in s.split()]
                    logger.debug("BINGO: {} matches {}".format(key,src_prop))
                    syscfg['msw scaling'][k_e_y] = corr_value
                else:
                    logger.debug("{} does not match {}".format(key,src_prop))
            if k_e_y not in  syscfg['msw scaling'].keys():
                raise KeyError("ERROR in {}: correction factor for beamline/radtype '{}' that does not match any source properties file in the beamlines directory {}".format(syscfg['sysconfig'],key,syscfg['beamlines']))


def get_user_roles(syscfg,sysprsr,a):
    user_roles=dict()
    user_aliases=dict()
    if 'user roles' not in sysprsr.sections():
        raise RuntimeError("Section 'user roles' is missing from system config file {}, please tell the admin (or fix it, if you are the admin).".format(syscfg['sysconfig']))
    user_roles = sysprsr['user roles']
    all_names=list()
    for k in user_roles.keys():
        names=k.replace(","," ").strip().split()
        username=names[0]
        for name in names:
            if name in all_names:
                raise RuntimeError("ERROR in 'user roles' section in {}: name {} is used multiple times.".format(syscfg['sysconfig'],name))
            user_aliases[name]=username
            #logger.debug("accepting username/alias {}".format(name))
        all_names += names
        role=user_roles[k]
        roles = ["admin","clinical","commissioning"]
        if role not in roles:
            raise RuntimeError("ERROR in 'user roles' section in {}: role '{}' for user '{}' is unknown, please use one of: {}.".format(syscfg['sysconfig'],role,username,roles))
        user_roles[username]=role
    syscfg["authorized users"] = user_aliases
    if a in user_aliases:
        username = user_aliases[a]
        syscfg["username"]=username
        syscfg["role"]=user_roles[username]
    else:
        raise RuntimeError("ERROR: user/alias '{}' is unknown, not listed in the 'user roles' section of the sysconfig file '{}', please specify one with the '-l' option and/or talk to the admin.".format(a,syscfg['sysconfig']))

def get_logging(syscfg,system_parser,want_logfile="default"):
    global logger
    if not bool(want_logfile):
        logging.basicConfig(level=syscfg['default logging level'])
        logger = logging.getLogger(__name__)
        return
    msg=""
    try:
        # try to get the logging directory before anything else
        logdir=system_parser['directories']['logging']
        if not os.path.isdir(logdir):
            raise IOError(f"logging dir '{logdir}' is not an existing directory?")
        msg = f"got logdir={logdir} from system config file"
        # TODO: maybe we should also check here that we can actually write something to this directory
    except Exception as e:
        msg = f"WARNING: failed to get a valid log directory from your system configuration: '{e}'."
        logdir='/tmp'
    if os.path.isabs(want_logfile):
        logger,logfilepath = get_dual_logging( level  = syscfg['default logging level'],
                                               daemon_file = want_logfile )
    else:
        logger,logfilepath = get_dual_logging( level  = syscfg['default logging level'],
                                               prefix = os.path.join(logdir, syscfg['username']+"_"+os.path.basename(syscfg['cmd'])))
    syscfg["log file path"]=logfilepath
    if logdir == '/tmp':
        logger.warn(msg)
    else:
        logger.debug(msg)
    return logger


def get_sysconfig(filepath=None,verbose=False,debug=False,username=None,want_logfile="default"):
    """
    This function parses the "system configuration file", checks the contents against
    trivial mistakes and then creates the ``system_configuration`` singleton object.
    :param filepath: file path of the system configuration file, use ``dirname(script)/../cfg/system.cfg`` by default
    :param verbose: boolean, write DEBUG level log information to standard output if True, INFO level if False.
    :param debug: boolean, clean up the bulky temporary data if False, leave it for debugging if True.
    :param username: who is running this, with what role?
    :param want_logfile: possible values are empty string (no logfile), absolute file path (implying no logs to stdout) or 'default' (generated log file path, some logs will be streamed to standard output)
    """

    # IDEAL directories: where is the currently running script installed?
    this_cmd = os.path.realpath(__file__) #sys.argv[0]
    impl_dir = os.path.dirname(this_cmd)
    ideal_dir = os.path.dirname(impl_dir)
    install_dir = os.path.dirname(ideal_dir)
    bin_dir = os.path.join(install_dir,"bin")
    cfg_dir = os.path.join(install_dir,"cfg")
    system_cfg_path = os.path.join(cfg_dir,'system.cfg')
    if not username:
        username = getpass.getuser()

    # the ``syscfg`` dictionary will be used to initialize the system_configuration singleton object
    syscfg = {"cmd":this_cmd,
              "default logging level": logging.DEBUG if verbose else logging.INFO,
              "IDEAL home":install_dir,
              "bindir":bin_dir,
              "username":username, # temporary
              "debug":debug,
              "config dir":cfg_dir, 
              "log file path": ''}

    syscfg['sysconfig'] = filepath if filepath else system_cfg_path
    # read and parse the contents of system config file
    if not os.path.exists(syscfg['sysconfig']):
        msg = "ERROR: system config file {} does not exist. If you just did a fresh install of IDEAL, please run the 'first_install.py' script first generate it, or alternatively use the '-s' option to specify the system configuration file that should be used.".format(syscfg['sysconfig'])
        raise RuntimeError(msg)
    system_parser = configparser.ConfigParser()
    system_parser.read(syscfg['sysconfig'])
    get_user_roles(syscfg,system_parser,username) # find out the 'real' user name and role
    syscfg["logdir"] = system_parser['directories']['logging']
    
    if want_logfile != "default":
        logger = get_logging(syscfg,system_parser,want_logfile)
#        log_filepath = syscfg["logdir"]+'/sysconfig_'+timestamp()
#        logger = create_logger('syslog',log_filepath)
    else:
        log_filepath = syscfg["logdir"]+'/sysconfig_'+timestamp()
        logger = create_logger('syslog',log_filepath)

    if filepath:
        logger.debug("You provided sysconfig={}".format(filepath))
    else:
        logger.debug("I will use the default sysconfig file {}".format(system_cfg_path))

    # now read all sections, check the settings, store everything in the syscfg dictionary
    get_basedirs(syscfg,system_parser,logger)
    get_commissioning_dirs(syscfg,logger)
    get_mc_stats_settings(syscfg,system_parser,logger)
    get_simulation_install(syscfg,system_parser,logger)
    get_condor_memory_req_fits(syscfg,system_parser,logger)
    get_phantoms(syscfg,logger)
    get_materials(syscfg,system_parser,logger)
    get_tmp_correction_factors(syscfg,system_parser,logger)
    get_msw_scaling(syscfg,system_parser,logger)

    # now create the singleton system configuration object (and return a reference)
    return system_configuration(syscfg)
    



# vim: set et softtabstop=4 sw=4 smartindent:
