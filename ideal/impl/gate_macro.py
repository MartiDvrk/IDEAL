# -----------------------------------------------------------------------------
#   Copyright (C): MedAustron GmbH, ACMIT Gmbh and Medical University Vienna
#   This software is distributed under the terms
#   of the GNU Lesser General  Public Licence (LGPL)
#   See LICENSE for further details
# -----------------------------------------------------------------------------

"""
This module defines the function which writes a Gate macro very similar to the
one that Alessio composed for his PSQA simulations.

Some commands which were executed by Alessio from an external macro file are
now copied here directly, in an effort to completely eliminate the use of
aliases. Writing all beam numbers and beamset names explicitly hopefully
facilitates later debugging.
"""

import os
import time
import logging
logger=logging.getLogger(__name__)
from impl.idc_details import MCStatType
from impl.system_configuration import system_configuration
import numpy as np

def roman_year():
    """
    Vanity function to print the year in Roman numerals.
    """
    y=int(time.strftime("%Y"))
    if y<2000 or y>=2100:
        raise ValueError("only 21st century")
    decade=["","X","XX","XXX","XL","L","LX","LXX","LXXX","XC"][(y//10)%10]
    year=["","I","II","III","IV","V","VI","VII","VIII","IX"][y%10]
    return "MM{}{}".format(decade,year)

def check(ct=True,**kwargs):
    """
    Function to check completeness and correctness of a long list of key word arguments for gate macro composition.
    """
    missing = list()
    undefined = list()
    mandatory = [ "beamset", "uid", "spotfile", "physicslist",
                  "isoC", "beamline", "beamnr", "beamname", "radtype",
                  "rsids", "rmids", "dose_nvoxels" ]
    if ct:
        mandatory += ["ct_mhd", "ct_bb", "mod_patient_angle", "gantry_angle",
                      "HU2mat", "HUmaterials", "dose_center", "dose_size"]
        # ct_mhd is the filepath of the resized & HU-overridden CT image to be used by Gate (path only; file does not yet exist when Gate macros are written)
        # ct_resized is the ITK image object that is resized (padded and cropped), to be used for geometry definitions
    else:
        mandatory += ["phantom"]
    #stattype = kwargs.get("mcstattype",-1)
    #if MCStatType.Xpct_unc_in_target == stattype:
    #    mandatory += ["unc_goal"]
    #elif MCStatType.Nminutes_per_job == stattype:
    #    mandatory += ["minutes_per_job"]
    for k in mandatory:
        if k not in kwargs.keys():
            missing.append(k)
    for k in kwargs.keys():
        if k not in mandatory:
            undefined.append(k)
    if bool(missing) or bool(undefined):
        raise RuntimeError("missing: {} undefined: {}".format(",".join(missing),",".join(undefined)))


def write_gate_macro_file(ct=True,**kwargs):
    """
    This function writes the main macro for PSQA calculations for a specific beam in a beamset.
    beamset: name of the beamset
    uid: DICOM identifier of the full plan/beamset
    user: string that identifies the person who is responsible for this calculation
    spotfile: text file 
    """
    logger.debug("going to write mac files, got {} keywords".format(len(kwargs.keys())))
    check(ct,**kwargs)
    logger.debug("basic checks went fine")
    syscfg = system_configuration.getInstance()
    beamline = kwargs["beamline"]
    kwargs["beamlinename"]=beamline.name
    kwargs["nmaxprimaries"]=2**31-1 # 2147483647, max number for a signed integer = max number of primaries that GATE can simulate
    radtype = kwargs["radtype"]
    isoC = kwargs["isoC"]
    # want Uncertainty?
    # wantU = (MCStatType.Xpct_unc_in_target == kwargs["mcstattype"])
    #coreMinutes = (MCStatType.Nminutes_per_job == kwargs["mcstattype"])
    #if coreMinutes:
    #    kwargs["coreSeconds"] = 60*int(kwargs["minutes_per_job"])
    #kwargs["wantU"] = str(wantU).lower()
    dose2w = True
    re_alanine = syscfg['write dicom re alanine']
    if ct:
        logger.debug("CT geometry for beamline {}".format(beamline.name))
        kwargs["geometry"] = "CT"
        ct_bb     = kwargs["ct_bb"]
        dose_center   = kwargs["dose_center"]
        dose_size     = kwargs["dose_size"]
        # in the dose actor, we need to give the position relative to the center of the CT image
        rot_box_size = 2.0001*np.max(np.abs(np.stack([ct_bb.mincorner-isoC,ct_bb.maxcorner-isoC])),axis=0)
        kwargs["xboxsize"] = float(rot_box_size[0])
        kwargs["yboxsize"] = float(rot_box_size[1])
        kwargs["zboxsize"] = float(rot_box_size[2])
        # for CT we always request "dose to water"
        dose2w = True
        kwargs["wantDose"]    = "false"
        kwargs["wantDose2W"]  = "true"
        kwargs["label"]="CT-{beamset}-B{beamnr}-{beamname}-{beamlinename}".format(**kwargs)
        #kwargs["massmhd"]=os.path.basename(kwargs["ct_mhd"]).replace(".mhd","_mass.mhd")
    else:
        # for a water box we request normal dose to material
        # for other materials (e.g. PMMA) we request "dose to water"
        phantom = kwargs["phantom"] # this is an object of class phantom_specs, see system_configuration.py
        dose2w = phantom.dose_to_water # this is a boolean
        kwargs["wantDose2W"]  = str(dose2w).lower()
        kwargs["wantDose"]    = str(not dose2w).lower()
        kwargs["geometry"]    = "PHANTOM"
        kwargs["phantom_macfile"] = os.path.join("data","phantoms",phantom.label+".mac")
        kwargs["phantom_name"]    = phantom.label
        # MFA 11/10/2022
        #kwargs["phantom_move_x"] = -1.*isoC[0]
        #kwargs["phantom_move_y"] = -1.*isoC[1]
        #kwargs["phantom_move_z"] = -1.*isoC[2]
        kwargs["label"]="PHANTOM-{phantom_name}-{beamset}-B{beamnr}-{beamname}".format(**kwargs)
        logger.debug("PHANTOM geometry for beamline {}".format(beamline.name))
    dosemhd="idc-{label}.mhd".format(**kwargs)
    dosedose="idc-{label}-".format(**kwargs)
    dosedose+="DoseToWater" if dose2w else "Dose"
    logger.debug("dose actor MHD output file name is {}".format(dosemhd))
    kwargs["dosemhd"]=dosemhd
    kwargs["dosedosemhd"]=dosedose+".mhd"
    kwargs["dosedoseraw"]=dosedose+".raw"
    kwargs["stop_on_script_every_n_seconds"] = syscfg["stop on script actor time interval [s]"]
    dose_nvoxels=kwargs["dose_nvoxels"]
    kwargs["dnx"]=dose_nvoxels[0]
    kwargs["dny"]=dose_nvoxels[1]
    kwargs["dnz"]=dose_nvoxels[2]
    kwargs["materialsdb"]=syscfg["materials database"]
    kwargs["fulldate"]=time.strftime("%A %d %B ")+roman_year()
    kwargs["hhmmss"]=time.strftime("%H:%M:%S")
    kwargs["yyyymmdd"]=time.strftime("%Y-%m-%d")
    kwargs["user"]=syscfg["username"]
    kwargs["isocx"]=isoC[0]
    kwargs["isocy"]=isoC[1]
    kwargs["isocz"]=isoC[2]
    kwargs["srcprops"] = os.path.basename(beamline.source_properties_file(radtype))
    if radtype.upper() != "PROTON":
        ion_details = radtype.split("_")
        if len(ion_details) != 4 or ion_details[0] != "ION":
            logger.error('Radiation type string should be either "PROTON" or of the form "ION_Z_A_Q"')
            raise RuntimeError("unsupported and/or ill-formatted radiation type {}".format(radtype))
        kwargs["ionZ"]=ion_details[1]
        kwargs["ionA"]=ion_details[2]
        kwargs["ionQ"]=ion_details[3]

    ########## HEADER ###########
    header = """
###############################################################################
# Gate-based IDC for PSQA on CT or PHANTOM geometries
#
###############################################################################
#
# Geometry          = {geometry}
# BeamSet           = {beamset}
# UID               = {uid}
# User              = {user}
# Beam Number       = {beamnr}
# Beam Name         = {beamname}
# Treatment Machine = {beamlinename}
# Spot file         = {spotfile}
# Date (full)       = {fulldate}
# Date              = {yyyymmdd}
# Time              = {hhmmss}
#
###############################################################################

/control/execute {{VISUMAC}}
""".format(**kwargs)
    logger.debug("defined header section")

    
    ########## MATERIALS ###########
    materials_section = """
#=====================================================
#= MATERIALS
/gate/geometry/setMaterialDatabase data/{materialsdb}
""".format(**kwargs)
    logger.debug("defined material section")

    ########## GEOMETRY ###########
    geometry_section = """
#=====================================================
#= GEOMETRY
#=====================================================

# WORLD
/gate/world/setMaterial Air
/gate/world/geometry/setXLength 5.0 m
/gate/world/geometry/setYLength 5.0 m
/gate/world/geometry/setZLength 5.0 m

# Make sure that G4_WATER is included in the geometry (and hence the dynamic database)
# so that "DoseToWater" works.
/gate/world/daughters/name water_droplet
/gate/world/daughters/insert box
/gate/water_droplet/setMaterial G4_WATER
/gate/water_droplet/geometry/setXLength 1 mm
/gate/water_droplet/geometry/setYLength 1 mm
/gate/water_droplet/geometry/setZLength 1 mm
/gate/water_droplet/placement/setTranslation -2499.5 -2499.5 -2499.5 mm
/gate/water_droplet/vis/setColor cyan
"""

    if beamline.beamline_details_mac_file:
        geometry_section += """
/control/execute mac/{}
""".format(os.path.basename(beamline.beamline_details_mac_file))
        logger.debug("defined beamline part of geometry section")

    for rs in kwargs["rsids"]:
        geometry_section += """
/control/execute mac/{}
""".format(os.path.basename(beamline.rs_details_mac_file(rs)))
        logger.debug("added mac file for RS={}".format(rs))
    for rm in kwargs["rmids"]:
        geometry_section += """
/control/execute mac/{}
""".format(os.path.basename(beamline.rm_details_mac_file(rm)))
        logger.debug("added mac file for rifi={}".format(rm))

    if ct:
        # TODO: correct rotation axis and order for gantry & patient angles
        # TranslateTheImageAtThisIsoCenter -> the image will be placed such that this isocenter is at position (0,0,0) of the mother volume
        geometry_section += """
# Patient virtual container: for couch rotation
/gate/world/daughters/name                      patient_box
/gate/world/daughters/insert                    box
/gate/patient_box/setMaterial Air
/gate/patient_box/geometry/setXLength {xboxsize} mm
/gate/patient_box/geometry/setYLength {yboxsize} mm
/gate/patient_box/geometry/setZLength {zboxsize} mm
/gate/patient_box/vis/setColor        yellow

# CT Volume
/gate/patient_box/daughters/name                      patient
/gate/patient_box/daughters/insert                    ImageNestedParametrisedVolume
/gate/geometry/setMaterialDatabase              {HUmaterials}
/gate/patient/geometry/setHUToMaterialFile      {HU2mat}
/gate/patient/geometry/setImage                 {ct_mhd}
/gate/patient/geometry/TranslateTheImageAtThisIsoCenter {isocx} {isocy} {isocz} mm
/gate/patient_box/placement/setRotationAxis 0 1 0
/gate/patient_box/placement/setRotationAngle {mod_patient_angle} deg
""".format(**kwargs)
        logger.debug("added HU material generator and voxelized image for CT")
    else:
        geometry_section += """
/control/alias phantom_name {phantom_name}
/control/execute {phantom_macfile}
""".format(**kwargs)
#/control/alias {phantom_name}_move_x {phantom_move_x}  #MFA 11/10/2022
#/control/alias {phantom_name}_move_y {phantom_move_y}
#/control/alias {phantom_name}_move_z {phantom_move_z}

    ########## PHYSICS ###########
    physics_section = """
#=====================================================
#= PHYSICS
/gate/physics/addPhysicsList {physicslist}
#=====================================================
#= GATE BASIC PARAMETERS
""".format(**kwargs)
    logger.debug("physics = {}".format(kwargs["physicslist"]))

    if ct:
        physics_section += """
/control/execute data/CT/ct-parameters.mac
"""
    else:
        physics_section += """
/control/execute data/phantoms/phantom-parameters.mac
"""

    ########## OUTPUT ###########
    output_section = """
#=====================================================
#= OUTPUTS
#=====================================================

#=====================================================
# Statistics actor
#=====================================================

/gate/actor/addActor SimulationStatisticActor stat
/gate/actor/stat/saveEveryNSeconds       120
/gate/actor/stat/save {{OUTPUTDIR}}/statActor-{label}.txt

#=====================================================
# Dose Actor
#=====================================================

/gate/actor/addActor    DoseActor               dose
/gate/actor/dose/save                     {{OUTPUTDIR}}/{dosemhd}
""".format(**kwargs)
    if ct:
        output_section += """
/gate/actor/dose/attachTo                 patient
""".format(**kwargs)
    else:
        output_section += """
/gate/actor/dose/attachTo                 {phantom_name}
""".format(**kwargs)
    output_section += """
/gate/actor/dose/setResolution            {dnx} {dny} {dnz}
/gate/actor/dose/stepHitType              random
/gate/actor/dose/enableEdep               false
/gate/actor/dose/enableUncertaintyEdep    false
/gate/actor/dose/enableDose               true
/gate/actor/dose/enableDoseToWater        {wantDose2W}
""".format(**kwargs)
# Add LET actor for alanine relative efficency calculation
    if re_alanine:
        alanine_data_path = '/usr/local/GATE/GateRTion-1.1/install/data/RE_Alanine/'
        kwargs['alanine_data_path'] = alanine_data_path
        re_alanine_name = dosedose + '_RE_alanine'
        kwargs['re_alanine_name'] = re_alanine_name
        output_section += """
#=====================================================
# LET Actor
#=====================================================
/gate/actor/addActor    LETActor               alanine_RE
/gate/actor/alanine_RE/save                     {{OUTPUTDIR}}/{re_alanine_name}.mhd
/gate/actor/alanine_RE/attachTo                 patient 
/gate/actor/alanine_RE/setResolution            {dnx} {dny} {dnz}
/gate/actor/alanine_RE/stepHitType              random
/gate/actor/alanine_RE/setType               RE
/gate/actor/alanine_RE/setOtherMaterial    G4_ALANINE
/gate/actor/alanine_RE/doParallelCalculation True
/gate/actor/alanine_RE/setDoseEfficiencyFileAndAtomicNumber {alanine_data_path}/RE_Z1.txt 1 
/gate/actor/alanine_RE/setDoseEfficiencyFileAndAtomicNumber {alanine_data_path}/RE_Z2.txt 2 
/gate/actor/alanine_RE/setDoseEfficiencyFileAndAtomicNumber {alanine_data_path}/RE_Z3.txt 3 
/gate/actor/alanine_RE/setDoseEfficiencyFileAndAtomicNumber {alanine_data_path}/RE_Z4.txt 4 
/gate/actor/alanine_RE/setDoseEfficiencyFileAndAtomicNumber {alanine_data_path}/RE_Z5.txt 5 
/gate/actor/alanine_RE/setDoseEfficiencyFileAndAtomicNumber {alanine_data_path}/RE_Z6.txt 6 
/gate/actor/alanine_RE/addFilter                   particleFilter
/gate/actor/alanine_RE/particleFilter/addParticleZ 1
/gate/actor/alanine_RE/particleFilter/addParticleZ 2
/gate/actor/alanine_RE/particleFilter/addParticleZ 3
/gate/actor/alanine_RE/particleFilter/addParticleZ 4
/gate/actor/alanine_RE/particleFilter/addParticleZ 5
/gate/actor/alanine_RE/particleFilter/addParticleZ 6

""".format(**kwargs)
    logger.debug("defined output section")

    noop_filepath="mac/stop_on_script_{label}_noop.mac".format(**kwargs)
    with open(noop_filepath,"w") as noop:
        noop.write("""/control/echo NOT adding any stop-on-script actor (QT test run)\n""")
    stop_on_script_filepath="mac/stop_on_script_{label}.mac".format(**kwargs)
    with open(stop_on_script_filepath,"w") as sosmac:
        sosmac.write( """
/control/echo "ADDING stop on script actor"
/gate/actor/addActor                      StopOnScriptActor stopWhenReady
/gate/actor/stopWhenReady/saveEveryNSeconds {stop_on_script_every_n_seconds}
/gate/actor/stopWhenReady/saveAllActors     true
/gate/actor/stopWhenReady/save              mac/check_the_flags_{label}.sh
""".format(**kwargs))
        with open("mac/check_the_flags_{label}.sh".format(**kwargs),"w") as chksh:
            chksh.write("""
workdir=$(dirname $(realpath tmp))
echo "workdir is $workdir"
stopflag="STOP_{dosedosemhd}"
if test -r "$workdir/$stopflag" ; then
    echo "found STOP flag $stopflag, exit with return value 1 to stop simulation"
    exit 1
fi
echo "did not find STOP flag $stopflag, returning 0 to continue simulation"
d=output.$clusterid.$procid
ls -ld $d
ls -lrt $d
ls -ld tmp/$d
ls -lrt tmp/$d
ls "$d/{dosedosemhd}" "$d/{dosedoseraw}" "$d/statActor-{label}.txt"
type locked_copy.py
locked_copy.py -v -d tmp/$d/ -l tmp/$d/{dosedosemhd}.lock "$d/{dosedosemhd}" "$d/{dosedoseraw}" "$d/statActor-{label}.txt"
exit 0
""".format(**kwargs))
    output_section += """
/control/strif {{RUNMAC}} == mac/run_all.mac {}
/control/strif {{RUNMAC}} == mac/run_qt.mac {}
""".format(stop_on_script_filepath,noop_filepath)
#    elif coreMinutes:
#        Nseconds_filepath="mac/stop_after_N_seconds_{label}.mac".format(**kwargs)
#        with open(Nseconds_filepath,"w") as Nsec_mac:
#            Nsec_mac.write( """
#/control/echo adding stop on script actor for runtime threshold
#/gate/actor/addActor                            StopOnScriptActor stopAfterNSeconds
#/gate/actor/stopAfterNSeconds/saveEveryNSeconds {coreSeconds}
#/gate/actor/stopAfterNSeconds/saveAllActors     true
#/gate/actor/stopAfterNSeconds/save              mac/stop_after_N_seconds_{label}.sh
#""".format(**kwargs))
#        with open("mac/stop_after_N_seconds_{label}.sh".format(**kwargs),"w") as chksh:
#            chksh.write("""
#echo "clusterid=$clusterid"
#echo "procid=$procid"
#echo "time is up!"
#date
#condor_q $clusterid.$procid -run -nobatch -allusers -global
#ddhhmmss=$(condor_q $clusterid.$procid -run -nobatch -allusers -global | tail -1 | awk '{{print $5}}')
#dd=$(echo $ddhhmmss | cut -f1 -d+ )
#hhmmss=$(echo $ddhhmmss | cut -f2 -d+ )
#hh=$(echo $hhmmss | cut -f1 -d: )
#mm=$(echo $hhmmss | cut -f2 -d: )
#ss=$(echo $hhmmss | cut -f3 -d: )
#hh=$[$hh+24*$dd]
#mm=$[$mm + 60*$hh]
#ss=$[$ss + 60*$mm + 30]
#mm2=$[$ss / 60]
#echo ran for mm=$mm mm2=$mm2 minutes
#exit 1
#echo this message will never be printed
#""".format(**kwargs))
#        output_section += """
#/control/strif {{RUNMAC}} == mac/run_all.mac {}
#/control/strif {{RUNMAC}} == mac/run_qt.mac {}
#""".format(Nseconds_filepath,noop_filepath)


    ########## INIT ###########
    init_section = """
#=====================================================
#= INITIALIZATION
#=====================================================
/gate/run/initialize
/gate/physics/print {{OUTPUTDIR}}/PHYSICS.txt
""".format(**kwargs)

    ########## SOURCE ###########
    source_section = """
#=====================================================
#= SOURCE
#=====================================================
/gate/source/addSource PBS TPSPencilBeam
/gate/source/PBS/setTestFlag false
"""
    if radtype.lower() == "proton":
        logger.debug("configuring protons")
        source_section += """
/gate/source/PBS/setParticleType proton
""".format(**kwargs)
    else:
        logger.debug("configuring ions (probably carbons)")
        source_section += """
/gate/source/PBS/setParticleType GenericIon
/gate/source/PBS/setIonProperties {ionZ} {ionA} {ionQ}
""".format(**kwargs)

    source_section += """
/gate/source/PBS/setPlan {spotfile}

# One beam at the time...
/gate/source/PBS/setAllowedFieldID {beamnr}

/gate/source/PBS/setSourceDescriptionFile data/{srcprops}
/gate/source/PBS/setSpotIntensityAsNbIons true
/gate/source/PBS/setBeamConvergenceXTheta {{{radtype}_CONVERGENCE_XTHETA}}
/gate/source/PBS/setBeamConvergenceYPhi {{{radtype}_CONVERGENCE_YPHI}}
""".format(**kwargs)

    logger.debug("defined source section")

    ########## RUN ###########
    run_section = """
#=====================================================
#= RUN
#=====================================================
/gate/random/setEngineName MersenneTwister
/control/execute {{RUNMAC}}

/gate/application/start
""".format(**kwargs)
    logger.debug("defined run section")


    ########## ALIAS ###########

    #main_fname=os.path.join("mac","main-{geometry}-{beamset}-{beamname}-{beamnr}-{beamlinename}.mac".format(**kwargs))
    main_fname=os.path.join("mac",dosemhd.replace(".mhd",".mac"))
    f=open(main_fname,"w")
    f.write(header)
    f.write(materials_section)
    f.write(geometry_section)
    f.write(physics_section)
    f.write(output_section)
    f.write(init_section)
    f.write(source_section)
    f.write(run_section)
    f.close()

    logger.debug("wrote main mac file: {}".format(main_fname))

    fname=os.path.join("mac","run_all.mac")
    f=open(fname,"w")
    f.write("""
/gate/random/setEngineSeed {{RNGSEED}}
/gate/application/setTotalNumberOfPrimaries {nmaxprimaries}
""".format(**kwargs))
    f.close()

    fname=os.path.join("mac","run_qt.mac")
    f=open(fname,"w")
    f.write("""
/gate/random/setEngineSeed 42
/gate/application/setTotalNumberOfPrimaries 100
""")
    f.close()

    fname=os.path.join("mac","visu.mac")
    f=open(fname,"w")
    f.write("""
/vis/enable
/vis/open OGLIQt
/vis/drawVolume
/vis/viewer/flush
/vis/scene/add/trajectories
/vis/scene/endOfEventAction accumulate
/vis/scene/add/axes            0 0 0 50 mm
/vis/scene/add/text            1 0 0 cm  20 0 0   X
/vis/scene/add/text            0 1 0 cm  20 0 0   Y
/vis/scene/add/text            0 0 1 cm  20 0 0   Z
/vis/viewer/set/viewpointThetaPhi 70 30
/vis/viewer/zoom 5
""")
    f.close()

    fname=os.path.join("mac","novisu.mac")
    f=open(fname,"w")
    f.write("""
/vis/disable
""")
    f.close()

    logger.debug("wrote auxiliary mac files")

    logger.debug("DONE: wrote all mac files")
    return main_fname,dosemhd

# vim: set et softtabstop=4 sw=4 smartindent:
