"""
This module serves to define the available Hounsfield lookup tables (HLUT) as defined
in the materials/HLUT/hlut.conf configuration file in the commissioning directory of
an IDEAL installation.

In the 'hlut.conf' configuration file, each section represents a supported CT
protocol and the data in the section should unambiguously describe for which
kinds of CTs the protocol is relevant and how to use it to define the
Gate/Geant4 material to use for each HU value in a given CT image with that
protocol.

There are two kinds of protocols: Schneider protocols which are intended for
both clinical and commissioning purposes.

For the Schneider type protocols the 'hlut.conf' file should specify two text
files, providing the "density curve" and the "material composition",
respectively. The density curve typically specifies the density for about a
dozen HU values between -1024 and +3000. The composition file specifies for
several HU intervals a mixture of elements. Gate/Geant4 has a method for
converting these two tables for a given density tolerance value into a "HU to
material" table, along with a database file that has the same format as the
"Gate material database" file (usually named "GateMaterials.db") which defines
lots of "interpolated materials", each with slightly different density, but the
composition exactly like one of the materials in the "composition" file.

The conversion of the Schneider density and composition tables into the
HU-to-materials table and the associated new database is taken care of by the
``gate_hlut_cache`` module, which saves the new table and database in a cache
directory such that they can be reused.

The protocols can either be directly selected through their name or through
a match with the metadata included in the CT DICOM files.

A match with DICOM metadata allows the "automatic detection" of the CT
protocol.  In ``the hlut.conf`` file, the commissioning physicist lists one or
more DICOM tags (by keyword or by hexadecimal tag) and which text strings the
DICOM values are expected to have for that CT protocol. It is the
responsibility of the commissioning physicist to define this in such a way that
for all data that will be used with this particular IDEAL installation, the
automatic selection will always result in one single match.  In case multiple
DICOM tags are given per CT protocol, the protocol with the most matches
"wins".

The user may also provide a name or keyword in the IDEAL command invocation
(e.g. the -c option from ``clidc.py``) which is then directly matched with the
CT protocol names.  First a case sensitive perfect match is attempted, if that
fails then a case insensitive partial match (the given keyword may be a
substring of the full CT protocoll name). In case the partial match results in
multiple results, an exception is thrown and the user should try again.
"""

from impl.system_configuration import system_configuration
from impl.gate_hlut_cache import generate_hlut_cache, hlut_cache_dir
import os
import configparser
import pydicom
import logging
import hashlib
logger=logging.getLogger(__name__)


class hlut:
    """
    This class serves to provide the HU to material conversion tables, based on
    the information from the HLUT configuration file.  This can be either a
    "Schneider" type conversion (Schneider tables, interpolated with a density
    tolerance value to a large list of interpolated materials) or a
    "commissioning" type conversion in which the HU to materials are directly
    coded by the commissioning physicist.
    """
    def __init__(self,name,prsr_section,hutol=None):
        syscfg = system_configuration.getInstance()
        self.name = name
        self.cache_dir = None
        self.dicom_match = dict()
        non_dicom_keys = list()
        if "density" in prsr_section and "composition" in prsr_section:
            # assuming this is a Schneider type HLUT
            logger.debug(f"{name} is a Schneider-type CT protocol")
            self.type = "Schneider"
            d = prsr_section["density"]
            c = prsr_section["composition"]
            self.density = os.path.join(syscfg["CT/density"],d)
            self.composition = os.path.join(syscfg["CT/composition"],c)
            if not os.path.exists(self.density):
                raise FileNotFoundError("For Schneider protocol '{name}' the density file '{d}' cannot be found.")
            if not os.path.exists(self.composition):
                raise FileNotFoundError("For Schneider protocol '{name}' the composition file '{c}' cannot be found.")
            self.hutol = syscfg['hu density tolerance [g/cm3]'] if hutol is None else float(hutol)
            non_dicom_keys += ["density","composition"]
        else:
            # assuming this is a commissioning type HLUT
            logger.debug(f"{name} is a CT protocol with direct HU-to-material tables (commissioning)")
            allowed_materials = syscfg['ct override list'].keys()
            self.type = "Commissioning"
            self.hu2mat_lines = list()
            self.density_lines = list()
            hulast=None
            for k in prsr_section.keys():
                if not "," in k:
                    continue
                try:
                    hufromstr,hutillstr = k.split(",")
                    hufrom = int(hufromstr)
                    if hulast is not None:
                        if hufrom != hulast:
                            raise ValueError(f"HU intervals are not consecutive for {name}: gap between HU={hulast} and HU={HUfrom}")
                    hutill = int(hutillstr)
                    if not (hutill>hufrom):
                        raise ValueError(f"wrong HU interval [{hufrom},{hutill}) for {name}: upper bound should be larger than lower bound")
                    material = str(prsr_section[k])
                    if material not in allowed_materials:
                        raise ValueError(f"material '{material}' not included in list of override materials in system configuration file '{syscfg['sysconfig']}'")
                    self.hu2mat_lines.append(f"{hufrom} {hutill} {material}")
                    density = syscfg['ct override list'][material]
                    self.density_lines.append(f"{hufrom} {hutill} {density}")
                    non_dicom_keys.append(k)
                except ValueError as ve:
                    logger.error(f"something wrong with Hounsfield configuration {name}, please fix: '{ve}'")
                    raise
            if len(self.density_lines)<2:
                raise ValueError(f"CT protocol definition {self.name} incomplete/incorrect: need either a density & composition file (Schneider) or at least two HU-interval to material lines.")
        for k,v in prsr_section.items():
            if k in non_dicom_keys:
                continue
            dk = str(k).replace(" ","")
            if dk in pydicom.datadict.keyword_dict.keys():
                m = str(v).strip()
                if m == "":
                    raise ValueError("DICOM match criterion cannot be empty or only white space")
                logger.debug(f"Got DICOM match line: {k} = {m}")
                self.dicom_match[dk] = m
            else:
                # k is not a DICOM keyword, not a Schneider table nor a HU interval
                raise KeyError(f"I do not know what to do with this hlut.conf line: '{k} = {v}'")
    def get_density_file(self):
        if self.type == "Schneider":
            return self.density
        elif self.type == "Commissioning":
            # instead of a density curve, we deliver a density table
            if not bool(self.cache_dir):
                self.get_hu2mat_files()
            hu3density = os.path.join(self.cache_dir,"hu3density.txt")
            if not os.path.exists(hu3density):
                with open(hu3density,"w") as fp:
                    fp.write("\n".join(self.density_lines))
                    fp.write("\n")
            return hu3density
        else:
            raise RuntimeError("OOPSIE: programming error")
    def get_hu2mat_files(self,hutol=None):
        if hutol is not None:
            self.hutol = hutol
        if self.type == "Schneider":
            self.cache_dir = hlut_cache_dir(self.density,self.composition,self.hutol,create=False)
            if self.cache_dir is None:
                # this will run Gate to run the db/hu2mat generation if necessary
                ok,self.cache_dir = generate_hlut_cache(self.density,self.composition,self.hutol)
                if not ok:
                    raise RuntimeError(f"{self.name} failed to create HU2material tables in cache directory {self.cache_dir}.")
            humatdb = os.path.join(self.cache_dir,'patient-HUmaterials.db')
            hu2mattxt = os.path.join(self.cache_dir,'patient-HU2mat.txt')
        elif self.type == "Commissioning":
            h4sh = hashlib.md5()
            for line in self.hu2mat_lines:
                h4sh.update(bytes(line,encoding='utf-8'))
            syscfg = system_configuration.getInstance()
            self.cache_dir = os.path.join(syscfg['CT/cache'],h4sh.hexdigest())
            os.makedirs(self.cache_dir,exist_ok=True)
            hu2mattxt = os.path.join(self.cache_dir,'commissioning-HU2mat.txt')
            humatdb = os.path.join(self.cache_dir,'commissioning-HUmaterials.db')
            if not os.path.exists(hu2mattxt):
                with open(hu2mattxt,"w") as hu2matfp:
                    # write hu2mat cache file
                    hu2matfp.write("\n".join(self.hu2mat_lines))
                    hu2matfp.write("\n")
            if not os.path.exists(humatdb):
                with open(humatdb,"w"):
                    # write empty file
                    pass
        else:
            raise RuntimeError("OOPSIE: programming error")
        assert(os.path.exists(humatdb))
        assert(os.path.exists(hu2mattxt))
        return hu2mattxt,humatdb
    def match(self,ct):
        """
        For now, check that all requirements match.
        :param ct: a DCM data set as read with pydicom.dcmread.
        """
        if len(self.dicom_match)==0:
            logger.debug(f"Skipping '{self.name}': no DICOM match items given, can only be selected manually.")
            return False
        logger.debug(f"Trying to match the '{self.name}' CT protocol")
        print(self.dicom_match)
        for k,m in self.dicom_match.items():
            if k in ct:
                ctk=str(ct[k].value)
                if ctk==m:
                    logger.debug(f"OK: DICOM tag {k} = {m} matches.")
                else:
                    logger.debug(f"Nope: DICOM tag {k} mismatch: {ctk}!={m}.")
                    return False
            else:
                logger.debug(f"Nope: DICOM tag {k} missing.")
                return False
        logger.debug(f"Match!")
        return True

class hlut_conf:
    """
    This is a 'singleton' class, only one HLUT configuration object is supposed
    to exist. The singleton HLUT configuration object behaves like a Python
    dictionary with limited functionality: no comparisons (since there is only
    one such object), no additions or removals, but the elements are not
    strictly 'const'. The HLUT configuration is initialized the first time its
    instance is acquired with the static `getInstance()` method.
    """
    __instance = None
    @staticmethod
    def getInstance(fname=None):
        """
        Get a reference to the one and only "HLUT configuration" instance, if
        it exists.  The 'fname' argument should ONLY be used by the unit tests.
        In normal use, the HLUT conf file should be defined by the system
        configuration.
        """
        if __class__.__instance is None:
            return __class__.__read_hlut_conf(fname)
        if fname is not None:
            # this should not happen
            # TODO: better throw a runtime error?
            logging.warn("Ignoring HLUT conf file '{fname}' (programming error?)")
        return __class__.__instance
    @staticmethod
    def __read_hlut_conf(fname=None):
        """
        Use the ``configparser`` module to read the config file, then try to
        convert each section into a HLUT definition.  The 'fname' argument
        should ONLY be different from None during unit tests.
        """
        hlut_parser = configparser.ConfigParser()
        hlut_parser.optionxform = str # keys (CT protocol names) should be case sensitive
        syscfg = system_configuration.getInstance()
        if fname is None:
            fname = os.path.join(syscfg['CT'], "hlut.conf" )
        if not os.path.exists(fname):
            raise RuntimeError(f"ERROR: cannot find the HLUT configuration file '{fname}'.")
        all_hluts=dict()
        with open(fname,"r") as fp:
            hlut_parser.read_file(fp)
            for s in hlut_parser.sections():
                try:
                    h=hlut(s,hlut_parser[s])
                    all_hluts[s]=h
                except Exception as e:
                    logger.error(f"OOPSIE skipping section '{s}' in hlut.conf because: '{e}'")
        return hlut_conf(all_hluts)
    def __init__(self,s=dict()):
        """
        Initialize the singleton object, just once.
        """
        if __class__.__instance is not None:
            raise RuntimeError("HLUT configuration initialized more than once?")
        self.__all_hluts = s
        __class__.__instance = self
    def __contains__(self,name):
        return self.__all_hluts.__contains__(name)
    def __iter__(self):
        return self.__all_hluts.__iter__()
    def __len__(self):
        return self.__all_hluts.__len__()
    def __getitem__(self,name):
        """
        Get a copy of any (available) system configuration item.
        """
        if name in self.__all_hluts:
            # return a reference
            return self.__all_hluts[name]
        raise KeyError(f"CT protocol name not found: '{name}'")
    def keys(self):
        return self.__all_hluts.keys()
    def values(self):
        return self.__all_hluts.values()
    def items(self):
        return self.__all_hluts.items()
    def hlut_match_keyword(self,kw):
        """
        This method compares a user-provided keyword (e.g. from the -c option for
        the ``clidc.py`` script) to the names of the CT protocols provided in
        the hlut.conf file.

        In case of success, it returns the matching CT protocol name.

        In case of failure, this function will throw a ``KeyError`` with some explaining
        text that should help the user understand what is going wrong.
        """
        if kw in self.__all_hluts:
            # exact match
            ctprotocolname = str(kw)
        else:
            allmatches = [ k for k in self.__all_hluts.keys() if kw.lower() in k.lower() ]
            if len(allmatches) == 0:
                ctlist='\n'.join(list(self.__all_hluts.keys()))
                raise KeyError(f"HLUT keyword '{kw}' does not seem to match any of these supported CT protocol names: {ctlist}")
            elif len(allmatches) > 1:
                ctlist='\n'.join(list(allmatches))
                raise KeyError(f"HLUT keyword '{kw}' matches more than one of the supported CT protocol names: {ctlist}")
            ctprotocolname = str(allmatches[0])
        return ctprotocolname
    def hlut_match_dicom(self,ctslice):
        """
        This method tries to find the HLUT for which the DICOM metadata matching
        rules uniquely matches user-provided CT file object (a data set returned by
        pydicom.dcmread, e.g. first CT slice in the series coming with the DICOM
        treatment plan file).

        In case of success, it returns the CT protocol name.

        In case of failure, this function will throw a ``KeyError`` with some explaining
        text that should help the user understand what is going wrong.
        """
        matches = [hpl for hpl,cthulhu in self.__all_hluts.items() if cthulhu.match(ctslice)]
        logger.debug(f"found {len(matches)} matches")
        if len(matches) == 0:
            raise KeyError("Failed to find any match for CT protocol, please check and fix 'hlut.conf'!")
        if len(matches)>1:
            matchtxt="\n".join(matches)
            raise KeyError(f"Failed to find unique match for CT protocol, please check and fix 'hlut.conf'! These protocols all match: {matchtxt}")
        return matches[0]


#######################################################################
# TESTING
#######################################################################
import unittest
import tempfile
from impl.system_configuration import get_sysconfig


class test_good_hlut_conf(unittest.TestCase):
    def setUp(self):
        # bunch of test CT files (each from a different imaginary image series...)
        syscfg = get_sysconfig(filepath='./cfg/unit_test_system.cfg',username="foobar_admin",want_logfile="")
        self.ct1 = pydicom.Dataset()
        self.ct1.SeriesDescription="Foo Bar"
        self.ct1.KVP="100"
        self.ct1.ConvolutionKernel="UB"
        self.ct2 = pydicom.Dataset()
        self.ct2.SeriesDescription="Bar Foo"
        self.ct2.KVP="110"
        self.ct2.ConvolutionKernel="UB"
        self.ct3 = pydicom.Dataset()
        self.ct3.SeriesDescription="foo barr"
        self.ct3.KVP="120"
        self.ct3.ConvolutionKernel="UB"
        self.ct4 = pydicom.Dataset()
        self.ct4.SeriesDescription="quux"
        self.ct4.KVP="210"
        self.ct4.ConvolutionKernel="BU"
        self.hlut_conf = tempfile.NamedTemporaryFile(prefix="test_hlut",suffix=".conf",delete=False)
        self.hlut_conf.write(b"""
[First test CT protocol]
density = Schneider2000DensitiesTable.txt
composition = Schneider2000MaterialsTable.txt
Series Description = Foo Bar

[Second test CT protocol]
-1024,-50 = G4_AIR
-50,5000 = G4_GRAPHITE
Series Description = Bar Foo
Convolution Kernel = UB
KVP = 110

[Third test CT protocol]
-1024,-50 = G4_AIR
-50,50 = G4_WATER
50,3000 = G4_STAINLESS-STEEL
Series Description = foo barr
""")
        self.hlut_conf.close()
        self.all_hluts = hlut_conf.getInstance(fname=self.hlut_conf.name)
    def tearDown(self):
        os.unlink(self.hlut_conf.name)
    def test_hlut_conf_dict_properties(self):
        self.assertEqual(3,len(self.all_hluts))
