"""
API for POP's configure
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import run_cmd_no_fail, expect

import glob, shutil
logger = logging.getLogger(__name__)

def configure_pop(case):

    exeroot = case.get_value("EXEROOT")
    srcroot = case.get_value("SRCROOT")
    caseroot = case.get_value("CASEROOT")

    # determine configure directory 
    configure_dir = os.path.join(srcroot,"components","pop","bld")
    if os.path.isfile(os.path.join(caseroot,"SourceMods","src.pop","configure")):
        configure_dir = os.path.join(caseroot,"SourceMods","src.pop")

    # create $CASEROOT/Buildconf/popconf if it does not exist
    popconf_dir = os.path.join(caseroot,"Buildconf","popconf")
    if not os.path.exists(popconf_dir):
        os.makedirs(popconf_dir)

    # run configure from $CASEROOT/Buildconf/popconf
    # running configure writes out the following xml variables in env_build.xml
    # POP_BLCKX, POP_BLCKY, POP_MXBLCKS, POP_DECOMPTYPE, POP_NX_BLOCKS, POP_NY_BLOCKS, POP_CPPDEFS

    command = os.path.join(configure_dir,"configure")
    cmd = "%s %s" %(command, caseroot)
    rc, out, err = run_cmd(cmd, from_dir=popconf_dir)

    expect(rc==0,"Command %s failed rc=%d\nout=%s\nerr=%s"%(command,rc,out,err))
    if out is not None:
        logger.info("     %s"%out)
    if err is not None:
        logger.info("     %s"%err)

    # Verify that config_cache.xml exists
    if not os.path.isfile(os.path.join(popconf_dir,"config_cache.xml")):
        expect(False, "config_cache.xml is missing after configure command")

    # since env_build.xml has been updated above - must update case with the new file
    case.read_xml()

    pop_cppdefs = case.get_value("POP_CPPDEFS")
    return pop_cppdefs
