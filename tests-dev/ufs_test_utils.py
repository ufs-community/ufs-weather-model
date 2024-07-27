import os
import sys
import glob
import yaml
import shutil
import subprocess

def sync_testscripts():
    """symlink sharable rt.sh test scripts
    """
    dst= os.getcwd()
    src= os.path.split(os.getcwd())[0]+'/tests'    
    for name in os.listdir(src):
        src_name= src +'/'+ name
        dst_name= dst +'/'+ name
        if not os.path.exists(dst_name):
            if "/compile.sh" in dst_name:
                shutil.copyfile(src_name, dst_name)
                subprocess.call(['chmod', '755', dst_name])
                with open(dst_name) as rfile:
                    buildsh = rfile.read().replace("${PATHTR}/tests/", "${PATHTR}/tests-dev/")
                    rfile.close()
                with open(dst_name, "w") as wfile:
                    wfile.write(buildsh)
                    wfile.close()
            else:
                os.symlink(src_name, dst_name)
                
def machine_check_off(machine_id, val):
    """Check turned-off machine from yaml configuration

    Args:
        machine_id (str): local machine name
        val (dic): build and test config dictionary list
    Returns:
        pass_machine: logical flag to pass local machine
    """
    pass_machine = True
    if 'turnoff' in val.keys():
        if machine_id in val['turnoff']:
            pass_machine = False
    if 'turnon' in val.keys():
        if not machine_id in val['turnon']:
            pass_machine = False
    return pass_machine

def delete_files(deletefiles):
    """Remove specified filepath

    Args:
        deletefiles (str): filepath to remove e.g. tests/rocoto.*
    """
    fileList = glob.glob(deletefiles, recursive=True)    
    for filePath in fileList:
        try:
            os.remove(filePath)
        except OSError:
            print("Error while deleting ",deletefiles)
        
def link_new_baselines():
    """Create symlinks for newly generated baselines.
    """ 
    USER = str(os.environ.get('USER'))
    MACHINE_ID = os.getenv('MACHINE_ID')        
    PATHRT     = os.getenv('PATHRT')
    with open("baseline_setup.yaml", 'r') as f:
        exp_config = yaml.load(f) #, Loader=yaml.FullLoader)
        base  = exp_config[MACHINE_ID]
        DISKNM= str(base['DISKNM'])
        STMP  = str(base['STMP'])
        PTMP  = str(base['PTMP'])
        path  = STMP+'/'+USER
        RTPWD = path + '/FV3_RT/REGRESSION_TEST'
        f.close()
    #--- capture user's NEW_BASELINE location ----
    logfile    = PATHRT+'/logs/RegressionTests_'+MACHINE_ID+'.log'
    with open(logfile,'r') as flog:
        logheads= flog.readlines()
        for line in logheads:   
            if "BASELINE DIRECTORY:" in line:
                NEW_BASELINE=line.split(" ")[1]
                break
        flog.close()
    #--- symlink verified baseline cases to users new baseline ---
    os.environ["RTPWD"] = RTPWD
    os.environ["NEW_BASELINE"] = NEW_BASELINE
    symlink_baselines = subprocess.Popen(['bash', '-c', '. ufs_test_utils.sh; link_new_baselines'])
    symlink_baselines.wait()
    
def get_testdep(casename,val):
    """Retrieve test case dependencies

    Args:
        casename (str): Test case name
        val (dict): Test case attributes e.g. val['compiler']

    Returns:
        dict: Test case and config for the specified dependency
    """    
    test_dep = None
    for test in val:
        case, config = get_testcase(test)
        if case == casename:
            test_dep = {case:config}
    return test_dep

def get_testcase(test):
    """Retrieve test case names and configs from given dict from pyaml

    Args:
        test (dict): dict retrieved from reading in yaml test file

    Returns:
        str, dict: test name and python dict of test configuration
    """
    case_name = None
    case_config = None
    for case, configs in test.items():
        case_name=case
        case_config=configs
    return case_name, case_config
    
def write_logfile(logfile, openmod, output="", subproc=""):
    """Append given output into log file

    Args:
        logfile (str): Log filename
        openmod (str): mode to open file in
        output (str): Content to append to log file. Defaults to "".
        subproc (str): Command to run within the shell. Defaults to "".
    """
    with open(logfile, openmod) as rtlog:
        if (not subproc == "") :
            subprocess.call(subproc, shell=True, stdout=rtlog)
        if (not output == "") :
            rtlog.writelines(output)
    rtlog.close()

def rrmdir(path):
    """Remove all files and directories in specified path.

    Args:
        path (str): File path to remove
    """ 
    for entry in os.scandir(path):
        if entry.is_dir():
            rrmdir(entry)
        else:
            os.remove(entry)
            os.rmdir(path)

