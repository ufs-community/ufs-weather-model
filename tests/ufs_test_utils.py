import os
import sys
import glob
import yaml
import subprocess

def delete_files(deletefiles):
    fileList = glob.glob(deletefiles, recursive=True)    
    for filePath in fileList:
        try:
            os.remove(filePath)
        except OSError:
            print("Error while deleting ",deletefiles)
        
def link_new_baselines():
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
    test_dep = None
    for test in val:
        case, config = get_testcase(test)
        if case == casename:
            test_dep = {case:config}
    return test_dep

def get_testcase(test):
    case_name = None
    case_config = None
    for case, configs in test.items():
        case_name=case
        case_config=configs
    return case_name, case_config
    
def write_logfile(logfile, openmod, output="", subproc=""):
    with open(logfile, openmod) as rtlog:
        if (not subproc == "") :
            subprocess.call(subproc, shell=True, stdout=rtlog)
        if (not output == "") :
            rtlog.writelines(output)
    rtlog.close()

def rrmdir(path):
    for entry in os.scandir(path):
        if entry.is_dir():
            rrmdir(entry)
        else:
            os.remove(entry)
            os.rmdir(path)

