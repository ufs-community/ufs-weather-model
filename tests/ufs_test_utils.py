import os
import sys
import subprocess

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
#def link_newbaseline():
    #CREATE_BASELINE
    #NEW_BASELINES_FILE
    ## IF -c AND -b; LINK VERIFIED BASELINES TO NEW_BASELINE
    #if CREATE_BASELINE == 'true' and not NEW_BASELINES_FILE == '':
    #    for dir in "${RTPWD}"/*/; do
    #    dir=${dir%*/}
    #[[ -d "${NEW_BASELINE}/${dir##*/}" ]] && continue
    #ln -s "${dir%*/}" "${NEW_BASELINE}/"
  #done
