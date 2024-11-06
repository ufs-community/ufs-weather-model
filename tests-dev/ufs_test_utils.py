import os
import sys
import re
import glob
import yaml
import shutil
import subprocess

def update_testyaml(input_list):
    """Generates temporary test YAML based on list of tests received

    Args:
        input_list (list): list of tests to run
    """
    UFS_TEST_YAML = "ufs_test.yaml" # default ufs_test.yaml
    new_yaml = {}
    yaml_item_count = None
    with open(UFS_TEST_YAML, 'r') as file_yaml:
        rt_yaml = yaml.load(file_yaml)#, Loader=yaml.FullLoader)
        for apps, jobs in rt_yaml.items():
            app_temp    = None
            build_temp  = None            
            for key, val in jobs.items():
                if (str(key) == 'build'):
                    #--- build information ---
                    build_val   = val
                    compiler_val= val['compiler']
                if (str(key) == 'tests'):
                    #--- serach for test cases given with -n or -b option ---
                    test_list   = []
                    temp_list   = []
                    app_temp    = None
                    build_temp  = None
                    test_temp   = None
                    test_temp_dep = None
                    for test in val:
                        case, config = get_testcase(test)
                        i=0
                        ilist= None
                        #--- search input_list test cases from ufs_test.yaml ---
                        for line in input_list:
                            case_check    = line.split(" ")[0]
                            compiler_check= line.split(" ")[1]
                            if case == case_check and compiler_val == compiler_check:
                                ilist=i
                                app_temp  = apps
                                build_temp= build_val
                                test_temp = {case:config}
                                temp_list.append(str(case))
                                if 'dependency' in config.keys():
                                    if not str(config['dependency']) in temp_list:
                                        test_temp_dep = get_testdep(str(config['dependency']),val)
                            i+=1
                        #--- pop input_list element if a test case is found ---
                        if not ilist is None:
                            input_list.pop(ilist)
                        #--- append test cases to new test list ---
                        if not test_temp_dep is None:
                            test_list.append(test_temp_dep)
                            test_temp_dep = None
                        if not test_temp is None:
                            test_list.append(test_temp)
                            test_temp = None
            if not app_temp is None:
                new_yaml[app_temp]={'build':build_temp,'tests':test_list}
            #--- check all search is done for input_list ---
            if len(input_list) == 0:
                break
        #--- dump into temporary test yaml file ---
        if len(new_yaml) > 0:
            yaml_item_count = len(new_yaml)
        try:
            yaml_item_count
        except NameError:
            print("*** Test cases given with runtime options -n or -b are not found in ufs_test.yaml! ***")
        else:
            with open(r'ufs_test_temp.yaml', 'w') as yaml_file:
                outputs = yaml.dump(new_yaml, yaml_file)
                yaml_file.close()
        file_yaml.close()

def update_testyaml_n():
    """Updates test YAML file for a single test specified in ``-n <test_name> <compiler>``
    """
    try:
        SRT_NAME     = str(os.getenv('SRT_NAME'))
        SRT_COMPILER = str(os.getenv('SRT_COMPILER'))
    except NameError:
        print("*** SRT_NAME or SRT_COMPILER are not given with runtime option -n! ***")
    input_list=[SRT_NAME+" "+SRT_COMPILER]
    update_testyaml(input_list)

def update_testyaml_b():
    """Updates test YAML file for tests specified in ``-b <file>``
    """
    NEW_BASELINES_FILE = str(os.getenv('NEW_BASELINES_FILE'))
    input_list=[]
    with open(NEW_BASELINES_FILE) as input_file:
        for line in input_file:
            line=line.strip('\n')
            line=line.strip()
            input_list.append(str(line))
        input_file.close()
    update_testyaml(input_list)

def string_clean(str_in):
    """Strips out RUN or COMPILE whitespace and separates with commas.

    Args:
        str_in (str): RUN or COMPILE line read in from ``rt.conf``

    Returns:
        str: Whitespace stripped and comma separated values
    """
    return "'"+("','".join(str_in.split()))+"'"

def parse_line(str_in):
    """Parses ``rt.conf`` line into list

    Args:
        str_in (str): RUN or COMPILE line from rt.conf

    Returns:
        build_attr: List of RUN or COMPILE test attributes
    """
    build_attr = " ".join(str_in.split()).split('|')
    build_attr = [attr.strip() for attr in build_attr]
    return build_attr

def create_yaml():
    """Parses default ``rt.conf`` into ``ufs_test.yaml``

    """
    with open('ufs_test.yaml', 'w') as yaml_file, open("rt.conf") as conf_file:
        for line in conf_file:
            line = line.strip()
            if not line:  # skip: line is blank
                continue
            if line.startswith("#"):  # skip: comment line
                continue
            if line.startswith("COMPILE"):  # COMPILE line
                build = parse_line(line)
                apps = build[1]
                compiler = f"'{build[2]}'"
                options = f"'{build[3]}'"
                machine = build[4]
                off_machine = None
                on_machine = None
                if (machine.find('-') != -1):
                    off_machine = machine.replace("-", "").strip()
                    off_machine = string_clean(off_machine)
                if (machine.find('+') != -1):
                    on_machine = machine.replace("+", "").strip()
                    on_machine = string_clean(on_machine)                
                yaml_file.write(f"{apps}_{build[2].strip()}:\n")
                yaml_file.write(f"  build: \n")
                yaml_file.write(f"    compiler: {compiler}\n")
                yaml_file.write(f"    option: {options}\n")
                if not (off_machine is None):
                    yaml_file.write(f"    turnoff: [{off_machine}]\n")
                if not (on_machine is None):
                    yaml_file.write(f"    turnon: [{on_machine}]\n")
                prev_line = 'COMPILE'
            if line.startswith("RUN"):  # RUN line
                build = parse_line(line)
                test = build[1]
                machine = build[2]
                baseline = f"'{build[3]}'"
                depend = build[4]
                if (machine.find('-') != -1):
                    off_machine = machine.replace("-", "").strip()
                    off_machine = string_clean(off_machine)
                if (machine.find('+') != -1):
                    on_machine = machine.replace("+", "").strip()
                    on_machine = string_clean(on_machine)
                tests = f"    - {test}: {{'project':['daily']"
                if baseline.isalnum():
                    tests += f",'baseline': {baseline}"
                if depend and depend.strip():
                    tests += f",'dependency':'{depend}'"
                if not (off_machine is None):
                    tests += f",'turnoff':[{off_machine}]"
                if not (on_machine is None):
                    tests += f",'turnon':[{on_machine}]"                    
                if prev_line == "COMPILE":
                    yaml_file.write("  tests: \n")
                yaml_file.write(tests+"}\n")
                prev_line = 'RUN'

    yaml_file.close(); conf_file.close()    

def sync_testscripts():
    """Symlinks sharable ``rt.sh`` test scripts
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
    """Checks turned-off machine from YAML configuration

    Args:
        machine_id (str): Local machine name
        val (dict): Build and test config dictionary list
    Returns:
        pass_machine: Logical flag to pass local machine
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
    """Removes specified filepath

    Args:
        deletefiles (str): filepath to remove, e.g., ``tests/rocoto.*``
    """
    fileList = glob.glob(deletefiles, recursive=True)    
    for filePath in fileList:
        try:
            os.remove(filePath)
        except OSError:
            print("Error while deleting ",deletefiles)
        
def link_new_baselines():
    """Creates symlinks for newly generated baselines
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
    """Retrieves test case dependencies

    Args:
        casename (str): Test case name
        val (dict): Test case attributes, e.g., val['compiler']

    Returns:
        test_dep: Dictionary with test case and configuration for the specified dependency
    """    
    test_dep = None
    for test in val:
        case, config = get_testcase(test)
        if case == casename:
            test_dep = {case:config}
    return test_dep

def get_testcase(test):
    """Retrieves test case names and configs from given dictionary from PyYAML

    Args:
        test (dict): Dictionary retrieved from reading in YAML test file

    Returns:
        case_name, case_config: Test name (str) and Python dictionary of test configuration
    """
    case_name = None
    case_config = None
    for case, configs in test.items():
        case_name=case
        case_config=configs
    return case_name, case_config
    
def write_logfile(logfile, openmod, output="", subproc=""):
    """Appends given output into log file

    Args:
        logfile (str): Log filename
        openmod (str): Mode to open file in
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
    """Removes all files and directories in specified path.

    Args:
        path (str): File path to remove
    """
    shutil.rmtree(path)
    #for entry in os.scandir(path):
    #    if entry.is_dir():
    #        rrmdir(entry)
    #    else:
    #        os.remove(entry)
    #        os.rmdir(path)

#if __name__ == "__main__":
#    create_yaml()
