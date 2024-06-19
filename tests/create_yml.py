import os
import re
import yaml
from ufs_test_utils import get_testcase, get_testdep

def update_testyaml(input_list):
    """Generate temporary test yaml based on list of tests received

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
    """Update test yaml file for a single test specified in -n <test_name> <compiler>
    """
    try:
        SRT_NAME     = str(os.getenv('SRT_NAME'))
        SRT_COMPILER = str(os.getenv('SRT_COMPILER'))
    except NameError:
        print("*** SRT_NAME or SRT_COMPILER are not given with runtime option -n! ***")
    input_list=[SRT_NAME+" "+SRT_COMPILER]
    update_testyaml(input_list)

def update_testyaml_b():
    """Update test yaml file for tests specified in -b <file>
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
    """Strip out RUN or COMPILE whitespace and separate with commas.

    Args:
        str_in (str): RUN or COMPILE line read in from rt.conf

    Returns:
        str: whitespace stripped and comma separated values
    """
    return "'"+("','".join(str_in.split()))+"'"

def parse_line(str_in):
    """Parse rt.conf line into list

    Args:
        str_in (str): RUN or COMPILE line from rt.conf

    Returns:
        list: list of RUN or COMPILE test attributes
    """
    build_attr = " ".join(str_in.split()).split('|')
    build_attr = [attr.strip() for attr in build_attr]
    return build_attr

if __name__ == "__main__":
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

