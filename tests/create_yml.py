import os
import re
import yaml
from ufs_test_utils import get_testcase

def update_testyaml():
    try:
        SRT_NAME     = str(os.getenv('SRT_NAME'))
        SRT_COMPILER = str(os.getenv('SRT_COMPILER'))
    except NameError:
        print("*** SRT_NAME or SRT_COMPILER are not given with runtime option -n! ***")
    UFS_TEST_YAML = "ufs_test.yaml"
    test_dep = None
    test_temp_dep = None
    with open(UFS_TEST_YAML, 'r') as file_yaml:
        rt_yaml = yaml.load(file_yaml)
        for apps, jobs in rt_yaml.items():
            for key, val in jobs.items():
                if (str(key) == 'build'):
                    #--- build information ---
                    build_val   = val
                    rt_compiler = val['compiler']
                if (str(key) == 'tests'):
                    #--- serach for test case given -n option ---
                    for test in val:
                        case, config = get_testcase(test)
                        print(case,SRT_NAME, rt_compiler, SRT_COMPILER)
                        if case == SRT_NAME and rt_compiler == SRT_COMPILER:
                            app_temp  = apps
                            build_temp= build_val
                            test_temp = {case:config}
                            if 'dependency' in config.keys():
                                test_dep = str(config['dependency'])
                    #--- check if dependency exists with serached test case ---
                    if test_dep is not None:
                        for test in val:
                            case, config = get_testcase(test)
                            if case == test_dep:
                                test_temp_dep = {case:config}
        file_yaml.close()
    #--- dump into temporary test yaml file ---
    try:
        app_temp
    except NameError:
        print("*** Test case given with runtime option -n is not found in ufs_test.yaml! ***")
    else:
        test_list = []
        if test_temp_dep is not None:
            test_list.append(test_temp_dep)
        test_list.append(test_temp)
        test_yaml = {app_temp:{'build':build_temp,'tests':test_list}}
        with open(r'ufs_test_temp.yaml', 'w') as yaml_file:
            outputs = yaml.dump(test_yaml, yaml_file)
            yaml_file.close()
            print(test_yaml)

def string_clean(str_in):
    str_in=str_in.replace("  "," ").strip()
    str_in=str_in.replace("  "," ").strip()
    str_in=str_in.replace("  "," ").strip()
    str_in=str_in.replace("  "," ").strip()
    str_out="'"+str_in.replace(" ","','")+"'"
    return str_out

if __name__ == "__main__":
    with open('ufs_test.yaml', 'w') as yaml_file, open("rt.conf") as conf_file:
        for line in conf_file:
            line = line.strip()
            if not line:  # skip: line is blank
                continue
            if line.startswith("#"):  # skip: comment line
                continue
            if line.startswith("COMPILE"):  # COMPILE line
                build    = line.split('|')
                apps     =     build[1].strip()
                compiler = "'"+build[2].strip()+"'"
                options  = "'"+build[3].strip()+"'"
                machine  =     build[4].strip()
                off_machine = None
                if (machine.find('-') != -1):
                    off_machine=machine.replace("-","").strip()
                    off_machine=string_clean(off_machine)
                yaml_file.write(apps+'_'+build[2].strip()+":"+ '\n')
                yaml_file.write("  build: "+ '\n')
                yaml_file.write("    compiler: "+compiler+ '\n')
                yaml_file.write("    option: "+options+ '\n')
                if not (off_machine is None):
                    yaml_file.write("    turnoff: ["+off_machine+"]"+ '\n')
                prev_line='COMPILE'
            if line.startswith("RUN"):  # RUN line
                build    = line.split('|')
                test     =     build[1].strip()
                machine  =     build[2].strip()
                baseline = "'"+build[3].strip()+"'"
                depend   =     build[4].strip()
                if (machine.find('-') != -1):
                    off_machine=machine.replace("-","").strip()
                    off_machine=string_clean(off_machine)
                tests = "    "+"- "+test+": {'project':['daily']"
                if baseline.isalnum(): tests = tests + ",'baseline': "+baseline
                if depend and depend.strip(): tests = tests + ",'dependency':'"+depend+"'"                
                if not (off_machine is None): tests = tests +",'turnoff':["+off_machine+"]"
                if prev_line == "COMPILE": yaml_file.write("  tests: "+ '\n')
                yaml_file.write(tests+"}"+ '\n')                
                prev_line='RUN'
        
    yaml_file.close(); conf_file.close()

