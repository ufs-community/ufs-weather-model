import yaml

def get_testcase(test):
    test_cases=[]
    for case, configs in test.items():
        case_name=case
        case_config=configs
        return case_name, case_config

if __name__ == "__main__":

    match_tests=[{'name':'cpld_control_gfsv17', 'compiler':'intel'}]
    
    with open("rt.yaml", 'r') as f:
        rt_yaml = yaml.load(f, Loader=yaml.FullLoader)

        match_apps={}
        for match_test in match_tests:
            match_case=match_test['name']
            match_compiler=match_test['compiler']
            for apps, jobs in rt_yaml.items():
                for key, val in jobs.items():
                    if (str(key) == 'build'):
                        build_val = val
                        build_compiler = val['compiler']
                    if (str(key) == 'tests'):
                        test_list=[]
                        for test in val:
                            case, config = get_testcase(test)
                            if (case==match_case and match_compiler==build_compiler):
                                test_list.append( {case:{'recurring':['pre-test']}} )
                if test_list:
                    match_apps.update({ apps: {'build': build_val, 'tests': test_list}})
                            
    with open("rt-test.yaml", 'w') as outfile:
        yaml.dump(match_apps, outfile, default_flow_style=False, sort_keys=False)
        #jkim
