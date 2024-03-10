import yaml

def get_testcase(test):
    test_cases=[]
    for case, configs in test.items():
        case_name=case
        case_config=configs
        return case_name, case_config

if __name__ == "__main__":
    with open("rt.yaml", 'r') as f:
        rt_yaml = yaml.load(f, Loader=yaml.FullLoader)
        for apps, jobs in rt_yaml.items():
            print(apps)
            for key, val in jobs.items():
                if (str(key) == 'build'):
                    if not ('turnoff' in val.keys()): print('   ',val['compiler'],val['option'])
                    if 'turnoff' in val.keys(): print('   ',val['compiler'],val['option'],'turnoff: ',val['turnoff'])
                if (str(key) == 'tests'):
                    for test in val:
                        case, config = get_testcase(test)
                        print('     ',case, config)

