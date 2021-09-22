#!/usr/bin/env python3
import json


def main():
    with open('ci.test', 'r') as setup_file:
        input_str = setup_file.read().splitlines()

    tests = []
    cases = []
    for i, e in enumerate(input_str):
        if i % 2 == 0:
            tests.append(e)
        else:
            cases.append([tests[(i-1)//2]+'_'+case for case in e.split()])

    bj = {'bld_set': [], 'include': []}
    tj = {'test_set': [], 'include': []}

    for i in range(len(tests)):
        for j in range(len(cases[i])):
            test = tests[i]
            case = cases[i][j]

            if any(e in case for e in ['thr', 'mpi', 'dcp', 'rst']):
                std = test+'_std'
                if not bj['bld_set'] or not any(std == e for e in bj['bld_set']):
                    bj['bld_set'].append(std)
                    bj['include'].append(
                        {'bld_set': std, 'name': test, 'case': std[-3:]})
                aj = {'test_set': case, 'name': test,
                      'case': case[-3:], 'artifact': std}
            else:
                bj['bld_set'].append(case)
                bj['include'].append(
                    {'bld_set': case, 'name': test, 'case': case[-3:]})
                aj = {'test_set': case, 'name': test,
                      'case': case[-3:], 'artifact': case}

            tj['test_set'].append(case)
            tj['include'].append(aj)

    print(json.dumps(bj), "|", json.dumps(tj))


if __name__ == '__main__':
    main()
