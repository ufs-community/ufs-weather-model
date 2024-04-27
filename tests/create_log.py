import os
import sys
import subprocess
import yaml
from datetime import datetime

def get_testcase(test):
    test_cases=[]
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
            
def finish_log():
    PATHRT     = os.getenv('PATHRT')
    MACHINE_ID = os.getenv('MACHINE_ID')
    REGRESSIONTEST_LOG = PATHRT+'/logs/RegressionTests_'+MACHINE_ID+'.log'
    filename   = REGRESSIONTEST_LOG

    run_logs= f"""
"""
    COMPILE_PASS = 0
    COMPILE_NR   = 0
    JOB_NR  = 0
    PASS_NR = 0
    FAIL_NR = 0
    with open("ufs_test.yaml", 'r') as f:
        rt_yaml = yaml.load(f)#, Loader=yaml.FullLoader)
        for apps, jobs in rt_yaml.items():
            for key, val in jobs.items():
                if (str(key) == 'build'):
                    COMPILE_NR += 1
                    if not ('turnoff' in val.keys()): print('   ',val['compiler'],val['option'])
                    if 'turnoff' in val.keys(): print('   ',val['compiler'],val['option'],'turnoff: ',val['turnoff'])
                    RT_COMPILER = val['compiler']
                    COMPILE_ID  = apps+'_'+RT_COMPILER
                    COMPILE_LOG = 'compile_'+COMPILE_ID+'.log'
                    COMPILE_START_CHECK = "+ date_strt='"
                    COMPILE_END_CHECK   = "+ date_end='"
                    COMPILE_ELAPSE_CHECK= "+ echo 'Elapsed time"
                    print('./logs/log_hera/'+COMPILE_LOG)
                    with open('./logs/log_hera/'+COMPILE_LOG) as f:
                        compilelog_file = f.readlines()
                        for line in compilelog_file:
                            if COMPILE_START_CHECK in line:
                                sdate = line.split(' ')[4].strip()
                                hh=sdate.split(":")[0]; mm=sdate.split(":")[1]; ss=sdate.split(":")[2]
                                sdate_sec=int(hh) * 3600 + int(mm) * 60 + float(ss)
                            if COMPILE_END_CHECK in line:
                                edate = line.split(' ')[4].strip()
                                hh=edate.split(":")[0]; mm=edate.split(":")[1]; ss=edate.split(":")[2]
                                edate_sec=int(hh) * 3600 + int(mm) * 60 + float(ss)
                            if COMPILE_ELAPSE_CHECK in line:
                                ctime = line.split(' ')[4].strip()
                                COMPILE_PASS += 1
                    f.close()
                    compile_log="PASS -- COMPILE "+COMPILE_ID+" ["+str(ctime)+"]\n"
                    run_logs += compile_log
                if (str(key) == 'tests'):
                    for test in val:
                        JOB_NR+=1
                        case, config = get_testcase(test)
                        TEST_NAME = case
                        TEST_ID   = TEST_NAME+'_'+RT_COMPILER
                        TEST_LOG  = 'rt_'+TEST_ID+'.log'
                        if 'dependency' in config.keys():
                            DEP_RUN = str(config['dependency'])+'_'+RT_COMPILER
                        else:
                            DEP_RUN = ""
                        print(TEST_LOG)
                        PASS_CHECK = 'Test '+TEST_ID+' PASS'
                        TIME_CHECK = 'The total amount of wall time'
                        MAXS_CHECK = 'The maximum resident set size (KB)'
                        pass_flag = False
                        with open('./logs/log_hera/'+TEST_LOG) as f:
                            if PASS_CHECK in f.read():
                                pass_flag = True
                        f.close()
                        with open('./logs/log_hera/'+TEST_LOG) as f:
                            if pass_flag :
                                rtlog_file = f.readlines()
                                for line in rtlog_file:
                                    if TIME_CHECK in line:
                                        timing = line.split('=')[1].strip()
                                        print(timing)
                                    if MAXS_CHECK in line:
                                        memsize= line.split('=')[1].strip()
                                        print(memsize)
                                test_log = 'PASS -- TEST '+TEST_ID+' ['+timing+']('+memsize+' MB)\n'
                                PASS_NR += 1
                            else:
                                test_log = 'FAIL -- TEST '+TEST_ID+'\n'
                                FAIL_NR += 1
                            run_logs += test_log
                        f.close()
    write_logfile(filename, "a", output=run_logs)

    TEST_START_TIME = os.getenv('TEST_START_TIME')
    TEST_END_TIME   = os.getenv('TEST_END_TIME')
    start_time      = datetime.strptime(TEST_START_TIME, "%Y%m%d %H:%M:%S")
    end_time        = datetime.strptime(TEST_END_TIME, "%Y%m%d %H:%M:%S")
    hours, remainder= divmod((end_time - start_time).total_seconds(), 3600)
    minutes, seconds= divmod(remainder, 60)
    hours = int(hours);    minutes=int(minutes);     seconds =int(seconds)
    hours = f"{hours:02}"; minutes= f"{minutes:02}"; seconds= f"{seconds:02}"
    elapsed_time = hours+'h:'+minutes+'m:'+seconds+'s'

    COMPILE_PASS = str(int(COMPILE_PASS))
    COMPILE_NR   = str(int(COMPILE_NR))
    JOB_NR       = str(int(JOB_NR))
    PASS_NR      = str(int(PASS_NR))
    FAIL_NR      = str(int(FAIL_NR))
    synop_log    = f"""
SYNOPSIS:
Starting Date/Time: {TEST_START_TIME}
Ending Date/Time: {TEST_END_TIME}
Total Time: {elapsed_time}
Compiles Completed: {COMPILE_PASS}/{COMPILE_NR}
Tests Completed: {PASS_NR}/{JOB_NR}

"""    
    write_logfile(filename, "a", output=synop_log)

    if (JOB_NR == PASS_NR):
        SUCCESS = "SUCCESS"
    else:
        SUCCESS = "FAILED"
        
    comment_log = f"""NOTES:
A file 'test_changes.list' was generated but is empty.
If you are using this log as a pull request verification, please commit 'test_changes.list'.

Result: {SUCCESS}

====END OF {MACHINE_ID} REGRESSION TESTING LOG====
"""
    write_logfile(filename, "a", output=comment_log)
    
if __name__ == '__main__':
    finish_log()
