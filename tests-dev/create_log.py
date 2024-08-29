import os
import sys
import subprocess
import yaml
from datetime import datetime
from ufs_test_utils import get_testcase, write_logfile, delete_files, machine_check_off

def finish_log():
    """Collect regression test results and generate log file.
    """
    UFS_TEST_YAML = str(os.getenv('UFS_TEST_YAML'))
    PATHRT     = os.getenv('PATHRT')
    MACHINE_ID = os.getenv('MACHINE_ID')
    REGRESSIONTEST_LOG = PATHRT+'/logs/RegressionTests_'+MACHINE_ID+'.log'
    filename   = REGRESSIONTEST_LOG
    KEEP_RUNDIR= str(os.getenv('KEEP_RUNDIR'))
    ROCOTO     = str(os.getenv('ROCOTO'))
    CREATE_BASELINE = str(os.getenv('CREATE_BASELINE'))
    COMPILE_ONLY = str(os.getenv('COMPILE_ONLY'))

    run_logs= f"""
"""
    COMPILE_PASS= 0
    COMPILE_NR  = 0
    JOB_NR = 0
    PASS_NR= 0
    FAIL_NR= 0
    failed_list= []
    test_changes_list= PATHRT+'/test_changes.list'
    with open(UFS_TEST_YAML, 'r') as f:
        rt_yaml = yaml.load(f, Loader=yaml.FullLoader)
        for apps, jobs in rt_yaml.items():
            for key, val in jobs.items():
                if (str(key) == 'build'):
                    machine_check = machine_check_off(MACHINE_ID, val)
                    PASS_TESTS = False
                    if machine_check:
                        COMPILE_NR += 1
                        RT_COMPILER = val['compiler']
                        COMPILE_ID  = apps
                        COMPILE_LOG = 'compile_'+COMPILE_ID+'.log'
                        COMPILE_LOG_TIME ='compile_'+COMPILE_ID+'_timestamp.txt'
                        with open('./logs/log_'+MACHINE_ID+'/'+COMPILE_LOG) as f:
                            if "[100%] Linking Fortran executable" in f.read():
                                COMPILE_PASS += 1
                                f.seek(0)
                                for line in f:
                                    if 'export RUNDIR_ROOT=' in line:
                                        RUNDIR_ROOT=line.split("=")[1]
                                        break
                                compile_err = RUNDIR_ROOT.strip('\n')+'/compile_'+COMPILE_ID+'/err'
                                with open(compile_err) as ferr:
                                    contents = ferr.read()
                                    count_warning = contents.count(": warning #")
                                    count_remarks = contents.count(": remark #")
                                    ferr.close()
                                warning_log = ""
                                if count_warning > 0:
                                    warning_log = "("+str(count_warning)+" warnings"
                                if count_remarks > 0:
                                    warning_log+= ","+str(count_remarks)+" remarks)"
                                flog = open('./logs/log_'+MACHINE_ID+'/'+COMPILE_LOG_TIME)
                                timing_data = flog.read()
                                first_line = timing_data.split('\n', 1)[0]
                                etime = int(first_line.split(",")[4].strip()) - int(first_line.split(",")[1].strip())
                                btime = int(first_line.split(",")[3].strip()) - int(first_line.split(",")[2].strip())
                                etime_min, etime_sec = divmod(int(etime), 60)
                                etime_min = f"{etime_min:02}"; etime_sec = f"{etime_sec:02}"
                                btime_min, btime_sec = divmod(int(btime), 60)
                                btime_min = f"{btime_min:02}"; btime_sec = f"{btime_sec:02}"
                                time_log = " ["+etime_min+':'+etime_sec+', '+btime_min+':'+btime_sec+"]"
                                flog.close()
                                compile_log = "PASS -- COMPILE "+COMPILE_ID+time_log+warning_log+"\n"
                            else:
                                compile_log = "FAIL -- COMPILE "+COMPILE_ID+"\n"                        
                            f.close()
                        run_logs += compile_log
                    else:
                        PASS_TESTS = True
                if (str(key) == 'tests' and COMPILE_ONLY == 'false' and not PASS_TESTS):
                    for test in val:
                        case, config = get_testcase(test)
                        machine_check = machine_check_off(MACHINE_ID, config)
                        if machine_check:
                            JOB_NR+=1
                            TEST_NAME = case
                            TEST_ID   = TEST_NAME+'_'+RT_COMPILER
                            TEST_LOG  = 'rt_'+TEST_ID+'.log'
                            TEST_LOG_TIME= 'run_'+TEST_ID+'_timestamp.txt'
                            if 'dependency' in config.keys():
                                DEP_RUN = str(config['dependency'])+'_'+RT_COMPILER
                            else:
                                DEP_RUN = ""
                            PASS_CHECK = 'Test '+TEST_ID+' PASS'
                            MAXS_CHECK = 'The maximum resident set size (KB)'
                            pass_flag = False
                            create_dep_flag = False
                            if (CREATE_BASELINE == 'true' and not DEP_RUN == ""):
                                create_dep_flag = True
                            if not create_dep_flag:
                                with open('./logs/log_'+MACHINE_ID+'/'+TEST_LOG) as f:
                                    if PASS_CHECK in f.read():
                                        pass_flag = True
                                        f.close()
                                if pass_flag:
                                    f = open('./logs/log_'+MACHINE_ID+'/'+TEST_LOG_TIME)
                                    timing_data = f.read()
                                    first_line = timing_data.split('\n', 1)[0]
                                    etime = str(int(first_line.split(",")[4].strip()) - int(first_line.split(",")[1].strip()))
                                    rtime = str(int(first_line.split(",")[3].strip()) - int(first_line.split(",")[2].strip()))
                                    etime_min, etime_sec = divmod(int(etime), 60)
                                    etime_min = f"{etime_min:02}"; etime_sec = f"{etime_sec:02}"
                                    rtime_min, rtime_sec = divmod(int(rtime), 60)
                                    rtime_min = f"{rtime_min:02}"; rtime_sec = f"{rtime_sec:02}"
                                    time_log = " ["+etime_min+':'+etime_sec+', '+rtime_min+':'+rtime_sec+"]"
                                    f.close()
                                with open('./logs/log_'+MACHINE_ID+'/'+TEST_LOG) as f:
                                    if pass_flag :
                                        rtlog_file = f.readlines()
                                        for line in rtlog_file:
                                            if MAXS_CHECK in line:
                                                memsize= line.split('=')[1].strip()
                                        test_log = 'PASS -- TEST '+TEST_ID+time_log+' ('+memsize+' MB)\n'
                                        PASS_NR += 1
                                    else:
                                        test_log = 'FAIL -- TEST '+TEST_ID+'\n'
                                        failed_list.append(TEST_NAME+' '+RT_COMPILER)
                                        FAIL_NR += 1
                                    run_logs += test_log
                                f.close()
                    run_logs += '\n'
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

    if (int(FAIL_NR) == 0):
        if os.path.isfile(test_changes_list):
            delete_files(test_changes_list)
        open(test_changes_list, 'a').close()
        SUCCESS = "SUCCESS"
        comment_log = f"""
NOTES:
A file test_changes.list was generated but is empty.
If you are using this log as a pull request verification, please commit test_changes.list.

Result: {SUCCESS}

====END OF {MACHINE_ID} REGRESSION TESTING LOG====
"""
        write_logfile(filename, "a", output=comment_log)
    else:
        with open(test_changes_list, 'w') as listfile:
            for line in failed_list:
                listfile.write(f"{line}\n")
            listfile.close()
        SUCCESS = "FAILED"
        comment_log = f"""
NOTES:
A file test_changes.list was generated with list of all failed tests.
You can use './rt.sh -c -b test_changes.list' to create baselines for the failed tests.
If you are using this log as a pull request verification, please commit test_changes.list.

Result: FAILURE

====END OF {MACHINE_ID} REGRESSION TESTING LOG====
"""
        write_logfile(filename, "a", output=comment_log)
   
    print("Performing Cleanup...")
    exefiles= PATHRT+'/fv3_*.*x*'; delete_files(exefiles)
    modfiles= PATHRT+'/modules.fv3_*'; delete_files(modfiles)
    modfiles= PATHRT+'modulefiles/modules.fv3_*'; delete_files(modfiles)
    tmpfiles= PATHRT+'/keep_tests.tmp'; delete_files(tmpfiles)
    if KEEP_RUNDIR == 'false':
        rundir = PATHRT+'/run_dir'
        os.unlink(rundir)
    if ROCOTO == 'true':
        rocotofiles=PATHRT+'/rocoto*'
        delete_files(rocotofiles)
        lockfiles=PATHRT+'/*_lock.db'
        delete_files(lockfiles)
    print("REGRESSION TEST RESULT: SUCCESS")    

#if __name__ == '__main__':

