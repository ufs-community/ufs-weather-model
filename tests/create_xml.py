import os
import subprocess
import yaml

def rocoto_create_entries(ROCOTO_XML):
    PATHRT=os.getenv('PATHRT')
    LOG_DIR=os.getenv('LOG_DIR')
    PATHTR=os.getenv('PATHTR')
    RTPWD=os.getenv('RTPWD')
    INPUTDATA_ROOT=os.getenv('INPUTDATA_ROOT')
    INPUTDATA_ROOT_WW3=os.getenv('INPUTDATA_ROOT_WW3')
    INPUTDATA_ROOT_BMIC=os.getenv('INPUTDATA_ROOT_BMIC')
    RUNDIR_ROOT=os.getenv('RUNDIR_ROOT')
    NEW_BASELINE=os.getenv('NEW_BASELINE')
    ROCOTO_SCHEDULER=os.getenv('ROCOTO_SCHEDULER')
    
    rocoto_entries = f"""<?xml version="1.0"?>
<!DOCTYPE workflow
[
  <!ENTITY PATHRT         "{PATHRT}">
  <!ENTITY LOG            "{LOG_DIR}">
  <!ENTITY PATHTR         "{PATHTR}">
  <!ENTITY RTPWD          "{RTPWD}">
  <!ENTITY INPUTDATA_ROOT "{INPUTDATA_ROOT}">
  <!ENTITY INPUTDATA_ROOT_WW3 "{INPUTDATA_ROOT_WW3}">
  <!ENTITY INPUTDATA_ROOT_BMIC "{INPUTDATA_ROOT_BMIC}">
  <!ENTITY RUNDIR_ROOT    "{RUNDIR_ROOT}">
  <!ENTITY NEW_BASELINE   "{NEW_BASELINE}">
]>
<workflow realtime="F" scheduler="{ROCOTO_SCHEDULER}" taskthrottle="20">
  <cycledef>197001010000 197001010000 01:00:00</cycledef>
  <log>&LOG;/workflow.log</log>    
"""
    with open(ROCOTO_XML,"w") as f:
        f.writelines(rocoto_entries)
    f.close()
    
def rocoto_create_compile_task(MACHINE_ID,COMPILE_ID,ROCOTO_COMPILE_MAXTRIES,MAKE_OPT,ACCNR,COMPILE_QUEUE,PARTITION,ROCOTO_XML):
    NATIVE=""
    BUILD_CORES="8"
    BUILD_WALLTIME="00:30:00"
    if ( MACHINE_ID == 'jet' ):
        BUILD_WALLTIME="02:00:00"
    if ( MACHINE_ID == 'hera'):
        BUILD_WALLTIME="01:00:00"
    if ( MACHINE_ID == 'orion'):
        BUILD_WALLTIME="01:00:00"
    if ( MACHINE_ID == 'hercules'):
        BUILD_WALLTIME="01:00:00"
    if ( MACHINE_ID == 's4' ):
        BUILD_WALLTIME="01:00:00"
    if ( MACHINE_ID == 'gaea' ):
        BUILD_WALLTIME="01:00:00"

    compile_task = f"""  <task name="compile_{COMPILE_ID}" maxtries="{ROCOTO_COMPILE_MAXTRIES}">
    <command>&PATHRT;/run_compile.sh &PATHRT; &RUNDIR_ROOT; "{MAKE_OPT}" {COMPILE_ID} 2>&amp;1 | tee &LOG;/compile_{COMPILE_ID}.log</\
command>
    <jobname>compile_{COMPILE_ID}</jobname>
    <account>{ACCNR}</account>
    <queue>{COMPILE_QUEUE}</queue>
"""
    
    if ( MACHINE_ID == 'gaea' ):
        compile_task+=f"""    <native>--clusters=es</native>
    <partition>eslogin_c5</partition>
"""
    if ( not PARTITION == "" or MACHINE_ID != "hera" ):
            compile_task+=f"""    <partition>{PARTITION}</partition>
"""
    compile_task+=f"""    <cores>{BUILD_CORES}</cores>
    <walltime>{BUILD_WALLTIME}</walltime>
    <join>&RUNDIR_ROOT;/compile_{COMPILE_ID}.log</join>
    {NATIVE}
  </task>
"""   
    with open(ROCOTO_XML,"a") as f:
        f.writelines(compile_task)
    f.close()

def write_metatask_begin(COMPILE_METATASK_NAME, filename):
    #filename='rocoto_workflow.xml'
    #COMPILE_METATASK_NAME=os.getenv('COMPILE_METATASK_NAME')
    metatask_name = f"""  <metatask name="compile_{COMPILE_METATASK_NAME}_tasks"><var name="zero">0</var>
"""
    with open(filename,"a") as f:
        f.writelines(metatask_name)
    f.close()

def write_metatask_end(filename):
    #filename='rocoto_workflow.xml'
    metatask_name = f"""  </metatask>
"""
    with open(filename,"a") as f:
        f.writelines(metatask_name)
    f.close()    
    
def write_compile_env():
    filename=str(os.getenv('RUNDIR_ROOT'))+"/compile_"+str(os.getenv('COMPILE_ID'))+".env"
    JOB_NR=os.getenv('JOB_NR')
    COMPILE_ID=os.getenv('COMPILE_ID')
    MACHINE_ID=os.getenv('MACHINE_ID')
    RT_COMPILER=os.getenv('RT_COMPILER')
    PATHRT=os.getenv('PATHRT')
    PATHTR=os.getenv('PATHTR')
    SCHEDULER=os.getenv('SCHEDULER')
    ACCNR=os.getenv('ACCNR')
    COMPILE_QUEUE=os.getenv('COMPILE_QUEUE')
    PARTITION=os.getenv('PARTITION')
    ROCOTO=os.getenv('ROCOTO')
    ECFLOW=os.getenv('ECFLOW')
    REGRESSIONTEST_LOG=os.getenv('REGRESSIONTEST_LOG')
    LOG_DIR=os.getenv('LOG_DIR')

    compile_envs = f"""export JOB_NR={JOB_NR}
export COMPILE_ID={COMPILE_ID}
export MACHINE_ID={MACHINE_ID}
export RT_COMPILER={RT_COMPILER}
export PATHRT={PATHRT}
export PATHTR={PATHTR}
export SCHEDULER={SCHEDULER}
export ACCNR={ACCNR}
export QUEUE={COMPILE_QUEUE}
export PARTITION={PARTITION}
export ROCOTO={ROCOTO}
export ECFLOW={ECFLOW}
export REGRESSIONTEST_LOG={REGRESSIONTEST_LOG}
export LOG_DIR={LOG_DIR}
"""

    with open(filename,"w+") as f:
        f.writelines(compile_envs)
    f.close()

def write_runtest_env():
    filename=str(os.getenv('RUNDIR_ROOT'))+"/run_test_"+str(os.getenv('TEST_ID'))+".env"
    
    JOB_NR=str(os.getenv('JOB_NR'))
    TEST_ID=str(os.getenv('TEST_ID'))
    MACHINE_ID=str(os.getenv('MACHINE_ID'))
    RT_COMPILER=str(os.getenv('RT_COMPILER'))
    RTPWD=str(os.getenv('RTPWD'))
    INPUTDATA_ROOT=str(os.getenv('INPUTDATA_ROOT'))
    INPUTDATA_ROOT_WW3=str(os.getenv('INPUTDATA_ROOT_WW3'))
    INPUTDATA_ROOT_BMIC=str(os.getenv('INPUTDATA_ROOT_BMIC'))
    PATHRT=str(os.getenv('PATHRT'))
    PATHTR=str(os.getenv('PATHTR'))
    NEW_BASELINE=str(os.getenv('NEW_BASELINE'))
    CREATE_BASELINE=str(os.getenv('CREATE_BASELINE'))
    RT_SUFFIX=str(os.getenv('RT_SUFFIX'))
    BL_SUFFIX=str(os.getenv('BL_SUFFIX'))
    SCHEDULER=str(os.getenv('SCHEDULER'))
    ACCNR=str(os.getenv('ACCNR'))
    QUEUE=str(os.getenv('QUEUE'))
    PARTITION=str(os.getenv('PARTITION'))
    ROCOTO=str(os.getenv('ROCOTO'))
    ECFLOW=str(os.getenv('ECFLOW'))
    REGRESSIONTEST_LOG=str(os.getenv('REGRESSIONTEST_LOG'))
    LOG_DIR=str(os.getenv('LOG_DIR'))
    DEP_RUN=str(os.getenv('DEP_RUN'))
    skip_check_results=str(os.getenv('skip_check_results'))
    delete_rundir=str(os.getenv('delete_rundir'))
    WLCLK=str(os.getenv('WLCLK'))
    MACHINE_ID=str(os.getenv('MACHINE_ID'))

    runtest_envs = f"""
export JOB_NR={JOB_NR}
export TEST_ID={TEST_ID}
export MACHINE_ID={MACHINE_ID}
export RT_COMPILER={RT_COMPILER}
export RTPWD={RTPWD}
export INPUTDATA_ROOT={INPUTDATA_ROOT}
export INPUTDATA_ROOT_WW3={INPUTDATA_ROOT_WW3}
export INPUTDATA_ROOT_BMIC={INPUTDATA_ROOT_BMIC}
export PATHRT={PATHRT}
export PATHTR={PATHTR}
export NEW_BASELINE={NEW_BASELINE}
export CREATE_BASELINE={CREATE_BASELINE}
export RT_SUFFIX={RT_SUFFIX}
export BL_SUFFIX={BL_SUFFIX}
export SCHEDULER={SCHEDULER}
export ACCNR={ACCNR}
export QUEUE={QUEUE}
export PARTITION={PARTITION}
export ROCOTO={ROCOTO}
export ECFLOW={ECFLOW}
export REGRESSIONTEST_LOG={REGRESSIONTEST_LOG}
export LOG_DIR={LOG_DIR}
export DEP_RUN={DEP_RUN}
export skip_check_results={skip_check_results}
export delete_rundir={delete_rundir}
export WLCLK={WLCLK}
"""
    if ( MACHINE_ID == 'jet' ):
         runtest_env+="export PATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/bin:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/bin:$PATH"
         runtest_env+="export PYTHONPATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/lib/python3.8/site-packages:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/lib/python3.8/site-packages"
     
    with open(filename,"w+") as f:
        f.writelines(runtest_envs)
    f.close()     
    
def compile_task():
    subprocess.Popen(['bash', '-c', '. ufs_test_utils.sh; set_compile_task'])

def run_task():
    subprocess.Popen(['bash', '-c', '. ufs_test_utils.sh; set_run_task'])

def get_testcase(test):
    test_cases=[]
    for case, configs in test.items():
        case_name=case
        case_config=configs
        return case_name, case_config

def main_loop():
    JOB_NR = 0
    in_metatask = False
    new_compile = False
    ROCOTO = True
    ROCOTO_XML=os.getenv('ROCOTO_XML')
    rocoto_create_entries(ROCOTO_XML)
    with open("rt.yaml", 'r') as f:
        rt_yaml = yaml.load(f)#, Loader=yaml.FullLoader)
        for apps, jobs in rt_yaml.items():
            print(apps)
            for key, val in jobs.items():
                if (str(key) == 'build'):
                    new_compile = True
                    if not ('turnoff' in val.keys()): print('   ',val['compiler'],val['option'])
                    if 'turnoff' in val.keys(): print('   ',val['compiler'],val['option'],'turnoff: ',val['turnoff'])
                    RT_COMPILER=val['compiler']
                    COMPILE_ID=apps+'_'+RT_COMPILER
                    MAKE_OPT=val['option']
                    os.environ["in_metatask"] = str(in_metatask)
                    os.environ["COMPILE_ID"] = str(COMPILE_ID)
                    os.environ["MAKE_OPT"] = str(MAKE_OPT)
                    MACHINE_ID=os.getenv('MACHINE_ID')
                    #ROCOTO_COMPILE_MAXTRIES=os.getenv('ROCOTO_COMPILE_MAXTRIES')
                    #if (ROCOTO_COMPILE_MAXTRIES == 'None'):
                    ROCOTO_COMPILE_MAXTRIES="3"
                    ACCNR=os.getenv('ACCNR')
                    COMPILE_QUEUE=os.getenv('COMPILE_QUEUE')
                    PARTITION=os.getenv('PARTITION')
                    write_compile_env()
                    rocoto_create_compile_task \
                        (MACHINE_ID,COMPILE_ID,ROCOTO_COMPILE_MAXTRIES,MAKE_OPT,ACCNR,COMPILE_QUEUE,PARTITION,ROCOTO_XML)
                if (str(key) == 'tests'):
                    JOB_NR+=1
                    if ( ROCOTO ):
                        write_metatask_begin(COMPILE_ID, ROCOTO_XML)
                        for test in val:
                            case, config = get_testcase(test)
                            TEST_NAME = case
                            TEST_ID   = TEST_NAME+'_'+RT_COMPILER
                            if 'dependency' in config.keys():
                                DEP_RUN = str(config['dependency'])+'_'+RT_COMPILER
                            else:
                                DEP_RUN = ""
                            RT_SUFFIX = ""
                            BL_SUFFIX = ""
                            ROCOTO_TEST_MAXTRIES="3"
                            os.environ["TEST_NAME"] = TEST_NAME
                            os.environ["DEP_RUN"] = DEP_RUN
                            os.environ["TEST_ID"] = TEST_ID
                            os.environ["RT_SUFFIX"] = RT_SUFFIX
                            os.environ["BL_SUFFIX"] = BL_SUFFIX
                            os.environ["RT_COMPILER"] = RT_COMPILER
                            os.environ["MACHINE_ID"] = MACHINE_ID
                            os.environ["ROCOTO_TEST_MAXTRIES"] = ROCOTO_TEST_MAXTRIES
                            print('     ',case, config)                            
                            rc_set_run_task =subprocess.Popen(['bash', '-c', '. ufs_test_utils.sh; set_run_task'])
                            rc_set_run_task.wait()   
                        write_metatask_end(ROCOTO_XML)
        print(rt_yaml)
    rocoto_close=f"""</workflow>
"""
    with open(ROCOTO_XML,"a") as f:
        f.writelines(rocoto_close)
    f.close()

#if __name__ == "__main__":

