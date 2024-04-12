import os
import subprocess

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
    
#if __name__ == "__main__":

