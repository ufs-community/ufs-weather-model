import os
import sys
import subprocess
import yaml
from ufs_test_utils import get_testcase, write_logfile, rrmdir, machine_check_off

def rocoto_create_entries(RTPWD,MACHINE_ID,INPUTDATA_ROOT,INPUTDATA_ROOT_WW3,INPUTDATA_ROOT_BMIC,RUNDIR_ROOT,NEW_BASELINE,ROCOTO_XML):
    """Generate header information for Rocoto xml file

    Args:
        RTPWD (str): Baseline directory
        MACHINE_ID (str): Machine ID i.e. Hera, Gaea, Jet, etc.
        INPUTDATA_ROOT (str): Input data directory
        INPUTDATA_ROOT_WW3 (str): WW3 input data directory
        INPUTDATA_ROOT_BMIC (str): BMIC input data directory
        RUNDIR_ROOT (str): Test run directory
        NEW_BASELINE (str): Directory for newly generated baselines
        ROCOTO_XML (str): Rocoto .xml filename to write to
    """
    PATHRT = os.getenv('PATHRT')
    LOG_DIR= PATHRT+'/logs/log_'+MACHINE_ID
    PATHTR, tail = os.path.split(PATHRT)
    ROCOTO_SCHEDULER = os.getenv('ROCOTO_SCHEDULER')
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
<workflow realtime="F" scheduler="{ROCOTO_SCHEDULER}" taskthrottle="10">
  <cycledef>197001010000 197001010000 01:00:00</cycledef>
  <log>&LOG;/workflow.log</log>    
"""
    with open(ROCOTO_XML,"w") as f:
        f.writelines(rocoto_entries)
    f.close()
    
def rocoto_create_compile_task(MACHINE_ID,COMPILE_ID,ROCOTO_COMPILE_MAXTRIES,MAKE_OPT,ACCNR,COMPILE_QUEUE,PARTITION,ROCOTO_XML):
    """Generate and append compile task into Rocoto xml file

    Args:
        MACHINE_ID (str): Machine ID i.e. Hera, Gaea, Jet, etc.
        COMPILE_ID (str): Compile identifier e.g. s2swa_intel
        ROCOTO_COMPILE_MAXTRIES (str): Max attempts for compile
        MAKE_OPT (str): Make build options
        ACCNR (str): Account to run the job with
        COMPILE_QUEUE (str): QOS i.e. batch, windfall, normal, etc.
        PARTITION (str): System partition i.e. xjet, c5
        ROCOTO_XML (str): Rocoto .xml filename to write to
    """
    NATIVE=""
    BUILD_CORES="8"
    BUILD_WALLTIME="00:30:00"
    if ( MACHINE_ID == 'jet' ):  BUILD_WALLTIME="02:00:00"
    if ( MACHINE_ID == 'hera'):  BUILD_WALLTIME="01:00:00"
    if ( MACHINE_ID == 'orion'): BUILD_WALLTIME="01:00:00"
    if ( MACHINE_ID == 'hercules'): BUILD_WALLTIME="01:00:00"
    if ( MACHINE_ID == 's4' ):   BUILD_WALLTIME="01:00:00"
    if ( MACHINE_ID == 'gaea' ): BUILD_WALLTIME="01:00:00"
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
        PARTITION= ""
    if ( PARTITION != "" and MACHINE_ID != "hera" ):
            compile_task+=f"""    <partition>{PARTITION}</partition>
"""
    compile_task+=f"""    <nodes>1:ppn={BUILD_CORES}</nodes>
    <walltime>{BUILD_WALLTIME}</walltime>
    <join>&RUNDIR_ROOT;/compile_{COMPILE_ID}.log</join>
    {NATIVE}
  </task>
"""   
    with open(ROCOTO_XML,"a") as f:
        f.writelines(compile_task)
    f.close()

def write_metatask_begin(COMPILE_METATASK_NAME, filename):
    """Write compile task metadata to Rocoto xml file

    Args:
        COMPILE_METATASK_NAME (str): Compile job name e.g. s2swa_intel
        filename (str): Rocoto xml filename to append to
    """
    metatask_name = f"""  <metatask name="compile_{COMPILE_METATASK_NAME}_tasks"><var name="zero">0</var>
"""
    with open(filename,"a") as f:
        f.writelines(metatask_name)
    f.close()

def write_metatask_end(filename):
    """Append closing metatask element to Rocoto xml

    Args:
        filename (str): Rocoto xml filename
    """
    metatask_name = f"""  </metatask>
"""
    with open(filename,"a") as f:
        f.writelines(metatask_name)
    f.close()    
    
def write_compile_env(SCHEDULER,PARTITION,JOB_NR,COMPILE_QUEUE,RUNDIR_ROOT):
    """Generate compile task .env file

    Args:
        SCHEDULER (str): Job scheduler e.g. pbs, slurm
        PARTITION (str): System partition i.e. xjet, c5
        JOB_NR (str): Job number
        COMPILE_QUEUE (str): QOS i.e. batch, windfall, normal, etc.
        RUNDIR_ROOT (str): Test run directory
    """
    filename   = RUNDIR_ROOT+"/compile_"+str(os.getenv('COMPILE_ID'))+".env"
    COMPILE_ID = os.getenv('COMPILE_ID')
    MACHINE_ID = os.getenv('MACHINE_ID')
    RT_COMPILER= os.getenv('RT_COMPILER')
    PATHRT     = os.getenv('PATHRT')
    PATHTR, tail = os.path.split(PATHRT)
    ACCNR      = os.getenv('ACCNR')
    ROCOTO     = os.getenv('ROCOTO')
    ECFLOW     = os.getenv('ECFLOW')
    REGRESSIONTEST_LOG = PATHRT+'/logs/RegressionTests_'+MACHINE_ID+'.log'
    LOG_DIR    = PATHRT+'/logs/log_'+MACHINE_ID
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
    """Generate run task .env file
    """
    filename   = str(os.getenv('RUNDIR_ROOT'))+"/run_test_"+str(os.getenv('TEST_ID'))+".env"
    JOB_NR     = str(os.getenv('JOB_NR'))
    TEST_ID    = str(os.getenv('TEST_ID'))
    MACHINE_ID = str(os.getenv('MACHINE_ID'))
    RT_COMPILER= str(os.getenv('RT_COMPILER'))
    RTPWD      = str(os.getenv('RTPWD'))
    INPUTDATA_ROOT     = str(os.getenv('INPUTDATA_ROOT'))
    INPUTDATA_ROOT_WW3 = str(os.getenv('INPUTDATA_ROOT_WW3'))
    INPUTDATA_ROOT_BMIC= str(os.getenv('INPUTDATA_ROOT_BMIC'))
    PATHRT = str(os.getenv('PATHRT'))
    PATHTR, tail    = os.path.split(PATHRT)
    NEW_BASELINE    = str(os.getenv('NEW_BASELINE'))
    CREATE_BASELINE =str(os.getenv('CREATE_BASELINE'))
    RT_SUFFIX = str(os.getenv('RT_SUFFIX'))
    BL_SUFFIX = str(os.getenv('BL_SUFFIX'))
    SCHEDULER = str(os.getenv('SCHEDULER'))
    ACCNR = str(os.getenv('ACCNR'))
    QUEUE = str(os.getenv('QUEUE'))
    PARTITION = str(os.getenv('PARTITION'))
    ROCOTO    = str(os.getenv('ROCOTO'))
    ECFLOW    = str(os.getenv('ECFLOW'))
    REGRESSIONTEST_LOG = PATHRT+'/logs/RegressionTests_'+MACHINE_ID+'.log'
    LOG_DIR = PATHRT+'/logs/log_'+MACHINE_ID
    DEP_RUN = str(os.getenv('DEP_RUN'))
    skip_check_results = str(os.getenv('skip_check_results'))
    delete_rundir = str(os.getenv('delete_rundir'))
    WLCLK         = str(os.getenv('WLCLK'))
    MACHINE_ID    = str(os.getenv('MACHINE_ID'))
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
export RTVERBOSE=false
"""
    if ( MACHINE_ID == 'jet' ):
         runtest_envs+="export PATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/bin:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/bin:$PATH"
         runtest_envs+="export PYTHONPATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/lib/python3.8/site-packages:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/lib/python3.8/site-packages"
     
    with open(filename,"w+") as f:
        f.writelines(runtest_envs)
    f.close()     

def make_loghead(ACCNR,MACHINE_ID,RUNDIR_ROOT,RTPWD,REGRESSIONTEST_LOG):
    """Generate log header information

    Args:
        ACCNR (str): Account to run the job with
        MACHINE_ID (str): Machine ID i.e. Hera, Gaea, Jet, etc.
        RUNDIR_ROOT (str): Test run directory
        RTPWD (str): Baseline directory
        REGRESSIONTEST_LOG (str): Regression Test log filename
    """
    filename   = REGRESSIONTEST_LOG
    TESTS_FILE = str(os.getenv('TESTS_FILE'))
    NEW_BASELINES_FILE = str(os.getenv('NEW_BASELINES_FILE'))
    CREATE_BASELINE    = str(os.getenv('CREATE_BASELINE'))
    DEFINE_CONF_FILE   = str(os.getenv('DEFINE_CONF_FILE'))
    RTPWD_NEW_BASELINE = str(os.getenv('RTPWD_NEW_BASELINE'))
    RUN_SINGLE_TEST = str(os.getenv('RUN_SINGLE_TEST'))
    COMPILE_ONLY    = str(os.getenv('COMPILE_ONLY'))
    delete_rundir   = str(os.getenv('delete_rundir'))
    skip_check_results = str(os.getenv('skip_check_results'))
    KEEP_RUNDIR = str(os.getenv('KEEP_RUNDIR'))
    ROCOTO      = str(os.getenv('ROCOTO'))
    ECFLOW      = str(os.getenv('ECFLOW'))
    RTVERBOSE   = str(os.getenv('RTVERBOSE'))
    SRT_NAME    = str(os.getenv('SRT_NAME'))
    SRT_COMPILER= str(os.getenv('SRT_COMPILER'))
    
    rtlog_head=f"""====START OF {MACHINE_ID} REGRESSION TESTING LOG====

UFSWM hash used in testing:
"""
    rtlog_submod=f"""
Submodule hashes used in testing:
"""
    write_logfile(filename, "w", output= rtlog_head)
    write_logfile(filename, "a", subproc="git rev-parse HEAD")
    write_logfile(filename, "a", output= rtlog_submod)
    write_logfile(filename, "a", subproc="git submodule status --recursive")

    with open(filename, 'r') as rtlog:
        filedata = rtlog.read(); rtlog.close()

    filedata = filedata.replace('../', '')
    write_logfile(filename, "w", output= filedata)

    info_note=f"""
NOTES:
[Times](Memory) are at the end of each compile/test in format [MM:SS](Size).
The first time is for the full script (prep+run+finalize).
The second time is specifically for the run phase.
Times/Memory will be empty for failed tests.

BASELINE DIRECTORY: {RTPWD}
COMPARISON DIRECTORY: {RUNDIR_ROOT}

UFS_TEST.SH OPTIONS USED:
"""
    write_logfile(filename, "a", output= info_note)

    write_logfile(filename, "a", output="* (-a) - HPC PROJECT ACCOUNT: "+ACCNR+"\n")
    if (not NEW_BASELINES_FILE == ""):
        write_logfile(filename, "a", output="* (-b) - NEW BASELINES FROM FILE: "+NEW_BASELINES_FILE+"\n")
    if (CREATE_BASELINE == "true"):
        write_logfile(filename, "a", output="* (-c) - CREATE NEW BASELINES"+"\n")
    if (DEFINE_CONF_FILE == "true"):
        write_logfile(filename, "a", output="* (-l) - USE CONFIG FILE: "+TESTS_FILE+"\n")
    if (RTPWD_NEW_BASELINE == "true"):
        write_logfile(filename, "a", output="* (-m) - COMPARE AGAINST CREATED BASELINES"+"\n")
    if (RUN_SINGLE_TEST == "true"):
        write_logfile(filename, "a", output="* (-n) - RUN SINGLE TEST: "+SRT_NAME+" "+SRT_COMPILER+"\n")
    if (COMPILE_ONLY == "true"):
        write_logfile(filename, "a", output="* (-o) - COMPILE ONLY, SKIP TESTS"+"\n")
    if (delete_rundir == "true"):
        write_logfile(filename, "a", output="* (-d) - DELETE RUN DIRECTORY"+"\n")
    if (skip_check_results == "true"):
        write_logfile(filename, "a", output="* (-w) - SKIP RESULTS CHECK"+"\n")
    if (KEEP_RUNDIR == "true"):
        write_logfile(filename, "a", output="* (-k) - KEEP RUN DIRECTORY"+"\n")
    if (ROCOTO == "true"):
        write_logfile(filename, "a", output="* (-r) - USE ROCOTO"+"\n")
    if (ECFLOW == "true"):
        write_logfile(filename, "a", output="* (-e) - USE ECFLOW"+"\n")
    if (RTVERBOSE == "true"):
        write_logfile(filename, "a", output="* (-v) - VERBOSE OUTPUT"+"\n")

def xml_loop():
    ACCNR      = str(os.getenv('ACCNR'))
    PATHRT     = str(os.getenv('PATHRT'))
    MACHINE_ID = str(os.getenv('MACHINE_ID'))
    RTPWD_NEW_BASELINE = str(os.getenv('RTPWD_NEW_BASELINE'))
    NEW_BASELINE       = str(os.getenv('NEW_BASELINE'))
    CREATE_BASELINE    = str(os.getenv('CREATE_BASELINE'))
    COMPILE_ONLY       = str(os.getenv('COMPILE_ONLY'))
    delete_rundir = str(os.getenv('delete_rundir'))
    if (delete_rundir == "true"): dependency_list= []
        
    with open('bl_date.conf', 'r') as bldate:
        bl_date = str(bldate.readline())
    BL_DATE = bl_date.split("=")[1].strip()
    with open("baseline_setup.yaml", 'r') as f:
        exp_config = yaml.load(f, Loader=yaml.FullLoader)
        base = exp_config[MACHINE_ID]
        USER = str(os.environ.get('USER')) #os.environ.get('USERNAME')) #os.getlogin()
        pid  = str(os.getpid())

        QUEUE         = str(base['QUEUE'])
        COMPILE_QUEUE = str(base['COMPILE_QUEUE'])
        PARTITION     = str(base['PARTITION'])
        if (PARTITION == "None"): PARTITION = ""
        dprefix       = str(base['dprefix']).replace("${USER}", str(USER))
        dprefix       = dprefix.replace("${ACCNR}", str(ACCNR))
        DISKNM        = str(base['DISKNM'])
        STMP          = str(base['STMP']).replace("${USER}", str(USER))
        STMP          = STMP.replace("${ACCNR}", str(ACCNR))
        PTMP          = str(base['PTMP']).replace("${USER}", str(USER))
        PTMP          = PTMP.replace("${ACCNR}", str(ACCNR))
        RUNDIR_ROOT   = str(base['RUNDIR_ROOT'])
        SCHEDULER     = str(base['SCHEDULER'])
        INPUTDATA_ROOT= str(base['INPUTDATA_ROOT'])
        INPUTDATA_ROOT_WW3 = str(base['INPUTDATA_ROOT_WW3'])
        INPUTDATA_ROOT_BMIC= str(base['INPUTDATA_ROOT_BMIC'])
        
    path = STMP+'/'+USER
    os.makedirs(path, exist_ok=True)
    NEW_BASELINE = path + '/FV3_RT/REGRESSION_TEST'
    if (RUNDIR_ROOT == "None"): RUNDIR_ROOT=PTMP+'/'+USER+'/FV3_RT/rt_'+pid  
    os.makedirs(RUNDIR_ROOT, exist_ok=True)
    if (os.path.islink(PATHRT+'/run_dir')): os.unlink(PATHRT+'/run_dir')
    if (os.path.isfile(PATHRT+'/run_dir')): os.remove(PATHRT+'/run_dir')
    if (os.path.isdir(PATHRT+'/run_dir')):  rrmdir(PATHRT+'/run_dir')
    print('Linking ',RUNDIR_ROOT,' to ',PATHRT,'/run_dir')
    os.symlink(RUNDIR_ROOT,PATHRT+'/run_dir')
    print('Run regression test in: ',RUNDIR_ROOT)
    LOG_DIR = PATHRT+'/logs/log_'+MACHINE_ID
    os.makedirs(LOG_DIR, exist_ok=True)
    
    if ( RTPWD_NEW_BASELINE == 'true' ):
        RTPWD = NEW_BASELINE
    else:
        RTPWD = DISKNM+'/NEMSfv3gfs/develop-'+BL_DATE

    if (CREATE_BASELINE == 'false'):
        if ( not os.path.isdir(RTPWD) ) :
            print("Baseline directory does not exist:",RTPWD)
            sys.exit("***Baseline directory trouble***")
    elif (len(os.listdir(RTPWD)) == 0):
        print("Baseline directory is empty:",RTPWD)
        sys.exit("***Baseline directory trouble***")
    else:
        if ( not os.path.isdir(NEW_BASELINE) ) :
            os.makedirs(NEW_BASELINE, exist_ok=True)
        else:
            rrmdir(NEW_BASELINE)
            os.makedirs(NEW_BASELINE, exist_ok=True)
        
    ROCOTO_TEST_MAXTRIES = "3"
    RTVERBOSE = False
    os.environ["MACHINE_ID"]  = MACHINE_ID
    os.environ["ROCOTO_TEST_MAXTRIES"] = ROCOTO_TEST_MAXTRIES
    os.environ["NEW_BASELINE"] = NEW_BASELINE
    os.environ["RUNDIR_ROOT"]  = RUNDIR_ROOT
    os.environ["QUEUE"]        = QUEUE
    os.environ["INPUTDATA_ROOT"]     = INPUTDATA_ROOT
    os.environ["INPUTDATA_ROOT_WW3"] = INPUTDATA_ROOT_WW3
    os.environ["INPUTDATA_ROOT_BMIC"]= INPUTDATA_ROOT_BMIC
    os.environ["PARTITION"] = PARTITION
    os.environ["SCHEDULER"] = SCHEDULER
    os.environ["RTPWD"]     = RTPWD
    os.environ["RTVERBOSE"] = str(RTVERBOSE)

    JOB_NR = 0
    ROCOTO = True
    ROCOTO_XML = os.getenv('ROCOTO_XML')
    rocoto_create_entries(RTPWD,MACHINE_ID,INPUTDATA_ROOT,INPUTDATA_ROOT_WW3,INPUTDATA_ROOT_BMIC,RUNDIR_ROOT,NEW_BASELINE,ROCOTO_XML)
    UFS_TEST_YAML = str(os.getenv('UFS_TEST_YAML'))
    with open(UFS_TEST_YAML, 'r') as f:
        rt_yaml = yaml.load(f, Loader=yaml.FullLoader)
        for apps, jobs in rt_yaml.items():
            for key, val in jobs.items():
                if (str(key) == 'build'):
                    machine_check = machine_check_off(MACHINE_ID, val)
                    PASS_TESTS = False
                    if machine_check:
                        RT_COMPILER = val['compiler']
                        COMPILE_ID  = apps
                        MAKE_OPT    = val['option']
                        os.environ["COMPILE_ID"]  = str(COMPILE_ID)
                        os.environ["MAKE_OPT"]    = str(MAKE_OPT)
                        ROCOTO_COMPILE_MAXTRIES = "3"
                        os.environ["RT_COMPILER"] = str(RT_COMPILER)
                        write_compile_env(SCHEDULER,PARTITION,str(JOB_NR),COMPILE_QUEUE,RUNDIR_ROOT)
                        rocoto_create_compile_task \
                            (MACHINE_ID,COMPILE_ID,ROCOTO_COMPILE_MAXTRIES,MAKE_OPT,ACCNR,COMPILE_QUEUE,PARTITION,ROCOTO_XML)
                    else:
                        PASS_TESTS = True
                if (str(key) == 'tests' and COMPILE_ONLY == 'false' and not PASS_TESTS):
                    JOB_NR+=1
                    if ( ROCOTO ):
                        write_metatask_begin(COMPILE_ID, ROCOTO_XML)
                        case_count=0
                        for test in val:
                            case, config = get_testcase(test)
                            machine_check = machine_check_off(MACHINE_ID, config)
                            if machine_check:
                                TEST_NAME = case
                                TEST_ID   = TEST_NAME+'_'+RT_COMPILER
                                if 'dependency' in config.keys():
                                    DEP_RUN = str(config['dependency'])+'_'+RT_COMPILER
                                    if (delete_rundir == "true"): dependency_list.append(config['dependency'])
                                else:
                                    DEP_RUN = ""
                                RT_SUFFIX = ""
                                BL_SUFFIX = ""
                                os.environ["TEST_NAME"] = TEST_NAME
                                os.environ["DEP_RUN"]   = DEP_RUN
                                os.environ["TEST_ID"]   = TEST_ID
                                os.environ["RT_SUFFIX"] = RT_SUFFIX
                                os.environ["BL_SUFFIX"] = BL_SUFFIX
                                os.environ["JOB_NR"]    = str(JOB_NR)
                                if (CREATE_BASELINE == 'true'):
                                    if (DEP_RUN == ""):
                                        rc_set_run_task = subprocess.Popen(['bash', '-c', '. ufs_test_utils.sh; set_run_task'])
                                        rc_set_run_task.wait()
                                        case_count+=1
                                else:
                                    rc_set_run_task = subprocess.Popen(['bash', '-c', '. ufs_test_utils.sh; set_run_task'])
                                    rc_set_run_task.wait()
                                    case_count+=1
                        if int(case_count) > 0:
                            write_metatask_end(ROCOTO_XML)
                        else:
                            lines = subprocess.check_output(['head', '-n', '-1', ROCOTO_XML])
                            lines = lines.decode('utf-8')
                            with open(ROCOTO_XML, 'w') as filetowrite:
                                filetowrite.write(lines)
                            filetowrite.close
    rocoto_close=f"""</workflow>
"""
    with open(ROCOTO_XML,"a") as f:
        f.writelines(rocoto_close)
    f.close()

    REGRESSIONTEST_LOG = PATHRT+'/logs/RegressionTests_'+MACHINE_ID+'.log'
    make_loghead(ACCNR,MACHINE_ID,RUNDIR_ROOT,RTPWD,REGRESSIONTEST_LOG)

    if (delete_rundir == "true" and len(dependency_list) > 0):
        with open('keep_tests.tmp', 'w+') as fdep:
            for i in dependency_list:
                fdep.write(str(i) + '\n')
            fdep.close()
            
#if __name__ == "__main__":6

