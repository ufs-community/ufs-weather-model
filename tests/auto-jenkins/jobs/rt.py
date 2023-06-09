# Imports
import datetime
import logging
import os
import re

#Logging filter to santize log output and ensure Github Access tokens are not leaked.
def loggingfilter(record):
    if re.search('ghp_(.*)\@', str(record.msg)):
       record.msg = re.sub(r'ghp_(.*)\@', r'ghp_****@', str(record.msg))
    elif re.search('github_pat_(.*)\@', str(record.msg)):
       record.msg = re.sub(r'github_pat_(.*)\@', r'github_pat_****@', str(record.msg))

    return True

def run(job_obj):
    logger = logging.getLogger('RT/RUN')
    workdir = set_directories(job_obj)
    branch, pr_repo_loc, repo_dir_str = clone_pr_repo(job_obj, workdir)
    run_regression_test(job_obj, pr_repo_loc)
    post_process(job_obj, pr_repo_loc, repo_dir_str, branch)


def set_directories(job_obj):
    logger = logging.getLogger('RT/SET_DIRECTORIES')
    if job_obj.machine == 'hera':
        workdir = '/scratch1/NCEPDEV/nems/role.epic/autort/pr'
    elif job_obj.machine == 'jet':
        workdir = '/lfs4/HFIP/hfv3gfs/role.epic/autort/pr'
    elif job_obj.machine == 'gaea':
        workdir = '/lustre/f2/pdata/ncep/role.epic/autort/pr'
    elif job_obj.machine == 'orion':
        workdir = '/work/noaa/epic-ps/role-epic-ps/autort/tests/auto/pr'
    elif job_obj.machine == 'cheyenne':
        workdir = '/glade/scratch/epicufsrt/autort/jenkins/autort/pr'
    else:
        print(f'Machine {job_obj.machine} is not supported for this job')
        raise KeyError

    logger.info(f'machine: {job_obj.machine}')
    logger.info(f'workdir: {workdir}')

    return workdir


def run_regression_test(job_obj, pr_repo_loc):
    logger = logging.getLogger('RT/RUN_REGRESSION_TEST')
    if job_obj.machine != 'hera':
        rt_command = [[f'cd tests && /bin/bash --login ./rt.sh -e', pr_repo_loc]]
    elif job_obj.machine == 'hera':
        rt_command = [[f'cd tests && /bin/bash --login ./rt.sh -r', pr_repo_loc]]
    job_obj.run_commands(logger, rt_command)


def remove_pr_data(job_obj, pr_repo_loc, repo_dir_str, rt_dir):
    logger = logging.getLogger('RT/REMOVE_PR_DATA')
    rm_command = [
                 [f'rm -rf {rt_dir}', pr_repo_loc],
                 [f'rm -rf {repo_dir_str}', pr_repo_loc]
                 ]
    job_obj.run_commands(logger, rm_command)


def clone_pr_repo(job_obj, workdir):
    ''' clone the GitHub pull request repo, via command line '''
    filename = ('jenkinsaccesstoken')
    f = open(filename)
    os.environ['ghapitoken'] = f.readline().strip('\n')

    logger = logging.getLogger('RT/CLONE_PR_REPO')
    logger.addFilter(loggingfilter)
    repo_name = job_obj.preq_dict['preq'].head.repo.name
    branch = job_obj.preq_dict['preq'].head.ref
    git_url = job_obj.preq_dict['preq'].head.repo.html_url.split('//')
    git_token = os.environ['ghapitoken']
    git_https_url = f'{git_url[0]}//{git_token}@{git_url[1]}'
    git_ssh_url = job_obj.preq_dict['preq'].head.repo.ssh_url
    logger.debug(f'GIT SSH_URL: {git_ssh_url}')
    logger.info('Starting repo clone')
    repo_dir_str = f'{workdir}/'\
                   f'{str(job_obj.preq_dict["preq"].id)}/'\
                   f'{datetime.datetime.now().strftime("%Y%m%d%H%M%S")}'
    pr_repo_loc = f'{repo_dir_str}/{repo_name}'
    job_obj.comment_text_append(f'[RT] Repo location: {pr_repo_loc}')
    create_repo_commands = [
        [f'mkdir -p "{repo_dir_str}"', os.getcwd()],
        [f'git clone -b {branch} {git_ssh_url}', repo_dir_str],
        ['git submodule update --init --recursive',
         f'{repo_dir_str}/{repo_name}'],
        [f'git remote add httpsorigin {git_https_url}',
        f'{repo_dir_str}/{repo_name}'],
        ['git config user.email "ecc.platform@noaa.gov"',
         f'{repo_dir_str}/{repo_name}'],
        ['git config user.name "epic-cicd-jenkins"',
         f'{repo_dir_str}/{repo_name}']
    ]

    job_obj.run_commands(logger, create_repo_commands)

    logger.info('Finished repo clone')
    return branch, pr_repo_loc, repo_dir_str


def post_process(job_obj, pr_repo_loc, repo_dir_str, branch):
    ''' This is the callback function associated with the "RT" command '''
    logger = logging.getLogger('RT/MOVE_RT_LOGS')
    rt_log = f'tests/logs/RegressionTests_{job_obj.machine}.log'
    filepath = f'{pr_repo_loc}/{rt_log}'
    rt_dir, logfile_pass = process_logfile(job_obj, filepath)
    if logfile_pass:
        #if job_obj.preq_dict['preq'].maintainer_can_modify:
        move_rt_commands = [
            [f'git pull --ff-only origin {branch}', pr_repo_loc],
            [f'git add {rt_log}', pr_repo_loc],
            [f'git commit -m "[AutoRT] {job_obj.machine} Job Completed.\n\n\n'
              'on-behalf-of @ufs-community <ecc.platform@noaa.gov>"',
             pr_repo_loc],
            ['sleep 10', pr_repo_loc],
            [f'git push httpsorigin {branch}', pr_repo_loc]
        ]
        job_obj.run_commands(logger, move_rt_commands)
    else:
        job_obj.comment_text_append(f'[RT] Log file shows failures.')
        job_obj.comment_text_append(f'[RT] Please obtain logs from {pr_repo_loc}')
        job_obj.preq_dict['preq'].create_issue_comment(job_obj.comment_text)


def process_logfile(job_obj, logfile):
    logger = logging.getLogger('RT/PROCESS_LOGFILE')
    rt_dir = []
    fail_string_list = ['Test', 'failed']
    if os.path.exists(logfile):
        with open(logfile) as f:
            for line in f:
                if all(x in line for x in fail_string_list):
                # if 'FAIL' in line and 'Test' in line:
                    job_obj.comment_text_append(f'[RT] Error: {line.rstrip(chr(10))}')
                elif 'working dir' in line and not rt_dir:
                    rt_dir = os.path.split(line.split()[-1])[0]
                elif 'SUCCESSFUL' in line:
                    return rt_dir, True
        job_obj.job_failed(logger, f'{job_obj.preq_dict["action"]}')
    else:
        logger.critical(f'Could not find {job_obj.machine}'
                        f'{job_obj.preq_dict["action"]} log')
        print(f'Could not find {job_obj.machine}'
              f'{job_obj.preq_dict["action"]} log')
        raise FileNotFoundError
