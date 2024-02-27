# Imports
import datetime
import logging
import os


def run(job_obj):
    logger = logging.getLogger('RT/RUN')
    workdir = set_directories(job_obj)
    branch, pr_repo_loc, repo_dir_str = clone_pr_repo(job_obj, workdir)
    run_regression_test(job_obj, pr_repo_loc)
    post_process(job_obj, pr_repo_loc, repo_dir_str, branch)


def set_directories(job_obj):
    logger = logging.getLogger('RT/SET_DIRECTORIES')
    if job_obj.machine == 'hera':
        workdir = '/scratch1/NCEPDEV/nems/emc.nemspara/autort/pr'
    elif job_obj.machine == 'jet':
        workdir = '/lfs4/HFIP/h-nems/emc.nemspara/autort/pr'
    elif job_obj.machine == 'gaea':
        workdir = '/lustre/f2/pdata/ncep/emc.nemspara/autort/pr'
    elif job_obj.machine == 'orion':
        workdir = '/work/noaa/nems/emc.nemspara/autort/pr'
    elif job_obj.machine == 'derecho':
        workdir = '/glade/scratch/dtcufsrt/autort/tests/auto/pr'
    else:
        print(f'Machine {job_obj.machine} is not supported for this job')
        raise KeyError

    logger.info(f'machine: {job_obj.machine}')
    logger.info(f'workdir: {workdir}')

    return workdir


def run_regression_test(job_obj, pr_repo_loc):
    logger = logging.getLogger('RT/RUN_REGRESSION_TEST')
    if job_obj.compiler == 'gnu':
        rt_command = [[f'export RT_COMPILER="{job_obj.compiler}" && cd tests '
                       '&& /bin/bash --login ./rt.sh -e -l rt_gnu.conf',
                       pr_repo_loc]]
    elif job_obj.compiler == 'intel':
        rt_command = [[f'export RT_COMPILER="{job_obj.compiler}" && cd tests '
                       '&& /bin/bash --login ./rt.sh -e', pr_repo_loc]]
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
    logger = logging.getLogger('RT/CLONE_PR_REPO')
    repo_name = job_obj.preq_dict['preq'].head.repo.name
    branch = job_obj.preq_dict['preq'].head.ref
    git_url = job_obj.preq_dict['preq'].head.repo.html_url.split('//')
    git_url = f'{git_url[0]}//${{ghapitoken}}@{git_url[1]}'
    logger.debug(f'GIT URL: {git_url}')
    logger.info('Starting repo clone')
    repo_dir_str = f'{workdir}/'\
                   f'{str(job_obj.preq_dict["preq"].id)}/'\
                   f'{datetime.datetime.now().strftime("%Y%m%d%H%M%S")}'
    pr_repo_loc = f'{repo_dir_str}/{repo_name}'
    job_obj.comment_text_append(f'Repo location: {pr_repo_loc}')
    create_repo_commands = [
        [f'mkdir -p "{repo_dir_str}"', os.getcwd()],
        [f'git clone -b {branch} {git_url}', repo_dir_str],
        ['git submodule update --init --recursive',
         f'{repo_dir_str}/{repo_name}'],
        ['git config user.email "brian.curtis@noaa.gov"',
         f'{repo_dir_str}/{repo_name}'],
        ['git config user.name "Brian Curtis"',
         f'{repo_dir_str}/{repo_name}']
    ]

    job_obj.run_commands(logger, create_repo_commands)

    logger.info('Finished repo clone')
    return branch, pr_repo_loc, repo_dir_str


def post_process(job_obj, pr_repo_loc, repo_dir_str, branch):
    ''' This is the callback function associated with the "RT" command '''
    logger = logging.getLogger('RT/MOVE_RT_LOGS')
    rt_log = f'tests/RegressionTests_{job_obj.machine}'\
             f'.{job_obj.compiler}.log'
    filepath = f'{pr_repo_loc}/{rt_log}'
    rt_dir, logfile_pass = process_logfile(job_obj, filepath)
    if logfile_pass:
        if job_obj.preq_dict['preq'].maintainer_can_modify:
            move_rt_commands = [
                [f'git pull --ff-only origin {branch}', pr_repo_loc],
                [f'git add {rt_log}', pr_repo_loc],
                [f'git commit -m "RT JOBS PASSED: {job_obj.machine}'
                 f'.{job_obj.compiler}. Log file uploaded.\n\n'
                  'on-behalf-of @ufs-community"',
                 pr_repo_loc],
                ['sleep 10', pr_repo_loc],
                [f'git push origin {branch}', pr_repo_loc]
            ]
            job_obj.run_commands(logger, move_rt_commands)
            remove_pr_data(job_obj, pr_repo_loc, repo_dir_str, rt_dir)
        else:
            job_obj.comment_text_append(f'Cannot upload {job_obj.machine}.'\
                                        f'{job_obj.compiler} RT Log'\
                                        'It is blocked by PR owner')
            job_obj.comment_text_append(f'Please obtain logs from {pr_repo_loc}')
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
                    job_obj.comment_text_append(f'{line.rstrip(chr(10))}')
                elif 'working dir' in line and not rt_dir:
                    rt_dir = os.path.split(line.split()[-1])[0]
                    job_obj.comment_text_append(f'Please manually delete: '
                                                f'{rt_dir}')
                elif 'SUCCESSFUL' in line:
                    return rt_dir, True
        job_obj.job_failed(logger, f'{job_obj.preq_dict["action"]}')
    else:
        logger.critical(f'Could not find {job_obj.machine}'
                        f'.{job_obj.compiler} '
                        f'{job_obj.preq_dict["action"]} log')
        print(f'Could not find {job_obj.machine}.{job_obj.compiler} '
              f'{job_obj.preq_dict["action"]} log')
        raise FileNotFoundError
