"""Automation of UFS Regression Testing

This script automates the process of UFS regression testing for code managers
at NOAA-EMC

This script should be started through rt_auto.sh so that env vars are set up
prior to start.
"""
from __future__ import print_function
from github import Github as gh
import argparse
import datetime
import subprocess
import re
import os
import logging

class GHInterface:
    '''
    This class stores information for communicating with GitHub
    ...

    Attributes
    ----------
    GHACCESSTOKEN : str
      API token to autheticate with GitHub
    client : pyGitHub communication object
      The connection to GitHub to make API requests
    '''
    def __init__(self):
        self.logger = logging.getLogger('GHINTERFACE')
        try:
            self.client = gh(os.getenv('ghapitoken'))
        except Exception as e:
            self.logger.critical(f'Exception is {e}')
            raise(e)

def parse_args_in():
    ''' Parse all input arguments coming from rt_auto.sh '''
    logger = logging.getLogger('PARSE_ARGS_IN')
    # Create Parse
    parser = argparse.ArgumentParser()

    # Setup Input Arguments
    choices = ['cheyenne', 'hera', 'orion', 'gaea', 'jet', 'wcoss_dell_p3']
    parser.add_argument('-m', '--machine', help='Machine name', required=True, choices=choices, type=str)
    parser.add_argument('-w', '--workdir', help='Working directory', required=True, type=str)

    # Get Arguments
    args = parser.parse_args()

    return args

def input_data(args):
    ''' Create dictionaries of data needed for processing UFS pull requests '''
    logger = logging.getLogger('INPUT_DATA')
    machine_dict = {
        'name': args.machine,
        'workdir': args.workdir
    }
    repo_list_dict = [{
        'name': 'ufs-weather-model',
        'address': 'ufs-community/ufs-weather-model',
        'base': 'develop'
    }]
    action_list_dict = [{
        'name': 'RT',
        'command': 'cd tests && /bin/bash --login ./rt.sh -e',
        'callback_fnc': 'move_rt_logs'
    }]

    return machine_dict, repo_list_dict, action_list_dict

def match_label_with_action(machine, actions, label):
    ''' Match the label that initiates a job with an action in the dict'''
    # <machine>-<compiler>-<test> i.e. hera-gnu-RT
    # RT = full regression test suite
    logger = logging.getLogger('MATCH_LABEL_WITH_ACTIONS')
    split_label = label.name.split('-')
    if len(split_label) != 3: return False, False #Make sure it has three parts
    if not re.match(split_label[0], machine['name']): return False, False #First check machine name matches
    compiler = split_label[1]
    if not str(compiler) in ["intel", "gnu"]: return False, False
    action_match = next((action for action in actions if re.match(action['name'], split_label[2])), False)
    action_match["command"] = f'export RT_COMPILER="{compiler}" && {action_match["command"]}'
    if split_label[2] == "RT" and compiler == "gnu":
        action_match["command"] = f'{action_match["command"]} -l rt_gnu.conf'
    return compiler, action_match


def get_preqs_with_actions(repos, machine, ghinterface_obj, actions):
    ''' Create list of dictionaries of a pull request and its machine label and action '''
    logger = logging.getLogger('GET_PREQS_WITH_ACTIONS')
    gh_preqs = [ghinterface_obj.client.get_repo(repo['address']).get_pulls(state='open', sort='created', base=repo['base']) for repo in repos]
    each_pr = [preq for gh_preq in gh_preqs for preq in gh_preq]
    preq_labels = [{'preq': pr, 'label': label} for pr in each_pr for label in pr.get_labels()]

    for i, pr_label in enumerate(preq_labels):
        compiler, match = match_label_with_action(machine, actions, pr_label['label'])
        if match:
            preq_labels[i]['action'] = match
            preq_labels[i]['compiler'] = compiler
        else:
            preq_labels[i] = False

    preq_dict = [x for x in preq_labels if x]

    return preq_dict

class Job:
    '''
    This class stores all information needed to run jobs on this machine.
    This class provides all methods needed to run all jobs.
    ...

    Attributes
    ----------
    preq_dict: dict
        Dictionary of all data that comes from the GitHub pull request
    ghinterface_obj: object
        An interface to GitHub setup through class GHInterface
    machine: dict
        Information about the machine the jobs will be running on
        provided by the bash script
    '''

    def __init__(self, preq_dict, ghinterface_obj, machine):
        self.logger = logging.getLogger('JOB')
        self.preq_dict = preq_dict
        self.ghinterface_obj = ghinterface_obj
        self.machine = machine

    def remove_pr_label(self):
        ''' Removes the pull request label that initiated the job run from PR '''
        self.logger.info(f'Removing Label: {self.preq_dict["label"]}')
        self.preq_dict['preq'].remove_from_labels(self.preq_dict['label'])

    def check_label_before_job_start(self):
        # LETS Check the label still exists before the start of the job in the
        # case of multiple jobs
        label_to_check = f'{self.machine["name"]}-{self.preq_dict["compiler"]}-{self.preq_dict["action"]["name"]}'
        labels = self.preq_dict['preq'].get_labels()
        label_match = next((label for label in labels if re.match(label.name, label_to_check)), False)

        return label_match


    def send_log_name_as_comment(self, log_filename):
        logger = logging.getLogger('JOB/SEND_LOG_NAME_AS_COMMENT')

        #Remove LAST MONTHS LOGS
        logger.info('Removing last months logs (if any)')
        last_month = datetime.date.today().replace(day=1) - datetime.timedelta(days=1)
        rm_command = [[f'rm rt_auto_*_{last_month.strftime("%Y%m")}*.log', os.getcwd()]]
        logger.info(f'Running "{rm_command}"')
        try:
            self.run_commands(rm_command)
        except Exception as e:
            logger.warning(f'"{rm_command}" failed with error:{e}')

        # Add log information to PR.
        comment_text = f'Log Name:{log_filename}\n'\
                       f'Log Location:{os.getcwd()}\n'\
                       'Logs are kept for one month'
        try:
            self.preq_dict['preq'].create_issue_comment(comment_text)
        except Exception as e:
            logger.warning('Creating comment with log location failed with:{e}')
        else:
            logger.info(f'{comment_text}')

    def run_commands(self, commands_with_cwd):
        logger = logging.getLogger('JOB/RUN_COMMANDS')
        for command, in_cwd in commands_with_cwd:
            logger.info(f'Running "{command}" in location "{in_cwd}"')
            try:
                output = subprocess.Popen(command, shell=True, cwd=in_cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                out, err = output.communicate()
                out = [] if not out else out.decode('utf8').split('\n')
                err = [] if not err else err.decode('utf8').split('\n')
            except Exception as e:
                logger.critical(e)
                [logger.critical(f'stdout: {item}') for item in out if not None]
                [logger.critical(f'stderr: {eitem}') for eitem in err if not None]
                assert(e)
            else:
                logger.info(f'Finished running: {command}')
                [logger.debug(f'stdout: {item}') for item in out if not None]
                [logger.debug(f'stderr: {eitem}') for eitem in err if not None]

    def remove_pr_dir(self):
        logger = logging.getLogger('JOB/REMOVE_PR_DIR')
        pr_dir_str = f'{self.machine["workdir"]}/{str(self.preq_dict["preq"].id)}'
        rm_command = [[f'rm -rf {pr_dir_str}', self.pr_repo_loc]]
        logger.info(f'Running "{rm_command}"')
        self.run_commands(rm_command)

    def clone_pr_repo(self):
        ''' clone the GitHub pull request repo, via command line '''
        logger = logging.getLogger('JOB/CLONE_PR_REPO')
        repo_name = self.preq_dict['preq'].head.repo.name
        self.branch = self.preq_dict['preq'].head.ref
        git_url = self.preq_dict['preq'].head.repo.html_url.split('//')
        git_url = f'{git_url[0]}//${{ghapitoken}}@{git_url[1]}'
        logger.info(f'GIT URL: {git_url}')
        logger.info('Starting repo clone')
        repo_dir_str = f'{self.machine["workdir"]}/{str(self.preq_dict["preq"].id)}/{datetime.datetime.now().strftime("%Y%m%d%H%M%S")}'

        create_repo_commands = [
            [f'mkdir -p "{repo_dir_str}"', self.machine['workdir']],
            [f'git clone -b {self.branch} {git_url}', repo_dir_str],
            [f'git submodule update --init --recursive', f'{repo_dir_str}/{repo_name}']
        ]

        self.run_commands(create_repo_commands)

        logger.info('Finished repo clone')
        self.pr_repo_loc = repo_dir_str+"/"+repo_name
        return self.pr_repo_loc

    def run_function(self):
        ''' Run the command associted with the label used to initiate this job '''
        logger = logging.getLogger('JOB/RUN_FUNCTION')
        compiler = self.preq_dict['compiler']
        logger.info(f'Compiler being used for command is {compiler}')
        command = self.preq_dict["action"]["command"]
        try:
            logger.info(f'Running: "{command}" in "{self.pr_repo_loc}"')
            output = subprocess.Popen(command, cwd=self.pr_repo_loc, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            out,err = output.communicate()
            out = [] if not out else out.decode('utf8').split('\n')
            err = [] if not err else err.decode('utf8').split('\n')
        except Exception as e:
            logger.critical(e)
            [logger.critical(f'stdout: {item}') for item in out if not None]
            [logger.critical(f'stderr: {eitem}') for eitem in err if not None]
            assert(e)
        else:
            if output.returncode != 0:
                comment_text = f'Script rt.sh failed \n'\
                               f'location: {self.pr_repo_loc} \n'\
                               f'machine: {self.machine["name"]} \n'\
                               f'compiler: {self.preq_dict["compiler"]}\n'\
                               f'STDOUT: {out} \n'\
                               f'STDERR: {err}'
                self.preq_dict['preq'].create_issue_comment(comment_text)
                logger.critical(f'{command} Failed')
                [logger.critical(f'stdout: {item}') for item in out if not None]
                [logger.critical(f'stderr: {eitem}') for eitem in err if not None]
            else:
                try:
                    logger.info(f'Attempting to run callback: {self.preq_dict["action"]["callback_fnc"]}')
                    getattr(self, self.preq_dict['action']['callback_fnc'])()
                except Exception as e:
                    logger.critical(f'Callback function {self.preq_dict["action"]["callback_fnc"]} failed with "{e}"')
                    goodexit = False
                    assert(e)
                else:
                    logger.info(f'Finished callback {self.preq_dict["action"]["callback_fnc"]}')
                    [logger.debug(f'stdout: {item}') for item in out if not None]
                    [logger.debug(f'stderr: {eitem}') for eitem in err if not None]

    # Add Callback Functions After Here
    def move_rt_logs(self):
        ''' This is the callback function associated with the "RT" command '''
        logger = logging.getLogger('JOB/MOVE_RT_LOGS')
        rt_log = f'tests/RegressionTests_{self.machine["name"]}.{self.preq_dict["compiler"]}.log'
        filepath = f'{self.pr_repo_loc}/{rt_log}'
        if os.path.exists(filepath):
            #check_for_success
            with open(filepath) as f:
                if 'SUCCESSFUL' in f.read():
                    move_rt_commands = [
                        [f'git pull --ff-only origin {self.branch}', self.pr_repo_loc],
                        [f'git add {rt_log}', self.pr_repo_loc],
                        [f'git commit -m "PASSED: {self.machine["name"]}.{self.preq_dict["compiler"]}. Log file uploaded. skip-ci"', self.pr_repo_loc],
                        ['sleep 10', self.pr_repo_loc],
                        [f'git push origin {self.branch}', self.pr_repo_loc]
                    ]
                    self.run_commands(move_rt_commands)
                else:
                    comment_text = f'REGRESSION TEST FAILED \n'\
                                   f'location: {self.pr_repo_loc} \n'\
                                   f'machine: {self.machine["name"]} \n'\
                                   f'compiler: {self.preq_dict["compiler"]}\n'\
                                   'Please make changes and add the following label back: '\
                                   f'{self.machine["name"]}-{self.preq_dict["compiler"]}-{self.preq_dict["action"]["name"]}'
                    self.preq_dict['preq'].create_issue_comment(comment_text)

        else:
            logger.critical(f'Could not find {self.machine["name"]}.{self.preq_dict["compiler"]} RT log')
            raise FileNotFoundError(f'Could not find {self.machine["name"]}.{self.preq_dict["compiler"]} RT log')

def main():

    # handle logging
    log_path = os.getcwd()
    log_filename = f'rt_auto_{datetime.datetime.now().strftime("%Y%m%d%H%M%S")}.log'
    # Please don't run the following on cron with level=logging.DEBUG
    # as it exposes the GH API Token
    # Only set it to DEBUG while debugging
    logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)
    logger = logging.getLogger('MAIN')
    logger.info('Starting Script')
    # handle input args
    logger.info('Parsing input args')
    args = parse_args_in()

    # get input data
    logger.info('Calling input_data().')
    machine, repos, actions = input_data(args)

    # setup interface with GitHub
    logger.info('Setting up GitHub interface.')
    ghinterface_obj = GHInterface()

    # get all pull requests from the GitHub object
    logger.info('Getting all pull requests, labels and actions applicable to this machine.')
    preq_dict = get_preqs_with_actions(repos, machine, ghinterface_obj, actions)
    # add Job objects and run them
    logger.info('Adding all jobs to an object list and running them.')
    jobs = [Job(pullreq, ghinterface_obj, machine) for pullreq in preq_dict]
    for job in jobs:
        logger.info(f'Starting Job: {job}')
        if job.check_label_before_job_start():
            try:
                logger.info('Calling remove_pr_label')
                job.remove_pr_label()
                logger.info('Calling clone_pr_repo')
                job.clone_pr_repo()
                logger.info('Calling run_function')
                job.run_function()
                logger.info('Calling remove_pr_dir')
                # job.remove_pr_dir()
                # logger.info('Calling send_log_name_as_comment')
                job.send_log_name_as_comment(log_filename)
            except Exception as e:
                logger.critical(e)
                assert(e)

    logger.info('Script Finished')


if __name__ == '__main__':
    main()
