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
    choices = ['hera.intel', 'orion.intel', 'gaea.intel', 'jet.intel', 'wcoss_dell_p3']
    parser.add_argument('-m', '--machine', help='Machine and Compiler combination', required=True, choices=choices, type=str)
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
    logger = logging.getLogger('MATCH_LABEL_WITH_ACTIONS')
    split_label = label.name.split('-')

    if len(split_label) != 3: return False
    if not re.match(split_label[0], 'Auto'): return False
    if not re.match(split_label[2], machine['name'].split('.')[0]): return False
    action_match = next((action for action in actions if re.match(action['name'], split_label[1])), False)

    return action_match


def get_preqs_with_actions(repos, machine, ghinterface_obj, actions):
    ''' Create list of dictionaries of a pull request and its machine label and action '''
    logger = logging.getLogger('GET_PREQS_WITH_ACTIONS')
    gh_preqs = [ghinterface_obj.client.get_repo(repo['address']).get_pulls(state='open', sort='created', base=repo['base']) for repo in repos]
    each_pr = [preq for gh_preq in gh_preqs for preq in gh_preq]
    preq_labels = [{'preq': pr, 'label': label} for pr in each_pr for label in pr.get_labels()]

    for i, pr_label in enumerate(preq_labels):
        match = match_label_with_action(machine, actions, pr_label['label'])
        if match:
            preq_labels[i]['action'] = match
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

    def add_pr_label(self):
        ''' adds the pull request label that initiated the job run to PR'''
        self.logger.info(f'Adding Label: {self.preq_dict["label"]}')
        self.preq_dict['preq'].add_to_labels(self.preq_dict['label'])

    def send_log_as_comment(self):
        with open('rt_auto.log', 'r') as f:
            filedata = f.read()
        self.preq_dict['preq'].create_issue_comment(filedata)

    def run_commands(self, commands_with_cwd):
        logger = logging.getLogger('RUN_COMMANDS')
        for command, in_cwd in commands_with_cwd:
            logger.info(f'Running "{command}" in location "{in_cwd}"')
            try:
                output = subprocess.Popen(command, shell=True, cwd=in_cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                out, err = output.communicate()
                out = [] if not out else out.decode('utf8').split('\n')
                err = [] if not err else err.decode('utf8').split('\n')
            except Exception as e:
                self.add_pr_label()
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
        rm_command = [f'rm -rf {pr_dir_str}', os.getcwd()]
        logger.info(f'Running "{rm_command}"')
        self.run_commands(rm_command)
        # try:
        #     output = subprocess.Popen(rm_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #     out, err = output.communicate()
        #     out = [] if not out else out.decode('utf8').split('\n')
        #     err = [] if not err else err.decode('utf8').split('\n')
        # except Exception as e:
        #     logger.warning('Removal of directory at end failed.')
        #     logger.warning(e)
        #     [logger.warning(f'stdout: {item}') for item in out if not None]
        #     [logger.warning(f'stderr: {eitem}') for eitem in err if not None]
        #     assert(e)
        # else:
        #     logger.info(f'Finished running: {rm_command}')
        #     [logger.debug(f'stdout: {item}') for item in out if not None]
        #     [logger.debug(f'stderr: {eitem}') for eitem in err if not None]

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
        # for command, in_cwd in create_repo_commands:
        #     logger.info(f'Running "{command}" in location "{in_cwd}"')
        #     try:
        #         output = subprocess.Popen(command, shell=True, cwd=in_cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #         out, err = output.communicate()
        #         out = [] if not out else out.decode('utf8').split('\n')
        #         err = [] if not err else err.decode('utf8').split('\n')
        #     except Exception as e:
        #         self.add_pr_label()
        #         logger.critical(e)
        #         [logger.critical(f'stdout: {item}') for item in out if not None]
        #         [logger.critical(f'stderr: {eitem}') for eitem in err if not None]
        #
        #         assert(e)
        #     else:
        #         logger.info(f'Finished running: {command}')
        #         [logger.info(f'stdout: {item}') for item in out if not None]
        #         [logger.info(f'stderr: {eitem}') for eitem in err if not None]

        logger.info('Finished repo clone')
        self.pr_repo_loc = repo_dir_str+"/"+repo_name
        return self.pr_repo_loc

    def runFunction(self):
        ''' Run the command associted with the label used to initiate this job '''
        logger = logging.getLogger('JOB/RUNFUNCTION')
        try:
            logger.info(f'Running: "{self.preq_dict["action"]["command"]}" in "{self.pr_repo_loc}"')
            output = subprocess.Popen(self.preq_dict['action']['command'], cwd=self.pr_repo_loc, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            out,err = output.communicate()
            out = [] if not out else out.decode('utf8').split('\n')
            err = [] if not err else err.decode('utf8').split('\n')
        except Exception as e:
            self.add_pr_label()
            logger.critical(e)
            [logger.critical(f'stdout: {item}') for item in out if not None]
            [logger.critical(f'stderr: {eitem}') for eitem in err if not None]
            assert(e)
        else:
            [logger.critical(f'stdout: {item}') for item in out if not None]
            [logger.critical(f'stderr: {eitem}') for eitem in err if not None]
            logger.debug(output.returncode)
            if output.returncode != 0:
                self.add_pr_label()
                logger.critical(f'{self.preq_dict["action"]["command"]} Failed')
                [logger.debug(f'stdout: {item}') for item in out if not None]
                [logger.debug(f'stderr: {eitem}') for eitem in err if not None]
            else:
                try:
                    logger.info(f'Attempting to run callback: {self.preq_dict["action"]["callback_fnc"]}')
                    getattr(self, self.preq_dict['action']['callback_fnc'])()
                except Exception as e:
                    self.add_pr_label()
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
        rt_log = f'tests/RegressionTests_{self.machine["name"]}.log'
        filepath = f'{self.pr_repo_loc}/{rt_log}'
        rm_filepath = '/'.join((self.pr_repo_loc.split('/'))[:-1])
        if os.path.exists(filepath):
            move_rt_commands = [
                [f'git add {rt_log}', self.pr_repo_loc],
                [f'git commit -m "Auto: Added Updated RT Log file: {rt_log}"', self.pr_repo_loc],
                [f'git pull --no-edit origin {self.branch}', self.pr_repo_loc],
                ['sleep 10', self.pr_repo_loc],
                [f'git push origin {self.branch}', self.pr_repo_loc]
            ]
            self.run_commands(move_rt_commands)
            # for command, in_cwd in move_rt_commands:
            #     try:
            #         logger.info(f'Attempting to run: {command}')
            #         output = subprocess.Popen(command, shell=True, cwd=in_cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            #         out, err = output.communicate()
            #         out = [] if not out else out.decode('utf8').split('\n')
            #         err = [] if not err else err.decode('utf8').split('\n')
            #     except Exception as e:
            #         self.add_pr_label()
            #         logger.critical(e)
            #         [logger.critical(f'stdout: {item}') for item in out if not None]
            #         [logger.critical(f'stderr: {eitem}') for eitem in err if not None]
            #         assert(e)
            #     else:
            #         logger.info(f'Finished command {command}')
            #         [logger.debug(f'stdout: {item}') for item in out if not None]
            #         [logger.debug(f'stderr: {eitem}') for eitem in err if not None]
        else:
            logger.critical('Could not find RT log')
            raise FileNotFoundError('Could not find RT log')

def main():

    # handle logging
    log_path = os.getcwd()
    log_filename = 'rt_auto.log'
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
        try:
            logger.info('Calling remove_pr_label')
            job.remove_pr_label()
            logger.info('Calling clone_pr_repo')
            job.clone_pr_repo()
            logger.info('Calling runFunction')
            job.runFunction()
            logger.info('Calling remove_pr_dir')
            job.remove_pr_dir()
            logger.info('Calling send_log_as_comment')
            job.send_log_as_comment()
        except Exception as e:
            job.add_pr_label()
            logger.critical(e)
            assert(e)

    logger.info('Script Finished')


if __name__ == '__main__':
    main()
