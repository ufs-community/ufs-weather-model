"""Automation of UFS Regression Testing

This script automates the process of UFS regression testing for code managers
at NOAA-EMC

The data set required for this code to operate properly is rt_auto.yml:
  - Provides all the repository information in which to search through
  - Provides all information about the machines this code should be run on
  - Provides all acceptable commands to be run on machines

This script should be started through rt_auto.sh so that env vars are set up
prior to start.
"""
from __future__ import print_function
from github import Github as gh
import argparse
import datetime
import time
import threading
import subprocess
import re
import yaml
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
    def __init__(self) -> None:
        self.logger = logging.getLogger("GHInterface")
        try:
            with open('accesstoken.txt', 'rb') as f:
                filedata = f.read()
        except FileNotFoundError as e:
            ERROR_MESSAGE = f'Please create a file "accesstoken.txt" that'\
                ' contains your GitHub API Token only.\n'\
                'Make sure to set permission so others can not read it (400)'
            self.logger.critical(ERROR_MESSAGE)
            raise FileNotFoundError(e)
        else:
            self.GHACCESSTOKEN = str(filedata)[2:-3]
            self.client = gh(self.GHACCESSTOKEN)
            try:
                self.logger.debug(f'GHUSERNAME is {self.client.get_user().login}')
            except Exception as e:
                self.logger.critical(f'Exception is {e}')
                raise Exception(e)

def process_pr(pullreq_obj, ghinterface_obj, machine_obj, functions):
    '''This method processes each PulLReq obj by checking if the
    label in the PR is approved and clones the repository

    '''
    logger = logging.getLogger(f'PR#{pullreq_obj.preq_obj.id}')
    logger.debug(f'Start')
    for prlabel in pullreq_obj.labels:
        if prlabel.is_approved(machine_obj, functions):
            try:
                logger.debug(f'Removing Label Auto-{prlabel.name}-{prlabel.machine}')
                pullreq_obj.preq_obj.remove_from_labels(f'Auto-{prlabel.name}-{prlabel.machine}')
                clone_pr_repo(pullreq_obj, ghinterface_obj, machine_obj)
                goodexit = runFunction(prlabel.action_obj, pullreq_obj)
            except KeyboardInterrupt:
                print(f'KEY BOAR DINTERRUPT')
                pullreq_obj.preq_obj.add_to_labels(f'Auto-{prlabel.name}-{prlabel.machine}')
                raise
            except Exception as e:
                logger.critical(f'ERROR RUNNING RT {prlabel.action_obj.command} with error: {e}')
                logger.debug(f'Adding back label Auto-{prlabel.name}-{prlabel.machine}')
                pullreq_obj.preq_obj.add_to_labels(f'Auto-{prlabel.name}-{prlabel.machine}')
                raise
            else:
                if goodexit == False:
                    logger.debug(f'runFunction exited with error')
                    pullreq_obj.preq_obj.add_to_labels(f'Auto-{prlabel.name}-{prlabel.machine}')
                    logger.debug(f'Adding back label Auto-{prlabel.name}-{prlabel.machine}')
    logger.debug(f'End')

def runFunction(action_obj, pullreq_obj):
    goodexit = None
    proc = subprocess.Popen(action_obj.command, shell=True, cwd=pullreq_obj.clone_dir)
    proc.wait()
    proc.poll()
    if proc.returncode == 0:
        try:
            globals()[action_obj.callback](pullreq_obj)
        except:
            goodexit = False
        else:
            goodexit = True
    return goodexit

def parse_args_in():
    # Create Parse
    parser = argparse.ArgumentParser()

    # Setup Input Arguments
    parser.add_argument("machine_name", help="Machine name in <machine>.<compiler> format", type=str)
    parser.add_argument("workdir", help="Working directory for the machine", type=str)

    # Get Arguments
    args = parser.parse_args()

    # Check incoming args for proper formatting and type
    if type(args.workdir) != str or type(args.machine_name) != str:
        raise TypeError('All arguments need to be of type str')
    if len(args.machine_name.split('.'))!=2:
        raise argparse.ArgumentTypeError("Please use <machine>.<compiler> format for machine_name")

    return args

def input_data(args):
    machine_dict = {
        'name': args.machine_name,
        'workdir': args.workdir
    }
    repo_list_dict = [{
        'name': 'ufs-weather-model',
        'address': 'ufs-community/ufs-weather-model',
        'base': 'develop'
    }]
    action_list_dict = [{
        'name': 'RT',
        'command': './tests/rt.sh -e',
        'callback_fnc': 'move_rt_logs'
    }]

    return machine_dict, repo_list_dict, action_list_dict

def match_label_with_action(machine, actions, label):

    split_label = label.name.split('-')

    if len(split_label) != 3: return False
    if not re.match(split_label[0], 'Auto'): return False
    if not re.match(split_label[2], machine['name'].split('.')[0]): return False
    action_match = next((action for action in actions if re.match(action['name'], split_label[1])), False)

    return action_match


def get_preqs_with_actions(repos, machine, ghinterface_obj, actions):

    gh_preqs = [ghinterface_obj.client.get_repo(repo['address']).get_pulls(state='open', sort='created', base=repo['base']) for repo in repos]
    each_pr = [preq for gh_preq in gh_preqs for preq in gh_preq]
    preq_labels = [{'preq': pr, 'label': label} for pr in each_pr for label in pr.get_labels()]

    for i, pr_label in enumerate(preq_labels):
        match = match_label_with_action(machine, actions, pr_label['label'])
        if match:
            preq_labels[i]['action'] = match
        else:
            preq_labels[i] = False

    preq_labels = [x for x in preq_labels if x]

    return preq_labels

class Job:

    def __init__(self, pullreq_obj, ghinterface_obj, machine):
        self.pullreq_obj = pullreq_obj
        self.ghinterface_obj = ghinterface_obj
        self.machine = machine

    def clone_pr_repo(self):

        logger = logging.getLogger("clone_pr_repo()")
        repo_name = self.pullreq_obj['preq'].head.repo.name
        self.branch = self.pullreq_obj['preq'].head.ref
        git_url = self.pullreq_obj['preq'].head.repo.html_url.split('//')
        git_url = f'{git_url[0]}//{self.ghinterface_obj.GHACCESSTOKEN}@{git_url[1]}'
        repo_dir_str = f'{self.machine["workdir"]}/{str(self.pullreq_obj["preq"].id)}/{datetime.datetime.now().strftime("%Y%m%d%H%M%S")}'

        create_repo_commands = [
            [f'mkdir -p "{repo_dir_str}"', self.machine['workdir']],
            [f'git clone -b {self.branch} {git_url}', repo_dir_str],
            [f'git submodule update --init --recursive', f'{repo_dir_str}/{repo_name}']
        ]

        for command, in_cwd in create_repo_commands:
            logger.info(f'Attempting to run: {command}')
            try:
                output = subprocess.check_output(command, shell=True, cwd=in_cwd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                logger.critical(e)
                logger.critical(f'STDOUT: {output.stdout}')
                logger.critical(f'STDERR: {output.stderr}')
                assert(e)
            else:
                logger.info(f'Finished running: {command}')

        logger.debug(f'Finished Cloning {git_url}')
        self.pr_repo_loc = repo_dir_str+"/"+repo_name
        return self.pr_repo_loc

    def runFunction(self):
        logger = logging.getLogger("runFunction()")
        try:
            output = subprocess.check_output(self.pullreq_obj['action']['command'], cwd=self.pr_repo_loc, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            logger.critical(e)
            logger.critical(f'STDOUT: {output.stdout}')
            logger.critical(f'STDERR: {output.stderr}')
            assert(e)
        else:
            try:
                globals()[self.pullreq_obj['action']['callback']](self)
            except Exception as e:
                logger.critical(f'Callback function {self.pullreq_obj["action"]["callback"]} failed')
                logger.critical(e)
                goodexit = False
            else:
                logger.info(f'Finished callback {self.pullreq_obj["action"]["callback"]}')

    # Callback Function After Here
    def move_rt_logs(self):
        logger = logging.getLogger("clone_pr_repo()")
        rt_log = 'tests/RegressionTests_'+self.machine.name+'.log'
        filepath = self.pr_repo_loc+'/'+rt_log
        rm_filepath = '/'.join((self.pr_repo_loc.split('/'))[:-1])
        if os.path.exists(filepath):
            move_rt_commands = [
                ['git add '+rt_log, self.pr_repo_loc],
                ['git commit -m "Auto: Added Updated RT Log file: '+rt_log+'"', self.pr_repo_loc],
                ['git pull --no-edit origin '+self.branch, self.pr_repo_loc],
                ['sleep 10', self.pr_repo_loc],
                ['git push origin '+self.branch, self.pr_repo_loc]
            ]
            for command, in_cwd in move_rt_commands:
                try:
                    output = subprocess.check_output(command, shell=True, cwd=in_cwd, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as e:
                    logger.critical(e)
                    logger.critical(f'STDOUT: {output.stdout}')
                    logger.critical(f'STDERR: {output.stderr}')
                    assert(e)
                else:
                    logger.info(f'Finished command {command}')
        else:
            logger.critical('Could not find RT log')
            raise FileNotFoundError('Could not find RT log')

def main():
    # handle input args
    args = parse_args_in()

    # handle logging
    logging.basicConfig(filename="rt_auto.log", filemode='w', level=logging.DEBUG)
    logger = logging.getLogger("main")

    # get input data
    machine, repos, actions = input_data(args)

    # setup interface with GitHub
    ghinterface_obj = GHInterface()

    # get all pull requests from the GitHub object
    full_preqs = get_preqs_with_actions(repos, machine, ghinterface_obj, actions)

    # add Job objects and run them
    jobs = [Job(pullreq_obj, ghinterface_obj, machine) for pullreq_obj in full_preqs]
    for job in jobs:
        job.clone_pr_repo()
        job.runFunction()


if __name__ == '__main__':
    main()
