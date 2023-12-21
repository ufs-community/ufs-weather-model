"""Automation of UFS Regression Testing

This script automates the process of UFS regression testing for code managers
at NOAA-EMC

This script should be started through rt_auto.sh so that env vars are set up
prior to start.
"""
from github import Github as gh
import datetime
import subprocess
import re
import os
import logging
import importlib


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

        filename = 'accesstoken'

        if os.path.exists(filename):
            if oct(os.stat(filename).st_mode)[-3:] != 600:
                with open(filename) as f:
                    os.environ['ghapitoken'] = f.readline().strip('\n')
            else:
                raise Exception('File permission needs to be "600" ')
        else:
            raise FileNotFoundError('Cannot find file "accesstoken"')

        try:
            self.client = gh(os.getenv('ghapitoken'))
        except Exception as e:
            self.logger.critical(f'Exception is {e}')
            raise(e)


def set_action_from_label(machine, actions, label):
    ''' Match the label that initiates a job with an action in the dict'''
    # <machine>-<compiler>-<test> i.e. hera-gnu-RT
    logger = logging.getLogger('MATCH_LABEL_WITH_ACTIONS')
    logger.info('Setting action from Label')
    split_label = label.name.split('-')
    # Make sure it has three parts
    if len(split_label) != 3:
        return False, False
    # Break the parts into their variables
    label_machine = split_label[0]
    label_compiler = split_label[1]
    label_action = split_label[2]
    # check machine name matches
    if not re.match(label_machine, machine):
        return False, False
    # Compiler must be intel or gnu
    if not str(label_compiler) in ["intel", "gnu"]:
        return False, False
    action_match = next((action for action in actions
                         if re.match(action, label_action)), False)

    logging.info(f'Compiler: {label_compiler}, Action: {action_match}')
    return label_compiler, action_match


def get_preqs_with_actions(repos, machine, ghinterface_obj, actions):
    ''' Create list of dictionaries of a pull request
        and its machine label and action '''
    logger = logging.getLogger('GET_PREQS_WITH_ACTIONS')
    logger.info('Getting Pull Requests with Actions')
    gh_preqs = [ghinterface_obj.client.get_repo(repo['address'])
                .get_pulls(state='open', sort='created', base=repo['base'])
                for repo in repos]
    each_pr = [preq for gh_preq in gh_preqs for preq in gh_preq]
    preq_labels = [{'preq': pr, 'label': label} for pr in each_pr
                   for label in pr.get_labels()]

    jobs = []
    # return_preq = []
    for pr_label in preq_labels:
        compiler, match = set_action_from_label(machine, actions,
                                                pr_label['label'])
        if match:
            pr_label['action'] = match
            # return_preq.append(pr_label.copy())
            jobs.append(Job(pr_label.copy(), ghinterface_obj, machine, compiler))

    return jobs


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

    def __init__(self, preq_dict, ghinterface_obj, machine, compiler):
        self.logger = logging.getLogger('JOB')
        self.preq_dict = preq_dict
        self.job_mod = importlib.import_module(
                       f'jobs.{self.preq_dict["action"].lower()}')
        self.ghinterface_obj = ghinterface_obj
        self.machine = machine
        self.compiler = compiler
        self.comment_text = ''
        self.failed_tests = []

    def comment_text_append(self, newtext):
        self.comment_text += f'{newtext}\n'

    def remove_pr_label(self):
        ''' Removes the PR label that initiated the job run from PR '''
        self.logger.info(f'Removing Label: {self.preq_dict["label"]}')
        self.preq_dict['preq'].remove_from_labels(self.preq_dict['label'])

    def check_label_before_job_start(self):
        # LETS Check the label still exists before the start of the job in the
        # case of multiple jobs
        label_to_check = f'{self.machine}'\
                         f'-{self.compiler}'\
                         f'-{self.preq_dict["action"]}'
        labels = self.preq_dict['preq'].get_labels()
        label_match = next((label for label in labels
                            if re.match(label.name, label_to_check)), False)

        return label_match

    def run_commands(self, logger, commands_with_cwd):
        for command, in_cwd in commands_with_cwd:
            logger.info(f'Running `{command}`')
            logger.info(f'in location "{in_cwd}"')
            try:
                output = subprocess.Popen(command, shell=True, cwd=in_cwd,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT)
            except Exception as e:
                self.job_failed(logger, 'subprocess.Popen')
            else:
                try:
                    out, err = output.communicate()
                    out = [] if not out else out.decode('utf8').split('\n')
                    logger.info(out)
                except Exception as e:
                    err = [] if not err else err.decode('utf8').split('\n')
                    self.job_failed(logger, f'Command {command}', exception=e,
                                    STDOUT=True, out=out, err=err)
                else:
                    logger.info(f'Finished running: {command}')

    def run(self):
        logger = logging.getLogger('JOB/RUN')
        logger.info(f'Starting Job: {self.preq_dict["label"]}')
        self.comment_text_append(newtext=f'Machine: {self.machine}')
        self.comment_text_append(f'Compiler: {self.compiler}')
        self.comment_text_append(f'Job: {self.preq_dict["action"]}')
        if self.check_label_before_job_start():
            try:
                logger.info('Calling remove_pr_label')
                self.remove_pr_label()
                logger.info('Calling Job to Run')
                self.job_mod.run(self)
            except Exception:
                self.job_failed(logger, 'run()')
                logger.info('Sending comment text')
                self.send_comment_text()
        else:
            logger.info(f'Cannot find label {self.preq_dict["label"]}')

    def send_comment_text(self):
        logger = logging.getLogger('JOB/SEND_COMMENT_TEXT')
        logger.info(f'Comment Text: {self.comment_text}')
        self.comment_text_append('Please make changes and add '
                                 'the following label back:')
        self.comment_text_append(f'{self.machine}'
                                 f'-{self.compiler}'
                                 f'-{self.preq_dict["action"]}')

        self.preq_dict['preq'].create_issue_comment(self.comment_text)

    def job_failed(self, logger, job_name, exception=Exception, STDOUT=False,
                   out=None, err=None):
        logger.critical(f'{job_name} FAILED. Exception:{exception}')

        if STDOUT:
            logger.critical(f'STDOUT: {[item for item in out if not None]}')
            logger.critical(f'STDERR: {[eitem for eitem in err if not None]}')

def setup_env():
    hostname = os.getenv('HOSTNAME')
    if bool(re.match(re.compile('hfe.+'), hostname)):
        machine = 'hera'
    elif bool(re.match(re.compile('hecflow.+'), hostname)):
        machine = 'hera'
    elif bool(re.match(re.compile('fe.+'), hostname)):
        machine = 'jet'
        os.environ['ACCNR'] = 'h-nems'
    elif bool(re.match(re.compile('gaea.+'), hostname)):
        machine = 'gaea'
        os.environ['ACCNR'] = 'nggps_emc'
    elif bool(re.match(re.compile('Orion-login.+'), hostname)):
        machine = 'orion'
    elif bool(re.match(re.compile('chadmin.+'), hostname)):
        machine = 'derecho'
        os.environ['ACCNR'] = 'P48503002'
    else:
        raise KeyError(f'Hostname: {hostname} does not match '\
                        'for a supported system. Exiting.')

    # Dictionary of GitHub repositories to check
    repo_dict = [{
        'name': 'ufs-weather-model',
        'address': 'ufs-community/ufs-weather-model',
        'base': 'develop'
    }]

    # Approved Actions
    action_list = ['RT', 'BL']

    return machine, repo_dict, action_list


def main():

    # handle logging
    log_filename = f'rt_auto_'\
                   f'{datetime.datetime.now().strftime("%Y%m%d%H%M%S")}.log'
    logging.basicConfig(filename=log_filename, filemode='w',
                        level=logging.INFO)
    logger = logging.getLogger('MAIN')
    logger.info('Starting Script')

    # setup environment
    logger.info('Getting the environment setup')
    machine, repos, actions = setup_env()

    # setup interface with GitHub
    logger.info('Setting up GitHub interface.')
    ghinterface_obj = GHInterface()

    # get all pull requests from the GitHub object
    # and turn them into Job objects
    logger.info('Getting all pull requests, '
                'labels and actions applicable to this machine.')
    jobs = get_preqs_with_actions(repos, machine,
                                       ghinterface_obj, actions)
    [job.run() for job in jobs]

    logger.info('Script Finished')


if __name__ == '__main__':
    main()
