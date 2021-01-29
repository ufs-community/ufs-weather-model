# HEADER
# Written by Brian Curtis
# Automation of RT tasks using github cli
from __future__ import print_function

from github import Github as gh
import datetime
import time
import pandas as pd
import socket
import threading
import subprocess
import re
import sys
import yaml
import os
import logging

class RTData:

    def __init__(self):
        self.logger = logging.getLogger("RTData")
        self.DB_FILENAME = 'rt_auto.yml'
        self.data = self.read_yaml_data()

    def read_yaml_data(self):
        self.logger.info(f'Reading in default YAML data')
        stream = open(self.DB_FILENAME, 'r')
        yaml_data = yaml.load(stream, Loader=yaml.SafeLoader)
        stream.close()
        self.logger.info(f'Finished reading in YAML data')
        return yaml_data

    def get_yaml_subset(self, keyin):
        self.logger.info(f'Getting YAML subset: {keyin}')
        try:
            yaml_subset = self.data[keyin]
        except:
            self.logger.critical(f'Unable to get yaml subset: {keyin}. Quitting')
            sys.exit()
        df = pd.DataFrame.from_dict(yaml_subset)
        self.logger.info(f'Finished getting YAML subset {keyin}')
        return df

class Machine:

    def __init__(self, rtdata_obj):
        self.logger = logging.getLogger("Machine")
        self.rtdata_obj = rtdata_obj
        self.get_machine_info()
        self.machineid = os.environ.get('MACHINE_ID')

    def get_machine_info(self):
        self.logger.info(f'Getting machine information')
        hostname = socket.gethostname()
        self.logger.debug(f'hostname is {hostname}')
        machines = self.rtdata_obj.get_yaml_subset('machines')
        regexmachines = machines['regexhostname'].values
        for i, mach in enumerate(regexmachines):
            if re.search(mach, hostname):
                self.logger.info(f'{hostname} is an approved machine')
                mach_db_info = machines.iloc[i]
                self.name = mach_db_info['name']
                self.regexhostname = mach_db_info['regexhostname']
                self.workdir = mach_db_info['workdir']
                self.baselinedir = mach_db_info['baselinedir']
                self.logger.debug(f'self.name: {self.name}\n\
                    self.regexhostname: {self.regexhostname}\n\
                    self.workdir: {self.workdir}\n\
                    self.baselinedir: {self.baselinedir}')
                self.logger.info(f'Finished getting machine information')
                return
            else:
                continue
        self.logger.critical(f'Hostname:{hostname} does not match approved list. Quitting')
        sys.exit()

class Function:

    def __init__(self, name, command, callback):
        self.logger = logging.getLogger("Function")
        self.name = name
        self.command = command
        self.callback = callback
        self.logger.debug(f'Function object created')
        self.logger.debug(f'self.name: {self.name}\n\
            self.command: {self.command}\n\
            self.callback: {self.callback}')

    def verify_name(self, comparable):
        self.logger.debug(f'Verifying name: {comparable}')
        if re.match(self.name.lower(), comparable.lower()):
            self.logger.debug(f'Match found')
            return True
        else:
            self.logger.debug(f'Not a match')
            return False

    def verify_command(self, comparable):
        self.logger.debug(f'Verifying command: {comparable}')
        if re.match(self.command.lower(), comparable.lower()):
            self.logger.debug(f'Match found')
            return True
        else:
            self.logger.debug(f'Not a match')
            return False

class GHInterface:

    def __init__(self):
        self.logger = logging.getLogger("GHInterface")
        self.get_access_token()
        self.client = gh(self.GHACCESSTOKEN)

    def get_access_token(self):
        if os.path.exists('accesstoken.txt'):
            f = open('accesstoken.txt', 'rb')
            filedata = f.read()
            f.close()
            filedata = str(filedata)[2:-1].split('\\n')

            self.GHACCESSTOKEN = filedata[0]
            self.GHUSERNAME = filedata[1]
            self.logger.debug(f'GHACCESSTOKEN is {self.GHACCESSTOKEN}')
            self.logger.debug(f'GHUSERNAME is {self.GHUSERNAME}')
        else:
            self.logger.critical(f'Please create a file "accesstoken.txt" that'\
                ' contains your GitHub API Token on the first line, and GitHub'\
                ' Username on the second line.\nMake sure to set permission'\
                ' so others can not read it (400)')
            sys.exit()

# REPO STUFF
class Repo:

    def __init__(self, name, address, base, machine_obj, ghinterface_obj):
        self.logger  = logging.getLogger("Repo")
        self.logger.debug(f'Creating a Repo object')
        self.name = name
        self.address = address
        self.base = base
        self.logger.debug(f'self.name: {self.name}\n\
            self.address: {self.address}\n\
            self.base {self.base}')
        self.machine_obj = machine_obj
        self.ghinterface_obj = ghinterface_obj
        self.get_GHrepo_object()
        self.get_repo_preqs()

    def get_GHrepo_object(self):
        self.logger.debug(f'Getting list of GitHub Repo objects')
        try:
            self.ghrepo = self.ghinterface_obj.client.get_repo(self.address)
        except Exception as e:
            self.logger.critical(f'Failed to get repo object with error {e}')
            sys.exit()

    def get_repo_preqs(self):
        self.logger.debug(f'Creating list of PullReq objects')
        self.pullreq_list = []
        try:
            preqs = self.ghrepo.get_pulls(state='open', sort='created', base=self.base)
        except Exception as e:
            self.logger.critical(f'Failed to get pull object with error {e}')
            sys.exit()

        for preq in preqs:
            self.pullreq_list.append(PullReq(self, preq, self.machine_obj))
        self.logger.debug(f'Finished getting list of GutHub Repo objects')

class PullReq:

    def __init__(self, repo_obj, preq_obj, machine_obj):
        self.logger = logging.getLogger("PullReq")
        self.preq_obj = preq_obj
        self.repo_obj = repo_obj
        self.machine_obj = machine_obj

        self.branch = self.preq_obj.head.ref
        self.repo_name = self.preq_obj.head.repo.name
        self.git_url = self.preq_obj.head.repo.html_url
        self.repo_dir_str = self.machine_obj.workdir+'/'+str(preq_obj.id)+'/'+datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        self.logger.debug(f'self.branch: {self.branch}\n\
            self.repo_name: {self.repo_name}\n\
            self.git_url: {self.git_url}\n\
            self.repo_dir_str: {self.repo_dir_str}')
        self.get_pr_labels()

    def get_pr_labels(self):
        self.logger.debug(f'Creating list of PRLabel objects')
        self.labels = []
        pr_labels = self.preq_obj.get_labels()
        for pr_label in pr_labels:
            split_pr_label = pr_label.name.split('-')
            if len(split_pr_label) == 3 and split_pr_label[0] == "Auto":
                self.labels.append(PRLabel(split_pr_label[1], split_pr_label[2]))
            else:
                continue
        self.logger.debug(f'Finished creating list of PRLabel objects')

    def add_clone_dir(self, clone_dir):
        self.logger.debug(f'Adding Clone Dir: {clone_dir}')
        self.clone_dir = clone_dir

class PRLabel:

    def __init__(self, name, machine):
        self.logger = logging.getLogger("PRLabel")
        self.name = name
        self.machine = machine
        self.logger.debug(f'self.name: {self.name}\n\
            self.machine: {self.machine} ')
        self.function_obj = None

    def add_function(self, function_obj):
        self.logger.debug("Adding function Object: {function_obj}")
        self.function_obj = function_obj

    def is_approved(self, machine_obj, functions):
        self.logger.debug("Checking if PRLabel is approved")
        if re.match(self.machine.lower(), machine_obj.name.lower()):
            for function_obj in functions:
                if function_obj.verify_name(self.name):
                    self.logger.debug(f'Approved label: "{self.name}-{self.machine}"')
                    self.add_function(function_obj)
                    return True
                else:
                    self.logger.warning(f'Label not approved. Name: {self.name} not on approved list.')
                    return False
        else:
            self.logger.warning(f'Label not approved. Machine: {self.machine} not on approved list.')
            return False

def get_approved_functions(rtdata_obj):
    logger = logging.getLogger("get_approved_functions()")
    logger.debug(f'Start of get_approved_functions()')
    function_data = rtdata_obj.get_yaml_subset('functions').\
        values.tolist()
    logger.debug("Sending function data off to Function objects")
    function_list = [Function(name,command,callback_fnc)
        for name,command,callback_fnc in function_data]
    logger.debug(f'End of get_approved_functions()')
    return function_list

def clone_pr_repo(pullreq_obj, ghinterface_obj, machine_obj):
    logger = logging.getLogger("clone_pr_repo()")
    git_url_w_login = pullreq_obj.git_url.split('//')
    git_url_w_login = git_url_w_login[0]+"//"+ghinterface_obj.GHUSERNAME+":"+ghinterface_obj.GHACCESSTOKEN+"@"+git_url_w_login[1]

    create_repo_commands = [
        ['mkdir -p "'+pullreq_obj.repo_dir_str+'"', machine_obj.workdir],
        ['git clone -b '+pullreq_obj.branch+' '+git_url_w_login, pullreq_obj.repo_dir_str],
        ['git submodule update --init --recursive', pullreq_obj.repo_dir_str+'/'+pullreq_obj.repo_name]
        # ['module use modulefiles/{} && module load fv3'.format(machine_obj.machineid), pullreq_obj.repo_dir_str+'/'+pullreq_obj.repo_name]
    ]

    for command, in_cwd in create_repo_commands:
        logger.debug(f'Attempting to run: {command}'.format(command))
        try:
            retcode = subprocess.Popen(command, shell=True, cwd=in_cwd)
            retcode.wait()
            if retcode.returncode==1:
                logger.critical('Error Occured:')
                logger.critical(f'Stdout: {retcode.stdout}')
                logger.critical(f'Stderr: {retcode.stderr}')
                sys.exit()
        except OSError as e:
            logger.critical(f'Execution failed with error: {e}')
            sys.exit()
        else:
            logger.debug(f'Finished Cloning {git_url_w_login}')

    pullreq_obj.add_clone_dir(pullreq_obj.repo_dir_str+"/"+pullreq_obj.repo_name)

def process_pr(pullreq_obj, ghinterface_obj, machine_obj, functions):
    logger = logging.getLogger(f'PR#{pullreq_obj.preq_obj.id}')
    logger.debug(f'Start')
    for prlabel in pullreq_obj.labels:
        if prlabel.is_approved(machine_obj, functions):
            logger.debug(f'Removing Label Auto-{prlabel.name}-{prlabel.machine}')
            pullreq_obj.preq_obj.remove_from_labels(f'Auto-{prlabel.name}-{prlabel.machine}')
            try:
                clone_pr_repo(pullreq_obj, ghinterface_obj, machine_obj)
            except Exception as e:
                logger.critical(f'Error cloning repo: {pullreq_obj.git_url}')
                sys.exit()
            try:
                goodexit = runFunction(prlabel.function_obj, pullreq_obj)
                # thread, badexit = send_to_thread(prlabel.function_obj, pullreq_obj)
            except Exception as e:
                logger.critical(f'ERROR RUNNING RT {prlabel.function_obj.command} with error: {e}')
                logger.debug(f'Adding back label Auto-{prlabel.name}-{prlabel.machine}')
                pullreq_obj.preq_obj.add_to_labels(f'Auto-{prlabel.name}-{prlabel.machine}')
                sys.exit()
            else:
                if goodexit == False:
                    logger.debug(f'runFunction exited with error')
                    pullreq_obj.preq_obj.add_to_labels(f'Auto-{prlabel.name}-{prlabel.machine}')
                    logger.debug(f'Adding back label Auto-{prlabel.name}-{prlabel.machine}')
    logger.debug(f'End')

def move_rt_logs(pullreq_obj):
    logger = logging.getLogger("move_rt_logs()")
    rt_log = 'tests/RegressionTests_'+pullreq_obj.machine_obj.machineid+'.log'
    logger.debug(f'rt_log is: {rt_log}')
    filepath = pullreq_obj.clone_dir+'/'+rt_log
    logger.debug(f'filepath is: {filepath}')
    rm_filepath = '/'.join((pullreq_obj.clone_dir.split('/'))[:-1])
    logger.debug(f'rm_filepath is: {rm_filepath}')
    if os.path.exists(filepath):
        move_rt_commands = [
            ['git add '+rt_log, pullreq_obj.clone_dir],
            ['git commit -m "Auto: Added Updated RT Log file: '+rt_log+'"', pullreq_obj.clone_dir],
            ['git pull --no-edit origin '+pullreq_obj.branch, pullreq_obj.clone_dir],
            ['sleep 10', pullreq_obj.clone_dir],
            ['git push origin '+pullreq_obj.branch, pullreq_obj.clone_dir]
        ]
        for command, in_cwd in move_rt_commands:
            logger.debug(f'Attempting to run: {command}')
            try:
                retcode = subprocess.Popen(command, shell=True, cwd=in_cwd)
                retcode.wait()
                if retcode.returncode==1:
                    logger.critical('Error Occured:')
                    logger.critical(f'Stdout: {retcode.stdout}')
                    logger.critical(f'Stderr: {retcode.stderr}')
                    sys.exit()
            except OSError as e:
                logger.critical(f'Execution failed with error: {e}')
                sys.exit()
            else:
                logger.debug(f'Funished running command: {command}')
    else:
        logger.critical('ERROR: Could not find RT log')
        sys.exit()

# def send_to_thread(function_obj, pullreq_obj):
#     badexit = False
#     def create_threaded_call(function_obj, pullreq_obj, badexit):
#
#         def runInThread(incallback, incommand, cwd_in, badexit):
#             print(f'cwd_in is:{cwd_in}')
#             proc = subprocess.Popen(incommand, shell=True, cwd=cwd_in)
#             proc.wait()
#             proc.poll()
#             if proc.returncode == 0:
#                 print(f'Process successful, running callback function')
#                 globals()[incallback](pullreq_obj)
#             else:
#                 print(f'badexit is {badexit}')
#                 print(f'Process failed, skipping callback function')
#                 badexit = True
#                 print(f'new badexit is {badexit}')
#             return
#
#         thread = threading.Thread(target=runInThread,
#                                   args=(function_obj.callback, function_obj.command, pullreq_obj.clone_dir, badexit))
#         thread.start()
#         print(f'ouside badexit is {badexit}')
#         return thread, badexit # returns immediately after the thread starts
#     thread, exit = create_threaded_call(function_obj, pullreq_obj, badexit)
#     return thread, exit

def runFunction(function_obj, pullreq_obj):
    logger = logging.getLogger("runFunction()")
    goodexit = None
    logger.debug(f'Running command: {function_obj.command}')
    proc = subprocess.Popen(function_obj.command, shell=True, cwd=pullreq_obj.clone_dir)
    proc.wait()
    proc.poll()
    if proc.returncode == 0:
        logger.debug(f'Command {function_obj.command} successful')
        try:
            logger.debug(f'Running callback function: {function_obj.callback}')
            globals()[function_obj.callback](pullreq_obj)
        except:
            logger.critical(f'Callback function {function_obj.callback} failed,\
                skipping callback function')
            goodexit = False
        else:
            logger.debug(f'Callback function {function_obj.callback} ran successfully')
            goodexit = True
    return goodexit

def get_approved_repos(rtdata_obj, machine_obj, ghinterface_obj):
    logger = logging.getLogger(f'get_approved_repos')
    logger.debug(f'Start')
    repo_data = rtdata_obj.get_yaml_subset('repository').values.tolist()
    logger.debug(f'Creating list of repository data Repo objects')
    repo_list = [Repo(name, address, base, machine_obj, ghinterface_obj)
        for name, address, base in repo_data]
    if not isinstance(repo_list, list):
        logger.warning(f'repo_list not list type')
        repo_list = list(repo_list)
    return repo_list

def delete_old_pullreq(repo_list, machine_obj):
    logger = logging.getLogger(f'delete_old_pullreq()')
    # Get all PR ID Nums
    pr_id_list = [str(single_pr.preq_obj.id) for single_repo in repo_list for single_pr in single_repo.pullreq_list]
    dir_list = next(os.walk(machine_obj.workdir))[1]
    not_in_both = [item for item in dir_list if not item in pr_id_list]

    for not_in in not_in_both:
        if os.path.isdir(machine_obj.workdir+'/'+not_in):
            command = 'rm -rf '+machine_obj.workdir+'/'+not_in
            logger.info(f'Attempting to run: {command}')
            try:
                retcode = subprocess.Popen(command, shell=True)
                retcode.wait()
                if retcode.returncode==1:
                    logger.critital('Error Occured:')
                    logger.critical(f'Stdout: {retcode.stdout}')
                    logger.critical(f'Stderr: {retcode.stderr}')
                    sys.exit()
            except OSError as e:
                logger.critical(f'Execution failed with error: {e}')
                sys.exit()
            else:
                if os.path.isdir(machine_obj.workdir+'/'+not_in):
                    logger.warning(f'Successful command but dir was not removed')
                else:
                    logger.info(f'Command {command} ran successfully')
        else:
            logger.warning(f'Somehow directory called for removal does not exist, skipping..')

def wait_for_daemon_threads(thread_list):
    logger = logging.getLogger(f'wait_for_daemon_threads()')
    logger.debug(f'Starting to wait for all daemon threads to finish before'\
        'exiting main()')
    last_thread_completed = False
    while not last_thread_completed:
        if len(thread_list) > 0:
            logger.debug("Checking Threads")
            for thread in thread_list:
                logger.debug(f'Checking Thread {thread.name}')
                if thread.is_alive():
                    logger.debug(f'Thread {thread.name} Alive, continuing')
                else:
                    logger.debug(f'Thread {thread.name} is NOT alive, removing')
                    thread_list.remove(thread)
            logger.debug("Sleeping for 10 seconds")
            time.sleep(10)
        else:
            logger.debug("Thread List Empty, Exiting")
            last_thread_completed = True
    logger.debug("All deamon threads finished")

def main():
    #SETUP LOGGER
    logging.basicConfig(filename="rt_auto.log", filemode='w', level=logging.DEBUG)
    logger = logging.getLogger("main")
    ghinterface_obj = GHInterface()
    rtdata_obj= RTData()
    machine_obj = Machine(rtdata_obj)
    functions_obj = get_approved_functions(rtdata_obj)
    repo_list = get_approved_repos(rtdata_obj, machine_obj, ghinterface_obj)
    delete_old_pullreq(repo_list, machine_obj)
    thread_list = []
    for single_repo in repo_list:
        logger.info(f'Processing {single_repo.address} repository')
        for single_pr in single_repo.pullreq_list:
            logger.info(f'Sending PR #{single_pr.preq_obj.id} processing off to thread')
            thread = threading.Thread(name=f'PR#{single_pr.preq_obj.id}',
                target=process_pr, args=(single_pr, ghinterface_obj,
                machine_obj, functions_obj), daemon=True)
            thread.start()
            thread_list.append(thread)
            # process_pr(single_pr, ghinterface_obj, machine_obj, functions_obj)
            logger.info(f'Finished Sending PR #{single_pr.preq_obj.id} processing off to thread')
    wait_for_daemon_threads(thread_list)

if __name__ == '__main__':
    main()
