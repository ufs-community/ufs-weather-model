# HEADER
# Written by Brian Curtis
# Automation of RT tasks using github cli

from github import Github as gh
import datetime
import pandas as pd
import itertools
import socket
import threading
import subprocess
import re
import sys
import yaml
import os

class RTData:

    def __init__(self):
        self.DB_FILENAME = 'rt_auto.yml'
        self.data = self.read_yaml_data()

    def read_yaml_data(self):
        stream = open(self.DB_FILENAME, 'r')
        yaml_data = yaml.load(stream, Loader=yaml.SafeLoader)
        stream.close()
        return yaml_data

    def get_yaml_subset(self, keyin):
        try:
            yaml_subset = self.data[keyin]
        except:
            sys.exit(f'Unable to get yaml subset: {keyin}. Quitting')
        df = pd.DataFrame.from_dict(yaml_subset)
        return df

class Machine:

    def __init__(self, RTData):
        self.RTData = RTData

    def get_machine_info(self):
        hostname = socket.gethostname()
        machines = self.RTData.get_yaml_subset('machines')
        regexmachines = machines['regexhostname'].values
        for i, mach in enumerate(regexmachines):
            if re.search(mach, hostname):
                print(f'{hostname} is an approved machine')
                mach_db_info = machines.iloc[i]
                self.name = mach_db_info['name']
                self.regexhostname = mach_db_info['regexhostname']
                self.workdir = mach_db_info['workdir']
                self.baselinedir = mach_db_info['baselinedir']
                break
            else:
                continue
        sys.exit(f'Hostname:{hostname} does not match approved list. Quitting')

class Repos:

    def __init__(self, ghinterface, rtdata):
        self.client = ghinterface.client
        self.rtdata = rtdata
        self.repo_info()

    class RepoInfo:
        def __init__(self, name, address, base):
            self.name = name
            self.address = address
            self.base = base

    def repo_info(self):
        repo_data = self.rtdata.get_yaml_subset('repository').\
            values.tolist()
        repo_list = [self.RepoInfo(name,address,base)
            for name,address,base in repo_data]
        print(f'type repo_list is {type(repo_list)}')
        if not isinstance(repo_list, list):
            repo_list = list(repo_list)
        self.repo_list = repo_list

class ActionInfo:

    def __init__(self, rtdata):
        self.rtdata = rtdata
        self.action_info()

    class ActionData:
        def __init__(self, name, command, callback_fnc):
            self.name = name
            self.command = command
            self.callback_fnc = callback_fnc

    def action_info(self):
        action_data = self.rtdata.get_yaml_subset('actions').\
            values.tolist()
        self.action_list = [self.ActionData(name,command,callback_fnc)
            for name,command,callback_fnc in action_data]

    def is_action_name(self, action_name_cmp):
        for match_action in self.action_list:
            if re.match(action_name_cmp.lower(), match_action.name.lower()):
                return True
            else:
                return False

    def get_action(self, action_name_cmp):
        for match_action in self.action_list:
            if re.match(action_name_cmp.lower(), match_action.name.lower()):
                return match_action
            else:
                return None

class GHInterface:

    def __init__(self, GHUSERNAME):
        self.GHUSERNAME = GHUSERNAME
        self.get_access_token()
        self.client = gh(self.GHACCESSTOKEN)

    def get_access_token(self):
        if os.path.exists('accesstoken.txt'):
            f = open('accesstoken.txt', 'rb')
            self.GHACCESSTOKEN = f.read()
            f.close()
            self.GHACCESSTOKEN = self.GHACCESSTOKEN.decode('utf-8').strip('\n')
        else:
            sys.exit('Please create a file "accesstoken.txt" that contains your\
                GitHub API Token.\nMake sure to set permissions so others can\
                not read it (400)')

class ProcessRepo:

    def __init__(self, machine, ghinterface, actions, repo):
        self.machine = machine
        self.ghinterface = ghinterface
        self.actions = actions
        self.repo = repo
        self.get_repo_object()
        self.get_repo_preqs()

    def get_repo_object(self):
        try:
            self.ghrepo = self.ghinterface.client.get_repo(self.repo.address)
        except Exception as e:
            print(f'Failed to get repo object with error {e}')
            self.ghrepo = None

    def get_repo_preqs(self):
        self.preq_list = []
        try:
            preqs = self.ghrepo.get_pulls(state='open', sort='created', base=self.repo.base)
        except Exception as e:
            print(f'Failed to get pull object with error {e}')
            preqs = None

        for preq in preqs:
            pullreq = PullReq(preq)
            sys.exit()
            pullreq.process()
            # self.preq_list.append(self.PullReq(preq))

class PullReq(ProcessRepo):

    def __init__(self, preq):
        super(PullReq, self).__init__()
        print(f'MACHINE NAMEMEMEME {super().machine.name}')
        self.preq = preq
        self.branch = self.preq.head.ref
        self.repo_name = self.preq.head.repo.name
        self.git_url = self.preq.head.repo.html_url
        self.repo_dir_str = super().machine.workdir+'/'+str(pr.id)+'/'+datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        self.get_labels()



    def get_labels(self):
        self.prlabels = []
        pr_labels = super().preq.get_labels()
        for pr_label in pr_labels:
            split_pr_label = pr_label.name.split('-')
            self.prlabels.append(PRLabel(split_pr_label[0], split_pr_label[1]))

    def clone_pr_repo(self):

        git_url_w_login = self.git_url.split('//')
        git_url_w_login = git_url_w_login[0]+"//"+super().ghinterface.GHUSERNAME+":"+super().ghinterface.GHACCESSTOKEN+"@"+git_url_w_login[1]

        create_repo_commands = [
            ['mkdir -p "'+self.repo_dir_str+'"', super().machine.workdir],
            ['git clone -b '+self.branch+' '+git_url_w_login, self.repo_dir_str],
            ['git submodule update --init --recursive', self.repo_dir_str+'/'+self.repo_name],
            ['module use modulefiles/{}.intel && module load fv3'.format(super().machine.name.lower()), self.repo_dir_str+'/'+self.repo_name]
        ]

        for command, in_cwd in create_repo_commands:
            print(f'Attempting to run: {command}'.format(command))
            try:
                retcode = subprocess.Popen(command, shell=True, cwd=in_cwd)
                retcode.wait()
                if retcode.returncode==1:
                    print('Error Occured:')
                    print(f'Stdout: {retcode.stdout}')
                    print(f'Stderr: {retcode.stderr}')
                    sys.exit()
            except OSError as e:
                print(f'Execution failed: {e}')
                sys.exit()

        self.clone_location = self.repo_dir_str+"/"+self.repo_name

    def process(self):
        for prlabel in self.prlabels:
            if prlabel.is_approved():
                self.approved_action = self.actions.get_action(prlabel.command)
            try:
                self.clone_pr_repo()
            except Exception as e:
                sys.exit(f'Error cloning repo: {self.git_url}')
            try:
                self.process_actions(self.approved_action.callback_fnc, self.approved_action.command, self.clone_location)
            except Exception as e:
                print(f'ERROR RUNNING RT {self.approved_action.command} with error: {e}')
                continue
            self.preq.remove_from_labels(f'{prlabel.command}-{super().machine.name}')

    def move_rt_logs(self):
        rt_log = 'tests/RegressionTests_'+super().machine.name+'.intel.log'
        filepath = self.clone_location+'/'+rt_log
        print(f'File path issssss {filepath}')
        if os.path.exists(filepath):
            print(f'Branch used is: {self.branch}')

            move_rt_commands = [
                ['git add '+rt_log, self.clone_location],
                ['git commit -m "Auto: Added Updated RT Log file: '+rt_log+'"', self.clone_location],
                ['git push origin '+branch, self.clone_location]
            ]
            for command, in_cwd in move_rt_commands:
                print(f'Attempting to run: {command}')
                try:
                    retcode = subprocess.Popen(command, shell=True, cwd=in_cwd)
                    retcode.wait()
                    if retcode.returncode==1:
                        print('Error Occured:')
                        print(f'Stdout: {retcode.stdout}')
                        print(f'Stderr: {retcode.stderr}')
                        sys.exit()
                except OSError as e:
                    print(f'Execution failed: {e}')
                    sys.exit()
        else:
            print('ERROR: Could not find RT log')
            sys.exit()

    def process_actions(callback_fnc, command):

        def create_threaded_call(callback_fnc, command, cwd_in):

            def runInThread(callback_fnc, command, cwd_in):
                proc = subprocess.Popen(command, shell=True, cwd=cwd_in)
                proc.wait()
                globals()[callback_fnc](self)
                return

            thread = threading.Thread(target=runInThread,
                                      args=(callback_fnc, command, cwd_in))
            thread.daemon=True
            thread.start()

            return thread # returns immediately after the thread starts

        print(f'{machine.name.upper()} is running command {command}')
        thread = create_threaded_call(callback_fnc, command, self.clone_location+'/tests')

class PRLabel(PullReq):
    def __init__(self, command, machine):
        self.lblcommand = command
        self.lblmachine = machine

    def is_approved(self):
        if re.match(machine.lower(), super().machine.name.lower()):
            if super().actions.is_action_name(command):
                print(f'Approved label "{command}-{machine}"')
                return True
            else:
                print(f'Label not approved. Command: {command} not on approved list.')
                return False
        else:
            print(f'Label not approved. Machine: {machine} does not match.')
            return False

def main():
    GHUSERNAME = 'BrianCurtis-NOAA'

    rtdata = RTData()
    machine = Machine(rtdata)
    actions = ActionInfo(rtdata)
    # repo_list = get_repo_info(rtdata)
    ghinterface = GHInterface(GHUSERNAME)
    repos = Repos(ghinterface, rtdata)
    print(f'type of repos is {type(repos)}')
    # repo_objs = repos.get_repo_objects()
    for single_repo in repos.repo_list:
        process = ProcessRepo(machine, ghinterface, actions, single_repo)

if __name__ == '__main__':
    main()
