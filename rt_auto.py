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
        self.get_machine_info()

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
                return
            else:
                continue
        sys.exit(f'Hostname:{hostname} does not match approved list. Quitting')

class Function:

    def __init__(self, name, command, callback):
        self.name = name
        self.command = command
        self.callback = callback

    def verify_name(self, comparable):
        if re.match(self.name.lower(), comparable.lower()):
            return True
        else:
            return False

    def verify_command(self, comparable):
        if re.match(self.command.lower(), comparable.lower()):
            return True
        else:
            return False

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
            sys.exit('Please create a file "accesstoken.txt" that contains your'\
                ' GitHub API Token.\nMake sure to set permissions so others can'\
                ' not read it (400)')

# REPO STUFF
class Repo:

    def __init__(self, name, address, base, machine, ghinterface):
        self.name = name
        self.address = address
        self.base = base
        self.machine = machine
        self.ghinterface = ghinterface
        self.get_GHrepo_object()
        self.get_repo_preqs()

    def get_GHrepo_object(self):
        try:
            self.ghrepo = self.ghinterface.client.get_repo(self.address)
        except Exception as e:
            print(f'Failed to get repo object with error {e}')
            self.ghrepo = None

    def get_repo_preqs(self):
        self.pullreq_list = []
        try:
            preqs = self.ghrepo.get_pulls(state='open', sort='created', base=self.base)
        except Exception as e:
            print(f'Failed to get pull object with error {e}')
            preqs = None

        for preq in preqs:
            self.pullreq_list.append(PullReq(self, preq, self.machine))

class PullReq:

    def __init__(self, Repo, preq, machine):
        self.preq = preq
        self.repo = Repo
        self.machine = machine

        self.branch = self.preq.head.ref
        self.repo_name = self.preq.head.repo.name
        self.git_url = self.preq.head.repo.html_url
        self.repo_dir_str = self.machine.workdir+'/'+str(preq.id)+'/'+datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        self.get_pr_labels()

    def get_pr_labels(self):
        self.labels = []
        pr_labels = self.preq.get_labels()
        for pr_label in pr_labels:
            split_pr_label = pr_label.name.split('-')
            self.labels.append(PRLabel(split_pr_label[0], split_pr_label[1]))

    def add_clone_dir(self, clone_dir):
        self.clone_dir = clone_dir

class PRLabel:

    def __init__(self, name, machine):
        self.name = name
        self.machine = machine
        self.function = None

    def add_function(self, function):
        self.function = function

    def is_approved(self, machine, functions):
        if re.match(self.machine.lower(), machine.name.lower()):
            for function in functions:
                if function.verify_name(self.name):
                    print(f'Approved label "{self.name}-{self.machine}"')
                    self.add_function(function)
                    return True
                else:
                    print(f'Label not approved. Name: {self.name} not on approved list.')
                    return False
        else:
            print(f'Label not approved. Machine: {self.machine} not on approved list.')
            return False

def get_approved_functions(rtdata):
    function_data = rtdata.get_yaml_subset('functions').\
        values.tolist()
    function_list = [Function(name,command,callback_fnc)
        for name,command,callback_fnc in function_data]
    return function_list

def clone_pr_repo(pullreq, ghinterface, machine):

    git_url_w_login = pullreq.git_url.split('//')
    git_url_w_login = git_url_w_login[0]+"//"+ghinterface.GHUSERNAME+":"+ghinterface.GHACCESSTOKEN+"@"+git_url_w_login[1]

    create_repo_commands = [
        ['mkdir -p "'+pullreq.repo_dir_str+'"', machine.workdir],
        ['git clone -b '+pullreq.branch+' '+git_url_w_login, pullreq.repo_dir_str],
        ['git submodule update --init --recursive', pullreq.repo_dir_str+'/'+pullreq.repo_name],
        ['module use modulefiles/{}.intel && module load fv3'.format(machine.name.lower()), pullreq.repo_dir_str+'/'+pullreq.repo_name]
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

    pullreq.add_clone_dir(pullreq.repo_dir_str+"/"+pullreq.repo_name)

def process_pr(pullreq, ghinterface, machine, functions):
    for prlabel in pullreq.labels:
        if prlabel.is_approved(machine, functions):
            try:
                clone_pr_repo(pullreq, ghinterface, machine)
            except Exception as e:
                sys.exit(f'Error cloning repo: {pullreq.git_url}')
            try:
                send_to_thread(prlabel.function, pullreq)
            except Exception as e:
                print(f'ERROR RUNNING RT {prlabel.function.command} with error: {e}')
                continue
            pullreq.preq.remove_from_labels(f'{prlabel.name}-{prlabel.machine}')

def move_rt_logs(PullReq):
    filepath = PullReq.clone_dir+'/'+rt_log
    if os.path.exists(filepath):

        move_rt_commands = [
            ['git add '+rt_log, PullReq.clone_dir],
            ['git commit -m "Auto: Added Updated RT Log file: '+rt_log+'"', PullReq.clone_dir],
            ['git push origin '+PullReq.branch, PullReq.clone_dir]
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

def send_to_thread(Function, PullReq):

    def create_threaded_call(Function, PullReq):

        def runInThread(incallback, incommand, cwd_in):
            print(f'cwd_in is:{cwd_in}')
            proc = subprocess.Popen(incommand, shell=True, cwd=cwd_in)
            proc.wait()
            globals()[incallback](PullReq)
            return

        thread = threading.Thread(target=runInThread,
                                  args=(Function.callback, Function.command, PullReq.clone_dir))
        # thread.daemon=True
        thread.start()

        return thread # returns immediately after the thread starts
    thread = create_threaded_call(Function, PullReq)

def get_approved_repos(rtdata, machine, ghinterface):
    repo_data = rtdata.get_yaml_subset('repository').values.tolist()
    repo_list = [Repo(name, address, base, machine, ghinterface)
        for name, address, base in repo_data]
    if not isinstance(repo_list, list):
        repo_list = list(repo_list)
    return repo_list

def main():
    GHUSERNAME = 'BrianCurtis-NOAA'
    ghinterface = GHInterface(GHUSERNAME)
    rtdata = RTData()
    machine = Machine(rtdata)
    functions = get_approved_functions(rtdata)
    repo_list = get_approved_repos(rtdata, machine, ghinterface)
    for single_repo in repo_list:
        for single_pr in single_repo.pullreq_list:
            process_pr(single_pr, ghinterface, machine, functions)
        # process = Process(machine, ghinterface, actions, single_repo)

if __name__ == '__main__':
    main()
