# HEADER
# Written by Brian Curtis
# Automation of RT tasks using github cli

from github import Github as gh
import datetime
import pandas as pd
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

    def __init__(self, rtdata_obj):
        self.rtdata_obj = rtdata_obj
        self.get_machine_info()
        self.machineid = os.environ.get('MACHINE_ID')

    def get_machine_info(self):
        hostname = socket.gethostname()
        machines = self.rtdata_obj.get_yaml_subset('machines')
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

    def __init__(self, name, address, base, machine_obj, ghinterface_obj):
        self.name = name
        self.address = address
        self.base = base
        self.machine_obj = machine_obj
        self.ghinterface_obj = ghinterface_obj
        self.get_GHrepo_object()
        self.get_repo_preqs()

    def get_GHrepo_object(self):
        try:
            self.ghrepo = self.ghinterface_obj.client.get_repo(self.address)
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
            self.pullreq_list.append(PullReq(self, preq, self.machine_obj))

class PullReq:

    def __init__(self, repo_obj, preq_obj, machine_obj):
        self.preq_obj = preq_obj
        self.repo_obj = repo_obj
        self.machine_obj = machine_obj

        self.branch = self.preq_obj.head.ref
        self.repo_name = self.preq_obj.head.repo.name
        self.git_url = self.preq_obj.head.repo.html_url
        self.repo_dir_str = self.machine_obj.workdir+'/'+str(preq_obj.id)+'/'+datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        self.get_pr_labels()

    def get_pr_labels(self):
        self.labels = []
        pr_labels = self.preq_obj.get_labels()
        for pr_label in pr_labels:
            split_pr_label = pr_label.name.split('-')
            if len(split_pr_label) == 3 and split_pr_label[0] == "Auto":
                self.labels.append(PRLabel(split_pr_label[1], split_pr_label[2]))
            else:
                continue

    def add_clone_dir(self, clone_dir):
        self.clone_dir = clone_dir

class PRLabel:

    def __init__(self, name, machine):
        self.name = name
        self.machine = machine
        self.function_obj = None

    def add_function(self, function_obj):
        self.function_obj = function_obj

    def is_approved(self, machine_obj, functions):
        if re.match(self.machine.lower(), machine_obj.name.lower()):
            for function_obj in functions:
                if function_obj.verify_name(self.name):
                    print(f'Approved label "{self.name}-{self.machine}"')
                    self.add_function(function_obj)
                    return True
                else:
                    print(f'Label not approved. Name: {self.name} not on approved list.')
                    return False
        else:
            print(f'Label not approved. Machine: {self.machine} not on approved list.')
            return False

def get_approved_functions(rtdata_obj):
    function_data = rtdata_obj.get_yaml_subset('functions').\
        values.tolist()
    function_list = [Function(name,command,callback_fnc)
        for name,command,callback_fnc in function_data]
    return function_list

def clone_pr_repo(pullreq_obj, ghinterface_obj, machine_obj):

    git_url_w_login = pullreq_obj.git_url.split('//')
    git_url_w_login = git_url_w_login[0]+"//"+ghinterface_obj.GHUSERNAME+":"+ghinterface_obj.GHACCESSTOKEN+"@"+git_url_w_login[1]

    create_repo_commands = [
        ['mkdir -p "'+pullreq_obj.repo_dir_str+'"', machine_obj.workdir],
        ['git clone -b '+pullreq_obj.branch+' '+git_url_w_login, pullreq_obj.repo_dir_str],
        ['git submodule update --init --recursive', pullreq_obj.repo_dir_str+'/'+pullreq_obj.repo_name],
        ['module use modulefiles/{}.intel && module load fv3'.format(machine_obj.name.lower()), pullreq_obj.repo_dir_str+'/'+pullreq_obj.repo_name]
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

    pullreq_obj.add_clone_dir(pullreq_obj.repo_dir_str+"/"+pullreq_obj.repo_name)

def process_pr(pullreq_obj, ghinterface_obj, machine_obj, functions):
    for prlabel in pullreq_obj.labels:
        if prlabel.is_approved(machine_obj, functions):
            try:
                clone_pr_repo(pullreq_obj, ghinterface_obj, machine_obj)
            except Exception as e:
                sys.exit(f'Error cloning repo: {pullreq_obj.git_url}')
            try:
                goodexit = runInThread(prlabel.function_obj, pullreq_obj)
                # thread, badexit = send_to_thread(prlabel.function_obj, pullreq_obj)
            except Exception as e:
                print(f'ERROR RUNNING RT {prlabel.function_obj.command} with error: {e}')
                continue
            else:
                if goodexit == True:
                    pullreq_obj.preq_obj.remove_from_labels(f'Auto-{prlabel.name}-{prlabel.machine}')

def move_rt_logs(pullreq_obj):
    rt_log = 'tests/RegressionTests_'+pullreq_obj.machine_obj.machineid+'.log'
    filepath = pullreq_obj.clone_dir+'/'+rt_log
    if os.path.exists(filepath):

        move_rt_commands = [
            ['git add '+rt_log, pullreq_obj.clone_dir],
            ['git commit -m "Auto: Added Updated RT Log file: '+rt_log+'"', pullreq_obj.clone_dir],
            ['git pull --commit origin '+pullreq_obj.branch, pullreq_obj.clone_dir],
            ['git push origin '+pullreq_obj.branch, pullreq_obj.clone_dir],
            ['rm -rf '+pullreq_obj.clone_dir, pullreq_obj.clone_dir]
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

def runInThread(function_obj, pullreq_obj):
    goodexit = None
    proc = subprocess.Popen(function_obj.command, shell=True, cwd=pullreq_obj.clone_dir)
    proc.wait()
    proc.poll()
    if proc.returncode == 0:
        print(f'Process successful, running callback function')
        globals()[function_obj.callback](pullreq_obj)
        goodexit = True
    else:
        print(f'Process failed, skipping callback function')
        goodexit = False
    return goodexit

def get_approved_repos(rtdata_obj, machine_obj, ghinterface_obj):
    repo_data = rtdata_obj.get_yaml_subset('repository').values.tolist()
    repo_list = [Repo(name, address, base, machine_obj, ghinterface_obj)
        for name, address, base in repo_data]
    if not isinstance(repo_list, list):
        repo_list = list(repo_list)
    return repo_list

def delete_old_pullreq(repo_list, machine_obj):
    # Get all PR ID Nums
    pr_id_list = [str(single_pr.preq_obj.id) for single_repo in repo_list for single_pr in single_repo.pullreq_list]
    dir_list = next(os.walk(machine_obj.workdir))[1]
    not_in_both = [item for item in dir_list if not item in pr_id_list]

    for not_in in not_in_both:
        if os.path.isdir(machine_obj.workdir+'/'+not_in):
            command = ['rm -rf '+machine_obj.workdir+'/'+not_in]
            print(f'Attempting to run: {command}')
            try:
                retcode = subprocess.Popen(command, shell=True)
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
                if os.path.isdir(machine_obj.workdir+'/'+not_in):
                    print(f'WARNING: Successful command but dir was not removed')
                else:
                    print(f'Command {command} ran successfully')
        else:
            print(f'Somehow directory called for removal does not exist, skipping..')

def main():
    GHUSERNAME = 'BrianCurtis-NOAA'
    ghinterface_obj = GHInterface(GHUSERNAME)
    rtdata_obj= RTData()
    machine_obj = Machine(rtdata_obj)
    functions_obj = get_approved_functions(rtdata_obj)
    repo_list = get_approved_repos(rtdata_obj, machine_obj, ghinterface_obj)
    delete_old_pullreq(repo_list, machine_obj)
    for single_repo in repo_list:
        print(f'Processing {single_repo.address} repository')
        for single_pr in single_repo.pullreq_list:
            print(f'Processing pull request #{single_pr.preq_obj.id}')
            process_pr(single_pr, ghinterface_obj, machine_obj, functions_obj)
            print(f'Finished processing pull request #{single_pr.preq_obj.id}')

if __name__ == '__main__':
    main()
