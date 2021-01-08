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

class machine_info():

    def __init__(self, name, regexhostname, workdir, baselinedir):
        self.name = name
        self.regexhostname = regexhostname
        self.workdir = workdir
        self.baselinedir = baselinedir

def get_access_token():
    # CREATE FILE "accesstoken.txt" add API Token and change perms appropriately
    # GITHUB RATE LIMIT FOR AUTHENTICATED USERS IS 5000 REQUESTS PER HOUR
    f = open("accesstoken.txt", 'rb')
    GHACCESSTOKEN = f.read()
    f.close()
    GHACCESSTOKEN = GHACCESSTOKEN.decode("utf-8").strip('\n')

    return GHACCESSTOKEN

def read_yaml_data(filename):
    stream = open(db_filename, 'r')
    yaml_data = yaml.load(stream, Loader=yaml.SafeLoader)
    stream.close()
    return yaml_data

def get_yaml_subset(in_yaml, keyin):
    try:
        yaml_subset = in_yaml[keyin]
    except:
        sys.exit("Unable to get yaml subset: {}. Quitting".format(keyin))
    df = pd.DataFrame.from_dict(yaml_subset)
    return df

def get_machine_info(static_data):
    hostname = socket.gethostname()
    machines = get_yaml_subset(static_data, "machines")
    A = machines['name'][machines.index[machines.regexhostname.str.match(hostname)].tolist()]
    if len(A) == 1:
        app_machine_name = str(list(A)[0])
        print("Approved Machine Found: {}".format(app_machine_name))
    else:
        sys.exit("Hostname {} does not match approved list. Quitting".format(hostname))

    mach_db_info = machines.loc[machines[machines['name'] == app_machine_name].index.values[0]]
    machine = machine_info(mach_db_info['name'], mach_db_info['regexhostname'],
                           mach_db_info['workdir'], mach_db_info['baselinedir'])
    return machine

def connect_to_github(GHACCESSTOKEN): # In case there's more needed later
    client = gh(GHACCESSTOKEN) #will work if None for public repos.
    return client

def process_pulls(pulls, repo):
    valid_actions = get_yaml_subset(static_data, "actions")
    for pr in pulls:
        labels = pr.get_labels()
        for label in labels:
            if label.name.split('-')[1].lower() == machine.name.lower():
                for action_name, action_command in valid_actions.values.tolist():
                    if label.name.split('-')[0].lower() == action_name.lower():
                        try:
                            pr_workdir = clone_pr_repo(pr)
                            pa_ret = process_actions(action_command, pr_workdir, pr)
                            if pa_ret == 0:
                                pr.remove_from_labels(label.name.split('-')[0]+"-"+label.name.split('-')[1])

                        except:
                            print("ERROR RUNNING RT {}".format(action_name))
                            continue

def clone_pr_repo(pr):
    branch = pr.head.ref
    repo_name = pr.head.repo.name
    git_url = pr.head.repo.html_url
    repo_dir_str = machine.workdir+"/"+str(pr.id)+"/"+datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    git_url_w_login = git_url.split('//')
    git_url_w_login = git_url_w_login[0]+"//"+GHUSERNAME+":"+GHACCESSTOKEN+"@"+git_url_w_login[1]

    create_repo_commands = [
        ["mkdir -p \""+repo_dir_str+"\"", machine.workdir],
        ["git clone -b "+branch+" "+git_url_w_login, repo_dir_str],
        ["git submodule update --init --recursive", repo_dir_str+"/"+repo_name]
    ]

    for command, in_cwd in create_repo_commands:
        print("Attempting to run: {}".format(command))
        try:
            retcode = subprocess.Popen(command, shell=True, cwd=in_cwd)
            retcode.wait()
            if retcode.returncode==1:
                print("Error Occured:")
                print("Stdout: {}".format(retcode.stdout))
                print("Stderr: {}".format(retcode.stderr))
                sys.exit()
        except OSError as e:
            print("Execution failed: {}".format(e))
            sys.exit()

    return repo_dir_str+"/"+repo_name

def create_threaded_call(callback, run_fnc):
    def runInThread(callback, run_fnc):
        proc = run_fnc
        proc.wait()
        callback()
        return

    thread = threading.Thread(target=runInThread,
                              args=(callback, run_fnc))
    thread.start()

    return thread # returns immediately after the thread starts

def move_rt_logs(pr, pr_workdir):
    rt_log = 'tests/RegressionTests_hera.intel.log'
    # rt_log = 'RegressionTests_'+machine.name+'.intel.log'
    filepath = pr_workdir+'/'+rt_log
    print("File path is {}".format(filepath))
    if os.path.exists(filepath):
        branch = pr.head.ref
        print("Branch used is: {}".format(branch))
        repo = pr.head.repo

        move_rt_commands = [
            ['echo "GOOSEBUMPS2" >> '+rt_log, pr_workdir],
            ['git add '+rt_log, pr_workdir],
            ['git commit -m "Auto: Added Updated RT Log file: '+rt_log+'"', pr_workdir],
            ['git push origin '+branch, pr_workdir]
        ]
        for command, in_cwd in move_rt_commands:
            print("Attempting to run: {}".format(command))
            try:
                retcode = subprocess.Popen(command, shell=True, cwd=in_cwd)
                retcode.wait()
                if retcode.returncode==1:
                    print("Error Occured:")
                    print("Stdout: {}".format(retcode.stdout))
                    print("Stderr: {}".format(retcode.stderr))
                    sys.exit()
            except OSError as e:
                print("Execution failed: {}".format(e))
                sys.exit()
    else:
        print("ERROR: Could not find RT log")
        sys.exit()

def process_actions(command, pr_workdir, pr):
        try:
            print("{} attempting command {}".format(machine.name.upper(), command))
            thethread = create_threaded_call(move_rt_logs(pr, pr_workdir), subprocess.Popen(command, shell=True, cwd=pr_workdir+"/tests"))
        except:
            print("{} failed command {}".format(machine.name.upper(), command))
            return 1
        return 0
        # pr.create_issue_comment("AUTOMATED COMMENT: Submitted Command '{}' in directory '{}' on machine '{}'".format(command, pr_workdir, machine.name))

# START OF MAIN
db_filename = 'rt_auto.yml'
machine = get_machine_info(read_yaml_data(db_filename))
GHUSERNAME = "BrianCurtis-NOAA"
GHACCESSTOKEN = get_access_token()

# Initial items
static_data = read_yaml_data('rt_auto.yml')
repos = get_yaml_subset(static_data, "repository")
client = connect_to_github(GHACCESSTOKEN)

# Maybe a process_repo(repos) function start?
for name,address,base in repos.values.tolist():
    repo = client.get_repo(address)
    pull_reqs = repo.get_pulls(state='open', sort='created', base=base)
    triggers = process_pulls(pull_reqs, repo)
