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
    regexmachines = machines['regexhostname'].values
    for i, mach in enumerate(regexmachines):
        if re.search(mach, hostname):
            print("{} is an approved machine".format(hostname))
            mach_db_info = machines.iloc[i]
            machine = machine_info(mach_db_info['name'], mach_db_info['regexhostname'],
                                   mach_db_info['workdir'], mach_db_info['baselinedir'])
            return machine
        else:
            continue
    sys.exit("Hostname {} does not match approved list. Quitting".format(hostname))

def connect_to_github(GHACCESSTOKEN): # In case there's more needed later
    client = gh(GHACCESSTOKEN) #will work if None for public repos.
    return client

def process_pulls(pulls, repo):
    valid_actions = get_yaml_subset(static_data, "actions")
    for pr in pulls:
        labels = pr.get_labels()
        for label in labels:
            if label.name.split('-')[1].lower() == machine.name.lower():
                for action_name, action_command, action_callback in valid_actions.values.tolist():
                    if label.name.split('-')[0].lower() == action_name.lower():
                        try:
                            pr_workdir = clone_pr_repo(pr)
                            pa_ret = process_actions(action_callback, action_command, pr_workdir, pr)
                            if pa_ret == 0:
                                pr.remove_from_labels(label.name.split('-')[0]+"-"+label.name.split('-')[1])

                        except Exception as e:
                            print("ERROR RUNNING RT {} with error: {}".format(action_name, e))
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
        ["git submodule update --init --recursive", repo_dir_str+"/"+repo_name],
        ["module use modulefiles/{}.intel && module load fv3".format(machine.name), repo_dir_str+"/"+repo_name]
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

def move_rt_logs(pr, pr_workdir):
    rt_log = 'RegressionTests_'+machine.name+'.intel.log'
    filepath = pr_workdir+'/tests/'+rt_log
    print("File path issssss {}".format(filepath))
    if os.path.exists(filepath):
        branch = pr.head.ref
        print("Branch used is: {}".format(branch))

        move_rt_commands = [
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

def process_actions(callback_fnc, command, pr_workdir, pr):

    def create_threaded_call(callback_fnc, command, cwd_in):

        def runInThread(callback_fnc, command, cwd_in):
            proc = subprocess.Popen(command, shell=True, cwd=cwd_in)
            proc.wait()
            globals()[callback_fnc](pr, pr_workdir)
            return

        thread = threading.Thread(target=runInThread,
                                  args=(callback_fnc, command, cwd_in))
        thread.start()

        return thread # returns immediately after the thread starts

    print("{} attempting command {}".format(machine.name.upper(), command))
    thread = create_threaded_call(callback_fnc, command, pr_workdir+'/tests')
    # create_threaded_call(callback_fnc, [1,2,3,4])
    # print("Thread is {}".format(thread))

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
