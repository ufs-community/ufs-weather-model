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

# WARNING WARNING WARNING: DO NOT COMMIT WITH THE GHACCESSTOKEN
# GHACCESSTOKEN = None
# GITHUB RATE LIMIT FOR AUTHENTICATED USERS IS 5000 REQUESTS PER HOUR
db_filename = 'rt_auto.yml'
yaml_data = None
client = None
machine = None

class machine_info():

    def __init__(self, name, regexhostname, workdir):
        self.name = name
        self.regexhostname = regexhostname
        self.workdir = workdir



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
    machine = machine_info(mach_db_info['name'], mach_db_info['regexhostname'], mach_db_info['workdir'])
    return machine

# THE FOLLOWING FUNCTION REQUIRES AUTHENTICATION TOKEN TO BE SET
def is_collaborator(repo, userin):
    collaborator_bool = False
    allowed_collaborators = repo.get_collaborators()
    for collaborator in allowed_collaborators:
        collab_perms = repo.get_collaborator_permission(collaborator)
        print("Collab {}, Perms {}".format(collaborator.login, collab_perms))
        if collaborator.login == userin && 'write' in collab_perms.split():
            collaborator_bool = True
            
    return collaborator_bool

def get_yaml_subset(in_yaml, keyin):
    try:
        yaml_subset = in_yaml[keyin]
    except:
        sys.exit("Unable to get yaml subset: {}. Quitting".format(keyin))
    df = pd.DataFrame.from_dict(yaml_subset)
    return df

def connect_to_github(): # In case there's more needed later
    global GHACCESSTOKEN
    global client
    client = gh(GHACCESSTOKEN) #will work if None for public repos.

def read_yaml_data(filename):
    stream = open(db_filename, 'r')
    yaml_data = yaml.load(stream, Loader=yaml.SafeLoader)
    stream.close()
    return yaml_data

# def process_triggers(pr, triggers, commentids):
def process_triggers(comments_with_triggers, pr, repo):
    global machine
    if not machine:
        sys.exit("No machine information set, please use \"get_machine_info()\"")
    lastread = False
    approved_triggers = []
    trigger_comment_ids = []
    for comment in comments_with_triggers:
        body_split = comment.body.split()
        for word in body_split:
            if re.match("\+\+.+\+\+", word):
                if is_collaborator(repo, comment.user.login):
                    split_trigger = str(word).strip("+")
                    split_trigger = split_trigger.split(":")
                    if split_trigger[0].lower() == "read":
                        if bool(lastread):
                            print("Found old 'READ' comment: deleting")
                            pr.get_issue_comment(int(lastread)).delete()
                        lastread = comment.id
                    else:
                        approved_triggers.append(split_trigger)
                        trigger_comment_ids.append(comment.id)

    actions = []
    for (approved_trigger, trigger_comment_id) in zip(approved_triggers, trigger_comment_ids):
        if approved_trigger[0].lower() == machine.name.lower() or approved_trigger[0].lower() == "all":
            if trigger_comment_id > lastread:
                actions.append(approved_trigger[1])
                print("SAVED ACTION {} for machine {}".format(approved_trigger[1], machine.name))

    if not actions:
        return lastread, None
    actions = set(actions)
    return lastread, actions

def get_comments_with_triggers(comments):
    comment_list = []
    for comment in comments:
        if [re.match("\+\+.+\+\+", incomment) for incomment in comment.body.split()]:
            comment_list.append(comment)

    if not comment_list:
        return None
    else:
        return comment_list

def process_pulls(pulls, repo):
    for pr in pulls:
        comments = pr.get_issue_comments()
        comments_with_triggers = get_comments_with_triggers(comments)
        if not comments_with_triggers:
            # LETS STOP HERE AND UPDATE THE 'READ' TRIGGER
            # ASSUMING THERE'S ZERO 'READ' COMMENTS, ADDING ONE AT END
            pr.create_issue_comment("++READ++")
            continue

        lastread, actions = process_triggers(comments_with_triggers, pr, repo)
        if bool(lastread):
            if lastread != max(comment.id for comment in comments):
                pr.get_issue_comment(int(lastread)).delete()
                print("Deleted previous READ comment")
                pr.create_issue_comment("++READ++")


        if bool(actions):
            print("Accepted actions: {}".format(actions))
            validate_actions(actions)
            pr_workdir = clone_pr_repo(pr)
            process_actions(actions, pr_workdir, pr)
        else:
            print("No Actions Found in this PR, skipping.")

def validate_actions(actions):
    approved_commands = []
    valid_actions = get_yaml_subset(static_data, "actions")
    for name,command in valid_actions.values.tolist():
        for action in actions:
            if name == action.lower():
                print("Action {} approved as valid".format(action.lower()))
                approved_commands.append(command)
    return approved_commands



def clone_pr_repo(pr):
    global machine
    branch = pr.head.ref
    repo_name = pr.head.repo.name
    git_url = pr.head.repo.html_url
    repo_dir_str = machine.workdir+"/"+str(pr.id)+"/"+datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    print("repo_dir_str is: {}".format(repo_dir_str))

    create_repo_commands = [
        ["mkdir -p \""+repo_dir_str+"\"", machine.workdir],
        ["git clone -b "+branch+" "+git_url, repo_dir_str],
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

    return repo_dir_str+"/"+repo_name

def create_threaded_call(callback, run_fnc):
    def runInThread(callback, args, kwargs):
        proc = run_fnc
        proc.wait()
        callback()
        return

    thread = threading.Thread(target=runInThread,
                              args=(callback, args, kwargs))
    thread.start()

    return thread # returns immediately after the thread starts

def thread_end(pr, command, pr_workdir):
    global machine
    pr.create_issue_comment("AUTOMATED COMMENT: Finished Processing Command '{}' in directory '{}' on machine '{}'".format(command, pr_workdir, machine.name))
    return

def process_actions(actions, pr_workdir, pr):
    global machine

    approved_commands = validate_actions(actions)
    if bool(approved_commands):
        for command in approved_commands:
            print("{} is processing command {}".format(machine.name.upper(), command))
            pr.create_issue_comment("AUTOMATED COMMENT: Started Processing Command '{}' in directory '{}' on machine '{}'".format(command, pr_workdir, machine.name))
            create_threaded_call(thread_end(pr, command, pr_workdir+"/tests"), subprocess.Popen(command, shell=True, cwd=pr_workdir+"/tests"))

# Initial items
static_data = read_yaml_data('rt_auto.yml')
machine = get_machine_info(static_data)
repos = get_yaml_subset(static_data, "repository")
connect_to_github()

# Maybe a process_repo(repos) function start?
for name,address,base in repos.values.tolist():
    repo = client.get_repo(address)
    pull_reqs = repo.get_pulls(state='open', sort='created', base=base)
    triggers = process_pulls(pull_reqs, repo)
