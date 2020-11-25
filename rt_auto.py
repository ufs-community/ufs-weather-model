# HEADER
# Written by Brian Curtis
# Automation of RT tasks using github cli

from github import Github as gh
import pandas as pd
import socket
import re
import sys
import yaml

GHACCESSTOKEN = None

def get_machine_name():
    hostname = socket.gethostname()
    machines = read_yaml_data("machines")
    df = pd.DataFrame.from_dict(machines)
    A = df['name'][df.index[df.regexhostname.str.match(hostname)].tolist()]
    if len(A) == 1:
        machine = list(A)[0]
        print("Approved Machine Found: {}".format(machine))
        return machine
    else:
        sys.exit("Hostname {} does not match approved list. Quitting".format(hostname))

def is_user_approved(userin):
    approvedusers = read_yaml_data('approvedusers')
    userlist = pd.DataFrame.from_dict(approvedusers)
    A = userlist['username'][userlist.index[userlist.username.str.match(userin)].tolist()]
    if len(A) == 1:
        user = list(A)[0]
        print("Approved Github Username Found: {}".format(user))
        return True
    else:
        print("User not found in approved list.")
        return False

def get_repos_list():
    repos = read_yaml_data("repository")
    df = pd.DataFrame.from_dict(repos)
    return df

def connect_to_github(): # In case there's more needed later
    global GHACCESSTOKEN
    client = gh(GHACCESSTOKEN) #will work if None for public repos.
    return client

def read_yaml_data(keyin):
    stream = open("rt_auto.yml", 'r')
    dictionary = yaml.load(stream, Loader=yaml.SafeLoader)
    try:
        dataout = dictionary[keyin]
    except KeyError as e:
        sys.exit("No data available with name {}. Quitting".format(keyin))
    return dataout

def clone_pr_repo():
    return


def check_for_trigger(comments):
    # CHECK FOR TRIGGER WORD IN COMMENTS, IF TRIGGERED, CHECK IF USER IS
    # APPROVED BEFORE ADDING TO LISTS
    trigger = []
    commentid = []
    for comment in comments:
        body = comment.body
        body_split = body.split()
        for word in body_split:
            if re.match("\+\+.+\+\+", word):
                print("FOUND MATCH {} at comment {}".format(word, comment.id))
                print("USER IS: {}".format(comment.user.login))
                approve_bool = is_user_approved(comment.user.login)
                if approve_bool:
                    print("User {} is approved, continuing".format(comment.user.login))
                    trigger.append(word)
                    commentid.append(comment.id)
                else:
                    print("User is NOT approved, skipping trigger")
                    continue
    if not trigger:
        return None, None
    else:
        return trigger, commentid
# Initial items
machine = get_machine_name() #ALWAYS RUN THIS FIRST
repos = get_repos_list()
client = connect_to_github()
for name,address,base in repos.values.tolist():
    repo = client.get_repo(address)
    pulls = repo.get_pulls(state='open', sort='created', base=base)
    for pr in pulls:
        comments = pr.get_issue_comments()
        trigger, commentid = check_for_trigger(comments)
        if trigger == None:
            continue
        else:
            print("Out Triggers, commentid: {}, {}".format(trigger, commentid))
        print(pr.head.repo.clone_url)
        print(pr.head.ref)





# CHECK IS USER ON APPROVED LIST

# CHECK IF COMMENT HAS ALREADY BEEN ACKNOWLEDGED

# CHECK IF COMMENT HAS ACTION TRIGGER: ++<machine>
