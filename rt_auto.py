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
        return user
    else:
        sys.exit("Username {} does not match approved list. Quitting".format(hostname))

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


# Initial items
machine = get_machine_name() #ALWAYS RUN THIS FIRST
repos = get_repos_list()
client = connect_to_github()
for name,address,base in repos.values.tolist():
    repo = client.get_repo(address)
    pulls = repo.get_pulls(state='open', sort='created', base='develop')
    for pr in pulls:
        print(pr.head.repo.default_branch)
        print(pr.head.repo.master_branch)
        print(pr.number)
        comments = pr.get_issue_comments()
        for comment in comments:
            comment_user = comment.user.login
            print("comment user: {}".format(comment_user))
            is_user_approved(comment_user)
            print(comment.body)




# CHECK IS USER ON APPROVED LIST

# CHECK IF COMMENT HAS ALREADY BEEN ACKNOWLEDGED

# CHECK IF COMMENT HAS ACTION TRIGGER: ++<machine>
