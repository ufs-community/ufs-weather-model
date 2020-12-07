# HEADER
# Written by Brian Curtis
# Automation of RT tasks using github cli

from github import Github as gh
import pandas as pd
import socket
import re
import sys
import yaml
import numpy as np
from os import path
import pdb

GHACCESSTOKEN = None
db_filename = 'rt_auto.yml'
yaml_data = None
client = None
machine = None

def get_machine_name(static_data):
    global machine
    hostname = socket.gethostname()
    machines = get_yaml_subset(static_data, "machines")
    A = machines['name'][machines.index[machines.regexhostname.str.match(hostname)].tolist()]
    if len(A) == 1:
        machine = str(list(A)[0])
        print("Approved Machine Found: {}".format(machine))
        return machine
    else:
        sys.exit("Hostname {} does not match approved list. Quitting".format(hostname))

# THE FOLLOWING FUNCTION REQUIRES AUTHENTICATION TOKEN TO BE SET
def is_collaborator(repo, userin):
    collaborator_bool = False
    allowed_collaborators = repo.get_collaborators()
    for collaborator in allowed_collaborators:
        if collaborator.login == userin:
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

def clone_pr_repo():
    return

def process_triggers(triggers, commentids):
    global machine
    lastread = False
    trigger_arr = []
    commentid_arr = []
    for (trigger, commentid) in zip(triggers, commentids):
        split_trigger = str(trigger).strip("+")
        split_trigger = split_trigger.split(":")
        if split_trigger[0].lower() == "read":
            lastread = commentid
        else:
            trigger_arr.append(split_trigger)
            commentid_arr.append(commentid)
    actions = []
    for (mach_act, req_commentid) in zip(trigger_arr, commentid_arr):
        if mach_act[0].lower() == machine.lower() or mach_act[0].lower() == "all":
            if req_commentid > lastread:
                actions.append(mach_act[1])
                print("SAVED ACTION: {}".format(mach_act[1]))
            else:
                print("COMMENTID {} WAS LOWER THAN LASTREAD {}, SKIPPING...".format(req_commentid, lastread))

    actions = set(actions)
    return(lastread, actions)

def process_pulls(repo, address):
    pulls = repo.get_pulls(state='open', sort='created', base=base)
    # mess_with = pd.DataFrame.from_dict(yaml_data)
    # pdb.set_trace()
    # address_data = get_yaml_subset(address)
    # lasttriggerid = address_data['lasttriggerid']
    for pr in pulls:
        trigger, commentid = get_triggers(repo, pr)
        if trigger == None:
            continue
        else:
            lastread, actions = process_triggers(trigger, commentid)
            print("Lastread: {}".format(lastread))
            if bool(lastread):
                pr.get_issue_comment(int(lastread)).delete()
                print("Deleted previous READ comment")
            print("Accepted actions: {}".format(actions))
    if bool(actions):
        print("Lets Process Actions")
        process_actions(actions)
    else:
        print("No Actions Found in this PR, skipping.")

def process_actions(actions):
    return

def get_triggers(repo, pr):
    # CHECK FOR TRIGGER WORD IN COMMENTS, IF TRIGGERED, CHECK IF USER IS
    # APPROVED BEFORE ADDING TO LISTS
    trigger = []
    commentid = []
    for comment in pr.get_issue_comments():
        print("Comment ID: {}".format(comment.id))
        # if(comment.id <= lasttriggerid):
        #     continue
        body_split = comment.body.split()
        for word in body_split:
            if re.match("\+\+.+\+\+", word):
                print("Found trigger {} at comment {}".format(word, comment.id))
                # approve_bool = is_user_approved(comment.user.login)
                if is_collaborator(repo, comment.user.login):
                    trigger.append(word)
                    commentid.append(comment.id)
                else:
                    print("User is NOT approved, skipping trigger")
                    continue
    if not trigger:
        return None, None
    else:
        pr.create_issue_comment("Automated Comment: ++READ++")
        return trigger, commentid

# Initial items
static_data = read_yaml_data('rt_auto.yml') #MUST DO THIS FIRST
machine = get_machine_name(static_data) #ALWAYS RUN THIS FIRST
repos = get_yaml_subset(static_data, "repository")
connect_to_github()
for name,address,base in repos.values.tolist():
    repo = client.get_repo(address)
    triggers = process_pulls(repo, address)



# CHECK IS USER ON APPROVED LIST

# CHECK IF COMMENT HAS ALREADY BEEN ACKNOWLEDGED

# CHECK IF COMMENT HAS ACTION TRIGGER: ++<machine>
