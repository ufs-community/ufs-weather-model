# HEADER
# Written by Brian Curtis
# Automation of RT tasks using github cli

#from github import Github as gh
import socket
import re
import sys

USERNAME = None
PASSWORD = None

def get_machine_name():
    hostname = socket.gethostname()
#    hfe##  # HERA
    if re.match()
#    v##X#  # VENUS
#    m##X#  # MARS
#    slogin# # SURGE
#    Orion-login-3.HPC.MsState.Edu #ORION
    print(re.match("hfe\d\d", hostname))
    print(re.match("N[a-z][a-z][a-z]", hostname))

def connect_to_github(): # In case there's more needed later
    global USERNAME
    global PASSWORD
    client = gh(USERNAME, PASSWORD) #essentially
    return client

# Initial items
get_machine_name()
sys.exit()
client = connect_to_github()
repo = client.get_repo("BrianCurtis-NOAA/ufs-weather-model")
print(repo.name)

pulls = repo.get_pulls(state='open', sort='created', base='develop')
for pr in pulls:
    print(pr.number)
    comments = pr.get_issue_comments()
    for comment in comments:
        print(comment.body)

# CHECK IS USER ON APPROVED LIST

# CHECK IF COMMENT HAS ALREADY BEEN ACKNOWLEDGED

# CHECK IF COMMENT HAS ACTION TRIGGER: ++<machine>
