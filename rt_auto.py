# HEADER
# Written by Brian Curtis
# Automation of RT tasks using github cli

from github import Github as gh

USERNAME = None
PASSWORD = None

def connect_to_github(): # In case there's more needed later
    global USERNAME
    global PASSWORD
    client = gh(USERNAME, PASSWORD) #essentially
    return client

# Initial items
client = connect_to_github()
repo = client.get_repo("BrianCurtis-NOAA/ufs-weather-model")
print(repo.name)

pulls = repo.get_pulls(state='open', sort='created', base='develop')
for pr in pulls:
    print(pr.number)
