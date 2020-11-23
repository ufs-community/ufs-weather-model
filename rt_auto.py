# HEADER
# Written by Brian Curtis
# Automation of RT tasks using github cli

from github import Github as gh

# Initial items
g = gh("05630912c8d01499a0135f760f01b96f0f6fe127")
user = g.get_user()
repo = g.get_repo("BrianCurtis-NOAA/ufs-weather-model")
print(repo.name)

pulls = repo.get_pulls(state='open', sort='created', base='develop')
for pr in pulls:
    print(pr.number)
