#!/usr/bin/env python3
# The following env variables are assumed available and valid:
# GITHUB_ACTOR, GITHUB_RUN_ID, TRIGGER_ID, TRIGGER_BR
import os
import re
import sys
import json

def check_skip(data):
  msg = data["head_commit"]["message"]
  if re.search("skip-ci", msg):
    return "yes"
  else:
    return "no"

def cancel_workflow(data):
  wfs=[x["id"] for x in data if x["head_repository"] is not None and
        re.search(os.environ["GITHUB_ACTOR"], x["head_repository"]["owner"]["login"]) and
        x["id"]!=int(os.environ["GITHUB_RUN_ID"]) and
        x["id"]!=int(os.environ["TRIGGER_ID"]) and
        x["head_branch"]==os.environ["TRIGGER_BR"] and
        (x["status"]=="queued" or x["status"]=="in_progress")]

  return wfs

def main():

  if sys.argv[1]=="check_skip":
    data = json.load(sys.stdin)["workflow_run"]
    ans = check_skip(data)
    print(ans)
  elif sys.argv[1]=="get_trigger_id":
    print(json.load(sys.stdin)["workflow_run"]["id"])
  elif sys.argv[1]=="get_trigger_br":
    print(json.load(sys.stdin)["workflow_run"]["head_branch"])
  elif sys.argv[1]=="cancel_workflow":
    data = json.load(sys.stdin)["workflow_runs"]
    wfs = cancel_workflow(data)
    if len(wfs)==0:
      print("")
    else:
      print(*wfs)
  else:
    print("ERROR")

if __name__ == "__main__": main()
