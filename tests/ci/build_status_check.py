#!/usr/bin/env python3

import re
import sys
import json
import time
from urllib.request import urlopen

def update_url_data(response):
  data = json.loads(response.read().decode())
  indices=[]
  for n in range(data["total_count"]):
    if re.search("Build", data["jobs"][n]["name"]):
      indices.append(n)

  if len(indices) == 0:
    sys.exit(1)

  return data, indices

def main():

  time.sleep(180)
  url = json.load(sys.stdin)["workflow_run"]["jobs_url"]

  status="not-completed"
  no_completed_jobs = 0

  while status != "completed":
    response = urlopen(url)
    data, indices = update_url_data(response)

    for i in indices:
      if data["jobs"][i]["status"] == "completed":
        no_completed_jobs += 1

    if no_completed_jobs == len(indices):
      status = "completed"
    else:
      no_completed_jobs = 0
      time.sleep(10)

  time.sleep(10)
  conclusion="failure"
  no_successful_jobs = 0
  for i in indices:
    if data["jobs"][i]["conclusion"] == "success":
      no_successful_jobs += 1

  if no_successful_jobs == len(indices):
    conclusion = "success"

  print(conclusion)

if __name__ == "__main__": main()
