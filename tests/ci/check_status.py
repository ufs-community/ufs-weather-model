#!/usr/bin/env python3
import os
import re
import sys
import json
import time
from urllib.request import urlopen, Request
from datetime import datetime, timedelta


def check_build(request):
    """Check if all build jobs are completed successfully.
    API endpoint: api.github.com/repos/{owner}/{repo}/actions/runs/{run_id}/jobs
    """
    all_completed = False
    while not all_completed:
        time.sleep(20)
        response = urlopen(request)
        data = json.loads(response.read().decode())["jobs"]
        ids = [x["id"] for x in data if re.search("Build", x["name"])]
        if len(ids) == 1 and next(re.search("matrix", x["name"]) for x in data if x["id"] in ids):
            break
        if len(ids) == 0:
            continue
        all_completed = all(
            [x["status"] == "completed" for x in data if x["id"] in ids])
    return all([x["conclusion"] == "success" for x in data if x["id"] in ids])


def check_completion(request, job_name):
    """Check if a job is completed successfully.
    API endpoint: api.github.com/repos/{owner}/{repo}/actions/runs/{run_id}/jobs
    """
    completed = False
    while not completed:
        time.sleep(60)
        response = urlopen(request)
        data = json.loads(response.read().decode())["jobs"]
        cid = next((x["id"]
                   for x in data if x["name"] == job_name), "not found")
        if cid == "not found":
            continue
        completed = next(
            x["status"] == "completed" for x in data if x["id"] == cid)
    return next(x["conclusion"] == "success" for x in data if x["id"] == cid)


def check_test(request):
    """Wait for a workflow run to be completed.
    API endpoint: api.github.com/repos/{owner}/{repo}/actions/runs/{run_id}
    """
    completed = False
    while not completed:
        time.sleep(20)
        response = urlopen(request)
        data = json.loads(response.read().decode())
        completed = (data["status"] == "completed")


def check_ec2(url, request, myid):
    """Wait for all previous workflow runs to finish using ec2 instances.
    API endpoint: api.github.com/repos/{owner}/{repo}/actions/runs
    """
    response = urlopen(request)
    data = json.loads(response.read().decode())["workflow_runs"]
    tformat = "%Y-%m-%dT%H:%M:%SZ"
    mytime = datetime.strptime(next(x["created_at"]
                               for x in data if x["id"] == myid), tformat)

    workflows = {}
    in_progress = []
    for x in data:
        oldtime = datetime.strptime(x["created_at"], tformat)
        dt = mytime - oldtime
        if x["name"] == "Helpers" and dt >= timedelta() and x["id"] != myid:
            request = Request(url+"/"+str(x["id"])+"/jobs")
            token = os.environ["AUTH"]
            request.add_header("Authorization", "token %s" % token)
            workflows[x["id"]] = request
            in_progress.append(x["id"])

    while True:
        print("in_progress: ", in_progress)
        if len(in_progress) == 0:
            break
        time.sleep(20)
        done = []
        for cid in reversed(in_progress):
            data = json.loads(urlopen(workflows[cid]).read().decode())["jobs"]
            start_status = next(
                (x["status"] for x in data if x["name"] == "Start runners"), "not found")
            stop_status = next(
                (x["status"] for x in data if x["name"] == "Stop runners"), "not found")
            if start_status == "not found" or stop_status == "not found":
                break
            if start_status == "completed" and stop_status == "completed":
                done.append(cid)
        print("done: ", done)
        if len(done) != 0:
            [workflows.pop(k) for k in done]
            [in_progress.remove(k) for k in done]


def main():
    url = sys.stdin.read()
    request = Request(url)
    try:
        token = os.environ["AUTH"]
        request.add_header("Authorization", "token %s" % token)
    except KeyError:
        pass

    if sys.argv[1] == "build":
        print("success") if check_build(request) else print("failure")
    elif sys.argv[1] == "completion":
        print("success") if check_completion(
            request, sys.argv[2]) else print("failure")
    elif sys.argv[1] == "ec2":
        check_ec2(url, request, int(sys.argv[2]))
    elif sys.argv[1] == "test":
        check_test(request)


if __name__ == "__main__":
    main()
