import pytest
import numpy as np
from scripts.APICall import APICall
from scripts.Log import Log
from scripts.get_data import *

def test_init_hercules_Log(herc_log):
   log = Log("HeRcuLes", ["e30900f9", "730df864", "ca94c232"], "36361d99")
   assert herc_log.machine == log.machine
   assert herc_log.text_per_log == log.text_per_log
   assert herc_log.pr_head_commit == log.pr_head_commit
   assert herc_log.repo_commits == log.repo_commits

def test_get_pr_head(ursa_log):
   assert ursa_log.get_pr_head() == "24j3m23m"

def test_get_repo_commits(ursa_log):
   assert ursa_log.get_repo_commits() == ["8ffeb8e3", "34khj5bh", "j34bnj3m"]

def test_fetch_log_text_w_no_commits(herc_log, caplog): 

   with pytest.raises(SystemExit) as e:
      herc_log.pr_head_commit = None
      herc_log._fetch_log_text(herc_log.pr_head_commit)
   assert e.value.code == None
   assert caplog.records[0].message == "An appropriate commit(s) was not provided."

def test_fetch_log_text_for_pr_head(log, most_recent_commits, hercules_log_texts_2882): 
   """Check that the log texts extracted by the API are the same as the hercules log texts that we expect."""
   log.pr_head_commit = most_recent_commits[0]
   log._fetch_log_text([log.pr_head_commit])
   
   assert log.text_per_log[log.pr_head_commit] == hercules_log_texts_2882['4760a41a2ba012b236361d99a0248b718a06921b']

def test_fetch_log_text_for_develop(log, most_recent_commits, hercules_log_texts_2882): 
   """Check that the log texts extracted by the API are the same as the hercules log texts that we expect."""
   log.repo_commits = most_recent_commits[1:6]
   log._fetch_log_text(log.repo_commits)

   for hash in most_recent_commits[1:6]:
      assert log.text_per_log[hash] == hercules_log_texts_2882[hash]

def test_get_log_text_for_pr_head(log, hercules_log_texts_2882):
   log.pr_head_commit = ['369cead91c98eb5c72da81ff78925250dad08903'] # Commit from HEAD of PR #2882
   log.get_log_text(log.pr_head_commit)
   assert log.text_per_log[log.pr_head_commit[0]] == hercules_log_texts_2882['4760a41a2ba012b236361d99a0248b718a06921b']

def test_get_log_text_for_3_commits(most_recent_commits, hercules_log_texts_2882):
   log = Log("Hercules", most_recent_commits[:3], "369cead9") # Info for PR #2882
   log.get_log_text(log.repo_commits)
   for text in log.text_per_log.values():
      assert text in hercules_log_texts_2882.values()

def test_get_instance_test_data(log_PR_2882, hercules_log_texts_2882, log_instance_results_2882_0):
   """From the log for PR 2882, extract test data. Compare it with the expected data to be sure it's the same.
   """
   tests_for_log_instance = log_PR_2882._get_instance_test_data(hercules_log_texts_2882[log_PR_2882.pr_head_commit])
   assert tests_for_log_instance == log_instance_results_2882_0

def test_get_current_pr_data(log_PR_2882, hercules_log_texts_2882, log_instance_results_2882_0):
   
   log_PR_2882.text_per_log = hercules_log_texts_2882
   pr_log_data = log_PR_2882.get_current_pr_data()
   assert pr_log_data == log_instance_results_2882_0

def test_no_current_pr_data(herc_log, caplog):
   """Check that the final log message (from get_current_pr_data() rather than the functions it calls) is correct."""
   
   with pytest.raises(SystemExit) as e:
      herc_log.get_current_pr_data()
   assert e.value.code == None
   assert caplog.records[4].message == "Cannot fetch current PR data."

def test_get_current_pr_runtime_data(log_PR_2882, hercules_current_pr_data):
   
   log_PR_2882.current_pr_runtime_data = log_PR_2882.get_current_pr_data()
   for test in hercules_current_pr_data:
      assert log_PR_2882.current_pr_runtime_data[test][0] == hercules_current_pr_data[test][0]

   # Could change current_pr_log_data to pandas DataFrame to access by column
def test_get_current_pr_mem_data(log_PR_2882, hercules_current_pr_data):
   
   log_PR_2882.current_pr_mem_data = log_PR_2882.get_current_pr_data()
   for test in hercules_current_pr_data:
      assert log_PR_2882.current_pr_mem_data[test][1] == hercules_current_pr_data[test][1]

def test_compare_runtimes(log_PR_2882, log_instance_results_2882_0, hercules_runtime_results, 
                          hercules_mean_std, hercules_sample_historical_log_data):
   
   expected = hercules_runtime_results["hercules"]
   log_PR_2882.current_pr_log_data = log_instance_results_2882_0
   log_PR_2882.runtime_stats = {test: values[0:2] for test, values in hercules_mean_std.items()}
   # Zip together hercules_sample_historical_log_data with each value's corresponding hash (pairing required for compare_results())
   log_PR_2882.historical_runtime_data = {
      test: dict(zip(log_PR_2882.repo_commits, data["runtime"]))
      for test, data in hercules_sample_historical_log_data.items()
   }
   log_PR_2882.compare_runtimes()

   for test, value in expected.items(): 
      assert log_PR_2882.runtime_results[test] == value

def test_compare_mem(log_PR_2882, log_instance_results_2882_0, sample_memory_results, hercules_mean_std, hercules_sample_historical_log_data):

   expected = sample_memory_results["hercules"]
   log_PR_2882.current_pr_log_data = log_instance_results_2882_0
   log_PR_2882.mem_stats = {test: values[2:] for test, values in hercules_mean_std.items()}
   # Zip together hercules_sample_historical_log_data with each value's corresponding hash (pairing required for compare_results())
   log_PR_2882.historical_mem_data = {
         test: dict(zip(log_PR_2882.repo_commits, data["memory"]))
         for test, data in hercules_sample_historical_log_data.items()
      }
   log_PR_2882.compare_mem()

   for test, value in expected.items(): 
      assert log_PR_2882.mem_results[test] == value
