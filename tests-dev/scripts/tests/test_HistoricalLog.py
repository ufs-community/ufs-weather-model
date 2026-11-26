import pytest
import numpy as np
from scripts.APICall import APICall
from scripts.Log import Log
from scripts.get_data import *

def test_init_historical_Log(historical_herc_log):
   log = HistoricalLog("HeRcuLes", ["e30900f9", "730df864", "ca94c232"])
   assert historical_herc_log.machine == log.machine
   assert historical_herc_log.text_per_log == log.text_per_log
   assert historical_herc_log.repo_commits == log.repo_commits

def test_get_repo_commits(historical_herc_log):
   assert historical_herc_log.get_repo_commits() == ["e30900f9", "730df864", "ca94c232"]

def test_fetch_log_text_w_no_commits(historical_log, caplog): 

   with pytest.raises(SystemExit) as e:
      historical_log.pr_head_commit = None
      historical_log._fetch_log_text(historical_log.pr_head_commit)
   assert e.value.code == None
   assert caplog.records[0].message == "An appropriate commit(s) was not provided."

def test_fetch_log_text_for_develop(historical_log, most_recent_commits, hercules_log_texts_2882): 
   """Check that the log texts extracted by the API are the same as the hercules log texts that we expect."""
   historical_log.repo_commits = most_recent_commits[1:6]
   historical_log._fetch_log_text(historical_log.repo_commits)

   for hash in most_recent_commits[1:6]:
      assert historical_log.text_per_log[hash] == hercules_log_texts_2882[hash]

def test_get_log_text_for_3_commits(most_recent_commits, hercules_log_texts_2882):
   """Check that HistoricalLog successfully gets the log text for the three most recent repo commits in hercules_log_texts_2882"""
   log = HistoricalLog("Hercules", most_recent_commits[1:4]) # Info for PR #2882
   log.get_log_text(log.repo_commits)
   for text in log.text_per_log.values():
      assert text in hercules_log_texts_2882.values()

def test_get_instance_test_data(historical_log_PR_2882, hercules_log_texts_2882, log_instance_results_2882_1):
   """From the log for PR 2882, extract test data. Compare it with the expected data to be sure it's the same.
   """
   tests_for_log_instance = historical_log_PR_2882._get_instance_test_data(hercules_log_texts_2882[historical_log_PR_2882.repo_commits[0]])
   assert tests_for_log_instance == log_instance_results_2882_1

def test_compile_historical_log_data(historical_log_PR_2882, hercules_log_texts_2882, hercules_sample_historical_log_data): 

   historical_log_PR_2882.text_per_log = {commit: text for commit, text in hercules_log_texts_2882.items() if commit != '4760a41a2ba012b236361d99a0248b718a06921b'}
   historical_log_PR_2882.compile_historical_log_data()
   
   for test in hercules_sample_historical_log_data:
      assert list(historical_log_PR_2882.historical_runtime_data[test].values()) == hercules_sample_historical_log_data[test]["runtime"]
      assert list(historical_log_PR_2882.historical_mem_data[test].values()) == hercules_sample_historical_log_data[test]["memory"]

def test_fetch_historical_data(historical_log_PR_2882, hercules_sample_historical_log_data):
   """
   Extract runtime/memory data for the authoritative repository's last several commits.
   Actual number determined by Manager.get_hashes()
   """
   historical_log_PR_2882.fetch_historical_data()

   for test in hercules_sample_historical_log_data:
      assert list(historical_log_PR_2882.historical_runtime_data[test].values()) == hercules_sample_historical_log_data[test]["runtime"]
      assert list(historical_log_PR_2882.historical_mem_data[test].values()) == hercules_sample_historical_log_data[test]["memory"]

def test_get_historical_runtime_data(historical_log_PR_2882, hercules_sample_historical_log_data):
   historical_log_PR_2882.historical_runtime_data = historical_log_PR_2882.get_historical_runtime_data()
   for test in hercules_sample_historical_log_data:
      assert list(historical_log_PR_2882.historical_runtime_data[test].values()) == hercules_sample_historical_log_data[test]["runtime"]
   
def test_get_historical_mem_data(historical_log_PR_2882, hercules_sample_historical_log_data):
   historical_log_PR_2882.historical_mem_data = historical_log_PR_2882.get_historical_mem_data()
   for test in hercules_sample_historical_log_data:
      assert list(historical_log_PR_2882.historical_mem_data[test].values()) == hercules_sample_historical_log_data[test]["memory"]

def test_calculate_runtime_stats(historical_log_PR_2882, hercules_sample_historical_log_data):
   """Check that test: [mean, stdev] == test: [mean, std]"""
      
   historical_log_PR_2882.runtime_stats = historical_log_PR_2882.get_runtime_stats()
   historical_log_PR_2882.calculate_runtime_stats()

   assert historical_log_PR_2882.runtime_stats != None
   for test in hercules_sample_historical_log_data:
      assert historical_log_PR_2882.runtime_stats[test][0] == round(np.mean(hercules_sample_historical_log_data[test]["runtime"]), 5)
      assert historical_log_PR_2882.runtime_stats[test][1] == round(np.std(hercules_sample_historical_log_data[test]["runtime"]), 5)

def test_get_runtime_stats(historical_log_PR_2882):
   historical_log_PR_2882.get_runtime_stats()
   assert historical_log_PR_2882.runtime_stats != None

def test_calculate_mem_stats(historical_log_PR_2882, hercules_sample_historical_log_data):
   historical_log_PR_2882.mem_stats = historical_log_PR_2882.get_mem_stats()
   historical_log_PR_2882.calculate_mem_stats()

   assert historical_log_PR_2882.mem_stats != None
   for test in hercules_sample_historical_log_data:
      assert historical_log_PR_2882.mem_stats[test][0] == round(np.mean(hercules_sample_historical_log_data[test]["memory"]), 5)
      assert historical_log_PR_2882.mem_stats[test][1] == round(np.std(hercules_sample_historical_log_data[test]["memory"]), 5)

def test_get_mem_stats(historical_log_PR_2882):
   historical_log_PR_2882.get_mem_stats()
   assert historical_log_PR_2882.mem_stats != None
