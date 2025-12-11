import numpy as np
import pytest
from pathlib import Path
from scripts.get_data import *

@pytest.mark.parametrize("endpoint", [
   f"commits?path=tests/logs/RegressionTests_ursa.log&per_page=1", #fetch_repo_commits_endpoint
   f"pulls/2882", #get_pr_head_endpoint
   f"contents/tests/logs/RegressionTests_ursa.log", #fetch_log_text_endpoint
   ])
@pytest.mark.parametrize("num_commits", [1, 5, 7])
def test_init_APICall(set_env_vars, monkeypatch, endpoint, num_commits):
   
   set_env_vars
   # Set token env var for duration of test only
   monkeypatch.setenv("GITHUB_TOKEN", "fake_github_pat_12BWCMCFZkhj35klj3h34kjh4kkjm3whe4nr")
   api_call = APICall(endpoint, num_commits)

   assert api_call.token == "fake_github_pat_12BWCMCFZkhj35klj3h34kjh4kkjm3whe4nr"
   assert api_call.base_url == "https://api.github.com/repos/ufs-community/ufs-weather-model"
   assert api_call.endpoint == endpoint
   assert api_call.url == f"https://api.github.com/repos/ufs-community/ufs-weather-model/{endpoint}"
   assert api_call.num_commits == num_commits
   assert api_call.header == {
         "Accept": "application/vnd.github.v3+json",
         "Authorization": f"Bearer fake_github_pat_12BWCMCFZkhj35klj3h34kjh4kkjm3whe4nr",
         "X-GitHub-Api-Version": "2022-11-28",
         "Accept": "application/vnd.github.raw"
      }

def test_init_hercules_Log(herc_log):
   assert herc_log.machine == "hercules"
   assert herc_log.text_per_log == []

def test_fetch_repo_commits(herc_log, set_env_vars, hercules_most_recent_commits):
   """Test the API call and it's ability to get the 10 most recent commits.
   Because the actual commits will change, only the length is checked. 
   Ability to extract the proper commit(s) is tested in test_get_pr_head().
   When running tests locally, create a GitHub token and set it as an environment variable. 
   Then, try one of the following methods to set the token:
   1. In the command line: 
      export GITHUB_TOKEN=fake_github_pat_12BWCMCFZkhj35klj3h34kjh4kkjm3whe4nr
   OR
   2. In this script, add the monkeypatch fixture to the arguments, uncomment the following line, 
      and add the actual token value:
         monkeypatch.setenv("GITHUB_TOKEN", "fake_github_pat_12BWCMCFZkhj35klj3h34kjh4kkjm3whe4nr")
      Remove this line of code before committing anything. 
   """
   set_env_vars
   herc_log._fetch_repo_commits(10)

   assert len(herc_log.repo_commits) == len(hercules_most_recent_commits)

def test_fetch_repo_commits_w_no_commits(herc_log, set_env_vars, monkeypatch, caplog):
   """Test the ability to handle errors when no commits are returned
   """
   # Need to mock case where no commits or fewer commits than expected are returned.
   set_env_vars
   monkeypatch.setenv("GITHUB_TOKEN", "fake_github_pat_12BWCMCFZkhj35klj3h34kjh4kkjm3whe4nr")
   herc_log._fetch_repo_commits()

   assert caplog.records[0].message == "API Call failed. The sha does not exist!"

def test_get_pr_head(herc_log, set_env_vars):
   """Test the API call and it's ability to get the PR 2882's head commit. 
   When running tests locally, create a GitHub token and set it as an environment variable 
   using one of the methods listed in test_fetch_repo_commits() above.
   """
   set_env_vars
   herc_log._get_pr_head()

   assert herc_log.pr_head_commit == ["369cead91c98eb5c72da81ff78925250dad08903"]

def test_fetch_log_text_w_no_commits(herc_log, caplog): 
   herc_log.pr_head_commit = None
   herc_log._fetch_log_text(herc_log.pr_head_commit)
   assert caplog.records[0].message == "An appropriate commit(s) was not provided. Call _get_pr_head() or _fetch_repo_commits() first."

def test_fetch_log_text_for_pr_head(herc_log, hercules_most_recent_commits, hercules_log_texts_2882): 
   """Check that the log texts extracted by the API are the same as the hercules log texts that we expect."""
   herc_log.pr_head_commit = hercules_most_recent_commits[0]
   #herc_log.repo_commits = hercules_most_recent_commits[1:]
   # Need to mock API call
   #herc_log._fetch_log_text(herc_log.repo_commits)
   herc_log._fetch_log_text(herc_log.pr_head_commit)

   #assert herc_log.text_per_log == hercules_log_texts_2882[1:]
   assert herc_log.text_per_log[0] == hercules_log_texts_2882[0]

def test_fetch_log_text_for_develop(herc_log, hercules_most_recent_commits, hercules_log_texts_2882): 
   """Check that the log texts extracted by the API are the same as the hercules log texts that we expect."""
   herc_log.repo_commits = hercules_most_recent_commits[1:]
   herc_log._fetch_log_text(herc_log.repo_commits)
   
   assert herc_log.text_per_log[1:] == hercules_log_texts_2882[1:]

def test_get_instance_test_data(herc_log, hercules_log_texts_2882, log_instance_results_2882_0):
   """From the log for PR 2882, extract test data. Compare it with the expected data to be sure it's the same.
   """
   tests_for_log_instance = herc_log._get_instance_test_data(hercules_log_texts_2882[0])
   assert tests_for_log_instance == log_instance_results_2882_0

      
def test_compile_historical_log_data(herc_log, hercules_log_texts_2882, hercules_sample_historical_log_data): 
   
   herc_log.text_per_log = hercules_log_texts_2882
   herc_log._compile_historical_log_data()
   
   # Are all items in the hercules_sample_historical_log_data in herc_log.historical_rt_mem_data? 
   for test in hercules_sample_historical_log_data:
      assert herc_log.historical_rt_mem_data[test] == hercules_sample_historical_log_data[test]
               
def test_calculate_stats(herc_log, hercules_sample_historical_log_data, hercules_mean_std):
   
   herc_log.historical_rt_mem_data = hercules_sample_historical_log_data
   herc_log.calculate_stats()

   for test in hercules_mean_std: 
      assert hercules_mean_std[test] == herc_log.test_stats[test]

def test_compare_results(herc_log, hercules_log_texts_2882, log_instance_results_2882_0, hercules_mean_std): 

   current_log = log_instance_results_2882_0
   herc_log.text_per_log = hercules_log_texts_2882
   herc_log.test_stats = hercules_mean_std
   herc_log.compare_results()

   for test in herc_log.test_stats: 
      hi_runtime = herc_log.test_stats[test][0] + herc_log.test_stats[test][1]
      hi_memory = herc_log.test_stats[test][2] + herc_log.test_stats[test][3]

      # Could improve test to check for correct warn vs. fail status
      if current_log[test][0] > hi_runtime:
         assert herc_log.runtime_results[test] != '✅'
      if current_log[test][1] > hi_memory:
         assert herc_log.memory_results[test] != '✅'

def test_create_json(stats_dict_snippet):
   
   path = Path('data')
   path.mkdir(exist_ok = True)
   create_json(stats_dict_snippet, 'stats')

   with open('test_file_stats.json', 'r') as test_stats_file, open ('data/stats.json', 'r') as new_json:
      test_file_content = test_stats_file.read()
      new_json_content = new_json.read()
   
   assert test_file_content == new_json_content


def test_load_json(stats_dict_snippet):
   machine = "orion"
   orion_snippet = load_json('test_file_stats.json')[machine]
   assert orion_snippet == stats_dict_snippet['orion']

def test_main_e2e_cached_stats(monkeypatch):
   """Test that main function runs to completion."""

   monkeypatch.setenv("MACHINES", "hercules")
   monkeypatch.setenv("TEST_STATS", "test_file_stats.json")
   exit_code = main()

   assert exit_code == 0

def test_main_e2e_no_cached_stats(monkeypatch):
   """Test that main function runs to completion."""

   monkeypatch.setenv("MACHINES", "hercules")
   exit_code = main()

   assert exit_code == 0