import requests
import sys
import re
import numpy as np
import pandas as pd
import logging
from .APICall import APICall

class Log():
   """A Regression Test log file."""
   
   def __init__(self, machine, repo_commits, pr_head_commit):
      """Create the log file object for a specific machine."""
      self.machine = machine.lower()
      self.pr_head_commit = pr_head_commit
      self.repo_commits = repo_commits
      self.current_pr_log_data = {}
      self.current_pr_runtime_data = {}
      self.current_pr_mem_data = {}
      self.text_per_log = {}
      # Runtime/memory data for each repo commit in self.repo_commits
      self.historical_runtime_data = {}
      self.historical_mem_data = {}
      # Runtime or memory mean & STDev for each test 
      self.runtime_stats = {}
      self.mem_stats = {}
      # Contains pass/warn/fail results for each test based on whether runtime or memory are too high
      self.runtime_results = {}
      self.memory_results = {}

   def get_pr_head(self):
      return self.pr_head_commit

   def get_repo_commits(self):
      return self.repo_commits

   def _fetch_log_text(self, list_of_commits): 
      """For each commit of a log, extract the log text. Store in a dictionary w/hash as key.
      Args:
         commits (list): list of commits for the repository - even the PR head commit is expected to be in list form. 
      """
      try:
         api_call = APICall(f"contents/tests/logs/RegressionTests_{self.machine}.log")
         for item_num in range(len(list_of_commits)): 
            url = api_call.url + (f"?ref={list_of_commits[item_num]}") #Could use a path join?
            r = requests.get(url, headers=api_call.header)
            if r.status_code != 200:
               logging.error(f"{r.status_code} - Commit {list_of_commits[item_num]} does not exist for this log. API call: {api_call.url}.")
               list_of_commits[item_num] = None
               sys.exit()
            else:
               self.text_per_log[list_of_commits[item_num]] = r.text
            
      except:
         logging.error("An appropriate commit(s) was not provided.")
         sys.exit()

   def get_log_text(self, list_of_commits):
      try:
         if list_of_commits == [self.pr_head_commit] and self.pr_head_commit[0] not in self.text_per_log:
            self._fetch_log_text(list_of_commits)
         if list_of_commits == self.repo_commits and self.repo_commits[0] not in self.text_per_log:
            self._fetch_log_text(list_of_commits)
         if not self.text_per_log:
            sys.exit(f"No log text has been fetched. Please check the provided commit(s): {list_of_commits}.")
         return self.text_per_log
      except:
         self._fetch_log_text(list_of_commits)

   def _get_instance_test_data(self, log_text):
      """For each instance of a log at a given commit, extract runtime and memory data from the log text
         Args:
            log_text: Log text for a given commit
         Returns: 
            tests_for_log_instance: A dictionary of tests (keys) with an array of total runtime and memory use as the value for each test
      """

      tests_for_log_instance = {}
      
      pattern = r"TEST \'(.*)\' \[\d+:\d+, (\d+):(\d+)\]\((\d+) MB\)"
      log_text = log_text.splitlines()
      
      for line in log_text:
         test_match = re.search(pattern, line)
         if test_match:
            test_name, hh, mm, mem = test_match.groups()
            total_minutes = int(hh) * 60 + int(mm)
            tests_for_log_instance[test_name] = [total_minutes, int(mem)]
      
      return tests_for_log_instance
      
   def compile_historical_log_data(self): # Could split for runtime, mem to make more maintainable
      """Create a dictionary of data with runtime and memory usage for each test over time. Structure:  
         historical_[runtime/mem]_test_data = {
            test: {hash: value, hash: value, ...}
         }
      """
   
      # Skip self.pr_head_commit because it is the log from the PR
      for hash in self.text_per_log:
         if hash != self.pr_head_commit:
            data = self._get_instance_test_data(self.text_per_log[hash])
            for test in data:
               try:
                  self.historical_runtime_data[test].update({hash: data[test][0]})
                  self.historical_mem_data[test].update({hash: data[test][1]})
               except KeyError: 
                  logging.info("Test key doesn't exist yet. Creating test key.")
                  self.historical_runtime_data[test] = {hash: data[test][0]}
                  self.historical_mem_data[test] = {hash: data[test][1]}
   
   def get_current_pr_data(self):
      """Extract runtime/memory data for the PR's most recent commit.
      Returns:
         pr_log_data: Dictionary of tests (keys) with an array of total runtime and memory use as the value for each test
      """

      try: 
         self.get_log_text([self.pr_head_commit]) # get/fetch_log_text expects a LIST of commit(s)
         self.current_pr_log_data = self._get_instance_test_data(self.text_per_log[self.pr_head_commit])
         return self.current_pr_log_data
      except:
         logging.error("Cannot fetch current PR data.")
         sys.exit()

   # Could change current_pr_log_data to pandas DataFrame to access by column
   def get_current_pr_runtime_data(self):
      if not self.current_pr_log_data:
         self.get_current_pr_data()
      for test in self.current_pr_log_data:
         self.current_pr_runtime_data[test] = self.current_pr_log_data[test][0]
      return self.current_pr_runtime_data

   # Could change current_pr_log_data to pandas DataFrame to access by column
   def get_current_pr_mem_data(self):
      if not self.current_pr_log_data: # pragma: no cover
         self.get_current_pr_data()
      for test in self.current_pr_log_data:
         self.current_pr_mem_data[test] = self.current_pr_log_data[test][1]
      return self.current_pr_mem_data

   def fetch_historical_data(self):
      """
      Extract runtime/memory data for the authoritative repository's last several commits.
      Actual number determined by Manager.get_hashes()
      """
      self.get_log_text(self.repo_commits)
      self.compile_historical_log_data()

   def get_historical_runtime_data(self):
      if not self.historical_runtime_data:
         self.fetch_historical_data()
      return self.historical_runtime_data
   
   def get_historical_mem_data(self):
      if not self.historical_mem_data: # pragma: no cover
         self.fetch_historical_data()
      return self.historical_mem_data
               
   def calculate_stats(self, data):
      """For each test, calculate the mean and standard deviation of memory and runtime.
      Args:
         data (dict): Dictionary structured {
            test1: {hash: value, hash: value, ...},
            test2: {hash: value, hash: value, ...},
         }
      """
      test_stats = {}
      
      for test in data:
         mean = round(np.mean(list(data[test].values())), 5)
         stdev = round(np.std(list(data[test].values())), 5)
         test_stats[test] = [mean, stdev]

      return test_stats

   def calculate_runtime_stats(self):
      if not self.historical_runtime_data:
         self.fetch_historical_data()
      self.runtime_stats = self.calculate_stats(self.historical_runtime_data)

   def get_runtime_stats(self):
      if not self.runtime_stats:
         self.calculate_runtime_stats()
      return self.runtime_stats 

   def calculate_mem_stats(self):
      if not self.historical_mem_data:
         self.fetch_historical_data()
      self.mem_stats = self.calculate_stats(self.historical_mem_data)

   def get_mem_stats(self):
      if not self.mem_stats:
         self.calculate_mem_stats()
      return self.mem_stats

   def _compare_results(self, category): 
      """Check results from previous three commits to determine whether the test runtime/memory usage 
      is within normal bounds."""

      results = {}

      if category == "runtime":
         stats = self.get_runtime_stats()
         data = self.get_historical_runtime_data()
         num = 0
      if category == "memory":
         stats = self.get_mem_stats()
         data = self.get_historical_mem_data()
         num = 1

      current_log = self.current_pr_log_data
      
      # Most recent commits: 
      last = self.repo_commits[0]
      second_to_last = self.repo_commits[1]
      third_to_last = self.repo_commits[2]
      recent_hashes = [last, second_to_last, third_to_last]

      previous_logs = {test: {hash: values[hash] for hash in recent_hashes if hash in values} 
                       for test, values in data.items()}

      

      for test in current_log:
         try:
            hi = stats[test][0] + (2 * stats[test][1])
            if current_log[test][num] > hi and previous_logs[test][last] > hi and previous_logs[test][second_to_last] > hi and previous_logs[test][third_to_last] > hi:
               results[test] = '❌'
            elif current_log[test][num] > hi:
               results[test] = '⚠️'
            else:
               results[test] = '✅'
         except KeyError:
            logging.info(f"{test} is new. No comparison data.")
            results[test] = 'New'

      return results

   def compare_runtimes(self):
      self.runtime_results = self._compare_results("runtime")

   def get_runtime_results(self):
      if not self.runtime_results:
         self.runtime_results = self._compare_results("runtime")
      return self.runtime_results 

   def compare_mem(self):
      self.mem_results = self._compare_results("memory")

   def get_mem_results(self):
      if not self.memory_results:
         self.memory_results = self._compare_results("memory")
      return self.memory_results
