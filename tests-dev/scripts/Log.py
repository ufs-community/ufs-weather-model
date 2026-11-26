import requests
import sys
import re
import logging
from .APICall import APICall

class Log():
   """A Regression Test log file."""
   
   def __init__(self, machine, repo_commits, pr_head_commit):
      """Initialize current-PR and historical data for one machine.

      Args:
         machine (str): Machine name whose regression log is evaluated.
         repo_commits (list[str]): Historical commits ordered newest to oldest.
         pr_head_commit (str): Commit identifier for the pull request head.

      Returns:
         None.
      """
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
      """Return the pull request head commit hash.

      Returns:
         str: Pull request head commit.
      """
      return self.pr_head_commit

   def get_repo_commits(self):
      """Return the historical repository commit hashes.

      Returns:
         list[str]: Commit hashes in configured order.
      """
      return self.repo_commits

   def _fetch_log_text(self, list_of_commits): 
      """For each commit of a log, extract the log text. Store in a dictionary w/hash as key.
      Args:
         list_of_commits (list[str]): Commit hashes to fetch; pass even a single PR-head commit 
            as a one-item list.

      Returns:
         None. Stores raw log text in ``text_per_log`` keyed by commit hash.
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
      """Fetch requested log text when needed and return loaded log contents.

      Args:
         list_of_commits (list[str]): Commit hashes to fetch.

      Returns:
         dict[str, str]: Mapping from commit hash to raw regression test log text.
      """
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
      """For each instance of a log at a given commit, extract runtime and memory data from the log text.
         Args:
            log_text (str): Raw regression test log for one commit.

         Returns: 
            dict[str, list[int]]: ``{test_name: [runtime_minutes, memory_megabytes]}``.
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
   
   def get_current_pr_data(self):
      """Extract runtime/memory data for the PR's most recent commit.

      Args:
         None.

      Returns:
         dict[str, list[int]]: ``{test_name: [runtime_minutes, memory_megabytes]}``.
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
      """Return the PR-head runtime for each test, extracting the log if needed.

      Returns:
         dict[str, int]: ``{test_name: runtime_minutes}``.
      """
      if not self.current_pr_log_data:
         self.get_current_pr_data()
      for test in self.current_pr_log_data:
         self.current_pr_runtime_data[test] = self.current_pr_log_data[test][0]
      return self.current_pr_runtime_data

   # Could change current_pr_log_data to pandas DataFrame to access by column
   def get_current_pr_mem_data(self):
      """Return the PR-head memory use for each test, extracting the log if needed.

      Returns:
         dict[str, int]: ``{test_name: memory_megabytes}``.
      """
      if not self.current_pr_log_data: # pragma: no cover
         self.get_current_pr_data()
      for test in self.current_pr_log_data:
         self.current_pr_mem_data[test] = self.current_pr_log_data[test][1]
      return self.current_pr_mem_data

   def _compare_results(self, category): 
      """Classify current results against historical-mean-plus-two-standard-deviation thresholds.

      Args:
         category (str): ``"runtime"`` or ``"memory"``; selects the corresponding
            measurements and statistics.

      Returns:
         dict[str, str]: ``{test_name: status}``, where status is a pass, warning,
         failure, or ``"New"`` when no comparison statistics exist.
      """

      results = {}

      if category == "runtime":
         stats = self.runtime_stats
         data = self.historical_runtime_data
         num = 0
      elif category == "memory":
         stats = self.mem_stats
         data = self.historical_mem_data
         num = 1
      else:
         logging.error(f"{category} does not exist!")
         sys.exit()

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
      """Calculate and store classifications for current PR runtimes.

      Returns:
         None. Stores results in ``runtime_results``.
      """
      self.runtime_results = self._compare_results("runtime")

   def get_runtime_results(self):
      """Return cached runtime classifications, computing them when absent.

      Returns:
         dict[str, str]: ``{test_name: status}`` using pass/warning/fail symbols
         or ``"New"`` for tests without historical statistics.
      """
      if not self.runtime_results:
         self.runtime_results = self._compare_results("runtime")
      return self.runtime_results 

   def compare_mem(self):
      """Calculate and store classifications (pass/warn/fail) for current PR memory usage.

      Returns:
         None. Stores results in ``mem_results``.
      """
      self.mem_results = self._compare_results("memory")

   def get_mem_results(self):
      """Return cached memory classifications (pass/warn/fail), computing them when absent.

      Returns:
         dict[str, str]: ``{test_name: status}`` using pass/warning/fail symbols
         or ``"New"`` for tests without historical statistics.
      """
      if not self.memory_results:
         self.memory_results = self._compare_results("memory")
      return self.memory_results
