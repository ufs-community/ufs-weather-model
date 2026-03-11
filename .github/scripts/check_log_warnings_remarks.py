import requests
from mdutils.mdutils import MdUtils
import os, sys
import json
import re
import logging

class APICall():
   """A GitHub API call"""

   def __init__(self, endpoint='', num_commits=1):
      self.token = os.environ.get('GITHUB_TOKEN')
      self.base_url = os.environ.get('BASE_URL')
      self.endpoint = endpoint
      self.url = f"{self.base_url}/{self.endpoint}" #Could use a path join?
      self.num_commits = num_commits
      self.header = {
         "Accept": "application/vnd.github.v3+json",
         "Authorization": f"Bearer {self.token}",
         "X-GitHub-Api-Version": "2022-11-28",
         "Accept": "application/vnd.github.raw"
      }

class Log():
   """A Regression Test log file."""
   
   def __init__(self, machine):
      """Create the log file object for a specific machine."""
      self.machine = machine.lower()
      self.text_per_log = []

   def call_API(self, endpoint):
      """Call the GitHub API to get information about the log file."""

      api_call = APICall(endpoint)
      response = requests.get(api_call.url, headers=api_call.header)
      if response.status_code != 200:
         logging.warning(response)
         print(response)
         sys.exit(1)
      response = json.loads(response.text)
      
      return response

   def _get_commits(self):
      """Get PR head and base commits. Structure of response: 
         response = [{"head": {"sha": "a1b2c3d..."}, "base": {"sha": "b2c3d4e..."}}]
         See GitHub documentation for https://docs.github.com/en/rest/commits/commits?apiVersion=2022-11-28#list-commits
      """
      response = self.call_API(f"pulls/{os.environ.get('PR_NUM')}")
      self.pr_head_commit = response['head']['sha']
      self.pr_base_commit = response['base']['sha']

   def _fetch_log_text(self, commit): 
      """For each commit of a log, extract the log text."""

      try:
         api_call = APICall(f"contents/tests/logs/RegressionTests_{self.machine}.log")

         url = api_call.url + (f"?ref={commit}") #Could use a path join?
         r = requests.get(url, headers=api_call.header)
         return r.text
      except:
         logging.error("An appropriate commit(s) was not provided. Call _get_commits() first.")

   def _get_test_data(self, log_instance):
      """For each instance of a log at a given commit, extract runtime and memory data from the log text
         Args:
            log_instance: Log text for a given commit
         Returns: 
            tests_for_log_instance: A dictionary of tests (keys) with a tuple of warnings and remarks as the value for each test
      """

      tests_for_log_instance = {}

      pattern = r"COMPILE \'(.*)\' \[\d+:\d+, \d+:\d+\] \( (\d+) warnings (\d+) remarks \)"
      log_instance = log_instance.splitlines()

      for line in log_instance:
         test_match = re.search(pattern, line)
         if test_match:
            test_name, warnings, remarks = test_match.groups()
            tests_for_log_instance[test_name] = (int(warnings), int(remarks))

      return tests_for_log_instance

   def _get_pr_data(self, commit):
      """Extract warnings/remarks data for a particular commit.
      Returns:
         log_data: A dictionary of tests as the key with a tuple of (warnings, remarks) as the value
      """
      try:
         log_text = self._fetch_log_text(commit)
         log_data = self._get_test_data(log_text)
         return log_data
      except:
         logging.error(f"No commit found for the ref {commit}")
         sys.exit(1)

   def compare_results(self, pr_log, base_log): 
      """Compare warnings/remarks for PR head and base commits to determine whether warnings/remarks have increased."""

      increases = {'warnings': {}, 'remarks': {}}

      for test in pr_log:
         # Check warnings
         if pr_log[test][0] > base_log[test][0]:
            increases['warnings'].update({test: pr_log[test][0] - base_log[test][0]})
         # Check remarks
         if pr_log[test][1] > base_log[test][1]:
            increases['remarks'].update({test: pr_log[test][1] - base_log[test][1]})
      
      return increases

def print_html_results(dict):
   """Print the comparison results in HTML."""
   
   pr_num = os.environ.get('PR_NUM')
   mdFile = MdUtils(file_name='summary.md', title=f'Increased Warnings/Remarks for PR #{pr_num}')

   for machine, results in dict.items():
      for category in results.keys():
         if results[category]:
            mdFile.write(f"\n<h3>{machine.upper()}</h3>\n")
            unordered_list = [f"**{category.title()}:**", []]
            for test, value in dict[machine][category].items():
               unordered_list[1].append(f"{test}: {value}")
            mdFile.new_list(unordered_list, marked_with='*')
   return mdFile.get_md_text()

def main():
   """For each machine, create a log object, get current PR data, and determine 
   which tests increase warnings and/or remarks on each machine.""" 

   machines = os.environ.get('MACHINES').split()

   # For each machine, tests where warnings and/or remarks increase
   increased_warnings_remarks = {}

   for machine in machines:
      log = Log(machine)
      log._get_commits()
      log.pr_log_data = log._get_pr_data(log.pr_head_commit)
      log.base_log_data = log._get_pr_data(log.pr_base_commit)

      increased_warnings_remarks[machine] = log.compare_results(log.pr_log_data, log.base_log_data)

   results = print_html_results(increased_warnings_remarks)

   if len(results) > 81: # Length of HTML header
      print(results)
      sys.exit(1)
   else:
      sys.exit(0)

if __name__ == "__main__": # pragma: no coverage

   main()
