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
         logging.error(f"{response}: API call failed for {api_call.url}")
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
         logging.error(f"No commit found for the ref {commit}")
         sys.exit(1)

   def _get_data(self, log_text, pattern, handler, to_clean=True):
      """Extract data on warnings, remarks, and failed compiles/tests
      Args:
         log_text: Text from a log file at a given commit
         pattern: regex pattern to match
         handler: function indicating how to handle the incoming data
         to_clean: Flag indicating whether the data needs to be cleaned (True) or not (False)
      """

      data = {}
      
      log_text = log_text.splitlines()

      for line in log_text:
         test_match = re.search(pattern, line)
         if test_match:
            handler(data, *test_match.groups())
         
      if to_clean:
         data = self._clean_data(data)
      
      return data

   def handle_warn_rmk(self, data_dict, test_name, warnings, remarks):
      data_dict.update({test_name: (warnings, remarks)})

   def handle_failures(self, data_dict, reason, test):
      data_dict.setdefault(reason, []).append(test)

   def _clean_data(self, test_data):
      """Convert None values to zeros in the test_data dictionary"""
      clean_data = {
         k: tuple(0 if v is None else int(v) for v in values) 
         for k, values in test_data.items()
      }
      return clean_data

   def compare_results(self, pr_log, base_log): 
      """Compare warnings/remarks for PR head and base commits to determine whether warnings/remarks have increased."""

      increases = {'warnings': {}, 'remarks': {}}

      for test in pr_log:
         if test not in base_log:
            logging.info(f"Skipped test {test}; nothing to compare against.")
            continue
         # Check warnings
         if pr_log[test][0] > base_log[test][0]:
            increases['warnings'].update({test: pr_log[test][0] - base_log[test][0]})
         # Check remarks
         if pr_log[test][1] > base_log[test][1]:
            increases['remarks'].update({test: pr_log[test][1] - base_log[test][1]})
      
      return increases

def print_html_results(data_dict, title, file_name):
   """Print results in HTML format.

   Args:
      data_dict: Dictionary with structure {machine: {category: tests, ...}, ...}
                 where values can be either:
                 - A dict with {test: count, ...} (for warnings/remarks)
                 - A list of tests (for failures)
      title: Title for the markdown file
      file_name: Name for the markdown file
   
   Returns:
      Formatted markdown string
   """
   mdFile = MdUtils(file_name=file_name, title=title)

   for machine, machine_data in data_dict.items():
      # Skip if there is no data for the machine (cases: (1) where all RTs pass on the machine or (2) there is no increase in warnings & remarks)
      if not machine_data:
         continue
      
      # Skip if all categories (warnings/remarks/failure reasons) are empty
      if all(not category for category in machine_data.values()):
         continue
      
      mdFile.write(f"\n<h3>{machine.upper()}</h3>\n")
      
      for category, tests in machine_data.items():
         # Skip printing info if there are no tests listed for a given category
         if not tests:
            continue
         
         unordered_list = [f"**{category.title()}:**", []]
         
         # Handle both dict (warnings/remarks) and list (failures)
         if isinstance(tests, dict):
            for test, count in tests.items():
               unordered_list[1].append(f"{test}: {count}")
         else:  # list
            for test in tests:
               unordered_list[1].append(test)
         
         mdFile.new_list(unordered_list, marked_with='*')
   
   return mdFile.get_md_text()

def main():
   """For each machine, create a log object, get current PR data, and determine 
   which tests increase warnings and/or remarks on each machine.""" 

   machines = os.environ.get('MACHINES').split()
   # Use non-capturing groups in pattern to indicate warnings/remarks may or may not be present.
   compile_pattern = r"COMPILE \'(.*)\' \[\d+:\d+, \d+:\d+\](?: \( (?:(\d+) warnings)?\s*(?:(\d+) remarks)? \))?"
   # failure_pattern = r"^(?:FAILED|SKIPPED): (?!UNABLE TO (?:COMPLETE COMPARISON|START TEST))(.+?) -- (?:TEST|COMPILE) '([^']+)"
   failure_pattern = r"^(?:FAILED|SKIPPED): (.+?) -- (?:TEST|COMPILE) '([^']+)"
      

   # For each machine, increased_warnings_remarks records where warnings and/or remarks increase
   increased_warnings_remarks = {}
   # For each machine, failures records which tests fail by reason
   failures = {}

   for machine in machines:
      log = Log(machine)
      log._get_commits()
      log.pr_log_text = log._fetch_log_text(log.pr_head_commit)
      log.pr_warn_rmk = log._get_data(log.pr_log_text, compile_pattern, log.handle_warn_rmk)
      log.pr_failures = log._get_data(log.pr_log_text, failure_pattern, log.handle_failures, to_clean=False)

      log.base_log_text = log._fetch_log_text(log.pr_base_commit)
      log.base_warn_rmk = log._get_data(log.base_log_text, compile_pattern, log.handle_warn_rmk)
      
      increased_warnings_remarks[machine] = log.compare_results(log.pr_warn_rmk, log.base_warn_rmk)
      failures[machine] = log.pr_failures

   pr_num = os.environ.get('PR_NUM')
   warn_rmk_results = print_html_results(increased_warnings_remarks, 
                                     f"Increased Warnings/Remarks for PR #{pr_num}",
                                     "warn_rmk.md")
   failure_results = print_html_results(failures, 
                                    f"Compile and Test Failures for PR #{pr_num}", 
                                    "failures.md")

   if len(warn_rmk_results.splitlines()) > 3: # HTML header is 3 lines long
      print(warn_rmk_results)
      if len(failure_results.splitlines()) > 3:
         print(failure_results)
      sys.exit(1)
   elif len(failure_results.splitlines()) > 3:
      print(failure_results)
      sys.exit(0)
   else:
      print(f"No increase in warnings or remarks. All RTs passed.")
      sys.exit(0)

if __name__ == "__main__": # pragma: no coverage

   main()
