import os, sys
import json
import logging
from .APICall import APICall

class Manager():
   """The Manager class contains utility functions used by the LogManager and PlotManager subclasses to 
   (1) create a JSON file from log data and (2) load data from a JSON file into a python dictionary. 
   """

   def __init__(self):
      self.repo_hashes = self.set_hashes()
      self.pr_head_commit = self.set_pr_head()
      self.machines = os.environ.get('MACHINES').lower().split()
      self.categories = ["runtime", "memory"]

   def set_hashes(self, num=30):
      """Retrieve the last "num" commit hashes from the repository.
      Args: 
         num (int): The number of commit hashes to be retrieved. 
      Returns:
         repo_hashes: list of commit hashes
      """
      try: 
         hashes = []
         api_call = APICall(f"commits?per_page={num}")
         response = api_call.call_API()
         response = api_call.load_json_from_api_call(response)
         
         for item in response:
            hashes.append(item['sha'][:8])
         return hashes
      except: # pragma: no cover
         logging.error(response['message'])
         sys.exit()

   def get_hashes(self):
      return self.repo_hashes
   
   def set_pr_head(self):
      """Get SHA for the HEAD of the PR. Structure of response: 
         response = [{"head": {"sha": "a1b2c3d..."}}]
         See GitHub documentation for https://docs.github.com/en/rest/commits/commits?apiVersion=2022-11-28#list-commits
      """
      try:
         api_call = APICall(f"pulls/{os.environ.get('PR_NUM')}")
         response = api_call.call_API()
         response = api_call.load_json_from_api_call(response)
         return response['head']['sha']
      except:
         logging.error(f"{response['status']} {response['message']}. URL: {api_call.url}")
         sys.exit()

   def get_pr_head(self): # pragma: no cover
      return self.pr_head_commit

   def set_machines(self, machine_list): # pragma: no cover
      self.machines = [m.lower() for m in machine_list]
      print(self.machines)

   # Utilities for file I/O & plotting
   def create_json(self, dictionary, file_name):
      """Create a json file, e.g., with statistics for each test on each machine"""

      with open(f"data/{file_name}.json", 'w') as fh:
         json.dump(dictionary, fh, indent=4)

   def load_json_from_file(self, file_path):
      """Convert JSON file to python dictionary."""
      with open(file_path, 'r', encoding='utf-8') as file:
         data = json.load(file)

      if not data: # pragma: no cover
         logging.error(f"No data retrieved.")

      return data
