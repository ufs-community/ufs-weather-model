import pytest
import os
import json
from pathlib import Path
from scripts.APICall import APICall
from scripts.get_data import *


def test_init_manager(set_env_vars):
   """Test initialization of the manager, including the API call and its ability to get PR 2882's head commit. 
   set_hashes() and set_pr_head_commit() are tested implicitly by initiating the manager, which calls set_hashes() 
   and set_pr_head_commit() in its __init__ function.
   When running tests locally, create a GitHub token and set it as an environment variable.
   """
   set_env_vars
   manager = Manager()
      
   assert manager.machines == os.environ.get('MACHINES').split()
   assert manager.categories == ["runtime", "memory"]
   assert len(manager.get_hashes()) == 30
   assert manager.pr_head_commit == "369cead91c98eb5c72da81ff78925250dad08903"

def test_set_pr_head_error(set_env_vars, monkeypatch, caplog):

   set_env_vars
   monkeypatch.setenv("PR_NUM", "-1")
   with pytest.raises(SystemExit) as error:
      manager = Manager()

   assert caplog.messages[0] == "404 Not Found. URL: https://api.github.com/repos/ufs-community/ufs-weather-model/pulls/-1"

def test_set_machines(monkeypatch):
   machines = "Larry Moe Curly"
   monkeypatch.setenv("MACHINES", machines)
   manager = Manager()
   assert manager.machines == ["larry", "moe", "curly"]

def test_create_json():
   """Create a json file, e.g., with statistics for each test on each machine"""
   
   manager = Manager()
   
   data_to_save = {
      "machine1": {
         "test1": [101, 5],
         "test2": [111, 0],
         "test3": [234, 4.56],
         "test4": [222, 9.67],
         "test5": [999, 7.6],
      },
      "machine2": {
         "test1": [95.5, 6.1],
         "test6": [43.5, 3.2],
         "test3": [8745, 55.7],
         "test4": [249, 10],
         "test7": [558, 17.8],
      },
      "machine3": {
         "test8": [0, 1.2],
         "test2": [117, 3.4],
         "test3": [8976, 131.9],
         "test4": [231, 10.6],
         "test5": [987654321, 2345],
         "test6": [43.5, 3.2],
         "test7": [567.34, 21.3],
      },
   }

   path = Path('data')
   path.mkdir(exist_ok = True)
   manager.create_json(data_to_save, "actual_json_data")

   with open('data/actual_json_data.json', 'r') as actual_data, open ('data/expected_json_data.json', 'r') as expected_data:
      actual_data_contents = actual_data.read()
      expected_data_contents = expected_data.read()

   assert actual_data_contents == expected_data_contents

def test_load_json_from_file(set_env_vars):
   """Convert JSON file to python dictionary."""
   set_env_vars
   manager = Manager()
   actual_data = manager.load_json_from_file("data/test_load_json_from_file.json")
   expected_data = {
      "cat": {
         "kitten1": [1, 2, 3],
         "kitten2": ["meow1", "meow2", "meow3"],
         "kitten3": [True, False],
         "kitten4": [222, 9.67],
         "kitten5": [999, 7.6]
      },
      "dog": {
         "puppy1": ["woof", "arf", "bark"],
         "puppy6": [1, 2],
         "puppy3": [8745, 55.7],
         "puppy4": [249, 10],
         "puppy7": [True, False, True],
         "puppy5": [0, 1.2]
      }
   }

   assert actual_data == expected_data
