import pytest
import os
import json
from pathlib import Path
from scripts.APICall import APICall
from scripts.LogManager import LogManager
from scripts.get_data import *

#@pytest.mark.parametrize("attribute", attributes)
def test_init_LogManager(set_env_vars):
   set_env_vars
   manager = LogManager()
   assert len(manager.repo_hashes) == 30
   assert len(manager.pr_head_commit) == 40
   assert "ursa" in manager.machines
   assert manager.categories == ['runtime', 'memory']
   dicts = ['historical_runtime', 'historical_mem', 'runtime_stats_by_machine', 'mem_stats_by_machine', 'runtime_results_by_machine', 'mem_results_by_machine', 'current_pr_runtime_data', 'current_pr_mem_data']
   for item in dicts:
      assert not getattr(manager, item)





def test_update_log_manager_w_cached_data(set_env_vars):
   set_env_vars
   manager = LogManager()

def test_add_current_pr_data(set_env_vars):
   set_env_vars
   manager = LogManager()

def test_manage_data(set_env_vars):
   set_env_vars
   manager = LogManager()

def test_manage_preexisting_data(set_env_vars):
   set_env_vars
   manager = LogManager()

def test_collect_new_log_data(set_env_vars):
   set_env_vars
   manager = LogManager()

def test_save_data(set_env_vars):
   set_env_vars
   manager = LogManager()
