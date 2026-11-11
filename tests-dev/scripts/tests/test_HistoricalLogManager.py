import pytest
import json
from pathlib import Path
from unittest.mock import Mock, call, patch
from scripts.APICall import APICall
from scripts.HistoricalLogManager import HistoricalLogManager
from scripts.get_data import *

@pytest.fixture
def historical_log_manager(monkeypatch):
   monkeypatch.setenv("MACHINES", "hercules ursa")
   monkeypatch.setattr(
      HistoricalLogManager, "set_hashes", lambda self: ["commit1", "commit2", "commit3", "commit4"]
   )
   monkeypatch.setattr(
      HistoricalLogManager, "set_pr_head", lambda self: "pr-head"
   )
   return HistoricalLogManager()

#@pytest.mark.parametrize("attribute", attributes)
def test_init_HistoricalLogManager(set_env_vars):
   """Test initialization of HistoricalLogManager."""
   set_env_vars
   manager = HistoricalLogManager()
   assert len(manager.repo_hashes) == 30
   assert "ursa" in manager.machines
   assert "hercules" in manager.machines
   assert manager.categories == ['runtime', 'memory']
   dicts = ['historical_runtime', 'historical_mem', 'runtime_stats_by_machine', 'mem_stats_by_machine', ]
   for item in dicts:
      assert not getattr(manager, item)

def test_manage_data(historical_log_manager):
   """Checks whether every configured machine gets a log with the right commits 
   and is passed to collect_new_log_data()."""

   logs = [Mock(spec_set=HistoricalLog), Mock(spec_set=HistoricalLog)]

   # Replace HistoricalLog() construction with pre-created mock logs.
   # Replace collect_new_log_data() with a mock that records calls without collecting data.
   with patch(
      "scripts.HistoricalLogManager.HistoricalLog",
      side_effect=logs,
   ) as make_log, patch.object(historical_log_manager, "collect_new_log_data") as collect:
      historical_log_manager.manage_data()

   # Check that when we call manage_data(), a log was requested (from make_log) for each machine with the correct commits.
   assert make_log.call_args_list == [
      call("hercules", historical_log_manager.repo_hashes),
      call("ursa", historical_log_manager.repo_hashes),
   ]
   # Check that when we call manage_data(), each mock log is passed to collect (instead of collect_new_log_data()).
   assert collect.call_args_list == [
      call(logs[0]),
      call(logs[1]),
   ]

def test_collect_new_log_data(historical_log_manager):
   log = HistoricalLog("hercules", historical_log_manager.repo_hashes)
   log.historical_runtime_data = {
      "test1": {"commit1": 1, "commit2": 7, "commit3": 9, "commit4": 15},
   }
   log.historical_mem_data = {
      "test1": {"commit1": 110, "commit2": 116, "commit3": 124, "commit4": 134},
   }

   historical_log_manager.collect_new_log_data(log)

   assert historical_log_manager.historical_runtime == {
      "hercules": {"test1": {"commit1": 1, "commit2": 7, "commit3": 9, "commit4": 15}},
   }
   assert historical_log_manager.historical_mem == {
      "hercules": {"test1": {"commit1": 110, "commit2": 116, "commit3": 124, "commit4": 134}}
   }
   assert historical_log_manager.runtime_stats_by_machine == {
      "hercules": {"test1": [8.0, 5.0]},
   }
   assert historical_log_manager.mem_stats_by_machine == {
      "hercules": {"test1": [121.0, 9.0]},
   }

def test_collect_new_log_data_w_real_data(set_env_vars):
   set_env_vars
   manager = HistoricalLogManager()
   log = HistoricalLog("ursa", manager.repo_hashes)
   manager.collect_new_log_data(log)
   for test_info in [manager.historical_runtime["ursa"], manager.historical_mem["ursa"], 
                manager.runtime_stats_by_machine["ursa"], manager.mem_stats_by_machine["ursa"]]:
      assert 'control_c48_gnu' in test_info

def test_save_data(set_env_vars):
   set_env_vars
   manager = HistoricalLogManager()
   assert False

def test_save_data(historical_log_manager, tmp_path, monkeypatch):
   monkeypatch.chdir(tmp_path)
   (tmp_path / "data").mkdir()

   historical_log_manager.runtime_stats_by_machine = {
      "hercules": {"test1": [12.0, 2.0]},
   }
   historical_log_manager.mem_stats_by_machine = {
      "hercules": {"test1": [110.0, 10.0]},
   }
   historical_log_manager.historical_runtime = {
      "hercules": {"test1": {"commit1": 10, "commit2": 14}},
   }
   historical_log_manager.historical_mem = {
      "hercules": {"test1": {"commit1": 100, "commit2": 120}},
   }

   historical_log_manager.save_data()

   expected = {
      "runtime_stats.json": historical_log_manager.runtime_stats_by_machine,
      "memory_stats.json": historical_log_manager.mem_stats_by_machine,
      "historical_runtime.json": historical_log_manager.historical_runtime,
      "historical_memory.json": historical_log_manager.historical_mem,
   }

   assert {path.name for path in (tmp_path / "data").iterdir()} == set(expected)
   for filename, contents in expected.items():
      actual = json.loads((tmp_path / "data" / filename).read_text())
      assert actual == contents
