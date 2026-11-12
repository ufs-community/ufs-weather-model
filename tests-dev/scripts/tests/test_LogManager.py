import json
from unittest.mock import Mock, call, patch

import pytest

from scripts.HistoricalLogManager import HistoricalLogManager
from scripts.Log import Log
from scripts.LogManager import LogManager


@pytest.fixture
def log_manager(monkeypatch):
   monkeypatch.setenv("MACHINES", "hercules ursa")
   monkeypatch.setattr(LogManager, "set_hashes", lambda self: ["commit1", "commit2"])
   monkeypatch.setattr(LogManager, "set_pr_head", lambda self: "pr-head")
   return LogManager()

@pytest.fixture
def cached_data(log_manager):
   return {
      "historical_runtime": {
         machine: {"test1": {"commit1": 10 + index, "commit2": 14 + index},
                   "test2": {"commit1": 20 + index, "commit2": 24 + index}}
         for index, machine in enumerate(log_manager.machines)
      },
      "historical_memory": {
         machine: {"test1": {"commit1": 100 + index, "commit2": 120 + index},
                   "test2": {"commit1": 200 + index, "commit2": 220 + index}}
         for index, machine in enumerate(log_manager.machines)
      },
      "runtime_stats": {
         machine: {"test1": [12.0 + index, 2.0],
                   "test2": [22.0 + index, 2.0]}
         for index, machine in enumerate(log_manager.machines)
      },
      "memory_stats": {
         machine: {"test1": [110.0 + index, 10.0],
                   "test2": [210.0 + index, 10.0]}
         for index, machine in enumerate(log_manager.machines)
      },
   }

def test_init_LogManager(log_manager):
   """Verify correct initialization of LogManager."""
   assert log_manager.repo_hashes == ["commit1", "commit2"]
   assert log_manager.pr_head_commit == "pr-head"
   assert log_manager.machines == ["hercules", "ursa"]
   assert log_manager.categories == ["runtime", "memory"]
   for attribute in [
      "historical_runtime", "historical_mem",
      "runtime_stats_by_machine", "mem_stats_by_machine",
      "runtime_results_by_machine", "mem_results_by_machine",
      "current_pr_runtime_data", "current_pr_mem_data",
   ]:
      assert getattr(log_manager, attribute) == {}

def test_update_log_manager_w_cached_data(log_manager, cached_data, tmp_path, monkeypatch):
   """Verify that the LogManager is correctly updated with cached_data."""

   # Write cached data to a temporary directory.
   monkeypatch.setenv("TEST_STATS", str(tmp_path))
   for file, contents in cached_data.items():
      (tmp_path / f"{file}.json").write_text(json.dumps(contents))

   log_manager.update_log_manager_w_cached_data()

   assert log_manager.historical_runtime == cached_data["historical_runtime"]
   assert log_manager.historical_mem == cached_data["historical_memory"]
   assert log_manager.runtime_stats_by_machine == cached_data["runtime_stats"]
   assert log_manager.mem_stats_by_machine == cached_data["memory_stats"]

def test_collect_historical_log_data(log_manager, cached_data):
   """Verify that historical data is collected and correctly transferred to the LogManager."""

   historical_manager = Mock(spec=HistoricalLogManager)
   # Expected data (from cached_data)
   historical_manager.historical_runtime = cached_data["historical_runtime"]
   historical_manager.historical_mem = cached_data["historical_memory"]
   historical_manager.runtime_stats_by_machine = cached_data["runtime_stats"]
   historical_manager.mem_stats_by_machine = cached_data["memory_stats"]

   # Replace HistoricalLogManager() construction with pre-created mock historical log manager.
   with patch("scripts.LogManager.HistoricalLogManager", return_value=historical_manager):
      log_manager.collect_historical_log_data()

   # Is retrieved data the same as the cached data?
   assert log_manager.historical_runtime == cached_data["historical_runtime"]
   assert log_manager.historical_mem == cached_data["historical_memory"]
   assert log_manager.runtime_stats_by_machine == cached_data["runtime_stats"]
   assert log_manager.mem_stats_by_machine == cached_data["memory_stats"]


@pytest.mark.parametrize("attribute", ["current_pr_runtime_data", "current_pr_mem_data"])
def test_add_current_pr_data(log_manager, attribute):
   """Verify the log_manager receives current PR values without losing other tests or machines."""
   data = getattr(log_manager, attribute)
   data.update({"hercules": {"existing": 10, "updated": 20}, "ursa": {"other": 30}})
   pr_data = {"updated": 25, "new": 40}

   log_manager.add_current_pr_data("hercules", data, pr_data)
   log_manager.add_current_pr_data("wcoss2", data, {"test1": 50})

   assert data == {
      "hercules": {"existing": 10, "updated": 25, "new": 40},
      "ursa": {"other": 30},
      "wcoss2": {"test1": 50},
   }
   assert pr_data == {"updated": 25, "new": 40}

@pytest.mark.parametrize("machine", ["hercules", "ursa", "HeRcuLes"])
def test_manage_historical_data(log_manager, cached_data, machine):
   """Verify that cached/historical data is correctly transferred from the LogManager to the Log."""

   log_manager.historical_runtime = cached_data["historical_runtime"]
   log_manager.historical_mem = cached_data["historical_memory"]
   log_manager.runtime_stats_by_machine = cached_data["runtime_stats"]
   log_manager.mem_stats_by_machine = cached_data["memory_stats"]
   log = Log(machine, log_manager.repo_hashes, log_manager.pr_head_commit)

   log_manager.manage_historical_data(log)

   assert log.historical_runtime_data == cached_data["historical_runtime"][machine.lower()]
   assert log.historical_mem_data == cached_data["historical_memory"][machine.lower()]
   assert log.runtime_stats == cached_data["runtime_stats"][machine.lower()]
   assert log.mem_stats == cached_data["memory_stats"][machine.lower()]

def test_manage_data(log_manager):
   """Verify each machine's log is processed and its PR values and results stored."""
   logs = [Mock(spec=Log), Mock(spec=Log)]

   # Fill logs with dummy values to return when get_current_pr_*_data() or get_*_results() is called.
   for index, log in enumerate(logs):
      log.get_current_pr_runtime_data.return_value = {"test1": 15 + index}
      log.get_current_pr_mem_data.return_value = {"test1": 130 + index}
      log.get_runtime_results.return_value = {"test1": f"runtime-result-{index}"}
      log.get_mem_results.return_value = {"test1": f"memory-result-{index}"}

   # When manage_data() is called, replace Log creation with mock logs and
   # replace manage_historical_data() call with mock call (manage_history)
   with patch(
      "scripts.LogManager.Log", side_effect=logs
   ) as make_log, patch.object(log_manager, "manage_historical_data") as manage_history:
      log_manager.manage_data()

   assert make_log.call_args_list == [
      call(machine, log_manager.repo_hashes, log_manager.pr_head_commit)
      for machine in log_manager.machines
   ]
   assert manage_history.call_args_list == [call(log) for log in logs]

   assert log_manager.current_pr_runtime_data == {
      "hercules": {"test1": 15}, "ursa": {"test1": 16},
   }
   assert log_manager.current_pr_mem_data == {
      "hercules": {"test1": 130}, "ursa": {"test1": 131},
   }
   assert log_manager.runtime_results_by_machine == {
      "hercules": {"test1": "runtime-result-0"},
      "ursa": {"test1": "runtime-result-1"},
   }
   assert log_manager.mem_results_by_machine == {
      "hercules": {"test1": "memory-result-0"},
      "ursa": {"test1": "memory-result-1"},
   }

def test_save_pr_data(log_manager, tmp_path, monkeypatch):
   """Verify PR data and comparison results are saved to the correct JSON files."""

   # Use a temporary working directory to keep output files out of the repository
   monkeypatch.chdir(tmp_path)
   (tmp_path / "data").mkdir()
   log_manager.current_pr_runtime_data = {"hercules": {"test1": 15}}
   log_manager.current_pr_mem_data = {"hercules": {"test1": 130}}
   log_manager.runtime_results_by_machine = {"hercules": {"test1": "runtime-result"}}
   log_manager.mem_results_by_machine = {"hercules": {"test1": "memory-result"}}

   log_manager.save_pr_data()

   expected = {
      "current_pr_runtime_data.json": log_manager.current_pr_runtime_data,
      "current_pr_memory_data.json": log_manager.current_pr_mem_data,
      "runtime_results.json": log_manager.runtime_results_by_machine,
      "memory_results.json": log_manager.mem_results_by_machine,
   }

   # Compares the set of filenames in the temporary data/ directory with the filenames expected by the test.
   assert {path.name for path in (tmp_path / "data").iterdir()} == set(expected)

   # Compares the contents of a file in the temporary data/ dir with the contents expected by the test.
   for filename, contents in expected.items():
      assert json.loads((tmp_path / "data" / filename).read_text()) == contents
