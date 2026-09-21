import pytest
from scripts.Manager import *
from scripts.get_data import *

def test_main_e2e_cached_stats(monkeypatch, set_env_vars):
   """Test that main function runs to completion."""

   set_env_vars
   monkeypatch.setenv("MACHINES", "hercules")
   monkeypatch.setenv("TEST_STATS", "data")
   exit_code = main()

   assert exit_code == 0

def test_main_e2e_no_cached_stats(monkeypatch, set_env_vars):
   """Test that main function runs to completion."""

   set_env_vars
   monkeypatch.setenv("MACHINES", "hercules")
   exit_code = main()

   assert exit_code == 0
