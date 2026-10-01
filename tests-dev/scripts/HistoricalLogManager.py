import os
from .Manager import *
from .HistoricalLog import *

"""Caches historical log data and statistics for use in the resource check workflows."""

class HistoricalLogManager(Manager):
   """Manages log objects and information common to all logs."""

   def __init__(self):
      super().__init__() # set hashes, machine list, categories

      # Contains runtime/memory data by machine for the last X number of commits (set in get_hashes() )
      self.historical_runtime = {}
      self.historical_mem = {}
      # Contains mean and standard deviation of runtime/memory for each test on each machine over the past X commmits
      self.runtime_stats_by_machine = {}
      self.mem_stats_by_machine = {}
      # Contains information on whether test runtime/memory was more than 2 standard deviations above the mean. 

   def manage_data(self):

      for machine in self.machines:
         print(f"Fetching historical data for {machine.upper()}.")
         log = HistoricalLog(machine, self.repo_hashes)

         self.collect_new_log_data(log, machine)

   def collect_new_log_data(self,log,machine):
      """Download and process log data for a given machine and update log with that information; calculate runtime/memory statistics.
      Args:
         log (Log):
         machine (string):
      """
      self.historical_runtime[machine] = log.get_historical_runtime_data()
      self.historical_mem[machine] = log.get_historical_mem_data()
      self.runtime_stats_by_machine[machine] = log.get_runtime_stats() # Add stats to save/cache later
      self.mem_stats_by_machine[machine] = log.get_mem_stats() # Add stats to save/cache later
   
   def save_data(self):
      # Cache statistics on mean/standard deviation
      self.create_json(self.runtime_stats_by_machine, "runtime_stats")
      self.create_json(self.mem_stats_by_machine, "memory_stats")
      
      # Cache a record of historical runtime & memory values
      # (to use in plotting job and subsequent workflow runs)
      self.create_json(self.historical_runtime, "historical_runtime")
      self.create_json(self.historical_mem, "historical_memory")
