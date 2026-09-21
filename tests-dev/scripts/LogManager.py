import os
from .Manager import *
from .Log import *

class LogManager(Manager):
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
      self.runtime_results_by_machine = {}
      self.mem_results_by_machine = {}
      self.current_pr_runtime_data = {}
      self.current_pr_mem_data = {}

   def update_log_manager_w_cached_data(self):
      self.historical_runtime = self.load_json_from_file(f"{os.environ.get('TEST_STATS')}/historical_runtime.json")
      self.historical_mem = self.load_json_from_file(f"{os.environ.get('TEST_STATS')}/historical_memory.json")
      # Contains mean and standard deviation of runtime/memory for each test on each machine
      self.runtime_stats_by_machine = self.load_json_from_file(f"{os.environ.get('TEST_STATS')}/runtime_stats.json")
      self.mem_stats_by_machine = self.load_json_from_file(f"{os.environ.get('TEST_STATS')}/memory_stats.json")

   def add_current_pr_data(self, machine, data, pr_data): 
      
      for test in list(pr_data.keys()):
         data.setdefault(machine, {}).update({test: pr_data[test]})
         
   def manage_data(self):

      for machine in self.machines:
         print(machine.upper())
         log = Log(machine, self.repo_hashes, self.pr_head_commit)

         if os.environ.get('TEST_STATS'):
            self.manage_preexisting_data(log, machine)
         else:
            self.collect_new_log_data(log, machine)
         
         # Add current PR data
         self.add_current_pr_data(machine, self.current_pr_runtime_data, log.get_current_pr_runtime_data())
         self.add_current_pr_data(machine, self.current_pr_mem_data, log.get_current_pr_mem_data())
         
         # Compare current results to historical values and save results (pass/warn/fail)
         self.runtime_results_by_machine[machine] = log.get_runtime_results()
         self.mem_results_by_machine[machine] = log.get_mem_results()

   def manage_preexisting_data(self,log,machine):
      """Populate the Log data with cached data and statistics on runtime/memory.
      Args:
         log (Log):
         machine (string): 
      """
      # Historical runtime/memory data for each repo commit in self.repo_commits
      log.historical_runtime_data = self.historical_runtime[machine]
      log.historical_mem_data = self.historical_mem[machine]
      # Runtime or memory mean & STDev for each test 
      log.runtime_stats = self.runtime_stats_by_machine[machine]
      log.mem_stats = self.mem_stats_by_machine[machine]

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
      # If the statistics on mean/standard deviation have NOT already been cached, create file to cache.
      if not os.environ.get('TEST_STATS'):
         self.create_json(self.runtime_stats_by_machine, "runtime_stats")
         self.create_json(self.mem_stats_by_machine, "memory_stats")
      
      # Create a record of historical runtime & memory values w/current PR data for caching 
      # (to use in plotting job and subsequent workflow runs)
      self.create_json(self.historical_runtime, "historical_runtime")
      self.create_json(self.historical_mem, "historical_memory")
      
      # Save current_pr_data for plotting task; do not mix with cached data that is reused w/newer PR commits
      self.create_json(self.current_pr_runtime_data, "current_pr_runtime_data")
      self.create_json(self.current_pr_mem_data, "current_pr_memory_data")

      # Create runtime/memory resource summaries to use in write_test_summary.py 
      self.create_json(self.runtime_results_by_machine, "runtime_results")
      self.create_json(self.mem_results_by_machine, "memory_results")
