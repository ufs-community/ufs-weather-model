import os
from .Manager import *
from .LogManager import *

"""Entry point for processing regression test log data and saving pull request results.
Uses cached data when TEST_STATS is set; otherwise, collects historical log data.
"""

def get_data():
   """Initialize a LogManager and process log data.

   Uses cached data if TEST_STATS environment variable is set, otherwise collects 
   historical log data. Manages the data and saves PR-related results.

   Returns:
      int: Zero after the manager has processed and saved data.
   """

   log_manager = LogManager()
      
   if os.environ.get('TEST_STATS'):
      log_manager.update_log_manager_w_cached_data()
   else:
      log_manager.collect_historical_log_data()
   log_manager.manage_data()
   log_manager.save_pr_data()

   return 0

if __name__ == "__main__": # pragma: no coverage

   get_data()
