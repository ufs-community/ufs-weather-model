import os
from .create_images import *
from .Manager import *
from .LogManager import *

"""This script contains a main() function that gets log information from GitHub using the APICall class 
and extracts data from the RegressionTest_<machine>.log files for each machine. 
"""

def main():
   """For each machine, create a log object, get current PR data, gather historical runtime/memory data, 
   and compare results to determine which test/machine combinations fall more than 2 standard deviations 
   above the historical mean for each test.""" 

   log_manager = LogManager()
   
   if os.environ.get('TEST_STATS'):
      log_manager.update_log_manager_w_cached_data()
   log_manager.manage_data()
   log_manager.save_data()

   return 0

if __name__ == "__main__": # pragma: no coverage

   main()
