from .Manager import *
from .HistoricalLogManager import *

"""This script contains a main() function that gets log information from GitHub using the APICall class 
and extracts data from the RegressionTest_<machine>.log files for each machine. 
"""

def main():
   """For each machine, create a log object, get current PR data, gather historical runtime/memory data, 
   and compare results to determine which test/machine combinations fall more than 2 standard deviations 
   above the historical mean for each test.""" 

   log_manager = HistoricalLogManager()
   log_manager.manage_data()
   log_manager.save_data()

   return 0

if __name__ == "__main__": # pragma: no coverage

   main()
