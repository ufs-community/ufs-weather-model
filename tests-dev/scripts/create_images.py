import os, sys
import logging
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from .Manager import *

class PlotManager(Manager):
   """Manages information common to all plots.
   Attributes: 
      category:
      data: {'machine': 
               {'test': 
                  {'hash1': runtime1, 'hash2': runtime2}
                  }
               }
            }  
   """

   def __init__(self, category):
      """Create a Plot Manager object that holds data related to all plots. 
      Args:
         data (dict): Runtime and memory data for each test and machine. Primary key is machine. Secondary key is test. 
         category (str): 'runtime' or 'memory'
      """
      super().__init__() # set hashes, machine list, categories
      self.category = category

   def get_test_names(self):
      """Create a set containing all test names by extracting the tests (keys) from the data_by_machine
      Returns:
         all_tests: Set of all test names
      """
      all_tests = set()
      for data_by_machine in self.current_pr_data.values():
         all_tests.update(data_by_machine.keys())
      
      return all_tests

   def organize_data_by_test(self):
      """Creates new runtime/memory dictionaries that use test name as key and have data for each machine 
      under each test. 
      Returns:
         metrics (dict): Runtime or memory data for each test and machine. Primary key is test. Secondary key is machine:
            {'datm_cdeps_ciceC_cfsr_intel': 
               {'acorn': {'cb16f329': 188, 'ead2c35f': 188, ...}
               'gaeac6': {'cb16f329': 119, 'ead2c35f': 120, ...}}
            }
      """

      tests = self.get_test_names()

      # Create a three-level deep dictionary where any key access at the first or second level that doesn't exist will automatically be created
      #metrics = defaultdict(lambda: dict) # --> Need lambda? Could just use dict as arg? Or is lambda necessary for multiple layers...
      metrics = {}

      for test in tests:
         for machine, test_data in self.historical_data.items():
            if test not in test_data:
               continue # No data to add
            else:
               metrics.setdefault(test, {}).update({machine: test_data[test]}) # Add historical test data
               metrics[test].setdefault(machine, {}).update({self.pr_head_commit: self.current_pr_data[machine].setdefault(test, None)}) # Add data for current PR
      
      return metrics
   
   def detect_statistical_anomalies(self, test_data):
      """Detect statistical anomalies, aka tests w/runtime or memory usage greater than 2 standard deviations above the mean.
      Args:
         test_data (list): Data that needs to be checked for statistical anomalies (may include None values)
      Returns:
         anomalies (list): A boolean list where True indicates the presence of an anomalous value.
      """

      # Filter out None values to compute stats on valid historical values 
      # (exclude PR head commit value from calculation using test_data[:-1])
      valid_data = [value for value in test_data[:-1] if value is not None]

      if not valid_data: # pragma: no cover
         return [False] * len(test_data)
      
      # Calculate threshold value (mean + 2*stdev) using only valid numbers (not Nones)
      mean = np.mean(valid_data, dtype=float)
      stdev = np.std(valid_data, dtype=float)
      threshold = mean + (2 * stdev)

      # Determine anomalies
      anomalies = [True if value is not None and value > threshold else False for value in test_data]

      return anomalies

   def rearrange_hashes(self):
      """Rearrange hashes for use in plotting function
      Returns:
         hashes (list): Commit metadata; by default, 30 most recent hashes from the repository plus 'PR Head'.
      """
      hashes = self.get_hashes().copy() # Default 30 hashes; change quantity in utility function for consistency
      hashes.insert(0, "PR Head")
      hashes.reverse()
      return hashes

   def generate_figure(self, test):
      """For each test, create a plot containing figure metadata (e.g., title, xticks, yticks)."""
      plt.figure(figsize=(14, 6), dpi=200)
      plt.title(f"{self.category} for {test}", fontsize=16)
      plt.xlabel("Commit Hash: oldest --> newest", fontsize=14)
      plt.ylabel(self.category, fontsize=14)
      plt.xticks(np.arange(len(self.hashes)), labels=self.hashes, rotation=45, fontsize=10)
      plt.yticks(fontsize=12)
      plt.grid(True, linestyle='--', alpha=0.5)

      return plt

   def add_test_metrics_by_machine(self, plt, test):
      """For each test, plot lines with data for each machine."""
      
      styles = ['o-', 's--', '^-', 'd:', 'x-.', 'v--', '*-', 'p:']
      # Add one line to the plot with data for each machine
      for i, machine in enumerate(self.metrics[test]):
         test_data = self.metrics[test][machine]

         # Map test data to commits, which are arranged in self.hashes from oldest to newest
         y = [test_data.get(self.pr_head_commit if hash == "PR Head" else hash, None) for hash in self.hashes]
         # Skip plotting if the machine has no valid data points
         if all(value is None for value in y):
            continue
         anomalies = self.detect_statistical_anomalies(y)
         
         x = self.hashes 
         
         plt.plot(x, y, styles[i % len(styles)], label=f"{machine}", linewidth=2, markersize=6)
         [plt.plot(x[idx], y[idx], 'ro', markersize=8) for idx, is_anomaly in enumerate(anomalies) if is_anomaly]
         
      plt.legend(fontsize=12)
      plt.tight_layout()

      return plt

   def save_plot_image(self, plt, test):
      
      png_path = f"plots/{test}_{self.category}.png"
      os.makedirs("plots", exist_ok=True)
      plt.savefig(png_path)
      plt.close()

   def load_data(self, file_path):
      """Loads historical or current PR data from JSON file.
      Args:
         file_path (str): Path to file
      """
      try:
         return self.load_json_from_file(file_path)
      except FileNotFoundError:
         logging.error(f"Could not load JSON file {file_path}.")
         sys.exit()

   def process_data(self):
      """
      Combines and manipulates the data to prepare it for plotting.
      """
      self.historical_data = self.load_data(f"{os.environ.get('PLOT_DATA')}/historical_{self.category}.json")
      self.current_pr_data = self.load_data(f"{os.environ.get('PLOT_DATA')}/current_pr_{self.category}_data.json")
      self.metrics = self.organize_data_by_test()
      self.hashes = self.rearrange_hashes()

   def plot_results(self):
      """
      Generates anomaly-highlighted plots.
      """
      
      # Create one plot per test
      for test in self.metrics:
         plt = self.generate_figure(test)
         plt = self.add_test_metrics_by_machine(plt, test)
         self.save_plot_image(plt, test)


def main():

   # Don't need to commit/push old plots cuz that has presumably already been done
   # Maybe need to make new plots w/most recent commit? 
   for category in ["runtime", "memory"]:
      plot_manager = PlotManager(category)
      plot_manager.process_data()
      plot_manager.plot_results()

   return plot_manager

if __name__ == "__main__":

   main()