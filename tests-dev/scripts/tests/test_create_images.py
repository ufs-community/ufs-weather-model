from pathlib import Path
import pytest
import requests
import numpy as np
from scripts.Manager import *
from scripts.create_images import PlotManager
from scripts.create_images import *

"""Note that the following methods do not have their own test but are run as part of another test:
   * process_data()
   * save_plot_image()
   * load_data()
"""

@pytest.mark.parametrize('category', ['runtime', 'memory'])
def test_initialize_PlotManager(category):
   """Check that PlotManager is initialized properly."""
   plot_manager = PlotManager(category)
   assert plot_manager.category == category

def test_get_runtime_test_names(set_env_vars, sample_runtime_results):
   """Compare test names extracted from sample_runtime_results via get_test_names() to the expected list of test names."""

   set_env_vars
   plot_manager = PlotManager('runtime')
   plot_manager.current_pr_data = sample_runtime_results
   actual_test_names = plot_manager.get_test_names()
   
   expected_test_names = set({"cpld_control_p8_mixedmode_intel","cpld_control_gefs_intel","cpld_restart_gefs_intel",
                          "cpld_dcp_gefs_intel","cpld_control_gfsv17_intel","cpld_control_gfsv17_iau_intel",
                          "cpld_restart_gfsv17_intel","cpld_restart_gfsv17_iau_intel","cpld_mpi_gfsv17_intel",
                          "cpld_control_sfs_intel","cpld_debug_gfsv17_intel","cpld_control_p8_intel",
                          "cpld_control_p8.v2.sfc_intel","cpld_restart_p8_intel","cpld_control_qr_p8_intel",
                          "cpld_restart_qr_p8_intel","cpld_2threads_p8_intel","cpld_decomp_p8_intel",
                          "cpld_mpi_p8_intel","cpld_control_gfsv17_intelllvm" },)

   assert expected_test_names == actual_test_names

def test_get_memory_test_names(set_env_vars, sample_memory_results):

   set_env_vars
   plot_manager = PlotManager('memory')
   plot_manager.current_pr_data = sample_memory_results
   actual_test_names = plot_manager.get_test_names()
   
   expected_test_names = set({"cpld_control_p8_mixedmode_intel", "cpld_control_gefs_intel", 
                           "cpld_control_noaero_p8_agrid_intel", "control_c48_intel", "control_p8_intel", 
                           "control_restart_p8_intel", "hrrr_control_intel", "atmaero_control_p8_intel", 
                           "regional_atmaq_intel", "hafs_regional_docn_intel", "datm_cdeps_control_cfsr_intel", 
                           "control_c48_gnu", "control_p8_gnu", "control_debug_p8_gnu", "hrrr_control_gnu", 
                           "datm_cdeps_control_cfsr_gnu", "cpld_restart_gefs_intel", "cpld_dcp_gefs_intel",
                           "cpld_control_gfsv17_intel", "cpld_control_gfsv17_iau_intel", "cpld_restart_gfsv17_intel", 
                           "cpld_restart_gfsv17_iau_intel", "cpld_mpi_gfsv17_intel", "cpld_control_sfs_intel", 
                           "cpld_debug_gfsv17_intel", "cpld_control_p8_intel", "cpld_control_p8.v2.sfc_intel",
                           "cpld_restart_p8_intel", "cpld_control_qr_p8_intel", "cpld_restart_qr_p8_intel", 
                           "cpld_2threads_p8_intel", "cpld_decomp_p8_intel", "cpld_mpi_p8_intel", 
                           "cpld_control_gfsv17_intelllvm" },)

   assert expected_test_names == actual_test_names

@pytest.mark.parametrize('category', ['runtime', 'memory'])
def test_organize_data_by_test(set_env_vars, test_data, data_by_test, category, current_pr_data):
   """Check that organize_data_by_test() creates new dictionaries that use test name as primary key instead of machine as primary key. 
   """
   set_env_vars
   plot_manager = PlotManager(category)
   plot_manager.historical_data = test_data
   plot_manager.current_pr_data = current_pr_data[category]
   # Need to add current PR data? 
   actual_data_by_test = plot_manager.organize_data_by_test()
   expected_data_by_test = data_by_test
   
   for test in expected_data_by_test:
      for machine in ['hercules', 'orion', 'ursa']:
         try:
            assert expected_data_by_test[test][category][machine] == actual_data_by_test[test][category][machine]
         except KeyError:
            continue

def test_detect_statistical_anomalies():

   data = [2091, 1195, 2699, 1896, 2098, 2712, 2249, 1620, 1938, 1132, 1978, 1215, 1523, 2257, 1852, 1184, 
           1541, 1803, 2004, 1962, 2030, 2680, 1306, 1471, 2292, 1740, 2831, 1746, 1255, 1668, 2258]

   plot_manager = PlotManager('memory')

   actual_anomalies = plot_manager.detect_statistical_anomalies(data)
   mean = 1865.6
   stdev = 474.02543
   expected_anomalies = []

   for num in data: 
      if num > mean + 2 * stdev:
         expected_anomalies.append(True)
      else:
         expected_anomalies.append(False)

   assert actual_anomalies == expected_anomalies

def test_rearrange_hashes(set_env_vars, hashes):
   set_env_vars
   plot_manager = PlotManager('runtime')
   plot_manager.repo_hashes = hashes['repo_hashes']

   plot_manager.hashes = plot_manager.rearrange_hashes()

   assert plot_manager.hashes == hashes['expected_hashes']

@pytest.mark.parametrize('category', ['runtime', 'memory'])
def test_generate_figure(set_env_vars, category, hashes):
   set_env_vars
   plot_manager = PlotManager(category)
   plot_manager.hashes = hashes
   plot = plot_manager.generate_figure("cpld_sample_test")
   
   assert plot.gca().get_title() == f"{category} for cpld_sample_test"
   assert plot.gca().get_ylabel() == category
   assert plot.gca().get_xlabel() == "Commit Hash: oldest --> newest"

@pytest.mark.parametrize('category', ['runtime', 'memory'])
@pytest.mark.parametrize('test', ['cpld_control_gfsv17_intel', 'cpld_control_gfsv17_intelllvm', 'rap_control_dyn64_phy32_intel', 
                   'datm_cdeps_control_cfsr_intel', 'control_gfs_mpas_gnu', 'cpld_control_gefs_intel', 
                   'cpld_restart_gefs_intel', 'cpld_control_sfs_intel', 'cpld_restart_sfs_intel', 
                   'regional_control_intel', 'rap_control_dyn32_phy32_gnu', 'control_p8_intel',])
def test_add_test_metrics_by_machine(set_env_vars, monkeypatch, test_data_subset, data_by_test, category, test, current_pr_data):
   """Test hash retrieval and metrics restructuring prior to plotting."""
   set_env_vars
   monkeypatch.setenv("MACHINES", "hercules orion ursa")
   plot_manager = PlotManager(category)
   
   #print(test_data_subset)
   plot_manager.historical_data = test_data_subset[category]
   plot_manager.current_pr_data = current_pr_data[category]
   plot_manager.pr_head_commit = 'PR Head' #'369cead91c98eb5c72da81ff78925250dad08903'
   plot_manager.metrics = plot_manager.organize_data_by_test()
   #oldest to newest 
   plot_manager.hashes = ['900ef4d3', '9827d29d', 'f6ca8234', '5d0ff3d6', 'ae18d62a', 'cdca1f0c', 'e202ff3e', '1e549b4e', '55bfaafe', '5eaf8148', 'PR Head']
   
   plot = plot_manager.generate_figure(test)
   plot = plot_manager.add_test_metrics_by_machine(plot, test)
   #plot_manager.save_plot_image(plot, test)

   lines = plot.gca().get_lines()
   
   for line in lines:
      if not line.get_label().startswith("_child"): #anomaly points are labeled "_child#"
         expected_data = data_by_test[category][test][line.get_label()]
         actual_data = line.get_ydata()
         assert np.array_equal(actual_data, expected_data)

   plot.close()

@pytest.mark.parametrize('category', ['runtime', 'memory'])
def test_load_data_fail(set_env_vars, category, caplog):

   set_env_vars
   plot_manager = PlotManager(category)

   with pytest.raises(SystemExit) as error:
      plot_manager.load_data("nonexistent/path")

   assert caplog.messages[0] == f"Could not load JSON file nonexistent/path."

@pytest.mark.run_manual
@pytest.mark.parametrize('category', ['runtime', 'memory'])
def test_plot_results(set_env_vars, metrics_subset, category):
   """
   Check that plotting runs error-free and generates expected files. 
   """

   # Delete files in plots dir first?
   set_env_vars
   plot_manager = PlotManager(category)
   plot_manager.metrics = metrics_subset[category]
   plot_manager.hashes = ['900ef4d3', '9827d29d', 'f6ca8234', '5d0ff3d6', 'ae18d62a', 'cdca1f0c', 'e202ff3e', '1e549b4e', '55bfaafe', '5eaf8148', '369cead91c98eb5c72da81ff78925250dad08903']
   plot_manager.plot_results()

   tests = ['cpld_control_gfsv17_intel', 'cpld_control_gfsv17_intelllvm', 'rap_control_dyn64_phy32_intel', 
                'datm_cdeps_control_cfsr_intel', 'control_gfs_mpas_gnu', 'cpld_control_gefs_intel', 
                'cpld_restart_gefs_intel', 'cpld_control_sfs_intel', 'cpld_restart_sfs_intel', 
                'regional_control_intel', 'rap_control_dyn32_phy32_gnu', 'control_p8_intel',]

   for test in tests:
      filepath = Path(f"plots/{test}_{category}.png")
      assert filepath.is_file(), f"File not found: {filepath}"

@pytest.mark.run_manual
def test_main(set_env_vars, monkeypatch):

   set_env_vars
   monkeypatch.setenv("PLOT_DATA", "/Users/gpetro/wm-warn/tests-dev/scripts/tests/data")
      
   plot_manager = main()

   assert plot_manager
   assert isinstance(plot_manager,PlotManager)

