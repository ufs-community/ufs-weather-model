from mdutils.mdutils import MdUtils
import pandas as pd
from scripts.write_test_summary import *
from scripts.write_test_summary import _count_passes_per_machine, _count_passes_per_test

def test_load_json(stats_dict_snippet):

   content = load_json('test_file_stats.json')
   assert stats_dict_snippet == content

def test_create_mdFile():

   mdFile = create_mdFile()
   assert mdFile.get_md_text() == "\nTest Summary for PR #2882\n=========================\n"
   assert mdFile.file_name == 'summary.md'

def test_build_content(sample_runtime_results, actual_passes_per_test, actual_passes_per_machine):

   os.environ["RUNTIME_RESULTS"] = "runtime_results.json"
   content = build_content("runtime").sort_index()

   # Create comparison DataFrame from fixtures
   sample_runtime_results["Passing"] = actual_passes_per_test
   actual_results = pd.DataFrame.from_dict(sample_runtime_results).fillna("N/A")
   actual_passes_per_machine = pd.DataFrame.from_dict(actual_passes_per_machine, orient='index', columns=["hercules","orion","ursa","Passing"])
   actual_results = pd.concat([actual_results,actual_passes_per_machine]).sort_index()

   assert content.equals(actual_results)

def test_write_content(sample_runtime_results_complete, failing_results_table, actual_passes_per_machine):
   """Compare the results of write_content() with a markdown table containing the expected results.
   """
   
   # Set up and test write_content() method
   mdFile = create_mdFile()
   os.environ["MACHINES"] = "hercules orion ursa"
   results = pd.DataFrame.from_dict(sample_runtime_results_complete).fillna("N/A").sort_index()
   results = pd.concat([results, pd.DataFrame.from_dict(actual_passes_per_machine, orient='index', columns=["hercules","orion","ursa","Passing"])])
   results = write_content(results, mdFile)

   # Create comparison markdown table with only failing results
   table_header = "\nTest Summary for PR #2882\n=========================\n\n" + \
                  "|Test|hercules|orion|ursa|Passing|\n" + "| :---: | :---: | :---: | :---: | :---: |\n|"
   table_contents =  table_header + failing_results_table + "\n\n\n</details>"

   assert results.get_md_text() == table_contents

def test_create_summary(failing_results_table):
   """Compare the results of create_summary() with a markdown string containing the expected results.
   """
   
   summary_file = create_summary(['runtime'])

   # Create comparison markdown table with only failing results
   table_header = "\nTest Summary for PR #2882\n=========================\n" + \
                  "<details><summary><h3>RUNTIME Results Summary</h3></summary>\n" + \
                  "\n\n\n\n<h4>Key:</h4>\n\n" + "&nbsp;&nbsp;&nbsp;&nbsp;✅ = NORMAL runtime: Runtime falls within two standard deviations of the mean.\n\n" + \
                  "&nbsp;&nbsp;&nbsp;&nbsp;⚠️ = Runtime WARNING: Runtime is greater than two standard deviations above the mean.\n\n" + \
                  "&nbsp;&nbsp;&nbsp;&nbsp;❌ = Runtime FAIL: For the past 2+ PRs, runtime has been greater than two standard deviations above the mean.\n\n" + \
                  "&nbsp;&nbsp;&nbsp;&nbsp;N/A = Test does not run on this machine.\n\n\n\n" + \
                  f"|Test|hercules|orion|ursa|Passing|\n" + "| :---: | :---: | :---: | :---: | :---: |\n|"
   
   table_contents = table_header + failing_results_table + "\n\n\n</details>"

   assert summary_file.get_md_text() == table_contents


def test_count_passes_per_machine(sample_runtime_results, actual_passes_per_machine):
   """Tests whether the calculated number of tests passing per machine is the same as the actual number of tests passing per machine."""
   
   # Set up dataframe with test results
   results = pd.DataFrame()

   for machine in sample_runtime_results.keys():
      machine_results = pd.DataFrame.from_dict(sample_runtime_results[machine], orient='index',columns=[machine])
      results = pd.merge(results, machine_results, left_index=True, right_index=True, how='outer').fillna("N/A")

   # Calculate passing tests per machine
   results = _count_passes_per_machine(results)
   actual_values = pd.DataFrame.from_dict(actual_passes_per_machine, orient='index', columns=["hercules","orion","ursa","Passing"])
   
   assert results.equals(actual_values)

def test_count_passes_per_test(sample_runtime_results, actual_passes_per_test):
   """Tests whether the calculated number of tests passing is the same as the actual number of tests passing."""
   
   # Set up dataframe with test results
   results = pd.DataFrame()

   for machine in sample_runtime_results.keys():
      machine_results = pd.DataFrame.from_dict(sample_runtime_results[machine], orient='index',columns=[machine])
      results = pd.merge(results, machine_results, left_index=True, right_index=True, how='outer').fillna("N/A")

   # Calculate passing tests
   results = _count_passes_per_test(results)['Passing']

   # Sort by index before comparing calculated and actual values for equality
   assert results.sort_index().equals(pd.Series(actual_passes_per_test, name='Passing').sort_index())
