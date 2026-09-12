import os
from mdutils.mdutils import MdUtils
import pandas as pd
from scripts.write_test_summary import *

def test_init_HTMLBuilder(set_env_vars):

   set_env_vars
   html_builder = HTMLBuilder()
   assert html_builder.mdFile.get_md_text() == "\nTest Summary for PR #2882\n=========================\n"
   assert html_builder.mdFile.file_name == 'summary.md'
   assert len(html_builder.repo_hashes) == 30
   assert html_builder.machines == os.environ.get("MACHINES").split()
   assert html_builder.categories == ["runtime", "memory"]

def test_set_mdFile(set_env_vars, monkeypatch):

   set_env_vars
   html_builder = HTMLBuilder()
   assert html_builder.mdFile.get_md_text() == "\nTest Summary for PR #2882\n=========================\n"
   monkeypatch.setenv("PR_NUM", "3207")
   html_builder.set_mdFile()
   assert html_builder.mdFile.get_md_text() == "\nTest Summary for PR #3207\n=========================\n"
   assert html_builder.mdFile.file_name == 'summary.md'
   
def test_get_mdFile(set_env_vars, monkeypatch):

   set_env_vars
   monkeypatch.setenv("PR_NUM", "3190")
   html_builder = HTMLBuilder()
   assert html_builder.get_mdFile().get_md_text() == "\nTest Summary for PR #3190\n=========================\n"
   assert html_builder.get_mdFile().file_name == 'summary.md'

def test_add_legend():
   
   html_builder = HTMLBuilder()
   html_builder.add_legend("memory")
   predicted_text = "\nTest Summary for PR #2882\n=========================\n\n\n" + "<h4>Key:</h4>\n\n" + \
   f"&nbsp;&nbsp;&nbsp;&nbsp;✅ = NORMAL memory: Memory falls within two standard deviations of the mean.\n\n" + \
   f"&nbsp;&nbsp;&nbsp;&nbsp;⚠️ = Memory WARNING: Memory is greater than two standard deviations above the mean.\n\n" + \
   f"&nbsp;&nbsp;&nbsp;&nbsp;❌ = Memory FAIL: For the past 3+ PRs, memory has been greater than two standard deviations above the mean.\n\n" + \
   f"&nbsp;&nbsp;&nbsp;&nbsp;N/A = Test does not run on this machine.\n\n\n"
   assert html_builder.get_mdFile().get_md_text() == predicted_text

def test_build_content(set_env_vars,monkeypatch,sample_runtime_results, actual_passes_per_test, actual_passes_per_machine):

   set_env_vars
   monkeypatch.setenv("RUNTIME_RESULTS", "data/runtime_results.json")
   html_builder = HTMLBuilder()
   content = html_builder.build_content("runtime").sort_index()

   # Create comparison DataFrame from fixtures
   sample_runtime_results["Passing"] = actual_passes_per_test
   actual_results = pd.DataFrame.from_dict(sample_runtime_results).fillna("N/A")
   actual_passes_per_machine = pd.DataFrame.from_dict(actual_passes_per_machine, orient='index', columns=["hercules","orion","ursa","Passing"])
   actual_results = pd.concat([actual_results,actual_passes_per_machine]).sort_index()

   assert content.equals(actual_results)

def test_write_content(monkeypatch, set_env_vars, sample_runtime_results_complete, failing_results_table, actual_passes_per_machine):
   """Compare the results of write_content() with a markdown table containing the expected results.
   """
   
   # Set up and test write_content() method
   set_env_vars
   monkeypatch.setenv("RUNTIME_RESULTS", "data/runtime_results.json")
   monkeypatch.setenv("MACHINES", "hercules orion ursa")
   html_builder = HTMLBuilder()
   results = pd.DataFrame.from_dict(sample_runtime_results_complete).fillna("N/A").sort_index()
   results = pd.concat([results, pd.DataFrame.from_dict(actual_passes_per_machine, orient='index', columns=["hercules","orion","ursa","Passing"])])
   html_builder.write_content(results, "runtime")

   # Create comparison markdown table with only failing results
   table_header = "\nTest Summary for PR #2882\n=========================\n\n" + \
                  "|Test|hercules|orion|ursa|Passing|\n" + "| :---: | :---: | :---: | :---: | :---: |\n|"
   table_contents =  table_header + failing_results_table + "\n\n\n</details>"

   assert html_builder.get_mdFile().get_md_text() == table_contents

def test_create_summary(set_env_vars, monkeypatch, failing_results_table):
   """Compare the results of create_summary() with a markdown string containing the expected results.
   """
   
   set_env_vars
   monkeypatch.setenv("RUNTIME_RESULTS", "data/sample_runtime_results.json")
   monkeypatch.setenv("MACHINES", "hercules orion ursa")

   html_builder = HTMLBuilder()
   html_builder.categories = ["runtime"] # Only create summary for runtime to simplify testing
   html_builder.create_summary()
   
   # Create comparison markdown table with only failing results
   table_header = "\nTest Summary for PR #2882\n=========================\n" + \
                  "<details><summary><h3>RUNTIME Results Summary</h3></summary>\n" + \
                  "\n\n\n\n<h4>Key:</h4>\n\n" + "&nbsp;&nbsp;&nbsp;&nbsp;✅ = NORMAL runtime: Runtime falls within two standard deviations of the mean.\n\n" + \
                  "&nbsp;&nbsp;&nbsp;&nbsp;⚠️ = Runtime WARNING: Runtime is greater than two standard deviations above the mean.\n\n" + \
                  "&nbsp;&nbsp;&nbsp;&nbsp;❌ = Runtime FAIL: For the past 3+ PRs, runtime has been greater than two standard deviations above the mean.\n\n" + \
                  "&nbsp;&nbsp;&nbsp;&nbsp;N/A = Test does not run on this machine.\n\n\n\n" + \
                  "|Test|hercules|orion|ursa|Passing|\n" + "| :---: | :---: | :---: | :---: | :---: |\n|"
   
   table_contents = table_header + failing_results_table + "\n\n\n</details>"

   print(f"ACTUAL RESULTS: {html_builder.get_mdFile().get_md_text()}")

   print(f"EXPECTED RESULTS: {table_contents}")

   assert html_builder.get_mdFile().get_md_text() == table_contents


def test_count_passes_per_machine(sample_runtime_results, actual_passes_per_machine):
   """Tests whether the calculated number of tests passing per machine is the same as the actual number of tests passing per machine."""
   
   html_builder = HTMLBuilder()

   # Set up dataframe with test results
   results = pd.DataFrame()

   for machine in sample_runtime_results.keys():
      machine_results = pd.DataFrame.from_dict(sample_runtime_results[machine], orient='index',columns=[machine])
      results = pd.merge(results, machine_results, left_index=True, right_index=True, how='outer').fillna("N/A")

   # Calculate passing tests per machine
   results = html_builder._count_passes_per_machine(results)
   actual_values = pd.DataFrame.from_dict(actual_passes_per_machine, orient='index', columns=["hercules","orion","ursa","Passing"])
   
   assert results.equals(actual_values)

def test_count_passes_per_test(sample_runtime_results, actual_passes_per_test):
   """Tests whether the calculated number of tests passing is the same as the actual number of tests passing."""
   
   html_builder = HTMLBuilder()

   # Set up dataframe with test results
   results = pd.DataFrame()

   for machine in sample_runtime_results.keys():
      machine_results = pd.DataFrame.from_dict(sample_runtime_results[machine], orient='index',columns=[machine])
      results = pd.merge(results, machine_results, left_index=True, right_index=True, how='outer').fillna("N/A")

   # Calculate passing tests
   results = html_builder._count_passes_per_test(results)['Passing']

   # Sort by index before comparing calculated and actual values for equality
   assert results.sort_index().equals(pd.Series(actual_passes_per_test, name='Passing').sort_index())
