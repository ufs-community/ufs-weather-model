import os
import json
import re
from mdutils.mdutils import MdUtils
import pandas as pd

def load_json(file_path):
   """Convert JSON file to python dictionary."""
   with open(file_path, 'r', encoding='utf-8') as file:
      data = json.load(file)

   return data

def create_mdFile():
   """Create a markdown file named summary.md with the PR# in the title."""
   pr_num = os.environ.get('PR_NUM')
   mdFile = MdUtils(file_name='summary.md', title=f'Test Summary for PR #{pr_num}')

   return mdFile

def build_content(category):
   """Load the runtime or memory results dictionary, convert to dataframe, and return the results
   Args: 
      category (str): "runtime" or "memory"
   Returns:
      results: DataFrame containing the runtime/memory testing results. Rows are tests and columns are machines.
   """

   contents = load_json(os.environ.get(f"{category.upper()}_RESULTS"))
   results = pd.DataFrame()
   
   for machine in contents:
      
      machine_results = pd.DataFrame.from_dict(contents[machine], orient='index', columns=[machine])
      results = pd.merge(results, machine_results, left_index=True, right_index=True, how='outer').fillna("N/A")

   results = _count_passes_per_test(results)
   results = pd.concat([results, _count_passes_per_machine(results)])
   
   return results

def write_content(data, mdFile):
   
   machines = os.environ.get('MACHINES').split()
   
   # Create contents list starting with header row
   contents = ["Test"] + machines + ["Passing"]

   # Create table starting with one row (header)
   rows = 1
   for index, row in data.iterrows():
      warn = '⚠️'
      fail = '❌'
      # If there is a warn or fail in the row, add the row to contents to be printed; also add summary row
      if (data.loc[index] == warn).any() or (data.loc[index] == fail).any() or (index == 'Platform Total (Passing):'):
         rows += 1
         contents.append(str(index))
         for item in row:
            contents.append(item)

   mdFile.new_table(columns=(len(machines) + 2), rows=rows, text_align='center', text=contents)
   mdFile.new_paragraph('\n')
   mdFile.write('</details>')

   return mdFile

def _count_passes_per_machine(data):
   """Counts number of passing tests on each machine and procudes a row with the totals.
   Args:
      data(DataFrame): Table of tests and pass/warn/fail status by machine
   Returns:
      machine_total(DataFrame): Number of tests passing per machine
   """

   # Counts for passing tests
   passing_tests_by_machine = data.eq('✅').sum(axis=0).astype(str) + '/' + data.ne('N/A').sum(axis=0).astype(str)
   for machine in passing_tests_by_machine.index:
      passing_tests_by_machine[machine] = f"**{machine.upper()}:** " + passing_tests_by_machine[machine] + " passing"
   passing_tests_by_machine.name = 'Platform Total (Passing):'
   # Set bottom right corner to empty string
   passing_tests_by_machine.loc['Passing'] = ''
   machine_total = pd.DataFrame(passing_tests_by_machine).T
   
   return machine_total

def _count_passes_per_test(data):
   """Counts number of platforms on which a given test passes and adds a column to the table.
   Args:
      data (DataFrame): DataFrame containing pass/warn/fail status for each test on each machine
   Returns:
      data: with an extra column listing pass rates for each test 
   """

   passing_tests = data.eq('✅').sum(axis=1).astype(str) + "/" + data.ne('N/A').sum(axis=1).astype(str)
   passing_tests.name = 'Passing'
   data = pd.merge(data, pd.DataFrame(passing_tests), left_index=True, right_index=True, how='inner')

   return data

def create_summary(categories):
   """Append a runtime or memory header and key and call write_contents() to write the runtime/memory table to the file.
   Args:
      categories (list): Test categories. Currently 'runtime' and 'memory'.
   Returns:
      mdFile: A markdown file
   """

   mdFile = create_mdFile()

   for category in categories: 
      # Create <details> section
      mdFile.write(f"<details><summary><h3>{category.upper()} Results Summary</h3></summary>")
      mdFile.new_paragraph('\n')
      # Add key to section
      mdFile.new_paragraph("<h4>Key:</h4>")
      mdFile.new_paragraph(f"&nbsp;&nbsp;&nbsp;&nbsp;✅ = NORMAL {category}: {category.title()} falls within two standard deviations of the mean.")
      mdFile.new_paragraph(f"&nbsp;&nbsp;&nbsp;&nbsp;⚠️ = {category.title()} WARNING: {category.title()} is greater than two standard deviations above the mean.")
      mdFile.new_paragraph(f"&nbsp;&nbsp;&nbsp;&nbsp;❌ = {category.title()} FAIL: For the past 2+ PRs, {category} has been greater than two standard deviations above the mean.")
      mdFile.new_paragraph(f"&nbsp;&nbsp;&nbsp;&nbsp;N/A = Test does not run on this machine.")
      mdFile.new_paragraph('\n')
      # Create a DataFrame w/the runtime/memory results content
      data = build_content(category)

      # Write the content to a file
      mdFile = write_content(data, mdFile)
   
   return mdFile
   
def main(): # pragma: no cover

   summary = create_summary(['runtime', 'memory'])
   print(summary.get_md_text())

if __name__ == "__main__": # pragma: no cover
   
   main()
