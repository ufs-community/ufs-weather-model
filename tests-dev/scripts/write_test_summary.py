import os
from mdutils.mdutils import MdUtils
import pandas as pd
from .create_images import *
from .Manager import *

class HTMLBuilder(Manager):

   def __init__(self):
      super().__init__() # set hashes, machine list, categories
      self.set_mdFile()

   def set_mdFile(self):
      """Create a markdown file named summary.md with the PR# in the title."""
      pr_num = os.environ.get('PR_NUM')
      self.mdFile = MdUtils(file_name='summary.md', title=f'Test Summary for PR #{pr_num}')

   def get_mdFile(self):
      return self.mdFile

   def add_legend(self, category):
      """Add HTML for a key/legend describing the 4 possible statuses (✅ ⚠️ ❌ N/A)."""
      self.mdFile.new_paragraph("<h4>Key:</h4>")
      self.mdFile.new_paragraph(f"&nbsp;&nbsp;&nbsp;&nbsp;✅ = NORMAL {category}: {category.title()} falls within two standard deviations of the mean.")
      self.mdFile.new_paragraph(f"&nbsp;&nbsp;&nbsp;&nbsp;⚠️ = {category.title()} WARNING: {category.title()} is greater than two standard deviations above the mean.")
      self.mdFile.new_paragraph(f"&nbsp;&nbsp;&nbsp;&nbsp;❌ = {category.title()} FAIL: For the past 3+ PRs, {category} has been greater than two standard deviations above the mean.")
      self.mdFile.new_paragraph(f"&nbsp;&nbsp;&nbsp;&nbsp;N/A = Test does not run on this machine.")
      self.mdFile.new_paragraph('\n')

   def build_content(self, category):
      """Load the runtime or memory results dictionary, convert to dataframe, and return the results
      Returns:
         results: DataFrame containing the runtime/memory testing results. Rows are tests and columns are machines.
      """

      contents = self.load_json_from_file(os.environ.get(f"{category.upper()}_RESULTS"))
      results = pd.DataFrame()

      for machine in contents:
         
         machine_results = pd.DataFrame.from_dict(contents[machine], orient='index', columns=[machine])
         results = pd.merge(results, machine_results, left_index=True, right_index=True, how='outer').fillna("N/A")

      results = self._count_passes_per_test(results)
      results = pd.concat([results, self._count_passes_per_machine(results)])

      return results

   def write_content(self,data,category):
      
      # Create contents list starting with header row
      contents = ["Test"] + self.machines + ["Passing"]

      # Create table starting with one row (header)
      rows = 1
      for index, row in data.iterrows():
         warn = '⚠️'
         fail = '❌'
         # If there is a warn or fail in the row, add the row to contents to be printed; also add summary row
         if (data.loc[index] == warn).any() or (data.loc[index] == fail).any() or (index == 'Platform Total (Passing):'):
            rows += 1
            if (index != 'Platform Total (Passing):'):
               img_link = f"[{str(index)}](https://raw.githubusercontent.com/wiki/NOAA-EPIC/ufs-weather-model/plots/{str(index)}_{category}.png)"
               contents.append(img_link)
            else: 
               contents.append(index)
            for item in row:
               contents.append(item)

      self.mdFile.new_table(columns=(len(self.machines) + 2), rows=rows, text_align='center', text=contents)
      self.mdFile.new_paragraph('\n')
      self.mdFile.write('</details>')

   def _count_passes_per_machine(self,data):
      """Counts number of passing tests on each machine and procudes a row with the totals.
      Args:
         data(DataFrame): Table of tests and pass/warn/fail status by machine
      Returns:
         machine_total(DataFrame): Number of tests passing per machine
      """

      # Counts for passing tests
      passing_tests_by_machine = round((data.eq('✅').sum(axis=0)/data.ne('N/A').sum(axis=0) * 100),1).astype(str)
      for machine in passing_tests_by_machine.index:
         passing_tests_by_machine[machine] = f"**{machine.upper()}:** " + passing_tests_by_machine[machine] + "% passing"
      passing_tests_by_machine.name = 'Platform Total (Passing):'
      # Set bottom right corner to empty string
      passing_tests_by_machine.loc['Passing'] = ''
      machine_total = pd.DataFrame(passing_tests_by_machine).T
      
      return machine_total

   def _count_passes_per_test(self,data):
      """Counts number of platforms on which a given test passes and adds a column to the table.
      Args:
         data (DataFrame): DataFrame containing pass/warn/fail status for each test on each machine
      Returns:
         data: with an extra column listing pass rates for each test 
      """

      passing_tests = (round((data.eq('✅').sum(axis=1) / data.ne('N/A').sum(axis=1) * 100),1))
      passing_tests = passing_tests.astype(str).add('%')
      passing_tests.name = 'Passing'
      data = pd.merge(data, pd.DataFrame(passing_tests.astype(str)), left_index=True, right_index=True, how='inner')

      return data

   def create_summary(self):
      """Append a runtime or memory header and key and call write_contents() to write the runtime/memory table to the file.
      Args:
         categories (list): Test categories. Currently 'runtime' and 'memory'.
      Returns:
         mdFile: A markdown file
      """

      for category in self.categories: 
         # Create <details> section with legend
         self.mdFile.write(f"<details><summary><h3>{category.upper()} Results Summary</h3></summary>")
         self.mdFile.new_paragraph('\n')
         self.add_legend(category)

         # Create a DataFrame w/the runtime/memory results content
         data = self.build_content(category)

         # Write the content to a file
         self.write_content(data,category)

def main(): # pragma: no cover

   html_builder = HTMLBuilder()
   html_builder.create_summary()
   print(html_builder.get_mdFile().get_md_text())

if __name__ == "__main__": # pragma: no cover
   
   main()
