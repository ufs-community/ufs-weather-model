from mdutils.mdutils import MdUtils
import argparse

def get_test_output(file_path):
   """Read in test output from file."""
   with open(file_path, 'r', encoding='utf-8') as file:
      data = file.read().split('>')

   return data

def create_mdFile(text):
   """Create a markdown file named test_summary.md."""
   mdFile = MdUtils(file_name='test_summary.md', title=f'Test Summary')
   for line in text:
      mdFile.new_paragraph(f"{line}")
      mdFile.new_paragraph(f" ")

   return mdFile

def main():
   parser = argparse.ArgumentParser(description="Generate test summary markdown.")
   parser.add_argument("--file", type=str, default="output.txt", help="Path to the test output file.")
   args = parser.parse_args()

   data = get_test_output(args.file)
   mdFile = create_mdFile(data)
   print(mdFile.get_md_text())

if __name__ == "__main__":

   main()
