import pytest
from mdutils.mdutils import MdUtils
from scripts.generate_message_text import *

def test_init_MessageManager(set_env_vars):

   set_env_vars
   message_manager = MessageManager()

   assert list(message_manager.contents.keys()) == ["runtime", "memory"]
   assert list(message_manager.results.keys()) == ["runtime", "memory"]
   assert message_manager.message_content == ""

def test_get_test_names(set_env_vars,sample_runtime_results_complete,sample_memory_results,message_results):

   set_env_vars
   manager = MessageManager()
   manager.contents["runtime"].update(sample_runtime_results_complete)
   manager.contents["memory"].update(sample_memory_results)

   for category in manager.categories:

      actual_test_names = manager.get_test_names(manager.contents[category])
      expected_test_names = message_results[category].keys()
      assert actual_test_names == expected_test_names

def test_organize_data_by_test(set_env_vars,sample_runtime_results_complete,sample_memory_results,message_results):

   set_env_vars
   manager = MessageManager()
   manager.contents["runtime"].update(sample_runtime_results_complete)
   manager.contents["memory"].update(sample_memory_results)

   expected_results = {"runtime": 
         {
                  "cpld_control_p8_mixedmode_intel": {"hercules": "\u2705", "orion": "\u2705","ursa": "\u26a0\ufe0f","Passing": "66.7%",},
                  "cpld_control_gefs_intel": {"hercules": "\u2705", "orion": "\u2705","ursa": "\u2705","Passing": "100.0%",},
                  "cpld_restart_gefs_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u2705","Passing": "100.0%",},
                  "cpld_dcp_gefs_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u2705","Passing": "100.0%",},
                  "cpld_control_gfsv17_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u2705","Passing": "100.0%",},
                  "cpld_control_gfsv17_iau_intel": {"hercules": "\u26a0\ufe0f","orion": "\u2705","ursa": "\u2705","Passing": "66.7%",},
                  "cpld_restart_gfsv17_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u2705","Passing": "100.0%",},
                  "cpld_restart_gfsv17_iau_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u26a0\ufe0f","Passing": "66.7%",},
                  "cpld_mpi_gfsv17_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u2705","Passing": "100.0%",},
                  "cpld_control_sfs_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u26a0\ufe0f","Passing": "66.7%",},
                  "cpld_debug_gfsv17_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u2705","Passing": "100.0%",},
                  "cpld_control_p8_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u274c","Passing": "66.7%",},
                  "cpld_control_p8.v2.sfc_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u26a0\ufe0f","Passing": "66.7%",},
                  "cpld_restart_p8_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u2705","Passing": "100.0%",},
                  "cpld_control_qr_p8_intel": {"hercules": "\u26a0\ufe0f","orion": "\u2705","ursa": "\u26a0\ufe0f","Passing": "33.3%",},
                  "cpld_restart_qr_p8_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u2705","Passing": "100.0%",},
                  "cpld_2threads_p8_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u26a0\ufe0f","Passing": "66.7%",},
                  "cpld_decomp_p8_intel": {"hercules": "\u26a0\ufe0f","orion": "\u2705","ursa": "\u26a0\ufe0f","Passing": "33.3%",},
                  "cpld_mpi_p8_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u2705","Passing": "100.0%",},
                  "cpld_control_gfsv17_intelllvm": {"ursa": "\u2705","Passing": "100.0%",}
         },
         "memory": {
                  "cpld_control_p8_mixedmode_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u26a0\ufe0f",},
                  "cpld_control_gefs_intel": {"hercules": "\u2705","orion": "\u2705","ursa": "\u2705",},
                  "cpld_restart_gefs_intel": {"orion": "\u2705","ursa": "\u2705",},
                  "cpld_dcp_gefs_intel": {"orion": "\u2705","ursa": "\u2705",},
                  "cpld_control_gfsv17_intel": {"orion": "\u2705","ursa": "\u2705",},
                  "cpld_control_gfsv17_iau_intel": {"orion": "\u2705","ursa": "\u2705",},
                  "cpld_restart_gfsv17_intel": {"orion": "\u2705","ursa": "\u2705",},
                  "cpld_restart_gfsv17_iau_intel": {"orion": "\u2705","ursa": "\u26a0\ufe0f",},
                  "cpld_control_noaero_p8_agrid_intel": {"hercules": "\u2705",},
                  "cpld_mpi_gfsv17_intel": {"orion": "\u2705","ursa": "\u2705",},
                  "cpld_control_sfs_intel": {"orion": "\u2705","ursa": "\u26a0\ufe0f",},
                  "cpld_debug_gfsv17_intel": {"orion": "\u2705","ursa": "\u2705",},
                  "cpld_control_p8_intel": {"orion": "\u2705","ursa": "\u274c",},
                  "cpld_control_p8.v2.sfc_intel": {"orion": "\u2705","ursa": "\u26a0\ufe0f",},
                  "cpld_restart_p8_intel": {"orion": "\u2705","ursa": "\u2705",},
                  "cpld_control_qr_p8_intel": {"orion": "\u2705","ursa": "\u26a0\ufe0f",},
                  "cpld_restart_qr_p8_intel": {"orion": "\u2705","ursa": "\u2705",},
                  "cpld_2threads_p8_intel": {"orion": "\u2705","ursa": "\u26a0\ufe0f",},
                  "cpld_decomp_p8_intel": {"orion": "\u2705","ursa": "\u26a0\ufe0f",},
                  "cpld_mpi_p8_intel": {"orion": "\u2705","ursa": "\u2705",},
                  "cpld_control_gfsv17_intelllvm": {"ursa": "\u2705",},
                  "control_c48_intel": {"hercules": "\u2705",},
                  "control_p8_intel": {"hercules": "\u2705",},
                  "control_restart_p8_intel": {"hercules": "\u2705",},
                  "hrrr_control_intel": {"hercules": "\u2705",},
                  "atmaero_control_p8_intel": {"hercules": "\u2705",},
                  "regional_atmaq_intel": {"hercules": "\u2705",},
                  "hafs_regional_docn_intel": {"hercules": "\u2705",},
                  "datm_cdeps_control_cfsr_intel": {"hercules": "\u2705",},
                  "control_c48_gnu": {"hercules": "\u2705",},
                  "control_p8_gnu": {"hercules": "\u2705",},
                  "control_debug_p8_gnu": {"hercules": "\u2705",},
                  "hrrr_control_gnu": {"hercules": "\u26a0\ufe0f",},
                  "datm_cdeps_control_cfsr_gnu": {"hercules": "\u2705",},
               },
      }

   for category in manager.categories:
      actual_results = manager.organize_data_by_test(manager.contents[category])
      assert actual_results == message_results[category]

   

@pytest.mark.parametrize('category', ['runtime', 'memory'])
def test_create_message_content(set_env_vars, category, sample_runtime_results_complete, sample_memory_results):
   """If there are tests whose runtime/memory results are 2 standard deviations above the mean for 3+ tests, 
   append this information to self.message_content.
   """

   set_env_vars
   manager=MessageManager()
   manager.contents["runtime"].update(sample_runtime_results_complete)
   manager.contents["memory"].update(sample_memory_results)

   manager.create_message_content(manager.contents[category], category)

   hi = f"\nFor the past three PRs, {category.upper()} has been greater than two standard deviations above the mean for the following tests: \n\n" + \
            "  * cpld_control_p8_intel: ursa\n"
   #hi_mem = "\nFor the past three PRs, MEMORY has been greater than two standard deviations above the mean for the following tests: \n\n" + \
            #"  * cpld_control_p8_intel: ursa\n"

   assert manager.message_content == hi
   


# For runtime, mem, both, and neither
def test_create_md_file(set_env_vars,message_content,expected_messages):

   set_env_vars
   manager=MessageManager()

   for key in message_content:
      manager.message_content = message_content[key] 
      manager.create_md_file()

      with open("pr_post.md") as fh:
         actual_contents = fh.read()
         assert actual_contents == expected_messages[key]

def test_main(set_env_vars,monkeypatch):

   set_env_vars
   monkeypatch.setenv("RUNTIME_RESULTS", "/Users/gpetro/wm-warn/tests-dev/scripts/tests/data/runtime_results.json")
   monkeypatch.setenv("MEMORY_RESULTS", "/Users/gpetro/wm-warn/tests-dev/scripts/tests/data/memory_results.json")
   mdFile = main()
   print(mdFile)

   assert mdFile
   assert isinstance(mdFile, MdUtils)

