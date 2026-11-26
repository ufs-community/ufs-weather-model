import pytest
from scripts.APICall import APICall

@pytest.mark.parametrize("endpoint", [
   f"commits?path=tests/logs/RegressionTests_ursa.log&per_page=1", #fetch_repo_commits_endpoint
   f"pulls/2882", #get_pr_head_endpoint
   f"contents/tests/logs/RegressionTests_ursa.log", #fetch_log_text_endpoint
   ])
@pytest.mark.parametrize("num_commits", [1, 5, 7])
def test_init_APICall(set_env_vars, monkeypatch, endpoint, num_commits):
   """
   Test initialization of the api_call.
   When running tests locally, create a GitHub token and set it as an environment variable.
   """
   
   set_env_vars
   # Set token env var for duration of test only
   monkeypatch.setenv("GITHUB_TOKEN", "fake_github_pat_12BWCMCFZkhj35klj3h34kjh4kkjm3whe4nr")
   api_call = APICall(endpoint, num_commits)

   assert api_call.token == "fake_github_pat_12BWCMCFZkhj35klj3h34kjh4kkjm3whe4nr"
   assert api_call.base_url == "https://api.github.com/repos/ufs-community/ufs-weather-model"
   assert api_call.endpoint == endpoint
   assert api_call.url == f"https://api.github.com/repos/ufs-community/ufs-weather-model/{endpoint}"
   assert api_call.num_commits == num_commits
   assert api_call.header == {
         "Accept": "application/vnd.github.v3+json",
         "Authorization": f"Bearer fake_github_pat_12BWCMCFZkhj35klj3h34kjh4kkjm3whe4nr",
         "X-GitHub-Api-Version": "2022-11-28",
         "Accept": "application/vnd.github.raw"
      }

def test_set_endpoint():
   """Initialize an APICall object with no endpoint; set the endpoint and check that it updated correctly."""
   api_call = APICall()
   assert api_call.endpoint == ''

   api_call.set_endpoint("commits")
   assert api_call.endpoint == 'commits'
   assert api_call.url == f"{api_call.base_url}/{api_call.endpoint}"

def test_call_API():
   """Test the API call's ability to fetch the contents of a log file.
   """
   api_call = APICall("contents/tests/logs/RegressionTests_orion.log")
   assert api_call.call_API().text.startswith(r"====START OF ORION REGRESSION TESTING LOG====")

def test_load_json_from_api_call():
   """Test the API call's ability to fetch and parse the last X number of commits for the repository.
   """
   api_call = APICall(num_commits=5)
   api_call.set_endpoint(f"commits?per_page={api_call.num_commits}")
   response = api_call.call_API()
   response = api_call.load_json_from_api_call(response)
   assert len(response) == 5
   for item in response:
      assert item['sha']
