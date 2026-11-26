import os
import requests
import json

class APICall():
   """A GitHub API call"""

   def __init__(self, endpoint='', num_commits=1):
      """Configure a GitHub API request.

      Args:
         endpoint (str): API path relative to ``BASE_URL``; defaults to an empty path.
         num_commits (int): Number of commits associated with the request.

      Returns:
         None.
      """
      self.token = os.environ.get('GITHUB_TOKEN')
      self.base_url = os.environ.get('BASE_URL')
      self.endpoint = endpoint
      self.url = f"{self.base_url}/{self.endpoint}" #Could use a path join?
      self.num_commits = num_commits
      self.header = {
         "Authorization": f"Bearer {self.token}",
         "X-GitHub-Api-Version": "2022-11-28",
         "Accept": "application/vnd.github.raw"
      }
   
   def set_endpoint(self, endpoint):
      """Set the API path and rebuild the request URL.

      Args:
         endpoint (str): API path relative to the configured base URL.

      Returns:
         None.
      """
      self.endpoint = endpoint
      self.url = f"{self.base_url}/{self.endpoint}"

   def call_API(self):
      """Send the configured GET request.

      Returns:
         requests.Response: The HTTP response returned by ``requests.get``.
      """
      response = requests.get(self.url, headers=self.header)
      
      return response
   
   def load_json_from_api_call(self, response):
      """Decode a JSON response body.

      Args:
         response (requests.Response): HTTP response whose ``text`` is JSON.

      Returns:
         object: The decoded JSON value; the precise information depends on the endpoint.
      """

      results = json.loads(response.text)

      return results
