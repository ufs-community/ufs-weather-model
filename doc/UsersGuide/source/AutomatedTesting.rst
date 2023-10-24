.. _AutomatedTesting:

*****************
Automated Testing
*****************

The UFS Weather Model repository on GitHub employs two types of automated testing:

   #. CI/CD (Continuous Integration/Continuous Development) testing on the cloud 
   #. AutoRT on NOAA R&D platforms

Both are application level tests and utilize the regression testing framework
discussed in :numref:`Section %s <UsingRegressionTest>`.

=====
CI/CD
=====

The UFS Weather Model (:term:`WM`) uses GitHub Actions (GHA), a GitHub-hosted continuous integration service, 
to perform CI/CD testing. Build jobs are done on GHA-provided virtual machines. Test jobs are 
performed on the Amazon Web Services (AWS) cloud platform using a number of EC2 instances. 
Builds and tests are carried out in a Docker container. The container includes a pre-installed version of the
:term:`HPC-Stack`, which includes all prerequisite libraries. Input data needed to run the tests 
are stored as a separate Docker container.

When a developer makes a pull request (PR) to the UFS WM repository, a code
manager may add the `run-ci` label, which triggers the CI/CD workflow. 
The CI/CD workflow then executes the following steps:

   #. A check is performed to make sure the UFS Weather Model and its first level
      subcomponents are up to date with the top of the ``develop`` branch.

   #. If the check is successful, build jobs are started on GHA-provided virtual machines
      by downloading the HPC-Stack Docker container stored in Docker Hub.

   #. Once all build jobs are successful, the created executable files are stored as
      artifacts in GHA.

   #. A number of AWS EC2 instances are started.

   #. Test jobs are started on AWS after downloading the HPC-Stack Docker container,
      the executable file from the build job, and the input-data Docker container.

   #. When all tests are complete, EC2 instances are stopped. Test results are reported
      on GitHub.


The GHA-related ``yaml`` scripts are located in the ``.github/workflows/`` directory.
``build_test.yml`` is the main workflow file, and ``aux.yml`` is an auxiliary
file responsible for (1) checking that the PR branch is up-to-date and 
(2) starting/stopping the EC2 instances. 

Other CI-related scrips are located in the ``tests/ci/`` directory. ``ci.sh`` is the main script that 
invokes Docker build and run. ``Dockerfile`` is used to build the UFS Weather Model. 
Other shell and python scripts help with various tasks. For example:

   * ``repo_check.sh`` checks that the PR branch is up-to-date.
   * ``check_status.py`` checks the status of EC2 instances.
   * ``setup.py`` and ``ci.test`` configure the test cases to execute in the CI/CD workflow.

.. COMMENT: It sounds like aux.yml and repo_check.sh do the same thing... What's the difference?

=======
Auto RT
=======

The Automated Regression Testing (AutoRT) system is a python program that automates the process 
of regression testing on NOAA HPC platforms. 
It contains the files in :numref:`Table %s <autoRT-files>` below:

.. _autoRT-files:
.. table:: *Files for Automated Regression Testing (AutoRT) system*

   +-------------------+-----------------------------------------------------+
   | **File Name**     | **Description**                                     |
   +===================+=====================================================+
   |  start_rt_auto.sh | Verifies HPC name, sets the python paths            |
   +-------------------+-----------------------------------------------------+
   |  rt_auto.py       | Python interface between the HPC and the github API |
   +-------------------+-----------------------------------------------------+
   |  jobs/bl.py       | Functions for the baseline job                      |
   +-------------------+-----------------------------------------------------+
   |  jobs/rt.py       | Functions for the regression test job               |
   +-------------------+-----------------------------------------------------+

-----------------
AutoRT Workflow
-----------------

On supported HPC systems, a :term:`cron job` runs the ``start_rt_auto.sh`` bash script every 15 minutes. 
This script checks the HPC name and sets certain python paths. Then, it runs ``rt_auto.py``, 
which uses the Github API (through pyGitHub) to check the labels on pull requests to 
``ufs-weather-model``. If a PR label matches the HPC name 
(e.g., hera-intel-RT or cheyenne-gnu-BL), the label provides the HPC  
with the compiler and job information to run a test or task on the machine. 
If no PR label matches HPC name, the script exits.

For example, a PR labeled ``gaea-intel-BL`` will be recognized by the HPC machine 'Gaea'. 
It will set the ``RT_COMPILER`` variable to 'intel' and run the baseline creation script (``bl.py``).
This script creats a job class that contains all information from the machine that the job will need to run.
That information is sent into the ``jobs/rt[bl].py`` script. 

``rt.py`` sets directories for storage, gets repo information, runs the regression test, and 
completes any required post processing.

.. code-block:: python3

   def run(job_obj):
      logger = logging.getLogger('RT/RUN')
      workdir = set_directories(job_obj)
      branch, pr_repo_loc, repo_dir_str = clone_pr_repo(job_obj, workdir)
      run_regression_test(job_obj, pr_repo_loc)
      post_process(job_obj, pr_repo_loc, repo_dir_str, branch)

``bl.py``: (similar to ``rt.py``) Adds functionality to create baselines before running regression testing.

.. code-block:: python3
   :emphasize-lines: 5,6,7

      def run(job_obj):
         logger = logging.getLogger('BL/RUN')
         workdir, rtbldir, blstore = set_directories(job_obj)
         pr_repo_loc, repo_dir_str = clone_pr_repo(job_obj, workdir)
         bldate = get_bl_date(job_obj, pr_repo_loc)
         bldir = f'{blstore}/develop-{bldate}/{job_obj.compiler.upper()}'
         bldirbool = check_for_bl_dir(bldir, job_obj)
         run_regression_test(job_obj, pr_repo_loc)
         post_process(job_obj, pr_repo_loc, repo_dir_str, rtbldir, bldir)
