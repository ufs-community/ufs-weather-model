.. _AutomatedTesting:

*****************
Automated Testing
*****************

The UFS Weather Model repository on GitHub employs two types of automated testing.
One is the CI/CD on cloud and the other is the AutoRT on NOAA R&D platforms.
Both are application level tests, and utilize the regression testing framework
discussed in :numref:`Section %s <UsingRegressionTest>`.

=====
CI/CD
=====

The following summarizes the CI/CD used in the UFS Weather Model:

* GitHub Actions (GHA), a GitHub-hosted continuous integration service is used.
* Build jobs are done on GHA-provided virtual machines.
* Test jobs are performed on Amazon cloud using a number of EC2 instances.
* Builds and tests are carried out using a Docker container.
* Docker container has prerequisite libraries installed via the hpc-stack.
* Input data needed to run tests are stored as a separate Docker container.


When a developer makes a pull request (PR) to the UFS Weather Model repository, and a code
manager subsequently adds the `run-ci` label, the CI/CD workflow is triggerd:

#. A check is performed to make sure the UFS Weather Model and its first level
   subcomponents are up to date with the top of develop branch.

#. If the check is successful, build jobs are started on GHA-provided virtual machines
   by downloading the hpc-stack Docker container stored in Docker Hub.

#. Once all build jobs are successful, the created executable files are stored as
   artifacts in GHA.

#. A number of AWS EC2 instances are started.

#. Test jobs are started on Amazon cloud by downloading the hpc-stack Docker container,
   the executable file from the build job, and the input-data Docker container.

#. When all tests are finished, EC2 instances are stopped. Test results are reported
   on GitHub.


The GHA-related yaml scripts are located in the ``.github/workflows/`` directory.
``build_test.yml`` is the main workflow file, and ``aux.yml`` is an auxiliary
file responsible for checking the up-to-dateness of the PR branch, and starting
and stopping the EC2 instances. Other CI-related scrips are located in the ``tests/ci/``
directory. ``ci.sh`` is the main script that invokes Docker build and run. ``Dockerfile``
is used to build UFS Weather Model. Other shell and python scripts help with various
tasks such as checking the up-to-dateness of the PR branch (``repo_check.sh``),
checking the status of EC2 instances (``check_status.py``), and configuring the test cases
to carry out in the CI/CD workflow (``setup.py`` and ``ci.test``).


=======
Auto RT
=======

The Automated Regression Testing (AutoRT) system:

* Automates the process of regression testing on NOAA HPC platforms.

* Written in python.

* Contains the following files:

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

~~~~~~~~~~~~~~~
AutoRT Workflow
~~~~~~~~~~~~~~~
* Cron-job on supported HPC systems runs start_rt_auto.sh bash script every
  15 minutes.

  * This script verifies the HPC name, and sets the python paths. Runs
    rt_auto.py.

* rt_auto.py: Uses the Github API (Through pyGitHub)

  * Checks the pull requests to ufs-community/ufs-weather-model for
    labels specific to the HPC name. If no match to HPC name, exits.
    (i.e. hera-intel-RT or cheyenne-gnu-BL)

  * If the HPC name matches the label in ufs-weather-model pull
    request, the label provides the HPC with the compiler and job to run on
    the machine.

    * For example the label gaea-intel-BL will be recognized by the HPC
      machine 'Gaea', set the RT_COMPILER variable to 'intel' and run the
      baseline creation script (bl.py).

  * Creates a Job class that contains all information from the machine
    that the job will need to run. That is sent into the jobs/rt[bl].py script.

* rt.py: Sets directories for storage, gets repo information, runs RT,
  post processes.

.. code-block:: python3

    def run(job_obj):
        logger = logging.getLogger('RT/RUN')
        workdir = set_directories(job_obj)
        branch, pr_repo_loc, repo_dir_str = clone_pr_repo(job_obj, workdir)
        run_regression_test(job_obj, pr_repo_loc)
        post_process(job_obj, pr_repo_loc, repo_dir_str, branch)

* bl.py: (similar to rt.py) Adds functionality to create baselines before
  running regression testing.

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
