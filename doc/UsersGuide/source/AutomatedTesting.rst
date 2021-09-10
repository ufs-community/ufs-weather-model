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

* GitHub Actions (GHA), a GitHub-hosted continuous integration service is used
* Build jobs are done on GHA-provided virtual machines
* Test jobs are performed on Amazon cloud using a number of EC2 instances
* Builds and tests are carried out using a Docker container
* Docker container has prerequisite libraries installed via the hpc-stack
* Input data needed to run tests are stored as a separate Docker container


When a developer makes a pull request (PR) to the UFS Weather Model repository, and a code
manager subsequently adds the `run-ci` label, the CI/CD workflow is triggerd:

#. A check is performed to make sure the UFS Weather Model and its first level
   subcomponents are up to date with the top of develop branch

#. If the check is successful, build jobs are started on GHA-provided virtual machines
   by downloading the hpc-stack Docker container stored in Docker Hub

#. Once all build jobs are successful, the created executable files are stored as
   artifacts in GHA

#. A number of AWS EC2 instances are started

#. Test jobs are started on Amazon cloud by downloading the hpc-stack Docker container,
   the executable file from the build job, and the input-data Docker container

#. When all tests are finished, EC2 instances are stopped. Test results are reported
   on GitHub


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
