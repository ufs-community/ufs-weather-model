.. _ContributingDevelopment:
  
*************************
Contributing Development
*************************

* Overview; big picture

* Opening issues

* Description of development workflow (fork, branch, commit, test, PR submission, review)

* What are the requirements for a code to be accepted?

* What platforms are supported

================================
 ufs-community/ufs-weather-model
================================

-----------------------
 Prerequisite libraries:
-----------------------

Below is the list of libraries that are required to compile the model and run regression tests:

+-------------+---------------------------------------------------------------------+---------+----------------+
| Library name| source	 			                                    |version  | env variables  |
+=============+=====================================================================+=========+================+
| bacio	      |https://github.com/NOAA-EMC/NCEPLIBS-bacio                           | 2.0.3   | $BACIO_LIB4    |
+-------------+---------------------------------------------------------------------+---------+----------------+
| sp	      |https://github.com/NOAA-EMC/NCEPLIBS-sp                              | 2.0.3   | $SP_LIBd       |
+-------------+---------------------------------------------------------------------+---------+----------------+
| w3nco	      |https://github.com/NOAA-EMC/NCEPLIBS-w3nco                           | 2.0.7   | $W3NCO_LIBd    |
+-------------+---------------------------------------------------------------------+---------+----------------+
| w3emc       |https://github.com/NOAA-EMC/NCEPLIBS-w3emc                           | 2.3.1   | $W3EMC_LIBd    |
+-------------+---------------------------------------------------------------------+---------+----------------+
| nemsio      |https://github.com/NOAA-EMC/NCEPLIBS-nemsio                          | 2.2.4   | $NEMSIO_INC    |
+             +                                                                     +         +                +
|	      |					                                    |         | $NEMSIO_LIB    |
+-------------+---------------------------------------------------------------------+---------+----------------+
| g2	      |https://github.com/NOAA-EMC/NCEPLIBS-g2	                            |3.1.1    | $G2_LIB4       |
+-------------+---------------------------------------------------------------------+---------+----------------+
| g2tmpl      |https://github.com/NOAA-EMC/NCEPLIBS-g2tmpl(required for inline post)|1.5.1    | $G2TMPL_LIB    |
+-------------+---------------------------------------------------------------------+---------+----------------+
| crtm	      |https://github.com/NOAA-EMC/EMC_crtm (required for inline post)      |2.2.6    | $CRTM_LIB      |
+-------------+---------------------------------------------------------------------+---------+----------------+
| post	      |	http://github.com/NOAA-EMC/EMC_post (required for inline post)	    | 8.0.0   | $POST_INC      |
+             +                                                                     +         +                +
|	      |					                                    |         | $POST_LIB      |
+-------------+---------------------------------------------------------------------+---------+----------------+

External libraries:


+-------------+--------------------------------------------------------------+---------+-------------------------+
| Library name| source	 			                             |version  | env variables           |
+=============+==============================================================+=========+=========================+
| 	      |                                                              |         |$ESMFMKFILE,             |
+             +                                                              +         +                         +
| 	      |                                                              |         |$ESMF_F90COMPILEPATHS,   |
+             +                                                              +         +                         +
|esmf         |https://www.earthsystemcog.org/projects/esmf/download_800     | 8.0.0   |$ESMF_F90ESMFLINKPATHS,  |
+             +                                                              +         +                         +
| 	      |                                                              |         |$(ESMF_F90ESMFLINKRPATHS)|
+             +                                                              +         +                         +
| 	      |                                                              |         |$(ESMF_F90ESMFLINKLIBS)  |
+-------------+--------------------------------------------------------------+---------+-------------------------+
| 	      |                                                              |  4.5.0  |                         |
+             +                                                              +         +                         +
| NetCDF      |https://www.unidata.ucar.edu/downloads/netcdf/index.jsp       |  4.6.1  |$NETCDF or($NETCDF_DIR)  |
+             +                                                              +         +                         +
| 	      |                                                              |  4.6.3  |                         |
+             +                                                              +         +                         +
| 	      |                                                              |  4.7.0  |                         |
+-------------+--------------------------------------------------------------+---------+-------------------------+
| 	      |                                                              |  1.8.14 |                         |
+             +                                                              +         +                         +
|HDF5         | https://www.hdfgroup.org/downloads/hdf5/                     |  1.8.16 | $HDF5                   |
+             +                                                              +         +                         +
| 	      |                                                              |  1.10.1 |                         |
+             +                                                              +         +                         +
| 	      |                                                              |  1.10.4 |                         |
+-------------+--------------------------------------------------------------+---------+-------------------------+
|jasper       | https://github.com/mdadams/jasper                            | 1.900.1 | $JASPER_LIB             |
+-------------+--------------------------------------------------------------+---------+-------------------------+
|zlib 	      | https://github.com/madler/zlib                               | 1.2.11  | $Z_LIB                  |
+-------------+--------------------------------------------------------------+---------+-------------------------+
|libpng       | http://www.libpng.org/pub/png/libpng.html                    | 1.2.44  | $PNG_LIB                |
+-------------+--------------------------------------------------------------+---------+-------------------------+


Compiler:


+--------------------+------------------------------------------------------------------+
|Compiler name       |    Branch name                                                   |
+====================+==================================================================+
| Intel              |   16.0.3.210,16.3.210,18.0.1.163, 18.0.3.222,18.0.5.274, 19.0.2  |
+--------------------+------------------------------------------------------------------+
|impi                |   2018.0.4, 2018.4.274,18.0.1, 2019.2.187                                                        |
+--------------------+------------------------------------------------------------------+


These libraries are available on several platforms including Hera, Cheyenne and WCOSS. Installing these libraries on jet and gaea is a work in progress.

In addition to these libraries, other software is required such as a compiler (Intel or GNU), MPI library (impi), HDF5, cmake, Python, Perl and Ruby to compile the model and run regression tests.

-----------------------
Making code changes using a forking workflow
-----------------------

If developers would like to make code changes, they need to make a personal fork, set up upstream remote (for merging with the original ufs-weather-model), and create a branch for ufs-weather-model and each of the subcomponent repositories they want to change. They can then make code changes, perform testing and commit the changes to the branch in their personal fork. It is suggested that they merge their branch with the develop branch of the original repositories periodically to get the latest updates and bug fixes.

If developers would like to get their code committed back to the original repository, it is suggested to follow the steps below:

      # Create an issue in the authoritative repository. For example to commit code changes to fv3atm, please go to https://github.com/NOAA-EMC/fv3atm, under NOAA-EMC/fv3atm and find the “Issues” tab next to the “Code” tab. Click on “Issues” and a new page will appear. On the right side of the page, there is a green “New issue” button. Clicking on that will lead to a new issue page. Fill out the title, comments to describe the code changes, and also please provide personal fork and branch information. Lastly, click on the “Submit new issue” button, so that the new issue is created.

      # When the development is mature, tests have been conducted, and the developer is satisfied with the results, create a pull request to commit the code changes.

      * Merge developer’s branch to the latest ufs-weather-model develop branch in authoritative repository. If changes are made in model sub-components, developers need to merge their branches to branches with the corresponding authoritative repository (or original repository for some components). For this, code management practices of the subcomponents need to be followed.

      * Regression tests need to pass on at least one supported platform.

      * For each component branch where developers make changes, developers need to go to their personal fork on GitHub and click on the “New pull request” button. When a new page “Compare changes” appears, developers will choose the branch in their fork with code changes to commit and the branch in upstream repository that the changes will be committed to. The code differences between the two branches will be displayed. Developers can review the differences and click on “submit pull request” to make the pull request. After code changes are committed to the component repository, developers will make pull requests to ufs-weather-model repository.

It is suggested that the developers inform all the related code managers as the hierarchy structure of the ufs-weather-model repository may require collaboration among the code managers.

-----------------------
Engaging in the code review process
-----------------------

When code managers receive a pull request to commit the code changes, it is recommended that they add at least two code reviewers to review the code. The reviewers will write comments about the code changes and give a recommendation as to whether the code changes can be committed. What kinds of code changes will be accepted in the repository is beyond the scope of this document; future ufs-weather-model code management may have detailed answer for that.

Reviewers may suggest some code changes during the review process. Developers need to respond to these comments in order to get code changes committed. If developers make further changes to their branch, reviewers need to check the code changes again. When both reviewers give recommendation to commit the code, code managers will merge the changes into the repository.

-----------------------
Conducting regression tests
-----------------------

Only developers that are running on a limited set of platforms (Hera, Cheyenne, WCOSS) can compile and run regression tests using the ufs-weather-model.

To run regression test using rt.sh

.. code-block:: console

   % cd ufs-weather-model/tests
   % ./rt.sh -f

Regression test log files (ufs-weather-model/tests/Compile_$(MACHINE_ID).log and ufs-weather-model/tests/RegressionTests_$(MACHINE_ID).log ) will be updated.

To create new baseline:

.. code-block:: console

   % cd ufs-weather-model/tests
   % ./rt.sh -f -c

      * To run regression test using NEMSCompsetRun

.. code-block:: console

   % cd ufs-weather-model
   % ./NEMS/NEMSCompsetRun -f 

Regression test log files (ufs-weather-model/log/$MACHINE_ID/* ) will be updated.

To create new baseline:

.. code-block:: console

   % cd ufs-weather-model
   % ./NEMS/NEMSCompsetRun--baseline fv3 --platform=${PLATFORM}


The value of ${PLATFORM} can be found in ufs-weather-model/compsets/platforms.input.

Developers need to commit the regression test log files to their branch before making pull request.

-----------------------
Compiling the code and running a test
-----------------------

Currently developers running on Hera, Cheyenne or WCOSS can compile and run tests using the ufs-weather-model. In the document below, ufs-weather-model directory points to a branch in the developer’s personal fork.

      * Compile and run a test using rt.sh

	i) compile the code

Developers can compile the ufs-weather-model using the compile script ufs_weather_model/tests/compile.sh. To compile the code:

.. code-block:: console

   % git clone --recursive https://github.com/ufs-community/ufs-weather-model
   % cd ufs-weather-model/tests
   % ./compile.sh PATHTR MACHINE_ID MAKE_OPT BUILD_NR


Where:

.. code-block:: console

   PATHTR: the full path of FV3 directory under ufs-weather-model
   MACHINE_ID: machine ID, e.g: wcoss_cray, hera.intel
   MAKE_OPT: compile options, default (‘’, empty string) is for 64-bit OpenMP non-hydrostatic build using AVX2, other options are:
   -       ‘DEBUG=Y’: turn on debug option
   -       ‘VERBOSE=Y’: turn on VERBOSE mode to get additional details on compile
   -       ‘OPENMP=Y’: use openmp
   -       ‘AVX2=Y’: use AVX2 in Intel Haswell for better performance.
   -       ‘HYDRO=Y’: hydrostatic mode
   -       ‘CCPP=Y’: using ccpp framework and physics
   -       ‘STATIC=Y’: when CCPP=Y, using static build instead of dynamic build
   -       ‘SUITE=xxx': CCPP physics suite, e.g. “FV3_GFS_2017_gfdlmp”
   BUILD_NR: the number of build (there might be several copies of the executable). The final executable would be fv3_${BUILD_NR}.exe


Example:

.. code-block:: console

   %  ./compile.sh /gpfs/hps/emc/global/noscrub/First.Last/ufs-weather-model/FV3 wcoss_cray '' 1

The executable generated from the compile.sh can be used in global workflow to run experiments.

	ii) run a regression test

Developers can run a single regression test from the regression test suite. The full list of the regression tests can be found at:

.. code-block:: console

   us-weather-model/tests/fv3_conf

To run a specific regression test:

.. code-block:: console

   %cd ufs_weather_model/tests
   % cp rt.conf rt.conf1

Edit the rt.conf1 file, keep the test developers intended to run and remove all the rest tests.

COMPILE | fv3                           | standard| wcoss_cray  |             |
RUN     | fv3_control                   | standard|             |             |


To compile and run the test, do

.. code-block:: console

   % ./rt.sh -l rt.conf1

The code will be compiled and run on wcoss_cray. A log directory will be shown at: ufs-weather-model/tests /log_$MACHINE_ID. If there are compile errors, please check file: compile_${BUILD_NR}.log under above log directory. If the code is successfully compiled, but the job failed, please go to above log directory to look for rt_${test_number}_${test_id}.log for details.

      * Compile and run a test using NEMSCompsetRun

	i) Compile the code

.. code-block:: console

   % cd ufs-weather-model
   % ./NEMS/NEMSAppBuilder app=coupledFV3_WW3

ii) Run a regression test

NEMSCompsetRun runs the same tests as rt.sh. The two share the same baseline results. To run a single test using NEMSCompsetRun:

.. code-block:: console

   % cd ufs-weather-model
   % ./NEMS/NEMSCompsetRun “{test_name}”

The test name can be found in compsets/all.input.
























