.. _ContributingDevelopment:
  
*************************
Contributing Development
*************************

The ufs-weather-model repository contains the model code and external links needed to build the Unified Forecast System (UFS) atmosphere model and associated components, including the WaveWatch III model. This weather model is used in several of the UFS applications, including the medium-range weather application, the short-range weather application, and the sub-seasonal to seasonal application.

------------------------
 Prerequisite libraries:
------------------------

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
|impi                |   2018.0.4, 2018.4.274,18.0.1, 2019.2.187                        |
+--------------------+------------------------------------------------------------------+


These libraries are available on several platforms including Hera, Cheyenne and WCOSS. Installing these libraries on jet and gaea is a work in progress. Users who need to build the libraries. [https://github.com/NOAA-EMC/NCEPLIBS/wiki/Cloning-and-Compiling-NCEPLIBS]

In addition to these libraries, other software is required such as a compiler (Intel or GNU), MPI library (impi), HDF5, cmake, Python, Perl and Ruby to compile the model and run regression tests.

---------------------------------------------
Making code changes using a forking workflow
---------------------------------------------

If developers would like to make code changes, they need to make a personal fork, set up upstream remote (for merging with the original ufs-weather-model), and create a branch for ufs-weather-model and each of the subcomponent repositories they want to change. They can then make code changes, perform testing and commit the changes to the branch in their personal fork. It is suggested that they merge their branch with the develop branch of the original repositories periodically to get the latest updates and bug fixes.

If developers would like to get their code committed back to the original repository, it is suggested to follow the steps below:

1. Create an issue in the authoritative repository. For example to commit code changes to fv3atm, please go to https://github.com/NOAA-EMC/fv3atm, under NOAA-EMC/fv3atm and find the “Issues” tab next to the “Code” tab. Click on “Issues” and a new page will appear. On the right side of the page, there is a green “New issue” button. Clicking on that will lead to a new issue page. Fill out the title, comments to describe the code changes, and also please provide personal fork and branch information. Lastly, click on the “Submit new issue” button, so that the new issue is created.

2. When the development is mature, tests have been conducted, and the developer is satisfied with the results, create a pull request to commit the code changes.

      * Merge developer’s branch to the latest ufs-weather-model develop branch in authoritative repository. If changes are made in model sub-components, developers need to merge their branches to branches with the corresponding authoritative repository (or original repository for some components). For this, code management practices of the subcomponents need to be followed.

      * Regression tests associated with the ufs-weather-model are available in Tier 1 and Tier 2 platforms as described in https://github.com/ufs-community/ufs-weather-model/wiki/Weather-Model-Platform-and-Compiler-Support. If the developer has access to these platforms, the developer should pass the regression test on at least one supported platform. If the developer does not have access to these platforms, this should be stated in the PR so the code manager(s) can conduct the tests.

      * For each component branch where developers make changes, developers need to go to their personal fork on GitHub and click on the “New pull request” button. When a new page “Compare changes” appears, developers will choose the branch in their fork with code changes to commit and the branch in upstream repository that the changes will be committed to. Also developers in the commit comment must add the github issue title and number created in 1) in the comment box. The code differences between the two branches will be displayed. Developers can review the differences and click on “submit pull request” to make the pull request. After code changes are committed to the component repository, developers will make pull requests to ufs-weather-model repository.

3. When PRs are created, the creator must temporarily modify .gitmodules to point to his/her fork and branch if updates are required for submodules.

4. Merging code from PRs with submodules requires coordination with the person making the PRs. From the "innermost" nested PR up to the top-level PR, the PRs need to be merged as-is. After each merge, the person creating the PRs has to update his/her local code to check out the merged version, revert the change to .gitmodules, and push this to GitHub to update the PR. And so on and so forth.

5. Checking out the code ufs_release_1.0 should always be as follows:

.. code-block:: console

   git clone https://github.com/ufs-community/ufs-weather-model
   cd ufs-weather-model
   git checkout ufs_release_1.00
   git submodule update --init --recursive

6. Checking out a PR with id ID for testing it should always be as follows:

.. code-block:: console

   git clone https://github.com/ufs-community/ufs-weather-model
   cd ufs-weather-model
   git fetch origin pull/ID/head:BRANCHNAME
   git checkout BRANCHNAME
   git submodule update --init --recursive

It is suggested that the developers inform all the related code managers as the hierarchy structure of the ufs-weather-model repository may require collaboration among the code managers.

-----------------------------------
Engaging in the code review process
-----------------------------------

When code managers receive a pull request to commit the code changes, it is recommended that they add at least two code reviewers to review the code and at least one of the reviewers has write permission. The reviewers will write comments about the code changes and give a recommendation as to whether the code changes can be committed. What kinds of code changes will be accepted in the repository is beyond the scope of this document; future ufs-weather-model code management may have detailed answer for that.

Reviewers may suggest some code changes during the review process. Developers need to respond to these comments in order to get code changes committed. If developers make further changes to their branch, reviewers need to check the code changes again. When both reviewers give recommendation to commit the code, code managers will merge the changes into the repository.

----------------------------
Conducting regression tests
----------------------------

Only developers that are running on a limited set of platforms (Hera, Cheyenne, WCOSS) can compile and run regression tests using the ufs-weather-model.

To run regression test using rt.sh

rt.sh is a bash shell file to run the RT and has the following options:

.. code-block:: console

   Usage: $0 -c <model> | -f | -s | -l <file> | -m | -r | -e | -h
   -c create new baseline results for <model>
   -f run full suite of regression tests
   -s run standard suite of regression tests
   -l run test specified in <file>
   -m compare against new baseline results
   -r use Rocoto workflow manager
   -e use ecFlow workflow manager
   -h display this help

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

