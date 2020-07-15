.. _ContributingDevelopment:

*************************
Contributing Development
*************************

The ufs-weather-model repository contains the model code and external links needed to build the Unified Forecast System (UFS) atmosphere model and associated components, including the WaveWatch III model. This weather model is used in several of the UFS applications, including the medium-range weather application, the short-range weather application, and the sub-seasonal to seasonal application.

---------------------------------------------
Making code changes using a forking workflow
---------------------------------------------

If developers would like to make code changes, they need to make a personal fork, set up upstream remote (for merging with the original ufs-weather-model), and create a branch for ufs-weather-model and each of the subcomponent repositories they want to change. They can then make code changes, perform testing and commit the changes to the branch in their personal fork. It is suggested that they update their branch by merging the develop branch with the develop branch of the original repositories periodically to get the latest updates and bug fixes.

If developers would like to get their code committed back to the original repository, it is suggested to follow the steps below:

1. Create an issue in the authoritative repository. For example to commit code changes to fv3atm, please go to https://github.com/NOAA-EMC/fv3atm, under NOAA-EMC/fv3atm and find the “Issues” tab next to the “Code” tab. Click on “Issues” and a new page will appear. On the right side of the page, there is a green “New issue” button. Clicking on that will lead to a new issue page. Fill out the title, comments to describe the code changes, and also please provide personal fork and branch information. Lastly, click on the “Submit new issue” button, so that the new issue is created.

2. When the development is mature, tests have been conducted, and the developer is satisfied with the results, create a pull request to commit the code changes.

      * Merge developer’s branch to the latest ufs-weather-model develop branch in authoritative repository. If changes are made in model sub-components, developers need to merge their branches to branches with the corresponding authoritative repository (or original repository for some components). For this, code management practices of the subcomponents need to be followed.

      * Regression tests associated with the ufs-weather-model are available on Tier 1 and Tier 2 platforms as described in https://github.com/ufs-community/ufs-weather-model/wiki/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers. If the developer has access to these platforms, the developer should pass the regression test on at least one supported platform. If the developer does not have access to these platforms, this should be stated in the PR so the code manager(s) can conduct the tests.

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

When code managers receive a pull request to commit the code changes, it is recommended that they add at least two code reviewers to review the code and at least one of the reviewers has write permission. The reviewers will write comments about the code changes and give a recommendation as to whether the code changes can be committed. What kinds of code changes will be accepted in the repository is beyond the scope of this document; future ufs-weather-model code management documents may have a detailed answer for that.

Reviewers may suggest some code changes during the review process. Developers need to respond to these comments in order to get code changes committed. If developers make further changes to their branch, reviewers need to check the code changes again. When both reviewers give recommendation to commit the code, code managers will merge the changes into the repository.

.. _ConductingRegTests:

----------------------------
Conducting regression tests
----------------------------

Only developers using Tier 1 and 2 platforms can run the ufs-weather-model regression tests. Other developers need to work with the code managers to assure completion of the regression tests.

To run regression test using rt.sh

rt.sh is a bash shell file to run the RT and has the following options:

.. code-block:: console

   Usage: ./rt.sh -c | -f | -s | -l <file> | -m | -k | -r | -e | -h
   -c create new baseline results for <model>
   -f run full suite of regression tests
   -s run standard suite of regression tests
   -l run test specified in <file>
   -m compare against new baseline results
   -k  keep run directory (automatically deleted otherwise if all tests pass)
   -r use Rocoto workflow manager
   -e use ecFlow workflow manager
   -h display this help

.. code-block:: console

   % cd ufs-weather-model/tests
   % ./rt.sh -f

This command can only be used on platforms that have been configured for regression testing (Tier 1 and Tier 2 platforms as described in https://github.com/ufs-community/ufs-weather-model/wiki/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers). For information on testing the CCPP code, or using alternate computational platforms, see the following sections.

This command and all others below produce log output in ./tests/log_machine.compiler. These log files contain information on the location of the run directories that can be used as templates for the user. Each rt*.conf contains one or more compile commands preceding a number of tests.

Regression test log files (ufs-weather-model/tests/Compile_$(MACHINE_ID).log and ufs-weather-model/tests/RegressionTests_$(MACHINE_ID).log ) will be updated.

If developers wish to contribute code that changes the results of the regression tests (because of updates to the physics, for example), it is useful to run rt.sh as described above to make sure that the test failures are as expected. It is then useful to establish a new personal baseline:

./rt.sh -l rt.conf -c # create own reg. test baseline

Once the personal baseline has been created, future runs of the RT should be compared against the personal baseline using the -m option.

./rt.sh -l rt.conf -m # compare against own baseline

To create new baseline:

.. code-block:: console

   % cd ufs-weather-model/tests
   % ./rt.sh -f -c

An alternative/complementary regression test system is using NEMSCompsetRun, which focuses more on coupled model configurations than testing features of the standalone ufs-weather-model. To run regression test using NEMSCompsetRun:

.. code-block:: console

   % cd ufs-weather-model
   % ./NEMS/NEMSCompsetRun -f

Regression test log files (ufs-weather-model/log/$MACHINE_ID/* ) will be updated.

To create new baseline:

.. code-block:: console

   % cd ufs-weather-model
   % ./NEMS/NEMSCompsetRun --baseline fv3 --platform=${PLATFORM}

The value of ${PLATFORM} can be found in ufs-weather-model/compsets/platforms.input.

Developers need to commit the regression test log files to their branch before making pull request.
