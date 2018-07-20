################################################################################
# Theia/Intel
################################################################################

# OFFICIAL BASELINE
# /scratch4/NCEPDEV/nems/noscrub/emc.nemspara/RT/NEMSfv3gfs/trunk-20180524

# OWN BASELINE
# /scratch4/NCEPDEV/stmp4/$USER/FV3_RT/REGRESSION_TEST_INTEL

# RUN DIRS - substitute $USER with your username
# /scratch4/NCEPDEV/stmp3/$USER/FV3_RT/rt_...

export NEMS_COMPILER=intel # optional, default is intel if not set
export ACCNR=...           # set account name/number for job scheduler

# from the top-level directory of NEMSfv3gfs
cd tests

# Test standard compile
./compile.sh $PWD/../FV3 theia.intel 2>&1 | tee log.compile

# Regression tests against official baseline
./rt.sh -f 2>&1 | tee rt_full.log

################################################################################
# Theia/GNU
################################################################################

# OFFICIAL BASELINE
# N/A

# OWN BASELINE
# /scratch4/NCEPDEV/stmp4/$USER/FV3_RT/REGRESSION_TEST_GNU

# RUN DIRS
# /scratch4/NCEPDEV/stmp3/$USER/FV3_RT/rt_...

export NEMS_COMPILER=gnu
export ACCNR=...         # set account name/number for job scheduler

# from the top-level directory of NEMSfv3gfs
cd tests

# Test standard compile
./compile.sh $PWD/../FV3 theia.gnu 2>&1 | tee log.compile

# Create new baseline
./rt.sh -f -l rt_gnu_pgi.conf -c fv3 2>&1 | tee rt_create_sub.log

# Regression tests against new baseline
./rt.sh -f -l rt_gnu_pgi.conf -m 2>&1 | tee rt_sub.log

################################################################################
# Theia/PGI
################################################################################

# OFFICIAL BASELINE
# N/A

# OWN BASELINE
# /scratch4/NCEPDEV/stmp4/$USER/FV3_RT/REGRESSION_TEST_PGI

# RUN DIRS
# /scratch4/NCEPDEV/stmp3/$USER/FV3_RT/rt_...

export NEMS_COMPILER=pgi
export ACCNR=...           # set account name/number for job scheduler

# from the top-level directory of NEMSfv3gfs
cd tests

# Test standard compile
./compile.sh $PWD/../FV3 theia.pgi 2>&1 | tee log.compile

# Create new baseline
./rt.sh -f -l rt_gnu_pgi.conf -c fv3 2>&1 | tee rt_create_sub.log

# Regression tests against new baseline
./rt.sh -f -l rt_gnu_pgi.conf -m 2>&1 | tee rt_sub.log
