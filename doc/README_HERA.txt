################################################################################
# Hera/Intel
################################################################################

# OFFICIAL BASELINE
# /scratch1/NCEPDEV/nems/emc.nemspara/RT/NEMSfv3gfs/trunk-20190824

# OWN BASELINE
# /scratch1/NCEPDEV/stmp4/$USER/FV3_RT

# RUN DIRS - substitute $USER with your username
#  /scratch1/NCEPDEV/stmp2/$USER/FV3_RT/rt_...

export NEMS_COMPILER=intel # optional, default is intel if not set
export ACCNR=...           # set account name/number for job scheduler

# from the top-level directory of NEMSfv3gfs
cd tests

# Test standard compile
./compile.sh $PWD/../FV3 hera.intel 2>&1 | tee log.compile

# Regression tests against official baseline
./rt.sh -f 2>&1 | tee rt_full.log

################################################################################
# Hera/GNU
################################################################################

# OFFICIAL BASELINE
# N/A

# OWN BASELINE
# /scratch1/NCEPDEV/stmp4/$USER/FV3_RT/REGRESSION_TEST_GNU

# RUN DIRS
# /scratch1/NCEPDEV/stmp2/$USER/FV3_RT/rt_...

export NEMS_COMPILER=gnu
export ACCNR=...         # set account name/number for job scheduler

# from the top-level directory of NEMSfv3gfs
cd tests

# Test standard compile
./compile.sh $PWD/../FV3 hera.gnu 2>&1 | tee log.compile

# Create new baseline
./rt.sh -f -l rt_gnu_pgi.conf -c fv3 2>&1 | tee rt_create_sub.log

# Regression tests against new baseline
./rt.sh -f -l rt_gnu_pgi.conf -m 2>&1 | tee rt_sub.log

################################################################################
# Hera/PGI
################################################################################

# OFFICIAL BASELINE
# N/A

# OWN BASELINE
# /scratch1/NCEPDEV/stmp4/$USER/FV3_RT/REGRESSION_TEST_PGI

# RUN DIRS
# /scratch1/NCEPDEV/stmp2/$USER/FV3_RT/rt_...

export NEMS_COMPILER=pgi
export ACCNR=...           # set account name/number for job scheduler

# from the top-level directory of NEMSfv3gfs
cd tests

# Test standard compile
./compile.sh $PWD/../FV3 hera.pgi 2>&1 | tee log.compile

# Create new baseline
./rt.sh -f -l rt_gnu_pgi.conf -c fv3 2>&1 | tee rt_create_sub.log

# Regression tests against new baseline
./rt.sh -f -l rt_gnu_pgi.conf -m 2>&1 | tee rt_sub.log
