################################################################################
# Cheyenne/Intel
################################################################################

# OFFICIAL BASELINE
# /glade/p/ral/jntp/GMTB/NEMSfv3gfs/RT/trunk-20190912/INTEL

# OWN BASELINE
# /glade/scratch/$USER/FV3_RT/REGRESSION_TEST_INTEL

# RUN DIRS
# /glade/scratch/$USER/FV3_RT/rt_...

export NEMS_COMPILER=intel # optional, default is intel if not set
export ACCNR=...           # set account name/number for job scheduler

# from the top-level directory of NEMSfv3gfs
cd tests

# Test standard compile
./compile.sh $PWD/../FV3 cheyenne.intel 2>&1 | tee log.compile

# Regression tests against official baseline
./rt.sh -f 2>&1 | tee rt.log

# Create new baseline
./rt.sh -f -c 2>&1 | tee rt_create.log

# Regression tests against new baseline
./rt.sh -f -m 2>&1 | tee rt_verify.log

################################################################################
# Cheyenne/GNU
################################################################################

# OFFICIAL BASELINE
# /glade/p/ral/jntp/GMTB/NEMSfv3gfs/RT/trunk-20190912/GNU

# OWN BASELINE
# /glade/scratch/$USER/FV3_RT/REGRESSION_TEST_GNU

# RUN DIRS
# /glade/scratch/$USER/FV3_RT/rt_...

export NEMS_COMPILER=gnu
export ACCNR=...           # set account name/number for job scheduler

# from the top-level directory of NEMSfv3gfs
cd tests

# Test standard compile
./compile.sh $PWD/../FV3 cheyenne.gnu 2>&1 | tee log.compile

# Regression tests against official baseline
./rt.sh -l rt_gnu.conf 2>&1 | tee rt_gnu.log

# Create new baseline
./rt.sh -l rt_gnu.conf -c 2>&1 | tee rt_gnu_create.log

# Regression tests against new baseline
./rt.sh -l rt_gnu.conf -m 2>&1 | tee rt_gnu_verify.log

################################################################################
# Cheyenne/PGI
################################################################################

no longer supported
