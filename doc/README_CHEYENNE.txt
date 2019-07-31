################################################################################
# Cheyenne/Intel
################################################################################

# OFFICIAL BASELINE
# /glade/p/ral/jntp/GMTB/NEMSfv3gfs/RT/trunk-20190315/INTEL

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
./rt.sh -f 2>&1 | tee rt_full.log

################################################################################
# Cheyenne/GNU
################################################################################

# OFFICIAL BASELINE
# /glade/p/ral/jntp/GMTB/NEMSfv3gfs/RT/trunk-20190315/GNU

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

# Create new baseline
./rt.sh -f -l rt_gnu_pgi.conf -c fv3 2>&1 | tee rt_create_sub.log

# Regression tests against new baseline
./rt.sh -f -l rt_gnu_pgi.conf -m 2>&1 | tee rt_sub.log

################################################################################
# Cheyenne/PGI
################################################################################

# OFFICIAL BASELINE
# /glade/p/ral/jntp/GMTB/NEMSfv3gfs/RT/trunk-20190315/PGI

# OWN BASELINE
# /glade/scratch/$USER/FV3_RT/REGRESSION_TEST_PGI

# RUN DIRS
# /glade/scratch/$USER/FV3_RT/rt_...

export NEMS_COMPILER=pgi
export ACCNR=...           # set account name/number for job scheduler

# from the top-level directory of NEMSfv3gfs
cd tests

# Test standard compile
./compile.sh $PWD/../FV3 cheyenne.pgi 2>&1 | tee log.compile

# Create new baseline
./rt.sh -f -l rt_gnu_pgi.conf -c fv3 2>&1 | tee rt_create_sub.log

# Regression tests against new baseline
./rt.sh -f -l rt_gnu_pgi.conf -m 2>&1 | tee rt_sub.log
