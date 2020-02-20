# 4.7.0
defs_state MIGRATE state>:complete flag:message state_change:455 modify_change:3 server_state:HALTED
edit ECF_MICRO '%' # server
edit ECF_HOME '/scratch2/NCEPDEV/stmp3/Minsuk.Ji/ufs-weather-model/tests/ecflow_run' # server
edit ECF_JOB_CMD '%ECF_JOB% 1> %ECF_JOBOUT% 2>&1' # server
edit ECF_KILL_CMD 'kill -15 %ECF_RID%' # server
edit ECF_STATUS_CMD 'ps --sid %ECF_RID% -f' # server
edit ECF_URL_CMD '${BROWSER:=firefox} -remote 'openURL(%ECF_URL_BASE%/%ECF_URL%)'' # server
edit ECF_URL_BASE 'https://software.ecmwf.int' # server
edit ECF_URL 'wiki/display/ECFLOW/Home' # server
edit ECF_LOG '/scratch2/NCEPDEV/stmp3/Minsuk.Ji/ufs-weather-model/tests/ecflow_run/hfe11.22035.ecf.log' # server
edit ECF_INTERVAL '60' # server
edit ECF_LISTS '/scratch2/NCEPDEV/stmp3/Minsuk.Ji/ufs-weather-model/tests/ecflow_run/ecf.lists' # server
edit ECF_CHECK '/scratch2/NCEPDEV/stmp3/Minsuk.Ji/ufs-weather-model/tests/ecflow_run/hfe11.22035.check' # server
edit ECF_CHECKOLD '/scratch2/NCEPDEV/stmp3/Minsuk.Ji/ufs-weather-model/tests/ecflow_run/hfe11.22035.check.b' # server
edit ECF_CHECKINTERVAL '120' # server
edit ECF_CHECKMODE 'CHECK_ON_TIME' # server
edit ECF_TRIES '2' # server
edit ECF_VERSION '4.7.0' # server
edit ECF_PORT '22035' # server
edit ECF_NODE '%ECF_HOST%' # server
edit ECF_HOST 'hfe11' # server
edit ECF_PID '14275' # server
history / MSG:[23:29:44 19.2.2020] --restart :Minsuk.JiMSG:[23:29:44 19.2.2020] --load=/scratch2/NCEPDEV/stmp3/Minsuk.Ji/ufs-weather-model/tests/ecflow_run/regtest.def  :Minsuk.JiMSG:[23:29:44 19.2.2020] --begin=regtest :Minsuk.JiMSG:[00:41:15 20.2.2020] --halt=yes :Minsuk.Ji
suite regtest #  begun:1 state:complete dur:01:10:16
edit ECF_HOME '/scratch2/NCEPDEV/stmp3/Minsuk.Ji/ufs-weather-model/tests/ecflow_run'
edit ECF_INCLUDE '/scratch2/NCEPDEV/stmp3/Minsuk.Ji/ufs-weather-model/tests/ecflow_run'
edit ECF_KILL_CMD 'kill -15 %ECF_RID% > %ECF_JOB%.kill 2>&1'
edit ECF_TRIES '1'
limit max_builds 6
limit max_jobs 30
label rundir_root "/scratch1/NCEPDEV/stmp2/Minsuk.Ji/FV3_RT/rt_12495"
calendar initTime:2020-Feb-19 23:29:44 suiteTime:2020-Feb-20 00:41:00 duration:01:11:16 initLocalTime:2020-Feb-19 23:29:44 lastTime:2020-Feb-20 00:41:00 calendarIncrement:00:01:00
task compile_1 # try:1 state:complete dur:00:06:16
inlimit max_builds
task fv3_control # try:1 state:complete dur:00:11:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_decomp # try:1 state:complete dur:00:10:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_2threads # try:1 state:complete dur:00:15:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_restart # try:1 state:complete dur:00:13:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_read_inc # try:1 state:complete dur:00:15:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_gfdlmp # try:1 state:complete dur:00:12:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_gfdlmprad_gwd # try:1 state:complete dur:00:13:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_gfdlmprad_noahmp # try:1 state:complete dur:00:12:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_thompson # try:1 state:complete dur:00:16:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_wsm6 # try:1 state:complete dur:00:14:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_wrtGauss_netcdf_esmf # try:1 state:complete dur:00:13:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_wrtGauss_netcdf # try:1 state:complete dur:00:11:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_wrtGauss_nemsio # try:1 state:complete dur:00:10:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_wrtGauss_nemsio_c192 # try:1 state:complete dur:00:10:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_stochy # try:1 state:complete dur:00:12:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_iau # try:1 state:complete dur:00:15:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_csawmgshoc # try:1 state:complete dur:00:14:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_csawmg # try:1 state:complete dur:00:13:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_rasmgshoc # try:1 state:complete dur:00:14:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_csawmg3shoc127 # try:1 state:complete dur:00:14:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_satmedmf # try:1 state:complete dur:00:11:16
trigger compile_1 == complete
inlimit max_jobs
task fv3_lheatstrg # try:1 state:complete dur:00:11:16
trigger compile_1 == complete
inlimit max_jobs
task compile_2 # try:1 state:complete dur:00:06:16
inlimit max_builds
task fv3_gfdlmprad # try:1 state:complete dur:00:10:16
trigger compile_2 == complete
inlimit max_jobs
task fv3_wrtGauss_nemsio_c768 # try:1 state:complete dur:01:10:16
trigger compile_2 == complete
inlimit max_jobs
task fv3_appbuild # try:1 state:complete dur:00:11:16
trigger compile_2 == complete
inlimit max_jobs
task compile_3 # try:1 state:complete dur:00:06:16
inlimit max_builds
task fv3_control_32bit # try:1 state:complete dur:00:10:16
trigger compile_3 == complete
inlimit max_jobs
task fv3_gfdlmprad_32bit_post # try:1 state:complete dur:00:14:16
trigger compile_3 == complete
inlimit max_jobs
task fv3_stretched # try:1 state:complete dur:00:17:16
trigger compile_3 == complete
inlimit max_jobs
task fv3_stretched_nest # try:1 state:complete dur:00:18:16
trigger compile_3 == complete
inlimit max_jobs
task fv3_stretched_nest_quilt # try:1 state:complete dur:00:19:16
trigger compile_3 == complete
inlimit max_jobs
task fv3_regional_control # try:1 state:complete dur:00:20:16
trigger compile_3 == complete
inlimit max_jobs
task fv3_regional_restart # try:1 state:complete dur:00:26:16
trigger compile_3 == complete and fv3_regional_control == complete
inlimit max_jobs
task fv3_regional_quilt # try:1 state:complete dur:00:20:16
trigger compile_3 == complete
inlimit max_jobs
task fv3_regional_c768 # try:1 state:complete dur:00:23:16
trigger compile_3 == complete
inlimit max_jobs
task compile_4 # try:1 state:complete dur:00:02:16
inlimit max_builds
task fv3_control_debug # try:1 state:complete dur:00:05:16
trigger compile_4 == complete
inlimit max_jobs
task fv3_stretched_nest_debug # try:1 state:complete dur:00:08:16
trigger compile_4 == complete
inlimit max_jobs
endsuite
