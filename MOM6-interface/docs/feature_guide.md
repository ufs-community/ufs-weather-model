# MOM6 Output Logging User Guide

The MOM6 output logging feature is designed to track the completion of
MOM6 history and restart files during a model run. This feature is
specific to UFS Weather Model (UWM) operational requirements and
configurations (eg specific output frequencies in hours) and may break
if used outside the scope of intended use.

The feature is enabled by adding a namelist to the model `input.nml`
which can be used to define the history output frequencies, file names
and types (snapshot or averages) which are to be tracked. The user is
responsible for ensuring that these settings match the contents of the
`diag_table` in use by the model. From the user-configuration, one or
more tracking Alarms are enabled at the desired frequencies. When an
Alarm rings, a filename is constructed referencing the model time and
the initial state of the file is recorded. Depending on the
characteristics of the file at creation, the criteria to declare a
file complete is defined.

The file state will be checked at each succeeding ModelAdvance until
the appropriate completion criteria is met. The criteria used are
either:

1. the unlimited dimension in the file is greater than zero,
2. when the unlimited dimension is greater than zero and the filesize
is larger than the initial size.

When a file is determined to be complete, a log file is recorded
containing the forecast hour, the valid time, the name of the output
file and the last completed restart file. The log file can then be
used by any related workflow
(e.g. [global-workflow](https://github.com/NOAA-EMC/global-workflow)).

To illustrate the concepts implemented in the output logging feature,
consider the following diagram:

@image html logging_diagram1.png "File Sequencing" width=50%

In this case, MOM6 is providing 6-hourly output; the model begins at
hour=00 and runs for 30 hours; the output alarm will ring every 6
hours. It is important to remember that alarms ring at the
ModelAdvance nextTime. If the model in this case is running with a 30
minute coupling frequency, the alarm at hour=12:00 will ring when the
modelAdvance currTime = 11:30 and the ModelAdvance nextTime = 12:00.

MOM6 history files for time-averaged output are timestamped at the
middle of the averaging window. The first MOM6 history output file
produced is the hour = 03 file, representing the average between hours
0 and 6. It will be created when the alarm rings at nextTime = 12:00
and be completed on the next timestep (when the currTime =
12:00). That file sequencing pattern holds until the model stop
time. At finalization, there are two final history files that are
written: the "pending" file (h = 21) as well as the final h = 03 file.

Restarts for MOM6 are written by MOM6 directly, not using FMS as for
the history files. Model restarts are written, complete, when the
ModelAdvance nextTime matches a requested hour. In the above diagram,
restarts written at hour = 12 will be when the ModelAdvance currTime =
11:30 and the ModelAdvance nextTime = 12:00.

The output log file written at forecast hour = 18 for this case
(`20111001.180000.mom6.06h`) would contain the following information:

```text
completed: mom6.06h
forecast hour:    18.000
valid time:     2011      10       1      18       0       0
last output: ./MOM6_OUTPUT/ocn_2011_10_01_09_00.nc
last restart:     2011      10       1      12       0       0
```

In the above case, the log file for both the h = 21 and h = 03 file
are both written at the same model hour (FH = 30 in this
case). Therefore, the log for the very last file is denoted with
`lstop` appended to the file name (`20111002.060000.mom6.lstop.06h`).

@note The above sample log file also highlights the essential feature
of the output sequencing. At FH = 18, the only output files which are
complete are those at hours 03 and 09; the history file for the
interval which *ends* at FH = 18 (the h = 15 file) has not yet been
written.

# Design and Function

The outputlog feature consists of two modules in
`config_src/drivers/nuopc_cap`:

* mom_cap_outputlog.F90
* mom_outputlog_methods.F90

and three type structures :

* mom_outputlog_methods::outputlog_config_type : the time-invariant
    configuration at each tracked frequency
* mom_outputlog_methods::outputlog_state_type : the time-evolving
    state at each tracked frequency
* mom_outputlog_methods::outputlog_modeltime_type: the model time
    state during each ModelAdvance

## Configuration

The output log feature is enabled with a `MOM_outputlog_nml` namelist
added to the `input.nml` which will populate the
outputlog_config_type. For example, the following values will request
tracking 6-hourly files, which (using the `diag_table`) have a
filename prefix of `ocn` and are defined as time-averaged values. No
debug information will be added to the standard output file.

```text
&MOM_outputlog_nml
  outputlog_fh = 6
  outputlog_fnameprefix = 'ocn'
  outputlog_treduce = 'average'
  outputlog_debug = .false.
/
```

At a minimum, the desired tracking frequency must be provided, which
will default to time-averaged files with the file prefix `ocn`. The
following namelist will be treated identically as that above

```text
&MOM_outputlog_nml
  outputlog_fh = 6
/
```

The following rules apply the namelist options for output logging:

### Logging Frequency

File tracking can be enabled for 1,3,6 or 24 hourly files
only. Multiple frequencies can be requested and the listed order is
immaterial. However, each logging frequency must be uniquely
defined. For example, 6-hourly average files and 6-hourly snapshot
files are not allowed but 6-hourly average and 3-hourly snapshot files
are.

@note Because it is intended for operational purposes, logging
frequency is intentionally set to log 3 and 6 hour frequencies
relative to the operational forecast windows. This means that a 6 hour
tracking frequency is always implemented as starting at one of hours
00,06,12 or 18. A similar rule applies to a 3 hour tracking frequency
(i.e., 00,03,06). Tracking at 24 hour intervals is for consecutive
24 hour periods. For example, hour 09 on day one through to hour 09 on
day 2.

### Filename Prefix

If a single frequency is requested, no filename prefix is required. A
default prefix of `ocn` will be used. It is the user's responsibility
that the default matches the specification in the
`diag_table`. Otherwise, the filename prefix must be set to the actual
filename prefix called for in the `diag_table`. If more than a single
frequency is requested, the user must provide filename prefixes for
all frequencies (again ensuring matches to the `diag_table`); these
prefixes must also be distinct from one another, or the namelist will
be rejected. The filename prefix can have a maximum length of 12 (not
including the trailing underscore, which will be appended).

### File Time Reduction

Either instantaneous (snapshot) files or time averaged files are
supported. These are specified by the namelist `treduce` settings of
`none` and `average`, respectfully. The file name format for time
averaged output is assumed to be timestamped with the midpoint of
the averaging window (i.e. 09 for the 6h-12h average). For snapshot
output, the timestamp will be the time of the snapshot.  The user is
responsible for ensuring that the diag_table in use matches these
definitions.

### Feature Debugging

When enabled in the namelist, debugging print statements will be
written to standard out. These statements track the state of the
feature tracking at each step through the ModelAdvance. Print
statements will be prepended with the routine name, for example
`MOM_cap:(track_freqn)`. Utilizing this feature for a case with
tracking of 6-hourly average output produces:

```text
MOM_cap:(track_freqn) ./MOM6_OUTPUT/ocn_2021_03_22_21_00.nc exists 2021-03-23T05:30:00  2021-03-23T06:00:00 not complete, chkflag  T     9415276      9415276    1
MOM_cap:(track_freqn) ./MOM6_OUTPUT/ocn_2021_03_22_21_00.nc exists 2021-03-23T06:00:00  2021-03-23T06:30:00     complete, chkflag  F     9415276     90532460    1
MOM_cap:(track_freqn) ./MOM6_OUTPUT/ocn_2021_03_22_21_00.nc exists 2021-03-23T06:30:00  2021-03-23T07:00:00     complete, chkflag  F     9415276     90532460    1
MOM_cap:(track_freqn) ./MOM6_OUTPUT/ocn_2021_03_22_21_00.nc exists 2021-03-23T07:00:00  2021-03-23T07:30:00     complete, chkflag  F     9415276     90532460    1
```

As the model advances, the feature scans for a specific filename. When
that file is found, the debug print will indicate that the file is
present at a given `modelAdvance` time pair (the currTime and the
nextTime). The completion state is given, as well as the status of the
`chkfile_nextAdvance` logical. The next two columns report the size of the
file when created and the current size of the file. The final column
reports the length of the unlimited dimension. Once the file is
determined complete, the checking flag flips to false and no further
inquires on the state of that particular file will be made.

### Alarm Initialization

An alarm is initialized at each desired tracking frequency. As noted
previously, alarms are set to ring at multiples of the tracking
frequency and initialized with a time-offset to ensure that they ring
on intervals associated with the operational forecast hours.

### IO-layout

When IO-layout is enabled, the root PE associated with the IO domain
is co-located with the root PE of the computation domain. Each
IO-domain will produce a history file for the domain; the file names
will be appended with the IO-domain number, for example
`.nc.0000`. The number of history files expected is obtained using the
internal MOM6 function `mpp_get_io_domain_layout`. For a single
IO-domain, no file name suffix is used; when IO-layout is in use, only
the root IO task associated with the `.nc.0000` file will be queried
for the model state.

## File Tracking Sequence

### File State at Creation

When an alarm rings for a frequency, the filename for the expected
history output is constructed and the initial state of the file, if
present, is obtained. Two characteristics of the file are obtained:
the length of the unlimited dimension and the initial size of the
file. For configurations using MOM6 and a data atmosphere (DATM), the
inital unlimited dimension when the file is created is found to be
0. For active atmosphere configurations, an initial unlimited
dimension length of 1 is observed. In these cases, file size will also
be used to determine file completion. This is enabled using the state
variable `use_filesize`.

Once the initial file state is obtained at ring time, a flag will be
set to check the file on the next ModelAdvance. The file will continue
to be checked on each advance until it is reported complete. In
practice, the file is completed on the next ModelAdvance after
creation.

### Determining File Completion

Depending on the state variable `use_filesize`, a file is determined
to be complete in one of two ways. If the initial unlimited dimension
was zero (a DATM case), file completion occurs when the file state
obtains an unlimited dimension of 1. Otherwise, the secondary criteria
of file size will be used. The criteria here is only that the file
size is greater than the initial size.

Once a file is determined complete, the check flag for the next
Advance is turned off and the state is updated to store the last
restart file written.

### Finalization

When the model completes in the ModelFinalization phase, two calls to
the tracking feature are required. Both occur after the IO has been
shut down. The first call checks the status of the penultimate file;
the file which would have completed at the next ModelAdvance
timestep. Since in this case, there will be no next timestep, this
'pending' file is completed by the IO shutdown. The final file, for
the last averaging window, is also completed during the IO shutdown.

### Restart Pairing

The history files for MOM6 are at specific cadences,
e.g. 6-hourly. The restart files, however, can be written at either
specific (but different) cadences or at specific forecast hours. There
is in general no simple *a priori* relationship between the writing of
restarts and the writing of history files. Therefore, the concept of
'restart pairing' is used to refer to tracking the last restart file
which was written at the time a history file is written. The paired
restart is written to the log file produced by the feature.

The following example is based on the gfsv17 IAU regression
test. History output for MOM6 is set to 6-hourly averages output and
restarts are written at specified hours. In this case, the model is
restarting with `FHROT = 3` and since IAU is active, MOM6 averaging
begins at `03-22-12`. The restarts are requested at `restart_interval:
6 24 45 78`.

The logfile `20210324.060000.mom6.06h` will contain:

```text
completed: mom6.06h
forecast hour:    48.000
valid time:     2021       3      24       6       0       0
last output: ./MOM6_OUTPUT/ocn_2021_03_23_21_00.nc
last restart:     2021       3      24       3       0       0
```

@image html logging_diagram2.png "Restart Pairing Diagram" width=80%


# Unit Testing

The output logging feature is covered by a suite of unit tests. The tests
consist of three supporting modules :

| Module Name | Purpose |
| :--- | :--- |
| nc_fixture_mod.F90 | creates mock netCDF files mimicking possible history and restart states |
| test_helpers.F90 | convenience module with re-used functions |
| test_utils.F90 | assertion and error utilities |

and five tests:

| Program Name | Purpose | Target Functions |
| :--- | :--- | :--- |
| test_outputlog_completion.F90 | tests that nc_fixture_mod creates files of the correct state | get_file_state, file_is_complete |
| test_outputlog_readnml.F90 | tests that readnml correctly identifies invalid namelist settings | readnml |
| test_outputlog_alarminit.F90 | tests alarm initialization and ring times against possible start times and frequencies | alarminit |
| test_outputlog_freqn.F90 | tests the orchestration between file creation and file completion | track_freqn |
| test_outputlog_restn.F90 | tests the completion check of single and multi-part restart files | track_restn |
