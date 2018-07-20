# GMTB CCPP

[GMTB](http://www.dtcenter.org/GMTB/html/) Common Community Physics Package
(CCPP), including the Interoperable Physics Driver (IPD).

| Branch  | Linux/MacOS Build | Coverage |
|---      |---                |---       |
| Master  | [![Build Status](https://travis-ci.org/NCAR/gmtb-ccpp.svg?branch=master)](https://travis-ci.org/NCAR/gmtb-ccpp)  | [![Coverage Status](https://codecov.io/github/NCAR/gmtb-ccpp/coverage.svg?branch=master)](https://codecov.io/github/NCAR/gmtb-ccpp) |
| Develop | [![Build Status](https://travis-ci.org/NCAR/gmtb-ccpp.svg?branch=develop)](https://travis-ci.org/NCAR/gmtb-ccpp) | [![Coverage Status](https://codecov.io/github/NCAR/gmtb-ccpp/coverage.svg?branch=develop)](https://codecov.io/github/NCAR/gmtb-ccpp?branch=develop) |


## Notes to Users
This repository contains the Common Community Physics Packages (CCPP) and the driver 
for the CCPP (which has been until recently referred to as the GMTB IPD). To avoid
ambiguity, the term "IPD" will not be used for code contained within this repository.

The repository for the CCPP and the CCPP driver contains sufficient code for standalone
testing of the CCPP. The CCPP repository may also be used in conjunction with the 
GMTB Single Column Model (SCM). Please see the [GMTB SCM+CCPP page](http://www.dtcenter.org/GMTB/gmtb_scm_ccpp_doc/)
for more information on combining the GMTB SCM and the CCPP.

This is the release v0.1.0 of the CCPP. As this is the initial release,
the CCPP only has infrastructure to support the neccesary functioning of
the anticipated package, without having actual (i.e. physically valid)
physical parameterization schemes included. The included physical
parameterization schemes inside of the CCPP are "stub" only. While the
schemes do have arguments similar to what traditional schemes require
(wind, surface temperature, physical constants), the schemes immediately
return after a message "I am in this scheme" has been output.

This repository for the CCPP and the CCPP driver contains tests to verify
proper running of the CCPP and driver. Detailed information on how to
include fully functioning physical parameterizations schemes will be
provided once examples of fully functioning schemes are part of the CCPP.

## Requirements

### Compilers
The CCPP uses both the C and Fortran compilers. Note, the
Fortran compiler must be 2008 compliant. There are a number of Fortran
2003 pieces, and a single convenience right now with Fortran 2008.

1. [GNU Compiler Collection](https://gcc.gnu.org/)
2. [Intel 16.0.2](https://software.intel.com/en-us/intel-compilers) and beyond work.
3. [PGI](http://www.pgroup.com/) compilers do **not** easily support C functions
   calling Fortran routines. The PGI compilers attach the Fortran module name as a
   prefix to the Fortran symbol. This **breaks** the method that the CCPP uses to
   identify which schemes to call.

### [Cmake](https://cmake.org)

The CCPP build system uses cmake.

### [LibXML2](http://xmlsoft.org/)

The suite definition is currently written in XML, LibXML2 is currently used to
parse these files.

## Building
It is recommend to do an out of source build. This is "cmake" terminology
for creating a separate directory where all of the built code (objects,
libraries, executables) exist.

1. Clone the repository.
```
git clone https://github.com/NCAR/ccpp-framework ccpp
```
2. Change into the repository clone
```
cd ccpp
```
3. Specify the compiler to use. For example the GNU compilers,
   when it is available as a module called `gcc`.
  * For sh or bash
```
ml gcc
export CC=gcc
export FC=gfortran
export CXX=g++
```
  * For csh or tcsh
```
ml gcc
setenv CC gcc
setenv FC gfortran
setenv CXX g++
```
4. Make a build directory and change into it.
```
mkdir build
cd build
```
5. Create the makefiles.
```
cmake ..
```
6. Build the CCPP library and test programs.
```
make
```

## Running Tests
There are a few test programs within the `ccpp/src/tests` directory.
These should be built when the CCPP library is compiled.

To run the tests you have to add the CCPP check scheme library (`libcheck.so`)
to your `LD_LIBRARY_PATH` (`DYLD_LIBRARY_PATH` for OS X).

For sh or bash:
```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/schemes/check/src/check-build/
```

For csh or tcsh:
```
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${cwd}/schemes/check/src/check-build/
```

Note that if CCPP was built as part of a build system, you might have to load
the compiler and set environment variables that were used by the build system.


Then issue the following within the build directory.
  * `make test`

All tests should pass, if not, please open an issue. The output should be
similar to:
~~~~{.sh}
Running tests...
Test project /home/tbrown/Sources/ccpp-framework/build
    Start 1: XML_1
1/8 Test #1: XML_1 ............................   Passed    0.02 sec
    Start 2: XML_2
2/8 Test #2: XML_2 ............................   Passed    0.01 sec
    Start 3: XML_3
3/8 Test #3: XML_3 ............................   Passed    0.01 sec
    Start 4: XML_4
4/8 Test #4: XML_4 ............................   Passed    0.01 sec
    Start 5: XML_5
5/8 Test #5: XML_5 ............................   Passed    0.00 sec
    Start 6: XML_6
6/8 Test #6: XML_6 ............................   Passed    0.00 sec
    Start 7: FIELDS
7/8 Test #7: FIELDS ...........................   Passed    0.00 sec
    Start 8: CHECK
8/8 Test #8: CHECK ............................   Passed    0.01 sec

100% tests passed, 0 tests failed out of 8


Total Test time (real) =   0.08 sec
~~~~

## Validating XML
A suite is defined in XML. There are several test suites defined within
the `ccpp/src/tests` directory (which are able to test the build and
installation of the standalone CCPP). In the `ccpp/examples` directory
there are the XML files that call physical parameterization schemes.
There is also the XML Schema Definition in
that directory too. To validate a new test suite, you can use
`xmllint`. For example to validate `suite_RAP.xml`:
~~~~{.sh}
cd ccpp/examples
xmllint --schema suite.xsd --noout suite_RAP.xml
suite_RAP.xml validates
~~~~

Within the `ccpp/src/tests` directory there is a Fortran file
`test_init_finalize.f90` which will get built into an executable program
when the CCPP library is built. This program only calls:
  * `ccpp_init()`
  * `ccpp_finalize()`

It is a program to check the suite XML validation within the CCPP
library. The following is an example of using it from within the
`build` directory.
~~~~{.sh}
src/tests/test_init_finalize my_suite.xml
~~~~

For this to work, the library that is referenced in the xml file
must be added to the LD_LIBRARY_PATH (see above). To test the
correct functionality of CCPP itself, the suite suite_EXAMPLE.xml
in gmtb-ccp/sr/examples can be used.

There are two general types of XML files for the CCPP. The first is the
definition file for a suite. This has been mapped out, is fairly short,
and examples exist. Below is `examples/suite_RAP.xml`

~~~~{.xml}
<?xml version="1.0" encoding="UTF-8"?>

<suite name="RAP">
  <ipd part="1">
    <subcycle loop="1">
      <scheme>RRTMGLW</scheme>
      <scheme>RRTMGSW</scheme>
      <scheme>MYNNSFC</scheme>
      <scheme>RUCLSM</scheme>
      <scheme>MYNNPBL</scheme>
      <scheme>GF</scheme>
    </subcycle>
  </ipd>
  <ipd part="2">
    <subcycle loop="1">
      <scheme>THOMPSONAERO</scheme>
    </subcycle>
  </ipd>
</suite>
~~~~

*  suite
  * This text string "name" attribute is compared to the user-selected
  physics suite option at run-time.
*  ipd part
  * To allow for the design of the interface between the dynamics and
physical parameterization schemes, this attribute clearly associates particular
packages with the dynamical sections. In this XML example, there are two "part"
sections, with the second part only containing the "THOMPSONAERO" microphysics
scheme.
  * Users should carefully construct the XML file to map the schemes into the
existing sections of the code that calls the physical parameterization schemes.
*  subcycle
  * This functionality is not fully enabled. It is expected to be utilized for
early testing, and is included in the initial release.
*  scheme
  * The scheme elements fully describe the calling sequence of the physical
    parameterization schemes within the model.
  * For each scheme, an XML file (the scheme definition file) needs to exist.
    For the initial release, this XML file has not yet been designed.

## Physics Schemes
All physics schemes are kept in the repository under the `schemes`
directory.

To add a new scheme one needs to

1. Add/Create the scheme within `schemes`. You should create a
   sub-directory under the `schemes` directory. You will need to
   add a [`ExternalProject_Add()`](https://cmake.org/cmake/help/latest/module/ExternalProject.html).
   call to the `schemes/CMakeLists.txt` file.
2. Create a `cap` subroutine. The CCPP will call your
   cap routine.

  1. The cap routine must be labelled "schemename_cap".

     For example, the dummy scheme has a cap called
     "dummy_cap". The requirements are that it is
    1. The scheme name is lowercase (the symbol is called from a C 
       function).
    2. "_cap" is appended.
    
  2. Map all the inputs for the cap from the `cdata` encapsulating
     type (this is of the `ccpp_t` type). The cap will extract the
     fields from the fields array with the `ccpp_field_get()`
     subroutine. 

An example of a scheme is `schemes/check/test.f90`. It has the cap
routine and the run routine. The run routine prints out that the
scheme has been entered.


## Usage
The CCPP must first be initialized, this is done by calling `ccpp_init()`.
Once initialized, all variables that will be required in a physics scheme
have to be added to the ccpp data object (of type `ccpp_t`). These variables
can later be retrieved in a physics schemes cap.

Example usage, in an atmosphere component:
~~~~{.f90}
type(ccpp_t), target :: cdata
character(len=128)   :: scheme_xml_filename
integer              :: ierr

ierr = 0

! Initialize the CCPP and load the physics scheme.
call ccpp_init(scheme_xml_filename, cdata, ierr)
if (ierr /= 0) then
    call exit(1)
end if

! Add surface temperature (variable surf_t).
call ccpp_field_add(cdata, 'surface_temperature', surf_t, ierr, 'K')
if (ierr /= 0) then
    call exit(1)
end if

! Call the first physics scheme
call ccpp_ipd_run(cdata%suite%ipds(1)%subcycles(1)%schemes(1), cdata, ierr)
if (ierr /= 0) then
    call exit(1)
end if
~~~~

Example usage, in a physics cap:
~~~~{.f90}
type(ccpp_t), pointer      :: cdata
real, pointer              :: surf_t(:)
integer                    :: ierr

call c_f_pointer(ptr, cdata)
call ccpp_field_get(cdata, 'surface_temperature', surf_t, ierr)
if (ierr /= 0) then
    call exit(1)
end if
~~~~

Note, the cap routine must
* Accept only one argument of type `type(c_ptr)`.
* Be marked as `bind(c)`.

## Documentation
The code is documented with [doxygen](www.doxygen.org/).
To generate the documentation you must have [doxygen](www.doxygen.org/)
and [graphviz](http://www.graphviz.org/) installed. Then execute:
```
make doc
```

## Code Coverage
The code can be built and run to indicate code coverage. In order to do
this, you must have GNU [gcov](https://gcc.gnu.org/onlinedocs/gcc/Gcov.html)
and [lcov](http://ltp.sourceforge.net/coverage/lcov.php) installed.
To generate the coverage:

1. Make sure you are using the GNU compilers.
2. Configure the build for coverage.
  * `cmake -DCMAKE_BUILD_TYPE=Coverage ..`
3. Build the CCPP.
  * `make`
4. Build the coverage report
  * `make coverage`
The coverage report will be in the `coverage` directory within the build.
