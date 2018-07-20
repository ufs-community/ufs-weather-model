#!/usr/bin/env python
#
# Script to generate a cap module and subroutines
# from a scheme xml file.
#

from __future__ import print_function
import os
import sys
import getopt
import xml.etree.ElementTree as ET

###############################################################################

STANDARD_VARIABLE_TYPES = [ 'character', 'integer', 'logical', 'real' ]

CCPP_ERROR_FLAG_VARIABLE = 'error_flag'

###############################################################################

class Var(object):

    def __init__(self, **kwargs):
        self._standard_name = None
        self._long_name     = None
        self._units         = None
        self._local_name    = None
        self._type          = None
        self._rank          = None
        self._container     = None
        self._kind          = None
        self._intent        = None
        self._optional      = None
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)

    @property
    def standard_name(self):
        '''Get the name of the variable.'''
        return self._standard_name

    @standard_name.setter
    def standard_name(self, value):
        self._standard_name = value

    @property
    def long_name(self):
        '''Get the name of the variable.'''
        return self._long_name

    @long_name.setter
    def long_name(self, value):
        self._long_name = value

    @property
    def units(self):
        '''Get the units of the variable.'''
        return self._units

    @units.setter
    def units(self, value):
        self._units = value

    @property
    def local_name(self):
        '''Get the local variable name of the variable.'''
        return self._local_name

    @local_name.setter
    def local_name(self, value):
        self._local_name = value

    @property
    def type(self):
        '''Get the type of the variable.'''
        return self._type

    @type.setter
    def type(self, value):
        if value == 'character' and self._rank:
            raise Exception('Arrays of Fortran strings not implemented in CCPP')
        self._type = value

    @property
    def rank(self):
        '''Get the rank of the variable.'''
        return self._rank

    @rank.setter
    def rank(self, value):
        if not isinstance(value, int):
            raise TypeError('Invalid type for variable property rank, must be integer')
        elif self.type == 'character' and value > 0:
            raise Exception('Arrays of Fortran strings not implemented in CCPP')
        if (value == 0):
            self._rank = ''
        else:
            self._rank = '('+ ','.join([':'] * value) +')'

    @property
    def kind(self):
        '''Get the kind of the variable.'''
        return self._kind

    @kind.setter
    def kind(self, value):
        self._kind = value

    @property
    def intent(self):
        '''Get the intent of the variable.'''
        return self._intent

    @intent.setter
    def intent(self, value):
        if not value in ['none', 'in', 'out', 'inout']:
            raise ValueError('Invalid value {0} for variable property intent'.format(value))
        self._intent = value

    @property
    def optional(self):
        '''Get the optional of the variable.'''
        return self._optional

    @optional.setter
    def optional(self, value):
        if not value in ['T', 'F']:
            raise ValueError('Invalid value {0} for variable property optional'.format(value))
        self._optional = value

    @property
    def container(self):
        '''Get the container of the variable.'''
        return self._container

    @container.setter
    def container(self, value):
        self._container = value

    def compatible(self, other):
        # We accept character(len=*) as compatible with character(len=INTEGER_VALUE)
        if self.type == 'character':
            if (self.kind == 'len=*' and other.kind.startswith('len=')) or \
                   (self.kind.startswith('len=') and other.kind == 'len=*'):
                return self.standard_name == other.standard_name \
                    and self.units == other.units \
                    and self.type == other.type \
                    and self.rank == other.rank
        return self.standard_name == other.standard_name \
            and self.units == other.units \
            and self.type == other.type \
            and self.kind == other.kind \
            and self.rank == other.rank

    def print_def(self):
        '''Print the definition line for the variable.'''
        if self.type in STANDARD_VARIABLE_TYPES:
            if self.kind:
                str = "{s.type}({s._kind}), pointer :: {s.local_name}{s.rank}"
            else:
                str = "{s.type}, pointer :: {s.local_name}{s.rank}"
        else:
            if self.kind:
                error_message = "Generating variable definition statements for derived types with" + \
                                " kind attributes not implemented; variable: {0}".format(self.standard_name)
                raise Exception(error_message)
            else:
                str = "type({s.type}), pointer     :: {s.local_name}{s.rank}"
        return str.format(s=self)

    def print_get(self):
        '''Print the data retrieval line for the variable. Depends on the type and of variable'''
        if self.type in STANDARD_VARIABLE_TYPES:
            str='''
        call ccpp_field_get(cdata, '{s.standard_name}', {s.local_name}, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve {s.standard_name} from CCPP data structure')
            return
        end if'''
        # Derived-type variables, scalar
        elif self.rank == '':
            str='''
        call ccpp_field_get(cdata, '{s.standard_name}', cptr, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve {s.standard_name} from CCPP data structure')
            return
        end if
        call c_f_pointer(cptr, {s.local_name})'''
        # Derived-type variables, array
        else:
            str='''
        call ccpp_field_get(cdata, '{s.standard_name}', cptr, ierr, dims=cdims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve {s.standard_name} from CCPP data structure')
            return
        end if
        call c_f_pointer(cptr, {s.local_name}, cdims)
        deallocate(cdims)'''
        return str.format(s=self)

    def print_add(self, ccpp_data_structure):
        '''Print the data addition line for the variable. Depends on the type of variable.
        Since the name of the ccpp data structure is not known, this needs to be filled later.
        In case of errors a message is printed to screen; using 'return' statements as above
        for ccpp_field_get is not possible, since the ccpp_field_add statements may be placed
        inside OpenMP parallel regions.'''
        # Standard-type variables, scalar and array
        if self.type in STANDARD_VARIABLE_TYPES:
            str='''
            call ccpp_field_add({ccpp_data_structure}, '{s.standard_name}', {s.target}, ierr, '{s.units}')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "{s.standard_name}" to CCPP data structure')
            end if'''
        # Derived-type variables, scalar
        elif self.rank == '':
            str='''
            call ccpp_field_add({ccpp_data_structure}, '{s.standard_name}', '', c_loc({s.target}), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "{s.standard_name}" to CCPP data structure')
            end if'''
        # Derived-type variables, array
        else:
            str='''
            call ccpp_field_add({ccpp_data_structure}, '{s.standard_name}', '', c_loc({s.target}), rank=size(shape({s.target})), dims=shape({s.target}), ierr=ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "{s.standard_name}" to CCPP data structure')
            end if'''
        return str.format(ccpp_data_structure=ccpp_data_structure, s=self)

    def print_debug(self):
        '''Print the data retrieval line for the variable.'''
        str='''Contents of {s} (* = mandatory for compatibility):
        standard_name = {s.standard_name} *
        long_name     = {s.long_name}
        units         = {s.units} *
        local_name    = {s.local_name}
        type          = {s.type} *
        rank          = {s.rank} *
        kind          = {s.kind} *
        intent        = {s.intent}
        optional      = {s.optional}
        container     = {s.container}'''
        return str.format(s=self)
    
    @classmethod
    def from_table(cls, columns, data):
        var = cls()
        var.standard_name = data[columns.index('standard_name')]
        var.long_name     = data[columns.index('long_name')]
        var.units         = data[columns.index('units')]
        var.local_name    = data[columns.index('local_name')]
        var.rank          = int(data[columns.index('rank')])
        var.type          = data[columns.index('type')]
        var.kind          = data[columns.index('kind')]
        var.intent        = data[columns.index('intent')]
        var.optional      = data[columns.index('optional')]
        return var

    def to_xml(self, element):
        element.set('name', self._standard_name)
        sub_element = ET.SubElement(element, 'standard_name')
        sub_element.text = self._standard_name
        sub_element = ET.SubElement(element, 'long_name')
        sub_element.text = self._long_name
        sub_element = ET.SubElement(element, 'units')
        sub_element.text = self._units
        sub_element = ET.SubElement(element, 'local_name')
        sub_element.text = self._local_name
        sub_element = ET.SubElement(element, 'type')
        sub_element.text = self._type
        sub_element = ET.SubElement(element, 'rank')
        sub_element.text = self._rank
        sub_element = ET.SubElement(element, 'intent')
        sub_element.text = self._intent
        sub_element = ET.SubElement(element, 'optional')
        sub_element.text = self._optional
        sub_element = ET.SubElement(element, 'container')
        sub_element.text = self._container
        return element

###############################################################################
class Cap(object):

    header='''
!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Auto-generated cap module for the {module} scheme
!!
!
module {module}_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: {module}, &
                      only: {subroutines}
    ! Other modules required, e.g. type definitions
    {module_use}

    implicit none

    private
    public :: {subroutine_caps}

    contains

'''

    sub='''
    function {subroutine}_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
{var_defs}

        ierr = 0

        call c_f_pointer(ptr, cdata)

{var_gets}

        call {subroutine}({args})
        {ierr_assign}

    end function {subroutine}_cap
'''

    def __init__(self, **kwargs):
        self._filename = 'sys.stdout'
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)

    def write(self, module, module_use, data):
        if (self.filename is not sys.stdout):
            f = open(self.filename, 'w')
        else:
            f = sys.stdout

        subs = ','.join(["{0}".format(s) for s in data.keys()])
        sub_caps = ','.join(["{0}_cap".format(s) for s in data.keys()])

        f.write(Cap.header.format(module = module,
                                  module_use = module_use,
                                  subroutines = subs,
                                  subroutine_caps = sub_caps))

        for sub in data.keys():
            var_defs = "\n".join([" "*8 + x.print_def() for x in data[sub]])
            var_gets = "\n".join([x.print_get() for x in data[sub]])
            # Split args so that lines don't exceed 260 characters (for PGI)
            args = ''
            length = 0
            for x in data[sub]:
                arg = "{0}={0},".format(x.local_name)
                args += arg
                length += len(arg)
                if length > 70 and not x == data[sub][-1]:
                    args += ' &\n                  '
                    length = 0
            args = args.rstrip(',')
            # If CCPP_ERROR_FLAG_VARIABLE is present, assign to ierr
            ierr_assign = ''
            for x in data[sub]:
                if x.standard_name == CCPP_ERROR_FLAG_VARIABLE:
                    ierr_assign = 'ierr={x.local_name}'.format(x=x)
                    break
            # Write to scheme cap
            f.write(Cap.sub.format(subroutine=sub,
                                   var_defs=var_defs,
                                   var_gets=var_gets,
                                   args=args,
                                   ierr_assign=ierr_assign))
        f.write("end module {module}_cap\n".format(module = module))

        if (f is not sys.stdout):
            f.close()

    @property
    def filename(self):
        '''Get the filename of write the output to.'''
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = value

class CapsMakefile(object):

    header='''
# All CCPP caps are defined here.
#
# This file is auto-generated using ccpp_prebuild.py
# at compile time, do not edit manually.
#
CAPS_F90 ='''

    def __init__(self, **kwargs):
        self._filename = 'sys.stdout'
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)

    def write(self, caps):
        if (self.filename is not sys.stdout):
            f = open(self.filename, 'w')
        else:
            f = sys.stdout

        contents = self.header
        for cap in caps:
            contents += ' \\\n\t   {0}'.format(cap)
        f.write(contents)

        if (f is not sys.stdout):
            f.close()

    @property
    def filename(self):
        '''Get the filename of write the output to.'''
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = value

class CapsCMakefile(object):

    header='''
# All CCPP caps are defined here.
#
# This file is auto-generated using ccpp_prebuild.py
# at compile time, do not edit manually.
#
set(CAPS
'''
    footer=''')
'''

    def __init__(self, **kwargs):
        self._filename = 'sys.stdout'
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)

    def write(self, schemes):
        if (self.filename is not sys.stdout):
            f = open(self.filename, 'w')
        else:
            f = sys.stdout

        contents = self.header
        for scheme in schemes:
            contents += '      {0}\n'.format(scheme)
        contents += self.footer
        f.write(contents)

        if (f is not sys.stdout):
            f.close()

    @property
    def filename(self):
        '''Get the filename of write the output to.'''
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = value

class SchemesMakefile(object):

    header='''
# All CCPP schemes are defined here.
#
# This file is auto-generated using ccpp_prebuild.py
# at compile time, do not edit manually.
#
SCHEMES_F =

SCHEMES_F90 =

SCHEMES_f =

SCHEMES_f90 ='''

    def __init__(self, **kwargs):
        self._filename = 'sys.stdout'
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)

    def write(self, schemes):
        if (self.filename is not sys.stdout):
            f = open(self.filename, 'w')
        else:
            f = sys.stdout

        contents = self.header
        schemes_F = 'SCHEMES_F ='
        schemes_F90 = 'SCHEMES_F90 ='
        schemes_f = 'SCHEMES_f ='
        schemes_f90 = 'SCHEMES_f90 ='
        for scheme in schemes:
            if scheme.endswith('.F'):
                schemes_F += ' \\\n\t   {0}'.format(scheme)
            elif scheme.endswith('.F90'):
                schemes_F90 += ' \\\n\t   {0}'.format(scheme)
            elif scheme.endswith('.f'):
                schemes_f += ' \\\n\t   {0}'.format(scheme)
            elif scheme.endswith('.f90'):
                schemes_f90 += ' \\\n\t   {0}'.format(scheme)
        contents = contents.replace('SCHEMES_F =', schemes_F)
        contents = contents.replace('SCHEMES_F90 =', schemes_F90)
        contents = contents.replace('SCHEMES_f =', schemes_f)
        contents = contents.replace('SCHEMES_f90 =', schemes_f90)
        f.write(contents)

        if (f is not sys.stdout):
            f.close()

    @property
    def filename(self):
        '''Get the filename of write the output to.'''
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = value

class SchemesCMakefile(object):

    header='''
# All CCPP schemes are defined here.
#
# This file is auto-generated using ccpp_prebuild.py
# at compile time, do not edit manually.
#
set(SCHEMES
'''
    footer=''')
'''

    def __init__(self, **kwargs):
        self._filename = 'sys.stdout'
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)

    def write(self, schemes):
        if (self.filename is not sys.stdout):
            f = open(self.filename, 'w')
        else:
            f = sys.stdout

        contents = self.header
        for scheme in schemes:
            contents += '      {0}\n'.format(scheme)
        contents += self.footer
        f.write(contents)

        if (f is not sys.stdout):
            f.close()

    @property
    def filename(self):
        '''Get the filename of write the output to.'''
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = value

###############################################################################
if __name__ == "__main__":
    main()
