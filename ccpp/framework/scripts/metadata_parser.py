#!/usr/bin/env python

import collections
import logging
import subprocess
from xml.etree import ElementTree as ET

from common import indent, encode_container
from mkcap import Var

# The argument tables for schemes and variable definitions should have the following format:
# !! \section arg_table_SubroutineName (e.g. SubroutineName = SchemeName_run) OR \section arg_table_DerivedTypeName OR \section arg_table_ModuleName
# !! | local_name     | standard_name                                         | long_name                                | units   | rank | type      |    kind   | intent | optional |
# !! |----------------|-------------------------------------------------------|------------------------------------------|---------|------|-----------|-----------|--------|----------|
# !! | im             | horizontal_loop_extent                                | horizontal loop extent, start at 1       | index   |    0 | integer   |           | in     | F        |
# !! | ix             | horizontal_dimension                                  | horizontal dimension                     | index   |    0 | integer   |           | in     | F        |
# !! | ...            | ...                                                   |                                          |         |      |           |           |        |          |
# !! | errmsg         | error_message                                         | error message for error handling in CCPP | none    |    0 | character |           | out    | F        |
# !! | ierr           | error_flag                                            | error flag for error handling in CCPP    | none    |    0 | integer   |           | out    | F        |
# !!
# Notes on the input format:
# - if the argument table starts a new doxygen section, it should start with !> \section instead of !! \section
# - the "\section arg_table_{SubroutineName,DerivedTypeName,ModuleName}" command denotes the start of the table
#      - SubroutineName must match the name of the subroutine that the argument table describes
#      - DerivedTypeName must match the name of the derived type that the argument table describes
#      - ModuleName must match the name of the module whose variables the argument table describes
# - the table must be placed immediately before the subroutine / derived data type,
#   or immediately before the module variables (but within the module structure)
# - each line of the table must begin with the doxygen-delimiter '!!'
# - table headers are the first row, the second row must have the |---|-----| format
# - after the last row of the table, there must be a blank doxygen line (only '!!') to denote the end of the table
# - for variable type definitions and module variables, the intent and optional columns must be set to 'none' and 'F'
# - each argument table (and its subroutine) must accept the following two arguments for error handling:
#      - character(len=512), intent(out) :: errmsg
#          - errmsg must be initialized as '' and contains the error message in case an error occurs
#      - integer, intent(out) :: ierr
#          - ierr must be initialized as 0 and set to >1 in case of errors
# Output: This routine converts the argument tables for all subroutines / typedefs / module variables into an XML file
# suitable to be used with mkcap.py (which generates the fortran code for the scheme cap)
# - the script generates a separate file for each module within the given files


# Items in this dictionary are used for checking valid entries in metadata tables. For columsn with no keys/keys
# commented out, no check is performed. This is the case for 'type' and 'kind' right now, since models use their
# own derived data types and kind types.
VALID_ITEMS = {
    'header' : ['local_name', 'standard_name', 'long_name', 'units', 'rank', 'type', 'kind', 'intent', 'optional'],
    #'type' : ['character', 'integer', 'real', ...],
    #'kind' : ['default', 'kind_phys', ...],
    'intent' : ['none', 'in', 'out', 'inout'],
    'optional' : ['T', 'F'],
    }

# Mandatory variables that every scheme needs to have
CCPP_MANDATORY_VARIABLES = {
    'error_message' : Var(local_name    = 'errmsg',
                          standard_name = 'error_message',
                          long_name     = 'error message for error handling in CCPP',
                          units         = 'none',
                          type          = 'character',
                          rank          = '',
                          kind          = 'len=*',
                          intent        = 'out',
                          optional      = 'F'
                          ),
    'error_flag' : Var(local_name    = 'ierr',
                       standard_name = 'error_flag',
                       long_name     = 'error flag for error handling in CCPP',
                       units         = 'flag',
                       type          = 'integer',
                       rank          = '',
                       kind          = '',
                       intent        = 'out',
                       optional      = 'F'
                       ),
    }

def merge_metadata_dicts(x, y):
    """Merges two metadata dictionaries. For each list of elements
    (variables = class Var in mkcap.py) in one dictionary, we know
    that all entries are compatible. If one or more elements exist
    in both x and y, we therefore have to test compatibility of
    one of the items in each dictionary only."""
    z = {}
    x_keys = sorted(x.keys())
    y_keys = sorted(y.keys())
    z_keys = sorted(list(set(x_keys + y_keys)))
    for key in z_keys:
        z[key] = {}
        if key in x_keys and key in y_keys:
            # We know that all entries within each dictionary are comptable;
            # we need to test compatibility of one of the items in each only.
            if not x[key][0].compatible(y[key][0]):
                raise Exception('Incompatible entries in metadata for variable {0}:\n'.format(key) +\
                                '    {0}\n'.format(x[key][0].print_debug()) +\
                                'vs. {0}'.format(y[key][0].print_debug()))
            z[key] = x[key] + y[key]
        elif key in x_keys:
            z[key] = x[key]
        elif key in y_keys:
            z[key] = y[key]
    return z


def parse_variable_tables(filename):
    """Parses metadata tables on the host model side that define the available variables.
    Metadata tables can refer to variables inside a module or as part of a derived
    datatype, which itself is defined inside a module (depending on the location of the
    metadata table). Each variable (standard_name) can exist only once, i.e. each entry
    (list of variables) in the metadata dictionary contains only one element
    (variable = instance of class Var defined in mkcap.py)"""
    # Set debug to true if logging level is debug
    debug = logging.getLogger().getEffectiveLevel() == logging.DEBUG

    # Final metadata container for all variables in file
    metadata = collections.OrderedDict()

    # Registry of modules and derived data types in file
    registry = collections.OrderedDict()

    # Read all lines of the file at once
    with (open(filename, 'r')) as file:
        file_lines = file.readlines()

    lines = []
    buffer = ''
    for i in xrange(len(file_lines)):
        line = file_lines[i].rstrip('\n').strip()
        # Skip empty lines
        if line == '' or line == '&':
            continue
        # Remove line continuations: concatenate with following lines
        if line.endswith('&'):
            buffer += file_lines[i].rstrip('\n').replace('&', ' ')
            continue
        # Write out line with buffer and reset buffer
        lines.append(buffer + file_lines[i].rstrip('\n').replace('&', ' '))
        buffer = ''
    del file_lines

    # Find all modules within the file, and save the start and end lines
    module_lines = {}
    line_counter = 0
    for line in lines:
        words = line.split()
        if len(words) > 1 and words[0].lower() == 'module' and not words[1].lower() == 'procedure':
            module_name = words[1].strip()
            if module_name in registry.keys():
                raise Exception('Duplicate module name {0}'.format(module_name))
            registry[module_name] = {}
            module_lines[module_name] = { 'startline' : line_counter }
        elif len(words) > 1 and words[0].lower() == 'end' and words[1].lower() == 'module':
            try:
                test_module_name = words[2]
            except IndexError:
                logging.warning('Encountered closing statement "end module" without module name; assume module_name is {0}'.format(module_name))
                test_module_name = module_name
            if not module_name == test_module_name:
                raise Exception('Module names in opening/closing statement do not match: {0} vs {1}'.format(module_name, test_module_name))
            module_lines[module_name]['endline'] = line_counter
        line_counter += 1

    # Parse each module in the file separately
    for module_name in registry.keys():
        startline = module_lines[module_name]['startline']
        endline = module_lines[module_name]['endline']
        line_counter = 0
        in_type = False
        for line in lines[startline:endline]:
            current_line_number = startline + line_counter
            words = line.split()
            for j in range(len(words)):
                # Check for the word 'type', that it is the first word in the line,
                # and that a name exists afterwards. It is assumed that definitions
                # (not usage) of derived types cannot be nested - reasonable for Fortran.
                if words[j].lower() == 'type' and j == 0 and len(words) > 1 and not '(' in words[j+1]:
                    if in_type:
                        raise Exception('Nested definitions of derived types not supported')
                    in_type = True
                    type_name = words[j+1].split('(')[0].strip()
                    if type_name in registry[module_name].keys():
                        raise Exception('Duplicate derived type name {0} in module {1}'.format(
                                                                       type_name, module_name))
                    registry[module_name][type_name] = [current_line_number]
                elif words[j].lower() == 'type' and j == 1 and words[j-1].lower() == 'end':
                    if not in_type:
                        raise Exception('Encountered "end_type" without corresponding "type" statement')
                    try:
                        test_type_name = words[j+1]
                    except IndexError:
                        logging.warning('Encountered closing statement "end type" without type name; assume type_name is {0}'.format(type_name))
                        test_type_name = type_name
                    if not type_name == test_type_name:
                        raise Exception('Type names in opening/closing statement do not match: {0} vs {1}'.format(type_name, test_type_name))
                    in_type = False
                    registry[module_name][type_name].append(current_line_number)
            line_counter += 1
        logging.debug('Parsing file {0} with registry {1}'.format(filename, registry))

        # Variables can either be defined at module-level or in derived types - alongside with their tables
        line_counter = 0
        in_table = False
        in_type = False
        for line in lines[startline:endline]:
            current_line_number = startline + line_counter

            # Check for beginning of new table
            words = line.split()
            # This is case sensitive
            if len(words) > 2 and words[0] in ['!!', '!>'] and '\section' in words[1] and 'arg_table_' in words[2]:
                if in_table:
                    raise Exception('Encountered table start for table {0} while still in table {1}'.format(words[2].replace('arg_table_',''), table_name))
                table_name = words[2].replace('arg_table_','')
                if not (table_name == module_name or table_name in registry[module_name].keys()):
                    raise Exception('Encountered table with name {0} without corresponding module or type name'.format(table_name))
                in_table = True
                header_line_number = current_line_number + 1
                line_counter += 1
                continue
            elif (words[0].startswith('!!') or words[0].startswith('!>')) and '\section' in words[0]:
                raise Exception("Malformatted table found in {0} / {1} / {2}".format(filename, module_name, table_name))
            # If an argument table is found, parse it
            if in_table:
                words = line.split('|')
                # Separate the table headers
                if current_line_number == header_line_number:
                    # Check for blank table
                    if len(words) <= 1:
                        logging.debug('Skipping blank table {0}'.format(table_name))
                        in_table = False
                        line_counter += 1
                        continue
                    table_header = [x.strip() for x in words[1:-1]]
                    # Check that only valid table headers are used
                    for item in table_header:
                        if not item in VALID_ITEMS['header']:
                            raise Exception('Invalid column header {0} in argument table {1}'.format(item, table_name))
                    # Locate mandatory column 'standard_name'
                    try:
                        standard_name_index = table_header.index('standard_name')
                    except ValueError:
                        raise Exception('Mandatory column standard_name not found in argument table {0}'.format(table_name))
                    line_counter += 1
                    continue
                elif current_line_number == header_line_number + 1:
                    # Skip over separator line
                    line_counter += 1
                    continue
                else:
                    if len(words) == 1:
                        # End of table
                        if words[0].strip() == '!!':
                            in_table = False
                            line_counter += 1
                            continue
                        else:
                            raise Exception('Encountered invalid line "{0}" in argument table {1}'.format(line, table_name))
                    else:
                        var_items = [x.strip() for x in words[1:-1]]
                        if not len(var_items) == len(table_header):
                            raise Exception('Error parsing variable entry "{0}" in argument table {1}'.format(var_items, table_name))
                        var_name = var_items[standard_name_index]
                        # Skip variables without a standard_name (i.e. empty cell in column standard_name)
                        if var_name:
                            var = Var.from_table(table_header,var_items)
                            if table_name == module_name:
                                container = encode_container(module_name)
                            else:
                                container = encode_container(module_name, table_name)
                            var.container = container
                            # Check for incompatible definitions with CCPP mandatory variables
                            if var_name in CCPP_MANDATORY_VARIABLES.keys() and not CCPP_MANDATORY_VARIABLES[var_name].compatible(var):
                                raise Exception('Entry for variable {0}'.format(var_name) + \
                                                ' in argument table {0}'.format(table_name) +\
                                                ' is incompatible with mandatory variable:\n' +\
                                                '    existing: {0}\n'.format(CCPP_MANDATORY_VARIABLES[var_name].print_debug()) +\
                                                '     vs. new: {0}'.format(var.print_debug()))
                            # Add variable to metadata dictionary
                            if not var_name in metadata.keys():
                                metadata[var_name] = [var]
                            else:
                                for existing_var in metadata[var_name]:
                                    if not existing_var.compatible(var):
                                        raise Exception('New entry for variable {0}'.format(var_name) + \
                                                        ' in argument table {0}'.format(table_name) +\
                                                        ' is incompatible with existing entry:\n' +\
                                                        '    existing: {0}\n'.format(existing_var.print_debug()) +\
                                                        '     vs. new: {0}'.format(var.print_debug()))

                                metadata[var_name].append(var)
                        #else:
                        #    logging.debug('Skipping variable entry "{0}" without a standard_name'.format(var_items))

            line_counter += 1

        # Informative output to screen
        if debug and len(metadata.keys()) > 0:
            for module_name in registry.keys():
                logging.debug('Module name: {0}'.format(module_name))
                container = encode_container(module_name)
                vars_in_module = []
                for var_name in metadata.keys():
                    for var in metadata[var_name]:
                        if var.container == container:
                            vars_in_module.append(var_name)
                logging.debug('Module variables: {0}'.format(', '.join(vars_in_module)))
                for type_name in registry[module_name].keys():
                    container = encode_container(module_name, type_name)
                    vars_in_type = []
                    for var_name in metadata.keys():
                        for var in metadata[var_name]:
                            if var.container == container:
                                vars_in_type.append(var_name)
                    logging.debug('Variables in derived type {0}: {1}'.format(type_name, ', '.join(vars_in_type)))

        if debug and len(metadata.keys()) > 0:
            # Write out the XML for debugging purposes
            top = ET.Element('definition')
            top.set('module', module_name)
            container = encode_container(module_name)
            for var_name in metadata.keys():
                for var in metadata[var_name]:
                    if var.container == container:
                        sub_var = var.to_xml(ET.SubElement(top, 'variable'))
            for type_name in registry[module_name].keys():
                container = encode_container(module_name, type_name)
                sub_type = ET.SubElement(top, 'type')
                sub_type.set('type_name', type_name)
                for var_name in metadata.keys():
                    for var in metadata[var_name]:
                        if var.container == container:
                            sub_var = var.to_xml(ET.SubElement(sub_type, 'variable'))
            indent(top)
            tree = ET.ElementTree(top)
            xmlfile = module_name + '.xml'
            tree.write(xmlfile, xml_declaration=True, encoding='utf-8', method="xml")
            logging.info('Parsed variable definition tables in module {0}; output => {1}'.format(module_name, xmlfile))
        elif len(metadata.keys()) > 0:
            logging.info('Parsed variable definition tables in module {0}'.format(module_name))

    return metadata


def parse_scheme_tables(filename):
    """Parses metadata tables for a physics scheme that requests/requires variables as
    input arguments. Metadata tables can only describe variables required by a subroutine
    'subroutine_name' of scheme 'scheme_name' inside a module 'module_name'. Each variable
    (standard_name) can exist only once, i.e. each entry (list of variables) in the metadata
    dictionary  contains only one element (variable = instance of class Var defined in
    mkcap.py). The metadata dictionaries of the individual schemes are merged afterwards
    (called from ccpp_prebuild.py) using merge_metadata_dicts, where multiple instances
    of variables are compared for compatibility and collected in a list (entry in the 
    merged metadata dictionary). The merged metadata dictionary of all schemes (which
    contains only compatible variable instances in the list referred to by standard_name)
    is then compared to the unique definition in the metadata dictionary of the variables
    provided by the host model using compare_metadata in ccpp_prebuild.py."""

    # Set debug to true if logging level is debug
    debug = logging.getLogger().getEffectiveLevel() == logging.DEBUG

    # Valid suffices for physics scheme routines
    subroutine_suffices = [ 'init', 'run', 'finalize']

    # Final metadata container for all variables in file
    metadata = collections.OrderedDict()

    # Registry of modules and derived data types in file
    #registry = {}
    registry = collections.OrderedDict()

    # Argument lists of each subroutine in the file
    arguments = collections.OrderedDict()

    # Read all lines of the file at once
    with (open(filename, 'r')) as file:
        file_lines = file.readlines()

    lines = []
    original_line_numbers = []
    buffer = ''
    for i in xrange(len(file_lines)):
        line = file_lines[i].rstrip('\n').strip()
        # Skip empty lines
        if line == '' or line == '&':
            continue
        # Remove line continuations: concatenate with following lines
        if line.endswith('&'):
            buffer += file_lines[i].rstrip('\n').replace('&', ' ')
            continue
        # Write out line with buffer and reset buffer
        lines.append(buffer + file_lines[i].rstrip('\n').replace('&', ' '))
        original_line_numbers.append(i+1)
        buffer = ''
    del file_lines

    # Find all modules within the file, and save the start and end lines
    module_lines = {}
    line_counter = 0
    for line in lines:
        words = line.split()
        if len(words) > 1 and words[0].lower() == 'module' and not words[1].lower() == 'procedure':
            module_name = words[1].strip()
            if module_name in registry.keys():
                raise Exception('Duplicate module name {0}'.format(module_name))
            registry[module_name] = {}
            module_lines[module_name] = { 'startline' : line_counter }
        elif len(words) > 1 and words[0].lower() == 'end' and words[1].lower() == 'module':
            try:
                test_module_name = words[2]
            except IndexError:
                logging.warning('Warning, encountered closing statement "end module" without module name; assume module_name is {0}'.format(module_name))
                test_module_name = module_name
            if not module_name == test_module_name:
                raise Exception('Module names in opening/closing statement do not match: {0} vs {1}'.format(module_name, test_module_name))
            module_lines[module_name]['endline'] = line_counter
        line_counter += 1

    # Parse each module in the file separately
    for module_name in registry.keys():
        startline = module_lines[module_name]['startline']
        endline = module_lines[module_name]['endline']
        line_counter = 0
        in_subroutine = False
        for line in lines[startline:endline]:
            current_line_number = startline + line_counter
            words = line.split()
            for j in range(len(words)):
                # Check for the word 'subroutine', that it is the first word in the line,
                # and that a name exists afterwards. Nested subroutines are ignored.
                if words[j].lower() == 'subroutine' and j == 0 and len(words) > 1:
                    if in_subroutine:
                        logging.debug('Warning, ignoring nested subroutine in module {0} and subroutine {1}'.format(module_name, subroutine_name))
                        continue
                    subroutine_name = words[j+1].split('(')[0].strip()
                    # Consider the last substring separated by a '_' of the subroutine name as a 'postfix'
                    if subroutine_name.find('_') >= 0:
                        subroutine_suffix = subroutine_name.split('_')[-1]
                        if subroutine_suffix in subroutine_suffices:
                            scheme_name = subroutine_name[0:subroutine_name.rfind('_')]
                            if not scheme_name == module_name:
                                raise Exception('Scheme name differs from module name: module_name="{0}" vs. scheme_name="{1}"'.format(
                                                                                                             module_name, scheme_name))
                            if not scheme_name in registry[module_name].keys():
                                registry[module_name][scheme_name] = {}
                            if subroutine_name in registry[module_name][scheme_name].keys():
                                raise Exception('Duplicate subroutine name {0} in module {1}'.format(
                                                                       subroutine_name, module_name))
                            registry[module_name][scheme_name][subroutine_name] = [current_line_number]
                            in_subroutine = True
                elif words[j].lower() == 'subroutine' and j == 1 and words[j-1].lower() == 'end':
                    try:
                        test_subroutine_name = words[j+1]
                    except IndexError:
                        logging.warning('Warning, encountered closing statement "end subroutine" without subroutine name; ' +\
                                        ' assume subroutine_name is {0}'.format(subroutine_name))
                        test_subroutine_name = subroutine_name
                    if in_subroutine and subroutine_name == test_subroutine_name:
                        in_subroutine = False
                        registry[module_name][scheme_name][subroutine_name].append(current_line_number)
                # Avoid problems by enforcing end statements to carry a descriptor (subroutine, module, ...)
                elif in_subroutine and len(words) == 1 and words[0].lower() == 'end':
                    raise Exception('Encountered closing statement "end" without descriptor (subroutine, module, ...): ' +\
                                    'line {0}="{1}" in file {2}'.format(original_line_numbers[current_line_number], line, filename))
            line_counter += 1

        # Check that for each registered subroutine the start and end lines were found
        for scheme_name in registry[module_name].keys():
            for subroutine_name in registry[module_name][scheme_name].keys():
                if not len(registry[module_name][scheme_name][subroutine_name]) == 2:
                    raise Exception('Error parsing start and end lines for subroutine {0} in module {1}'.format(subroutine_name, module_name))
        logging.debug('Parsing file {0} with registry {1}'.format(filename, registry))

        for scheme_name in registry[module_name].keys():
            for subroutine_name in registry[module_name][scheme_name].keys():
                # Record the order of variables in the call list to each subroutine in a list
                if not module_name in arguments.keys():
                    arguments[module_name] = {}
                if not scheme_name in arguments[module_name].keys():
                    arguments[module_name][scheme_name] = {}
                if not subroutine_name in arguments[module_name][scheme_name].keys():
                    arguments[module_name][scheme_name][subroutine_name] = []
                # Find the argument table corresponding to each subroutine by searching
                # "upward" from the subroutine definition line for the "arg_table_SubroutineName" section
                table_found = False
                header_line_number = None
                for line_number in range(registry[module_name][scheme_name][subroutine_name][0], -1, -1):
                    line = lines[line_number]
                    words = line.split()
                    for word in words:
                        if (len(words) > 2 and words[0] in ['!!', '!>'] and '\section' in words[1] and 'arg_table_{0}'.format(subroutine_name) in words[2]):
                            table_found = True
                            header_line_number = line_number + 1
                            table_name = subroutine_name
                            break
                        else:
                            for word in words:
                                if 'arg_table_{0}'.format(subroutine_name) in word:
                                    raise Exception("Malformatted table found in {0} / {1} / {2} / {3}".format(filename, module_name, scheme_name, subroutine_name))
                    if table_found:
                        break
                # If an argument table is found, parse it
                if table_found:
                    # Separate the table headers
                    table_header = lines[header_line_number].split('|')
                    # Check for blank table
                    if len(table_header) <= 1:
                        logging.debug('Skipping blank table {0}'.format(table_name))
                        table_found = False
                        continue
                    # Extract table header
                    table_header = [x.strip() for x in table_header[1:-1]]
                    # Check that only valid table headers are used
                    for item in table_header:
                        if not item in VALID_ITEMS['header']:
                            raise Exception('Invalid column header {0} in argument table {1}'.format(item, table_name))
                    # Locate mandatory column 'standard_name'
                    try:
                        standard_name_index = table_header.index('standard_name')
                    except ValueError:
                        raise Exception('Mandatory column standard_name not found in argument table {0}'.format(table_name))
                    # Get all of the variable information in table
                    end_of_table = False
                    line_number = header_line_number + 2
                    while not end_of_table:
                        line = lines[line_number]
                        words = line.split('|')
                        if len(words) == 1:
                            if words[0] == '!!':
                                end_of_table = True
                            else:
                                raise Exception('Encountered invalid line "{0}" in argument table {1}'.format(line, table_name))
                        else:
                            var_items = [x.strip() for x in words[1:-1]]
                            if not len(var_items) == len(table_header):
                                raise Exception('Error parsing variable entry "{0}" in argument table {1}'.format(var_items, table_name))
                            var_name = var_items[standard_name_index]
                            # Column standard_name cannot be left blank in scheme_tables
                            if not var_name:
                                raise Exception('Encountered line "{0}" without standard name in argument table {1}'.format(line, table_name))
                            # Add standard_name to argument list for this subroutine
                            arguments[module_name][scheme_name][subroutine_name].append(var_name)
                            var = Var.from_table(table_header,var_items)
                            # Check for incompatible definitions with CCPP mandatory variables
                            if var_name in CCPP_MANDATORY_VARIABLES.keys() and not CCPP_MANDATORY_VARIABLES[var_name].compatible(var):
                                raise Exception('Entry for variable {0}'.format(var_name) + \
                                                ' in argument table of subroutine {0}'.format(subroutine_name) +\
                                                ' is incompatible with mandatory variable:\n' +\
                                                '    existing: {0}\n'.format(CCPP_MANDATORY_VARIABLES[var_name].print_debug()) +\
                                                '     vs. new: {0}'.format(var.print_debug()))
                            # Record the location of this variable: module, scheme, table
                            container = encode_container(module_name, scheme_name, table_name)
                            var.container = container
                            # Add variable to metadata dictionary
                            if not var_name in metadata.keys():
                                metadata[var_name] = [var]
                            else:
                                for existing_var in metadata[var_name]:
                                    if not existing_var.compatible(var):
                                        raise Exception('New entry for variable {0}'.format(var_name) + \
                                                        ' in argument table of subroutine {0}'.format(subroutine_name) +\
                                                        ' is incompatible with existing entry:\n' +\
                                                        '    existing: {0}\n'.format(existing_var.print_debug()) +\
                                                        '     vs. new: {0}'.format(var.print_debug()))
                                metadata[var_name].append(var)

                        line_number += 1

                    # After parsing entire metadata table for the subroutine, check that all mandatory CCPP variables are present
                    for var_name in CCPP_MANDATORY_VARIABLES.keys():
                        if not var_name in arguments[module_name][scheme_name][subroutine_name]:
                            raise Exception('Mandatory CCPP variable {0} not declared in metadata table of subroutine {1}'.format(
                                                                                                       var_name, subroutine_name))

        # For CCPP-compliant files (i.e. files with metadata tables, perform additional checks)
        if len(metadata.keys()) > 0:
            # Check that all subroutine "root" names in the current module are equal to scheme_name
            # and that there are exactly three subroutines for scheme X: X_init, X_run, X_finalize
            message = ''
            abort = False
            for scheme_name in registry[module_name].keys():
                # Pre-generate error message
                message += 'Check that all subroutines in module {0} have the same root name:\n'.format(module_name)
                message += '    i.e. scheme_A_init, scheme_A_run, scheme_A_finalize\n'
                message += 'Here is a list of the subroutine names for scheme {0}:\n'.format(scheme_name)
                message += '{0}\n\n'.format(', '.join(sorted(registry[module_name][scheme_name].keys())))
                if (not len(registry[module_name][scheme_name].keys()) == 3):
                    logging.exception(message)
                    abort = True
                else:
                    for suffix in subroutine_suffices:
                        subroutine_name = '{0}_{1}'.format(scheme_name, suffix)
                        if not subroutine_name in registry[module_name][scheme_name].keys():
                            logging.exception(message)
                            abort = True
            if abort:
                raise Exception(message)

        # Debugging output to screen and to XML
        if debug and len(metadata.keys()) > 0:
            # To screen
            logging.debug('Module name: {0}'.format(module_name))
            for scheme_name in registry[module_name].keys():
                logging.debug('Scheme name: {0}'.format(scheme_name))
                for subroutine_name in registry[module_name][scheme_name].keys():
                    container = encode_container(module_name, scheme_name, subroutine_name)
                    vars_in_subroutine = []
                    for var_name in metadata.keys():
                        for var in metadata[var_name]:
                            if var.container == container:
                                vars_in_subroutine.append(var_name)
                    logging.debug('Variables in subroutine {0}: {1}'.format(subroutine_name, ', '.join(vars_in_subroutine)))
            # To XML
            for scheme_name in registry[module_name].keys():
                top = ET.Element('scheme')
                top.set('module', scheme_name)
                for subroutine_name in registry[module_name][scheme_name].keys():
                    sub_sub = ET.SubElement(top, 'subroutine')
                    sub_sub.set('name', subroutine_name)
                    container = encode_container(module_name, scheme_name, subroutine_name)
                    # Variable output in order of calling arguments
                    for var_name in arguments[module_name][scheme_name][subroutine_name]:
                        for var in metadata[var_name]:
                            if var.container == container:
                                sub_var = var.to_xml(ET.SubElement(sub_sub, 'variable'))
                indent(top)
                tree = ET.ElementTree(top)
                xmlfile = scheme_name + '.xml'
                tree.write(xmlfile, xml_declaration=True, encoding='utf-8', method="xml")
                logging.info('Parsed tables in scheme {0}; output => {1}'.format(scheme_name, xmlfile))
        # Standard output to screen
        elif len(metadata.keys()) > 0:
            for scheme_name in registry[module_name].keys():
                logging.info('Parsed tables in scheme {0}'.format(scheme_name))

    # End of loop over all module_names

    return (metadata, arguments)
