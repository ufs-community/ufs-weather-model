#!/usr/bin/env python

# Standard modules
import argparse
import collections
import logging
import os
import sys

# DH* TODO
# CONSISTENCY CHECK BETWEEN OPTIONAL ARGUMENTS IN THE METADATA TABLE AND IN
# THE ACTUAL ARGUMENT LIST / FORTRAN VARIABLE DECLARATIONS (RANKS, TYPE, INTENT).
# *DH

# Local modules
from common import encode_container, execute
from metadata_parser import merge_metadata_dicts, parse_scheme_tables, parse_variable_tables
from mkcap import Cap, CapsMakefile, CapsCMakefile, SchemesMakefile, SchemesCMakefile
from mkdoc import metadata_to_html, metadata_to_latex

###############################################################################
# User definitions                                                            #
###############################################################################

# List of configured host models
HOST_MODELS = ["FV3v0", "FV3v1", "SCM"]

###############################################################################
# Set up the command line argument parser and other global variables          #
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--model', action='store', choices=HOST_MODELS, help='host model (case-sensitive)', required=True)
parser.add_argument('--debug', action='store_true', help='enable debugging output', default=False)

# BASEDIR is the current directory where this script is executed
BASEDIR = os.getcwd()

###############################################################################
# Functions and subroutines                                                   #
###############################################################################

def parse_arguments():
    """Parse command line arguments."""
    success = True
    args = parser.parse_args()
    host_model = args.model
    debug = args.debug
    return (success, host_model, debug)

def import_config(host_model):
    """Import the configuration file for a given host model"""
    success = True
    config = {}

    # Import the host-model specific CCPP prebuild config
    ccpp_prebuild_config_name = "ccpp_prebuild_config_{0}".format(host_model)
    ccpp_prebuild_config = __import__(ccpp_prebuild_config_name)

    # Definitions in host-model dependent CCPP prebuild config script
    config['variable_definition_files'] = ccpp_prebuild_config.VARIABLE_DEFINITION_FILES
    config['scheme_files']              = ccpp_prebuild_config.SCHEME_FILES
    config['scheme_files_dependencies'] = ccpp_prebuild_config.SCHEME_FILES_DEPENDENCIES
    config['schemes_makefile']          = ccpp_prebuild_config.SCHEMES_MAKEFILE
    config['schemes_cmakefile']         = ccpp_prebuild_config.SCHEMES_CMAKEFILE
    config['target_files']              = ccpp_prebuild_config.TARGET_FILES
    config['caps_makefile']             = ccpp_prebuild_config.CAPS_MAKEFILE
    config['caps_cmakefile']            = ccpp_prebuild_config.CAPS_CMAKEFILE
    config['caps_dir']                  = ccpp_prebuild_config.CAPS_DIR
    config['optional_arguments']        = ccpp_prebuild_config.OPTIONAL_ARGUMENTS
    config['module_include_file']       = ccpp_prebuild_config.MODULE_INCLUDE_FILE
    config['fields_include_file']       = ccpp_prebuild_config.FIELDS_INCLUDE_FILE
    config['html_vartable_file']        = ccpp_prebuild_config.HTML_VARTABLE_FILE
    config['latex_vartable_file']       = ccpp_prebuild_config.LATEX_VARTABLE_FILE

    # Template code in host-model dependent CCPP prebuild config script
    config['module_use_template_host_cap']   = ccpp_prebuild_config.MODULE_USE_TEMPLATE_HOST_CAP
    config['ccpp_data_structure']            = ccpp_prebuild_config.CCPP_DATA_STRUCTURE
    config['module_use_template_scheme_cap'] = ccpp_prebuild_config.MODULE_USE_TEMPLATE_SCHEME_CAP

    return(success, config)

def setup_logging(debug):
    """Sets up the logging module and logging level."""
    success = True
    if debug:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(format='%(levelname)s: %(message)s', level=level)
    if debug:
        logging.info('Logging level set to DEBUG')
    else:
        logging.info('Logging level set to INFO')
    return success

def gather_variable_definitions(variable_definition_files):
    """Scan all Fortran source files with variable definitions on the host model side."""
    logging.info('Parsing metadata tables for variables provided by host model ...')
    success = True
    metadata_define = {}
    for variable_definition_file in variable_definition_files:
        (filedir, filename) = os.path.split(variable_definition_file)
        # Change to directory of variable_definition_file and parse it
        os.chdir(os.path.join(BASEDIR,filedir))
        metadata = parse_variable_tables(filename)
        metadata_define = merge_metadata_dicts(metadata_define, metadata)
        # Return to BASEDIR
        os.chdir(BASEDIR)
    return (success, metadata_define)

def collect_physics_subroutines(scheme_files):
    """Scan all Fortran source files in scheme_files for subroutines with argument tables."""
    logging.info('Parsing metadata tables in physics scheme files ...')
    success = True
    # Parse all scheme files
    metadata_request = {}
    arguments_request = {}
    for scheme_file in scheme_files:
        (scheme_filepath, scheme_filename) = os.path.split(os.path.abspath(scheme_file))
        # Change to directory where scheme_file lives
        os.chdir(scheme_filepath)
        (metadata, arguments) = parse_scheme_tables(scheme_filename)
        metadata_request = merge_metadata_dicts(metadata_request, metadata)
        arguments_request.update(arguments)
        os.chdir(BASEDIR)
    # Return to BASEDIR
    os.chdir(BASEDIR)
    return (success, metadata_request, arguments_request)

def check_optional_arguments(metadata, arguments, optional_arguments):
    """Check if for each subroutine with optional arguments, an entry exists in the
    optional_arguments dictionary. This is required to generate the caps correctly
    and to assess whether the variable is required from the host model. Optional
    arguments that are not requested by the model as specified in the dictionary
    optional_arguments are deleted from the list of requested data individually
    for each subroutine."""
    logging.info('Checking optional arguments in physics schemes ...')
    success = True
    for var_name in sorted(metadata.keys()):
        # The notation metadata[var_name][:] is a convenient way to make a copy
        # of the metadata[var_name] list, which allows removing items as we go
        for var in metadata[var_name][:]:
            if var.optional in ['t', 'T']:
                for item in var.container.split(' '):
                    subitems = item.split('_')
                    if subitems[0] == 'MODULE':
                        module_name = '_'.join(subitems[1:])
                    elif subitems[0] == 'SCHEME':
                        scheme_name = '_'.join(subitems[1:])
                    elif subitems[0] == 'SUBROUTINE':
                        subroutine_name = '_'.join(subitems[1:])
                    else:
                        success = False
                        logging.error('Invalid identifier {0} in container value {1} of requested variable {2}'.format(
                                                                                 subitems[0], var.container, var_name))
                if not module_name in optional_arguments.keys() or not \
                        subroutine_name in optional_arguments[module_name].keys():
                    success = False
                    logging.error('No entry found in optional_arguments dictionary for optional argument ' + \
                                  '{0} to subroutine {1} in module {2}'.format(var_name, subroutine_name, module_name))
                if type(optional_arguments[module_name][subroutine_name]) is list:
                    if var_name in optional_arguments[module_name][subroutine_name]:
                        logging.debug('Optional argument {0} to subroutine {1} in module {2} is required, keep in list'.format(
                                                                                       var_name, subroutine_name, module_name))
                    else:
                        logging.debug('Optional argument {0} to subroutine {1} in module {2} is not required, remove from list'.format(
                                                                                               var_name, subroutine_name, module_name))
                        # Remove this var instance from list of var instances for this var_name
                        metadata[var_name].remove(var)
                        # Remove var_name from list of calling arguments for that subroutine
                        arguments[module_name][scheme_name][subroutine_name].remove(var_name)
                elif optional_arguments[module_name][subroutine_name] == 'all':
                    logging.debug('optional argument {0} to subroutine {1} in module {2} is required, keep in list'.format(
                                                                                   var_name, subroutine_name, module_name))

        # If metadata[var_name] is now empty, i.e. the variable is not
        # requested at all by the model, remove the entry from metadata
        if not metadata[var_name]:
            del metadata[var_name]

    return (success, metadata, arguments)

def compare_metadata(metadata_define, metadata_request):
    """Compare the requested metadata to the defined one. For each requested entry, a
    single (i.e. non-ambiguous entry) must be present in the defined entries. All optional
    arguments that are still in the list of required variables for a scheme are needed,
    since they were checked in the routine check_optional_arguments beforehand."""

    logging.info('Comparing metadata for requested and provided variables ...')
    success = True
    modules = []
    metadata = {}
    for var_name in sorted(metadata_request.keys()):
        # Check that variable is provided by the model
        if not var_name in metadata_define.keys():
            requested_by = ' & '.join(var.container for var in metadata_request[var_name])
            success = False
            logging.error('Variable {0} requested by {1} not provided by the model'.format(var_name, requested_by))
            continue
        # Check that an unambiguous target exists for this variable
        if len(metadata_define[var_name]) > 1:
            success = False
            requested_by = ' & '.join(var.container for var in metadata_request[var_name])
            provided_by = ' & '.join(var.container for var in metadata_define[var_name])
            error_message = '  error, variable {0} requested by {1} cannot be identified unambiguously.'.format(var_name, requested_by) +\
                            ' Multiple definitions in {0}'.format(provided_by)
            logging.error(error_message)
            continue
        # Check that the variable properties are compatible between the model and the schemes
        if not metadata_request[var_name][0].compatible(metadata_define[var_name][0]):
            success = False
            error_message = '  incompatible entries in metadata for variable {0}:\n'.format(var_name) +\
                            '    provided:  {0}\n'.format(metadata_define[var_name][0].print_debug()) +\
                            '    requested: {0}'.format(metadata_request[var_name][0].print_debug())
            logging.error(error_message)
            continue
        # Construct the actual target variable and list of modules to use from the information in 'container'
        var = metadata_define[var_name][0]
        target = ''
        for item in var.container.split(' '):
            subitems = item.split('_')
            if subitems[0] == 'MODULE':
                modules.append('_'.join(subitems[1:]))
            elif subitems[0] == 'TYPE':
                pass
            else:
                logging.error('Unknown identifier {0} in container value of defined variable {1}'.format(subitems[0], var_name))
        target += var.local_name
        # Copy the length kind from the variable definition to update len=* in the variable requests
        if var.type == 'character':
            kind = var.kind
        metadata[var_name] = metadata_request[var_name]
        # Set target and kind (if applicable)
        for var in metadata[var_name]:
            var.target = target
            logging.debug('Requested variable {0} in {1} matched to target {2} in module {3}'.format(var_name, var.container, target, modules[-1]))
            # Update len=* for character variables
            if var.type == 'character' and var.kind == 'len=*':
                logging.debug('Update kind information for requested variable {0} in {1} from {2} to {3}'.format(var_name, var.container, var.kind, kind))
                var.kind = kind

    # Remove duplicated from list of modules
    modules = sorted(list(set(modules)))
    return (success, modules, metadata)

def create_module_use_statements(modules, module_use_template_host_cap):
    """Create Fortran module use statements to be included in the host cap.
    The template module_use_template_host_cap must include the required
    modules for error handling of the ccpp_field_add statments."""
    logging.info('Generating module use statements ...')
    success = True
    module_use_statements = module_use_template_host_cap
    cnt = 1
    for module in modules:
        module_use_statements += 'use {0}\n'.format(module)
        cnt += 1
    logging.info('Generated module use statements for {0} module(s)'.format(cnt))
    return (success, module_use_statements)

def create_ccpp_field_add_statements(metadata, ccpp_data_structure):
    """Create Fortran code to add host model variables to the cdata
    structure. The metadata container may contain multiple entries
    of a variable with the same standard_name, but for different
    "callers" (i.e. subroutines using it) with identical or
    different local_name. We only need to add it once to
    the add_field statement, since the target (i.e. the
    original variable defined by the model) is the same."""
    logging.info('Generating ccpp_field_add statements ...')
    success = True
    ccpp_field_add_statements = ''
    cnt = 0
    for var_name in sorted(metadata.keys()):
        # Add variable with var_name = standard_name once
        var = metadata[var_name][0]
        ccpp_field_add_statements += var.print_add(ccpp_data_structure)
        cnt += 1
    logging.info('Generated ccpp_field_add statements for {0} variable(s)'.format(cnt))
    return (success, ccpp_field_add_statements)

def generate_include_files(module_use_statements, ccpp_field_add_statements,
                           target_files, module_include_file, fields_include_file):
    """Generate include files for modules and field-add statements for host model cap."""
    logging.info('Generating include files for host model cap {0} ...'.format(', '.join(target_files)))
    success = True
    target_dirs = []
    for target_file in target_files:
        target_dirs.append(os.path.split(target_file)[0])
    target_dirs = sorted(list(set(target_dirs)))
    for target_dir in target_dirs:
        # module use statements
        includefile = os.path.join(target_dir, module_include_file)
        logging.info('Generated module-use include file {0}'.format(includefile))
        with open(includefile, "w") as f:
            f.write(module_use_statements)
        # ccpp_field_add statements
        includefile = os.path.join(target_dir, fields_include_file)
        logging.info('Generated fields-add include file {0}'.format(includefile))
        with open(includefile, "w") as f:
            f.write(ccpp_field_add_statements)
    return success

def generate_scheme_caps(metadata, arguments, caps_dir, module_use_template_scheme_cap):
    """Generate scheme caps for all schemes parsed."""
    success = True
    # Change to physics directory
    os.chdir(caps_dir)
    # List of filenames of scheme caps
    scheme_caps = []
    for module_name in arguments.keys():
        for scheme_name in arguments[module_name].keys():
            for subroutine_name in arguments[module_name][scheme_name].keys():
                # Skip subroutines without argument table or with empty argument table
                if not arguments[module_name][scheme_name][subroutine_name]:
                    continue
            # Create cap
            cap = Cap()
            cap.filename = "{0}_cap.F90".format(scheme_name)
            scheme_caps.append(cap.filename)
            # Parse all subroutines and their arguments to generate the cap
            capdata = collections.OrderedDict()
            for subroutine_name in arguments[module_name][scheme_name].keys():
                capdata[subroutine_name] = []
                for var_name in arguments[module_name][scheme_name][subroutine_name]:
                    container = encode_container(module_name, scheme_name, subroutine_name)
                    for var in metadata[var_name]:
                        if var.container == container:
                            capdata[subroutine_name].append(var)
                            break
            # If required (not at the moment), add module use statements to module_use_template_scheme_cap
            module_use_statement = module_use_template_scheme_cap
            # Write cap
            cap.write(module_name, module_use_statement, capdata)
    #
    os.chdir(BASEDIR)
    return (success, scheme_caps)

def generate_schemes_makefile(schemes, schemes_makefile, schemes_cmakefile):
    """Generate makefile/cmakefile snippets for all schemes."""
    logging.info('Generating schemes makefile/cmakefile snippet ...')
    success = True
    makefile = SchemesMakefile()
    makefile.filename = schemes_makefile
    cmakefile = SchemesCMakefile()
    cmakefile.filename = schemes_cmakefile
    # Adjust relative file path to schemes from caps makefile
    schemes_with_path = []
    schemes_makefile_dir = os.path.split(os.path.abspath(schemes_makefile))[0]
    for scheme in schemes:
        (scheme_filepath, scheme_filename) = os.path.split(os.path.abspath(scheme))
        relative_path = './{0}'.format(os.path.relpath(scheme_filepath, schemes_makefile_dir))
        schemes_with_path.append(os.path.join(relative_path, scheme_filename))
    makefile.write(schemes_with_path)
    cmakefile.write(schemes_with_path)
    logging.info('Added {0} schemes to makefile/cmakefile {1}/{2}'.format(
           len(schemes_with_path), makefile.filename, cmakefile.filename))
    return success

def generate_caps_makefile(caps, caps_makefile, caps_cmakefile, caps_dir):
    """Generate makefile/cmakefile snippets for all caps."""
    logging.info('Generating caps makefile/cmakefile snippet ...')
    success = True
    makefile = CapsMakefile()
    makefile.filename = caps_makefile
    cmakefile = CapsCMakefile()
    cmakefile.filename = caps_cmakefile
    # Adjust relative file path to schemes from caps makefile
    caps_makefile_dir = os.path.split(os.path.abspath(caps_makefile))[0]
    relative_path = './{0}'.format(os.path.relpath(caps_dir, caps_makefile_dir))
    caps_with_path = [ os.path.join(relative_path, cap) for cap in caps]
    makefile.write(caps_with_path)
    cmakefile.write(caps_with_path)
    logging.info('Added {0} auto-generated caps to makefile/cmakefile {1}/{2}'.format(
                          len(caps_with_path), makefile.filename, cmakefile.filename))
    return success

def main():
    """Main routine that handles the CCPP prebuild for different host models."""
    # Parse command line arguments
    (success, host_model, debug) = parse_arguments()
    if not success:
        raise Exception('Call to parse_arguments failed.')

    success = setup_logging(debug)
    if not success:
        raise Exception('Call to setup_logging failed.')

    (success, config) = import_config(host_model)
    if not success:
        raise Exception('Call to import_config failed.')

    # Variables defined by the host model
    (success, metadata_define) = gather_variable_definitions(config['variable_definition_files'])
    if not success:
        raise Exception('Call to gather_variable_definitions failed.')

    # Create an HTML table with all variables provided by the model
    success = metadata_to_html(metadata_define, host_model, config['html_vartable_file'])
    if not success:
        raise Exception('Call to metadata_to_html failed.')

    # Variables requested by the CCPP physics schemes
    (success, metadata_request, arguments_request) = collect_physics_subroutines(config['scheme_files'])
    if not success:
        raise Exception('Call to collect_physics_subroutines failed.')

    # Process optional arguments based on configuration in above dictionary optional_arguments
    (success, metadata_request, arguments_request) = check_optional_arguments(metadata_request,arguments_request,
                                                                              config['optional_arguments'])
    if not success:
        raise Exception('Call to check_optional_arguments failed.')

    # Create a LaTeX table with all variables requested by the pool of physics and/or provided by the host model
    success = metadata_to_latex(metadata_define, metadata_request, host_model, config['latex_vartable_file'])
    if not success:
        raise Exception('Call to metadata_to_latex failed.')

    # Check requested against defined arguments to generate metadata (list/dict of variables for CCPP)
    (success, modules, metadata) = compare_metadata(metadata_define, metadata_request)
    if not success:
        raise Exception('Call to compare_metadata failed.')

    # Crate module use statements to inject into the host model cap
    (success, module_use_statements) = create_module_use_statements(modules, config['module_use_template_host_cap'])
    if not success:
        raise Exception('Call to create_module_use_statements failed.')

    # Crate ccpp_fiels_add statements to inject into the host model cap
    (success, ccpp_field_add_statements) = create_ccpp_field_add_statements(metadata, config['ccpp_data_structure'])
    if not success:
        raise Exception('Call to create_ccpp_field_add_statements failed.')

    # Generate include files for module_use_statements and ccpp_field_add_statements
    success = generate_include_files(module_use_statements, ccpp_field_add_statements, config['target_files'],
                                     config['module_include_file'], config['fields_include_file'])
    if not success:
        raise Exception('Call to generate_include_files failed.')

    # Generate scheme caps
    (success, scheme_caps) = generate_scheme_caps(metadata_request, arguments_request, config['caps_dir'],
                                                  config['module_use_template_scheme_cap'])
    if not success:
        raise Exception('Call to generate_scheme_caps failed.')

    # Add filenames of schemes to makefile - add dependencies for schemes
    success = generate_schemes_makefile(config['scheme_files_dependencies'] + config['scheme_files'],
                                             config['schemes_makefile'], config['schemes_cmakefile'])
    if not success:
        raise Exception('Call to generate_schemes_makefile failed.')

    # Add filenames of scheme caps to makefile
    success = generate_caps_makefile(scheme_caps, config['caps_makefile'], config['caps_cmakefile'], config['caps_dir'])
    if not success:
        raise Exception('Call to generate_caps_makefile failed.')

    logging.info('CCPP prebuild step completed successfully.')

if __name__ == '__main__':
    main()
