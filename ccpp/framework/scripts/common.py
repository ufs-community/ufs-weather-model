#!/usr/bin/env python

import logging
import subprocess


def execute(cmd, abort = True):
    """Runs a local command in a shell. Waits for completion and
    returns status, stdout and stderr. If abort = True, abort in
    case an error occurs during the execution of the command."""

    # Set debug to true if logging level is debug
    debug = logging.getLogger().getEffectiveLevel() == logging.DEBUG

    logging.debug('Executing "{0}"'.format(cmd))
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE,
                         stderr = subprocess.PIPE, shell = True)
    (stdout, stderr) = p.communicate()
    status = p.returncode
    if debug:
        message = 'Execution of "{0}" returned with exit code {1}\n'.format(cmd, status)
        message += '    stdout: "{0}"\n'.format(stdout.rstrip('\n'))
        message += '    stderr: "{0}"'.format(stderr.rstrip('\n'))
        logging.debug(message)
    if not status == 0:
        message = 'Execution of command {0} failed, exit code {1}\n'.format(cmd, status)
        message += '    stdout: "{0}"\n'.format(stdout.rstrip('\n'))
        message += '    stderr: "{0}"'.format(stderr.rstrip('\n'))
        if abort:
            raise Exception(message)
        else:
            logging.error(message)
    return (status, stdout.rstrip('\n'), stderr.rstrip('\n'))


def indent(elem, level=0):
    """Subroutine for writing "pretty" XML; copied from
    http://effbot.org/zone/element-lib.htm#prettyprint"""
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
          elem.tail = i


def encode_container(*args):
    """Encodes a container, i.e. the location of a metadata table for CCPP.
    Currently, there are three possibilities with different numbers of input
    arguments: module, module+typedef, module+scheme+subroutine."""
    if len(args)==3:
        container = 'MODULE_{0} SCHEME_{1} SUBROUTINE_{2}'.format(*args)
    elif len(args)==2:
        container = 'MODULE_{0} TYPE_{1}'.format(*args)
    elif len(args)==1:
        container = 'MODULE_{0}'.format(*args)
    else:
        raise Exception("encode_container not implemented for {0} arguments".format(len(args)))
    return container

def decode_container(container):
    """Decodes a container, i.e. the description of a location of a metadata table
    for CCPP. Currently, there are three possibilities with different numbers of
    input arguments: module, module+typedef, module+scheme+subroutine."""
    items = container.split(' ')
    if not len(items) in [1, 2, 3]:
        raise Exception("decode_container not implemented for {0} items".format(len(items)))
    for i in xrange(len(items)):
        items[i] = items[i][:items[i].find('_')] + ' ' + items[i][items[i].find('_')+1:]
    return ' '.join(items)

def escape_tex(text):
    """Substitutes characters for generating LaTeX sources files from Python."""
    return text.replace(
                '%', '\%').replace(
                '_', '\_')

