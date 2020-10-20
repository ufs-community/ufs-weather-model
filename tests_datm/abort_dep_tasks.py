#!/usr/bin/env python
import ecflow
import re

# this script will work ONLY for standalone nmmb regression test ecflow workflow

class DefsTraverser:

    def __init__(self, defs, ci):
        assert (isinstance(defs, ecflow.Defs)),"Expected ecflow.Defs as first argument"
        assert (isinstance(ci, ecflow.Client)),"Expected ecflow.Client as second argument"
        self.__defs = defs
        self.__ci = ci
        self.__suite = None

    def force_abort(self):
        for suite in self.__defs.suites:
            self.__suite = suite
            self.__walk_node(suite)

    def __walk_node(self, node_container):
        for node in node_container.nodes:
            if isinstance(node, ecflow.Task):
                self.__check_task(node)
            else:
                self.__walk_node(node)

    def __check_task(self, node):
        trigger_expr = node.get_trigger()
        if trigger_expr:
            tasks = re.findall( r'(\S*) ==', trigger_expr.get_expression())
            for t in tasks:
                task = self.__defs.find_abs_node( self.__suite.get_abs_node_path() + "/" + t)
                if task.get_state() == ecflow.State.aborted:
                    if node.get_state() != ecflow.State.aborted:
                        print "Will force aborted state for task", node.get_abs_node_path()
                        self.__ci.force_state(node.get_abs_node_path(), ecflow.State.aborted)

try:
    # Create the client. This will read the default environment variables
    ci = ecflow.Client()

    # Get the node tree suite definition as stored in the server
    # The definition is retrieved and stored on the variable 'ci'
    ci.sync_local()

    # access the definition retrieved from the server
    server_defs = ci.get_defs()

    if server_defs == None :
        print "The server has no definition"
        exit(1)

    traverser = DefsTraverser(server_defs, ci)
    traverser.force_abort()

except RuntimeError, e:
    print "failed: " + str(e)
