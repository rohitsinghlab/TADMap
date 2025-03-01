#!/usr/bin/env python

import pathlib

PROJ_DIR=str(pathlib.Path.home())  #"../.."

DEBUG_MODE=False #True

def dbg_print(s, *args, **kwargs):
    global DEBUG_MODE
    if DEBUG_MODE or s[:5] != "Flag ":
        print(s, *args, **kwargs)


import logging

tadmap_loglevel = logging.WARNING #WARNING

def set_loglevel(l):
    """Set the loglevel to one of logging.{ERROR,WARNING,INFO,DEBUG}. Default = WARNING
"""
    
    global tadmap_loglevel
    print("Changing loglevel from %s to %s" % (tadmap_loglevel, l))
    tadmap_loglevel = l
    if tadmap_loglevel == logging.DEBUG:
        global DEBUG_MODE
        DEBUG_MODE = True
    

def tadmap_debug(*args, **kwargs):
    if tadmap_loglevel <= logging.DEBUG: print("DEBUG: ", *args, **kwargs)

def tadmap_info(*args, **kwargs):
    if tadmap_loglevel <= logging.INFO: print("INFO: ", *args, **kwargs)

def tadmap_warning(*args, **kwargs):
    if tadmap_loglevel <= logging.WARNING: print("WARNING: ", *args, **kwargs)

def tadmap_error(*args, **kwargs):
    if tadmap_loglevel <= logging.ERROR: print("ERROR: ", *args, **kwargs)

    
########## for maintenance ###################
# def noop(*args, **kwargs):
#     pass
#
# logging.info = print
# logging.debug = noop
##############################################
