##
#
# FILE: UnitTest.py
#
# AUTHORS: Alex Welton and Varun Ravishanker
#
# DESCRIPTION: Runs a series of unit tests (including benchmarking)
# and places results in a log file with the current date/time in the filename.
#

import time
import timeit
import unittest


# The "base" filename that current date/time info is appended to.
FILENAME_BASE = "TEST_LOG_"
FILENAME_EXT = ".txt"

# Return the constructed filename of the log file
def get_current_log_filename():

    # Get the date/time in nicely formatted output
    date_string = time.strftime("%m/%d/%y-%H:%M:%S")
    return FILENAME_BASE + date_string + FILENAME_EXT
