
""" Logging functions for the mkrefs

This module configures the logger for the jwst_reffiles package.
It also provides decorators to log the execution of modules. This code is
based on the JWQL package logging configuration code.

Authors
-------

    - Bryan Hilbert 2018
    - Catherine Martlin 2018 (JWQL version)
    - Alex Viana, 2013 (WFC3 QL Version)

Use
---

    To log the execution of a module, use:
    ::

        import os
        import logging

        from jwql.logging.logging_functions import configure_logging
        from jwql.logging.logging_functions import log_info
        from jwql.logging.logging_functions import log_fail

        @log_info
        @log_fail
        def my_main_function():
            pass

        if __name__ == '__main__':

            module = os.path.basename(__file__).replace('.py', '')
            configure_logging(module)

            my_main_function()

References
----------
    This code is adopted and updated from python routine
    ``logging_functions.py`` written by Alex Viana, 2013 for the WFC3
    Quicklook automation platform.
"""

import datetime
import getpass
import importlib
import logging
from logging import StreamHandler
import logging.config
from logging.handlers import RotatingFileHandler
import os
import socket
import sys
import time
import traceback

from functools import wraps

from jwst_reffiles.utils.permissions import set_permissions


def configure_logging(module, path='./', log_file_level='info', log_screen_level="info"):
    """Configure the log file with a standard logging format.

    Parameters
    ----------
    module : str
        The name of the module being logged.

    path : str
        Where to write the log if user-supplied path; default to working dir.

    log_file_level : str
        Minimum message level to route to the output log file.
        Allowed values: "debug", "info", "warning", "ciritical"

    log_screen_level : str
        Minimum message level to route to the screen.
        Allowed values: "debug", "info", "warning", "ciritical"

    Returns
    -------
    log_file : str
        Name and path of the output log file
    """

    # Determine log file location
    log_file = make_log_file(module, path=path)

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Create the Handler for logging data to a file
    logger_handler = RotatingFileHandler(log_file)
    logger_handler.setLevel(log_file_level.upper())

    print_to_screen = True
    # Create the Handler for logging data to console.
    if not os.path.exists(log_file):
        # if a log file exists, we don't set screen output
        console_handler = StreamHandler()
        console_handler.setLevel(log_screen_level.upper())
    else:
        print_to_screen = False

    # Create a Formatter for formatting the log messages
    logger_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                                         datefmt='%m/%d/%Y %H:%M:%S %p')

    # Add the Formatter to the Handler
    logger_handler.setFormatter(logger_formatter)

    if print_to_screen:
        console_handler.setFormatter(logger_formatter)

    # Add the Handler to the Logger
    root_logger.addHandler(logger_handler)

    if print_to_screen:
        root_logger.addHandler(console_handler)

    # Currently need to be on ST network for this to work??
    set_permissions(log_file, verbose=False)
    return log_file


def make_log_file(module, path='./'):
    """Create the log file name based on the module name.

    The name of the ``log_file`` is a combination of the name of the
    module being logged and the current datetime.

    Parameters
    ----------
    module : str
        The name of the module being logged.
    path : str
        Where to write the log if user-supplied path; default to
        working dir.

    Returns
    -------
    log_file : str
        The full path to where the log file will be written to.
    """

    timestamp = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M')
    filename = '{0}_{1}.log'.format(module, timestamp)
    log_file = os.path.join(path, filename)

    return log_file


def log_info(func):
    """Decorator to log useful system information.

    This function can be used as a decorator to log user environment
    and system information. Future packages we want to track can be
    added or removed as necessary.

    Parameters
    ----------
    func : func
        The function to decorate.

    Returns
    -------
    wrapped : func
        The wrapped function.
    """

    @wraps(func)
    def wrapped(*a, **kw):
        # For right now, skip this info, it clutters the output.
        print('HELLO wrapped')
    return wrapped

    def wrapped_old(*a, **kw):

        # Log environment information
        logging.info('User: ' + getpass.getuser())
        logging.info('System: ' + socket.gethostname())
        logging.info('Python Version: ' + sys.version.replace('\n', ''))
        logging.info('Python Executable Path: ' + sys.executable)

        # Log common module version information
        module_list = ['numpy', 'astropy', 'jwst']
        for module in module_list:
            try:
                mod = importlib.import_module(module)
                logging.info(module + ' Version: ' + mod.__version__)
                logging.info(module + ' Path: ' + mod.__path__[0])
            except ImportError as err:
                logging.warning(err)

        # Call the function and time it
        t1_cpu = time.clock()
        t1_time = time.time()
        func(*a, **kw)
        t2_cpu = time.clock()
        t2_time = time.time()

        # Log execution time
        hours_cpu, remainder_cpu = divmod(t2_cpu - t1_cpu, 60 * 60)
        minutes_cpu, seconds_cpu = divmod(remainder_cpu, 60)
        hours_time, remainder_time = divmod(t2_time - t1_time, 60 * 60)
        minutes_time, seconds_time = divmod(remainder_time, 60)
        logging.info('Elapsed Real Time: {0:.0f}:{1:.0f}:{2:f}'.format(hours_time, minutes_time,
                                                                       seconds_time))
        logging.info('Elapsed CPU Time: {0:.0f}:{1:.0f}:{2:f}'.format(hours_cpu, minutes_cpu,
                                                                      seconds_cpu))

    return wrapped


def log_fail(func):
    """Decorator to log crashes in the decorated code.

    Parameters
    ----------
    func : func
        The function to decorate.

    Returns
    -------
    wrapped : func
        The wrapped function.
    """

    @wraps(func)
    def wrapped(*a, **kw):

        try:

            # Run the function
            func(*a, **kw)
            logging.info('Completed Successfully')

        except Exception:
            logging.critical(traceback.format_exc())
            logging.critical('CRASHED')

    return wrapped
