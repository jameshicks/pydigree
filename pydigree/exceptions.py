"""
Exceptions raised by pydigree routines
"""


class IterationError(Exception):
    """
    The procedure has run out of iterations without finding a result
    """

    pass


class NotMeaningfulError(Exception):
    """
    The error raised for something that does not a meaningul result"
    """

    pass


class SimulationError(Exception):
    """
    Error raised when an error has occurred during a simulation
    """
    pass


class FileFormatError(Exception):
    """
    Error raised when a problem is encountered parsing a file
    """
    
    pass
