class IterationError(Exception):
    pass


class NotMeaningfulError(Exception):
    "The error raised for something that does not a meaningul result"
    pass


class SimulationError(Exception):
    "Error raised when an error has occurred during a simulation"
    pass


class FileFormatError(Exception):
    "Raised when a problem is encountered parsing a file"
    pass
