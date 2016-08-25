import gzip
import bz2


def smartopen(filename, mode='r'):
    """
    Seamlessly open gzipped and bzipped2 files. Use like regular open.

	Python's compression modules iterate over compressed files as bytes,
	not strings. Unless 'b' (binary) is specified in the mode, we add 't' (text)
	to the mode to force iteration as strings, operating consistently 
    """
    if 't' not in mode:
    	mode = mode + 't'
    if filename.endswith('.gz'):
        return gzip.open(filename, mode, compresslevel=5)
    elif filename.endswith('.bz2'):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)
