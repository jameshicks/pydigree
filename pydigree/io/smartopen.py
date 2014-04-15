import gzip
import bz2

def smartopen(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.GzipFile(filename, mode)
    elif filename.endswith('.bz2'):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)