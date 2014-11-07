import gzip
import bz2

def smartopen(filename, mode='r'):
    "Seamlessly open gzipped and bzipped2 files. Use like regular open"
    if filename.endswith('.gz'):
        return gzip.GzipFile(filename, mode, compresslevel=5)
    elif filename.endswith('.bz2'):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)
