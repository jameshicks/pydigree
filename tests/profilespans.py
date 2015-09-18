try:
    import line_profiler
except ImportError:
    print("No line profiler, skipping test.")
    import sys
    sys.exit(0)

from pydigree.common import spans
func = spans

test = [1] * 10000

profile = line_profiler.LineProfiler(spans)
profile.runcall(spans, test)
profile.print_stats()
