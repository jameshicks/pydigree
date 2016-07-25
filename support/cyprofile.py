from pydigree.common import spans
from pydigree.genotypes import Alleles
try:
    import line_profiler
except ImportError:
    print("No line profiler, skipping test.")
    import sys
    sys.exit(0)

test_data = Alleles([0] * 10000)
func = spans

for start in range(1, 10000, 200):
	test_data[start:(start+100)] = 1

profile = line_profiler.LineProfiler(func)
profile.runcall(func, test_data)
profile.print_stats()