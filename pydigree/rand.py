import random
import numpy as np

def choice(iter):
	l = len(iter)
	idx = np.random.randint(0, l)
	return iter[idx]

def sample_with_replacement(seq, n):
	l = len(seq)
	return [seq[np.random.randint(0,l)] for x in range(n)]