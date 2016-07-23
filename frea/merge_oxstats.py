import contextlib
import itertools
import sys

from .algorithms import *
from .formats import *

with contextlib.ExitStack() as stack:
    data = [stack.enter_context(oxstats_genotypes(*args)) for args in kwise(sys.argv[1:], 2)]
    merged = merge_oxstats([d for _, _, _, d in data])
    for row in merged:
        print(*row)
