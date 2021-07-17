
# import dolfin once and for all
from dolfin import *

# From Rick Teachy @ StackOverflow
def mpiprint(*args):
    if MPI.rank(MPI.comm_self) == 0:
        print(*args, flush=True)
