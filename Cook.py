from lib.Solvers.NonlinearSolvers import Newton, BFGS
from dolfin import *
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import sys
assert len(sys.argv) == 4, "Must give argument for number of refinements, FE degree and method (Newton|BFGS)"

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["mesh_partitioner"] = "ParMETIS"
method = sys.argv[3]  # Newton|BFGS

# Set options for KSP
if method == "Newton":
    ksp_type = "gmres"
    ksp_atol = 1e-8
    ksp_rtol = 1e-6
elif method in ("Newton-inexact", "BFGS-inexact"):
    ksp_type = "gmres"
    ksp_atol = 0.0
    ksp_rtol = 1e-2
elif method == "BFGS":
    ksp_type = "preonly"
    ksp_atol = 0.0
    ksp_rtol = 0.0
else:
    assert False, "method must be one of Newton|Newton-inexact|BFGS-inexact|BFGS"

PETSc.Options().setValue("-ksp_type", ksp_type)
PETSc.Options().setValue("-ksp_norm_type", "unpreconditioned")
PETSc.Options().setValue("-ksp_gmres_modifiedgramschmidt", None)
PETSc.Options().setValue("-ksp_gmres_restart", 1000)
PETSc.Options().setValue("-ksp_rtol", ksp_rtol)
PETSc.Options().setValue("-ksp_atol", ksp_atol)
PETSc.Options().setValue("-ksp_max_it", 200)

# HYPRE Options
# PETSc.Options().setValue("-pc_hypre_boomeramg_grid_sweeps_all", 1)
# PETSc.Options().setValue("-pc_hypre_boomeramg_P_max", 4)
# PETSc.Options().setValue("-pc_hypre_boomeramg_agg_nl", 1)
# PETSc.Options().setValue("-pc_hypre_boomeramg_agg_num_paths", 2)
# PETSc.Options().setValue("-pc_hypre_boomeramg_coarsen_type", "HMIS")
# PETSc.Options().setValue("-pc_hypre_boomeramg_interp_type", "ext+i")
# PETSc.Options().setValue("-pc_hypre_boomeramg_no_CF", "true")

# Solver params
method = sys.argv[3]  # Newton|BFGS
atol, rtol, maxit, alpha, verbose = 1e-8, 1e-6, 100, 1.0, True
# Set options for KSP
if method == "Newton":
    ksp_type = "gmres"
    ksp_atol = 1e-8
    ksp_rtol = 1e-6
elif method in ("Newton-inexact", "BFGS-inexact"):
    ksp_type = "gmres"
    ksp_atol = 0.0
    ksp_rtol = 1e-2
elif method == "BFGS":
    ksp_type = "preonly"
    ksp_atol = 0.0
    ksp_rtol = 0.0
else:
    assert False, "method must be one of Newton|Newton-inexact|BFGS-inexact|BFGS"

PETSc.Options().setValue("-ksp_type", ksp_type)
PETSc.Options().setValue("-ksp_norm_type", "unpreconditioned")
PETSc.Options().setValue("-ksp_gmres_modifiedgramschmidt", None)
PETSc.Options().setValue("-ksp_gmres_restart", 1000)
PETSc.Options().setValue("-ksp_rtol", ksp_rtol)
PETSc.Options().setValue("-ksp_atol", ksp_atol)
PETSc.Options().setValue("-ksp_max_it", 200)

anderson_depth = 0  # Used with BFGS

# Model params
N = int(sys.argv[1])
u_deg = int(sys.argv[2])
nx = N * 6
ny = N * 6
nz = N * 1
lx = 48e-3  # mm
ly = 44e-3  # mm
lz = 10e-3  # mm
w = 16e-3  # mm
tau = 1.995e6
mesh = BoxMesh(Point(0, 0, 0), Point(lx, ly, lz), nx, ny, nz)

# Convert box into Cook membrane
for i, coord in enumerate(mesh.coordinates()):
    x, y, z = coord
    mesh.coordinates()[i, 1] = (1 - (ly - w)/(ly*lx)*x)*y + ly * x / lx

V = VectorFunctionSpace(mesh, "Lagrange", u_deg)
if MPI.COMM_WORLD.rank == 0:
    print("Problem dofs = {}, dofs/core = {}".format(V.dim(), V.dim() / MPI.COMM_WORLD.size), flush=True)
    assert V.dim() / MPI.COMM_WORLD.size <= 250000, "We allow up to 250k dofs per core, sorry!"

# Mark boundary subdomians
c = Constant((0.0, 0.0, 0.0))

boundary_markers = MeshFunction("size_t", mesh, 2)
boundary_markers.set_all(0)


class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], lx)


class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)


bx0 = Right()
bx1 = Left()
RIGHT = 1
LEFT = 2
bx0.mark(boundary_markers, RIGHT)
bx1.mark(boundary_markers, LEFT)
bcl = DirichletBC(V, c, boundary_markers, LEFT)
bcs = [bcl]
ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)

FILE = File('images/cook.pvd')
# Problem constants
ap = 126  # 126kPa
bp = 252  # 252kPa
lp = 81512  # 81512kPa
# Neo-hookean params
#C1 = 0.4
#D1 = 2.5e-4

mu = Constant(80.194e6)
nu = Constant(0.3)
lmbd = 2*mu*nu / (1 - 2 * nu)
kk = lmbd + 0.666666 * mu
C1 = mu / 2

# Find coefficients of parabolic shear stress
tA = np.array([[ly**2, ly, 1], [(ly+w)**2, ly+w, 1], [(ly + 0.5 * w)**2, ly + 0.5 * w, 1]])
tb = np.array([0, 0, tau])
coeffs = np.linalg.solve(tA, tb)

# Define functions
du = TrialFunction(V)            # Incremental displacement
v = TestFunction(V)             # Test function
u = Function(V)                 # Displacement
u.rename("d", "d")
T = Expression(("0.0",  "tt", "0.0"), degree=1, tt=tau)  # Traction force on the boundary


# Kinematics
I = Identity(3)    # Identity tensor
F = I + grad(u)             # Deformation gradient
C, J = F.T*F, det(F)
H = J * inv(F).T

# Problem definition
psi = C1 * (J**(-2/3) * tr(C) - 3) + kk * (J**2 - 1 - 2 * ln(J))
Pi = psi*dx - dot(H * T, u)*ds(RIGHT)
Res = derivative(Pi, u, v)
Jac = derivative(Res, u, du)


def res_fun(res):
    assemble(Res, tensor=res)


def jac_fun(jac):
    assemble(Jac, tensor=jac)


jac = PETScMatrix()
jac.mat().setBlockSize(3)
assemble(Jac, tensor=jac)
for bb in bcs:
    bb.apply(jac)
if method in ("Newton", "Newton-inexact"):
    solver = Newton(jac.mat(), atol, rtol, maxit, alpha, verbose)
    status = solver.solve(u.vector(), res_fun, jac_fun, bcs, bs=3)
elif method in ("BFGS", "BFGS-inexact"):
    solver = BFGS(jac.mat(), atol, rtol, maxit, alpha, verbose,
                  anderson_depth, LM=True, LM_order=100)
    status = solver.solve(u.vector(), res_fun, bcs)
FILE << u
if MPI.COMM_WORLD.rank == 0 and not status:
    print("DIVERGED!")
