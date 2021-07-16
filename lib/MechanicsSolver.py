# -*- coding: utf-8 -*-

from lib.GLOBAL_VARIABLES import *
from lib.Solvers.NonlinearSolvers import BFGS, Newton


class MechanicsSolver:
    """
    """

    def __init__(self, name):
        self.name = name
        self.ys_n = self.ys_nn = None
        self.bcs = None

    def setup(self, mesh, physicalParametersSolid, Markers, OutputParameters):
        self.physicalParameters = physicalParametersSolid
        self.outputParameters = OutputParameters
        self.markers = Markers

        self.mesh = mesh
        self.V = VectorFunctionSpace(mesh, 'CG', self.physicalParameters.ys_degree)

        # Initialize solutions
        self.ys_n = Function(self.V, name="displacement")
        self.ys_nn = self.ys_n.copy(True)

        if OutputParameters.storeSolutions:
            self.ysList = [self.ys_n.copy(True)]
            self.usList = [self.us_n.copy(True)]
            self.ufList = [self.uf_n.copy(True)]
            self.pList = [self.p_n.copy(True)]
        if OutputParameters.exportSolutions:
            self.xdmf = XDMFFile("images/{}/{}.xdmf".format(self.name, OutputParameters.name))
            self.exportSolution(0)

        # Setup measures
        self.dx = Measure('dx', domain=mesh, metadata={'optimize': True})
        self.dx = self.dx(degree=5)
        self.ds = Measure('ds', domain=mesh, subdomain_data=Markers.markers,
                          metadata={'optimize': True})
        dsEmpty = ds(Markers.NONE)
        self.dSN = sum([self.ds(i) for i in Markers.neumannMarkers], dsEmpty)
        self.dSRob = sum([self.ds(i) for i in Markers.robinMarkers], dsEmpty)

    def setBC(self, *args):
        self.bcs = args  # tuple by construction

    def solveTimeStep(self, t, method="Newton", anderson_depth=0, LBFGS_order=20):

        Sol = Function(self.V)
        Sol.assign(self.ys_n)
        F = self.generateForm(Sol, t)
        Jac = derivative(F, Sol)
        jac = PETScMatrix()
        assemble(Jac, tensor=jac)
        jac.mat().setBlockSize(3)
        if method == "Newton":
            def res_fun(res):
                assemble(F, tensor=res)

            def jac_fun(jac):
                assemble(Jac, tensor=jac)

            atol = 1e-8
            rtol = 1e-6
            maxit = 100
            alpha = 1.0
            solver = Newton(jac.mat(), atol, rtol, maxit, alpha, verbose=self.outputParameters.verbose)
            solver.solve(Sol.vector(), res_fun, jac_fun, bcs=None, bs=3)

        else:  # method == "BFGS"

            atol = 1e-8
            rtol = 1e-6
            maxit = 100
            alpha = 1.0
            bfgs = BFGS(jac.mat(), atol, rtol, maxit, alpha,
                        verbose=self.outputParameters.verbose, anderson_depth=anderson_depth, LM=True, LM_order=LBFGS_order)

            def res_fun(res):
                assemble(F, tensor=res)
            bfgs.solve(Sol.vector(), res_fun)
        self.ys_nn.assign(self.ys_n)
        self.ys_n.assign(Sol)

    def solve(self, printEvery=1):

        mpiprint("----- Solving problem")
        Niter = int(self.physicalParameters.Niter)

        from time import time
        for n in range(Niter):
            t = self.physicalParameters.t0+self.physicalParameters.dt*n
            current_t = time()
            mpiprint("----- Solving t={:.5f}".format(t))
            self.solveTimeStep(t)
            if n % self.outputParameters.saveEvery == 0 and self.outputParameters.exportSolutions:
                self.exportSolution(t+self.physicalParameters.dt)
            if n % printEvery == 0:
                mpiprint("----- Solved t={:.5f} in {:.2f}s".format(t, time() - current_t))
        if self.outputParameters.exportSolutions:
            self.xdmf.close()

    def generateForm(self, Sol, t):

        ys = Sol
        ws = TestFunction(self.V)
        dt = Constant(self.physicalParameters.dt)
        ys_n = self.ys_n
        ys_nn = self.ys_nn

        # Auxiliary variables
        idt = 1/dt
        us = idt * (ys - ys_n)
        us_n = idt * (ys_n - ys_nn)

        F = Identity(3) + grad(ys)
        F = variable(F)
        J = det(F)

        # Usyk,. mc Culloch 2002
        Cg = .88e3   # [Pa]
        bf = 8       # [-]
        bs = 6       # [-]
        bn = 3       # [-]
        bfs = 12      # [-]
        bfn = 3       # [-]
        bsn = 3       # [-]
        k = 5e4

        E = 0.5*(F.T*F - Identity(3))
        f0, s0, n0 = self.getFibers()
        Eff, Efs, Efn = inner(E*f0, f0), inner(E*f0, s0), inner(E*f0, n0)
        Esf, Ess, Esn = inner(E*s0, f0), inner(E*s0, s0), inner(E*s0, n0)
        Enf, Ens, Enn = inner(E*n0, f0), inner(E*n0, s0), inner(E*n0, n0)

        Q = Constant(bf) * Eff**2 \
            + Constant(bs) * Ess**2 \
            + Constant(bn) * Enn**2 \
            + Constant(bfs) * 2.0 * Efs**2 \
            + Constant(bfn) * 2.0 * Efn**2 \
            + Constant(bsn) * 2.0 * Esn**2
        WP = 0.5*Constant(Cg)*(exp(Q)-1)
        WV = Constant(k)/2*(J-1)*ln(J)

        W = WP + WV

        P = diff(W, F)

        # Setup forms

        # Pfaller et al.
        k_perp = Constant(2e5)  # [Pa/m]
        c_perp = Constant(5e3)  # [Pa*s/m]

        n = FacetNormal(self.mesh)
        ts = Constant((0, 0, 0))
        rhos = Constant(self.physicalParameters.rhos)
        ts_robin = -outer(n, n)*(k_perp*ys + c_perp*us) - (Identity(3) - outer(n, n))*k_perp/10*ys
        Fs = (inner(P, grad(ws)) + rhos/dt*dot(us - us_n, ws))*self.dx \
            - dot(ts, ws)*self.dSN \
            - dot(ts_robin, ws)*self.dSRob

        return Fs + inner(self.getActive(F, t), grad(ws))*self.dx

    def getActive(self, F, t):
        f0 = self.f0
        C = F.T*F
        I4f = dot(f0, C*f0)
        T_wave = 0.6
        def time_stim(_t): return sin(2*DOLFIN_PI*t)**2
        def time_stim(_t): return sin(2*DOLFIN_PI*t/T_wave) if t <= T_wave/2 else 0
        Ta = Constant(self.physicalParameters.AS*time_stim(t))
        Pa = Ta*outer(F*f0, f0)/sqrt(I4f)
        return Pa

    def setFibers(self, f0, s0, n0):
        self.f0 = f0
        self.n0 = n0
        self.s0 = s0

    def getFibers(self):
        return self.f0, self.s0, self.n0

    def exportSolution(self, t):
        self.xdmf.write(self.ys_n, t)

    def V(self):
        return self.V
