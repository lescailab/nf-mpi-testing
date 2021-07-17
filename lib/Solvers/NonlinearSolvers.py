import numpy as np
import dolfin as df
from mpi4py import MPI
from petsc4py import PETSc
from time import perf_counter as time
from lib.AndersonAcceleration import AndersonAcceleration


class Newton:

    def __init__(self, jac0, atol, rtol, maxit, alpha, verbose=False):
        self.atol = atol
        self.rtol = rtol
        self.maxit = maxit
        self.alpha = alpha
        self.verbose = verbose
        ksp = PETSc.KSP().create()
        ksp.setType('gmres')
        ksp.setOperators(jac0)
        pc = ksp.getPC()
        pc.setType('hypre')
        pc.setFromOptions()
        ksp.setTolerances(1e-8, 1e-10, 1e20, 1000)
        ksp.setFromOptions()
        self.ksp = ksp
        self.do_fieldsplit = False
        self.sub_allocated = False

    def apply_bcs(self, obj, bcs):
        if bcs:
            for bb in bcs:
                bb.apply(obj)

    def set_mechanics_fieldsplit_pc(self, V0, V1):
        dofmap_s = V0.dofmap().dofs()
        dofmap_p = V1.dofmap().dofs()
        is_s = PETSc.IS().createGeneral(dofmap_s)
        is_p = PETSc.IS().createGeneral(dofmap_p)

        pc = self.ksp.getPC()
        pc.setType('fieldsplit')
        pc.setFromOptions()
        pc.setFieldSplitIS((None, is_s), (None, is_p))
        pc.setUp()  # Must be called after set from options
        ksps = pc.getFieldSplitSubKSP()
        for ksp in ksps:
            ksp.setFromOptions()

    def solve(self, sol, res_fun, jac_fun, bcs=None, bs=None):

        if bcs:
            for bb in bcs:
                bb.apply(sol)
                bb.homogenize()
        t0 = time()
        err_abs = 1.0
        err_rel = 1.0
        it = 0
        du = sol.copy()
        du.zero()
        res = df.PETScVector()
        jac = df.PETScMatrix()

        res_fun(res)
        self.apply_bcs(res, bcs)
        jac_fun(jac)
        self.apply_bcs(jac, bcs)

        res0 = res.norm('l2')
        total_krylov = 0
        while err_abs > self.atol and err_rel > self.rtol and it < self.maxit:
            it += 1

            self.ksp.setOperators(jac.mat())
            self.ksp.setUp()
            self.ksp.solve(res.vec(), du.vec())
            krylov_its = self.ksp.getIterationNumber()
            # krylov_its = df.solve(jac, du, res, 'gmres', 'hypre_amg')
            # solve(A, du, b)
            sol.vec().axpy(-1, du.vec())  # we solve jac du = - res
            sol.apply("")
            res_fun(res)
            self.apply_bcs(res, bcs)
            jac_fun(jac)
            if bs:
                jac.mat().setBlockSize(bs)
            self.apply_bcs(jac, bcs)
            err_abs = res.norm('l2')
            err_rel = err_abs/res0
            if self.verbose and MPI.COMM_WORLD.rank == 0:
                print('it {}, err abs = {:.3e}  err rel = {:.3e} in {} GMRES iterations'.format(
                    it, err_abs, err_rel, krylov_its), flush=True)
            total_krylov += krylov_its
            if min(err_abs, err_rel) > 1e14 or np.isnan(err_abs) or np.isnan(err_rel):
                if MPI.COMM_WORLD.rank == 0:
                    print("\t Newton diverged")
                return False
        if(MPI.COMM_WORLD.rank == 0 and it < self.maxit):
            print("\t Newton converged in {} nonlinear its, {} total krylov its, {} avg krylov its, {:.3f}s".format(
                it, total_krylov, total_krylov/it, time() - t0), flush=True)
        if it == self.maxit:
            return False
        else:
            return True


class BFGS:
    def __init__(self, jac, atol, rtol, maxit, alpha, verbose=False, anderson_depth=0, LM=False, LM_order=10):
        self.atol = atol
        self.rtol = rtol
        self.maxit = maxit
        self.alpha = alpha
        self.verbose = verbose
        self.anderson_depth = anderson_depth
        self.LM = LM  # Low memory
        self.LM_order = LM_order  # Number of iterations to be used with LM
        ksp = PETSc.KSP().create()
        ksp.setType('preonly')
        ksp.setOperators(jac)
        pc = ksp.getPC()
        pc.setType('hypre')
        ksp.setTolerances(1e-6, 1e-8, 1e20, 1000)
        ksp.setFromOptions()
        pc.setFromOptions()
        self.ksp = ksp
        self.b = None
        self.sks = []
        self.yks = []
        self.rhoks = []

        def H0(vec, k):
            assert k == 0, "Preconditioner in BFGS must me called at last level in recursion always"
            # return 0.1*vec
            if not self.b:
                self.b = vec.copy()
            vec.copy(self.b)
            self.ksp.solve(self.b, vec)
            # return vec
        self.Hs = [H0]
        self.k = 0

    def set_mechanics_fieldsplit_pc(self, V0, V1):
        dofmap_s = V0.dofmap().dofs()
        dofmap_p = V1.dofmap().dofs()
        is_s = PETSc.IS().createGeneral(dofmap_s)
        is_p = PETSc.IS().createGeneral(dofmap_p)

        pc = self.ksp.getPC()
        pc.setType('fieldsplit')
        pc.setFromOptions()
        pc.setFieldSplitIS((None, is_s))
        pc.setFieldSplitIS((None, is_p))
        pc.setUp()  # Must be called after set from options
        ksps = pc.getFieldSplitSubKSP()
        for ksp in ksps:
            ksp.setFromOptions()

    def update(self, sk, yk):
        self.sks.append(sk)
        self.yks.append(yk)
        self.rhoks.append(1 / (yk.dot(sk)))

        if self.LM:
            # Truncate iterations
            if len(self.sks) > self.LM_order:
                self.sks.pop(0)
                self.yks.pop(0)
                self.rhoks.pop(0)
            assert len(self.sks) <= self.LM_order, "More previous iterations than allowed."
        else:
            # Add recursive functions
            def Hk(vec, k):
                rhok = self.rhoks[k-1]
                skdvec = self.sks[k-1].dot(vec)
                vec.axpy(-rhok*skdvec, self.yks[k-1])
                self.Hs[k-1](vec, k-1)
                vec.axpy(-rhok * self.yks[k-1].dot(vec), self.sks[k-1])
                vec.axpy(rhok * skdvec, self.sks[k-1])
            self.Hs.append(Hk)

    def apply(self, vec):
        if self.LM:
            iters = len(self.sks)  # Number of iterations
            self.alphas = np.zeros(iters)
            self.betas = np.zeros(iters)

            # Use two-loop recursion
            for ii in range(iters):
                jj = iters - ii - 1  # invert direction
                self.alphas[jj] = self.rhoks[jj] * vec.dot(self.sks[jj])
                vec.axpy(-self.alphas[jj], self.yks[jj])
            self.Hs[0](vec, 0)
            for ii in range(iters):
                self.betas[ii] = self.rhoks[ii] * self.yks[ii].dot(vec)
                vec.axpy(self.alphas[ii] - self.betas[ii], self.sks[ii])
        else:
            return self.Hs[-1](vec, len(self.Hs)-1)

    def apply_bcs(self, obj, bcs):
        if bcs:
            for bb in bcs:
                bb.apply(obj)

    def apply_bc_petsc(self, vec, bcs):
        if bcs:
            for bb in bcs:
                vals = bb.get_boundary_values()
                for dof in vals:
                    vec.setValue(dof, vals[dof])

    def solve(self, sol, res_fun, bcs=None):
        t0 = time()
        # If there are BCs, homogenize them after initializing solution
        if bcs:
            for bb in bcs:
                bb.apply(sol)
                bb.homogenize()
        soln = sol.vec().copy()  # resn is PETSc vector

        res = df.PETScVector()
        res_fun(res)  # Initial residual
        self.apply_bcs(res, bcs)
        resn = res.vec().copy()
        anderson = AndersonAcceleration(self.anderson_depth)
        error0 = res.norm('l2')
        error_rel = error_abs = 1
        it = 0
        total_krylov = 0
        while error_abs > self.atol and error_rel > self.rtol and it < self.maxit:
            self.apply(res.vec())
            krylov_its = self.ksp.getIterationNumber()

            do_LS = False  # TODO: Line search
            if do_LS:
                pass
            sol.vec().axpy(-self.alpha, res.vec())
            anderson.get_next_vector(sol.vec())
            sol.apply("")
            res_fun(res)
            self.apply_bcs(res, bcs)

            # They must be allocated as the iteration count is unkown
            self.update(sol.vec() - soln, res.vec() - resn)  # sk, yk

            sol.vec().copy(soln), res.vec().copy(resn)
            error_abs = res.norm('l2')
            error_rel = error_abs / error0
            if(MPI.COMM_WORLD.rank == 0 and self.verbose):
                print("it {}, err_abs={:.3e}, err_rel={:.3e} in {} GMRES iterations".format(
                    it, error_abs, error_rel, krylov_its), flush=True)
            it += 1
            total_krylov += krylov_its
            if min(error_abs, error_rel) > 1e14 or np.isnan(error_abs) or np.isnan(error_rel):
                if MPI.COMM_WORLD.rank == 0:
                    print("\t BFGS diverged")
                return False
        if(MPI.COMM_WORLD.rank == 0 and it < self.maxit):
            print("\t BFGS converged in {} nonlinear its, {} total krylov its, {} avg krylov its, {:.3f}s".format(
                it, total_krylov, total_krylov/it, time() - t0), flush=True)
        if it == self.maxit:
            return False
        else:
            return True
