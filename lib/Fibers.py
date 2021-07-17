# -*- coding: utf-8 -*-
# Code by: Francesco Regazzoni

import dolfin as df


def getH1projection(f, V):

    mesh = V.mesh()
    uv = df.TrialFunction(V)
    vv = df.TestFunction(V)
    A = df.dot(uv, vv) + df.inner(df.grad(uv), df.grad(vv))
    L = df.dot(f, vv) + df.inner(df.grad(f), df.grad(vv))
    sol = df.Function(V)
    df.solve(A * df.dx(mesh) == L * df.dx(mesh), sol)
    return sol


def getTransmuralCoordinate(mesh, boundary_markers, ENDO, EPI, degree=2):

    df.set_log_level(30)
    dx = df.dx(mesh)

    V = df.FunctionSpace(mesh, 'CG', degree)
    bc1 = df.DirichletBC(V, df.Constant(0.0), boundary_markers, ENDO)
    bc2 = df.DirichletBC(V, df.Constant(1.0), boundary_markers, EPI)

    phi = df.Function(V, name="phi_transmural")
    phi_trial = df.TrialFunction(V)
    psi = df.TestFunction(V)

    df.solve(df.dot(df.grad(phi_trial), df.grad(psi)) * dx ==
             df.Constant(0.0) * psi * dx, phi, [bc1, bc2])

    return phi


def getApicobasalCoordinate(mesh, boundary_markers, BASE, degree=2):

    df.set_log_level(30)
    dx = df.dx(mesh)

    V = df.FunctionSpace(mesh, 'CG', degree)

    phi_trial = df.TrialFunction(V)
    psi = df.TestFunction(V)

    a = df.dot(df.grad(phi_trial), df.grad(psi)) * df.dx
    L = df.Constant(0.0) * psi * df.dx


#    quotes = np.matmul(self.mesh.coordinates(),self.centerline)  # centerline??
    quotes = mesh.coordinates()[:, 2]
    min_quote = df.MPI.max(df.MPI.comm_world, quotes.max())
#    min_quote = 0.060184572583718024

    def apex(x):
        result = abs(x[2] - min_quote) < df.DOLFIN_EPS
        return result
    bcs = [df.DirichletBC(V, df.Constant(1.0), boundary_markers, BASE),
           df.DirichletBC(V, df.Constant(0.0), apex, method='pointwise')]

    phi = df.Function(V, name="phi_apicobasal")
    df.solve(a == L, phi, bcs)
    return phi

# TODO: works in parallel?
def generateFibers(mesh, boundary_markers, ENDO, EPI, BASE, output_dir=None):

    if df.MPI.rank(df.MPI.comm_world) == 0:
        print("generating fibers...", flush=True)
    theta_endo = 60.
    theta_epi = -60.

    degree = 1
    phi_transmural = getTransmuralCoordinate(
        mesh, boundary_markers, ENDO, EPI, degree=degree)

    def Max(a, b): return (a + b + abs(a - b)) / df.Constant(2.)
    def Min(a, b): return (a + b - abs(a - b)) / df.Constant(2.)

    W = df. VectorFunctionSpace(mesh, 'CG', degree)
    n = df.Function(W)  # W?

    # -1: analytical fibers
    # 0: SR (on fields) NB: f and s are not perfectly orthogonal
    # 1: SR (on dofs)
    # 2: BT (on dofs)
    alg_type = 2  # 0 faster, but f and s are not perfectly orthogonal
    if alg_type == -1:

        #            apex_x = 0.0469341
        #            apex_y = 1.34562

        #            f = df.project(df.Expression(('-(x[1]-apex_y)/sqrt(pow(x[0]-apex_x,2)+pow(x[1]-apex_y,2))',
        #                                                ' (x[0]-apex_x)/sqrt(pow(x[0]-apex_x,2)+pow(x[1]-apex_y,2))',
        #                                                '0'), degree = 2,apex_x=apex_x,apex_y=apex_y),self.W)
        #            s = df.project(df.Expression((' (x[0]-apex_x)/sqrt(pow(x[0]-apex_x,2)+pow(x[1]-apex_y,2))',
        #                                                ' (x[1]-apex_y)/sqrt(pow(x[0]-apex_x,2)+pow(x[1]-apex_y,2))',
        #                                                '0'), degree = 2,apex_x=apex_x,apex_y=apex_y),self.W)
        f = df.project(df.Expression(('1.0', '0.0', '0.0'), degree=degree), W)
        s = df.project(df.Expression(('0.0', '1.0', '0.0'), degree=degree), W)
        n = df.project(df.cross(f, s), W)

    elif alg_type == 0:
        s = df.grad(phi_transmural)
#            s = s / Max(1e-10, df.sqrt(df.inner(s,s)))
        s = s / df.sqrt(df.inner(s, s))
        k = df.Constant((.0, .0, -1.))  # TODO: move to option file
        kp_tilde = k - df.dot(k, s) * s
        kp = kp_tilde / df.sqrt(df.inner(kp_tilde, kp_tilde))
        f_tilde = df.cross(s, kp)
        f_tilde = f_tilde / df.sqrt(df.inner(f_tilde, f_tilde))
        theta = (theta_endo + (theta_epi - theta_endo)
                 * phi_transmural) * df.pi / 180.0
        f = f_tilde + df.sin(theta) * df.cross(s, f_tilde) + 2.2 * \
            (df.sin(theta * .5))**2 * df.cross(s, df.cross(s, f_tilde))
        f = - f / df.sqrt(df.inner(f, f))
#            n = df.cross(f,s)
#            n = n / df.sqrt(df.inner(n,n))
        s = df.project(s, W)
        f = df.project(f, W)
#            n = df.project(df.cross(f,s),self.W)

#            sx,sy,sz = df.split(df.project(s,self.W))
#            sx = sx.vector().get_local()
#            sy = sy.vector().get_local()
#            sz = sz.vector().get_local()
#
#            sx,sy,sz = df.split(df.project(s,self.W))
#            sx = df.project(sx,self.V).vector().get_local()
#            sy = df.project(sy,self.V).vector().get_local()
#            sz = df.project(sz,self.V).vector().get_local()
#            s = np.concatenate((sx[:,None],sy[:,None],sz[:,None]),axis=1)

    elif alg_type == 1 or alg_type == 2:

        import numpy as np
        ndof_local = phi_transmural.vector().get_local().size
        s_vec = np.empty((ndof_local, 3))
        s_vec_tot = df.project(df.grad(phi_transmural), W).vector().get_local()
#            s_vec[:,0] = s_vec_tot[self.W.sub(0).dofmap().dofs()]
#            s_vec[:,1] = s_vec_tot[self.W.sub(1).dofmap().dofs()]
#            s_vec[:,2] = s_vec_tot[self.W.sub(2).dofmap().dofs()]
        s_vec[:, 0] = s_vec_tot[0::3]
        s_vec[:, 1] = s_vec_tot[1::3]
        s_vec[:, 2] = s_vec_tot[2::3]
#        print("differenza %s" % np.linalg.norm(s-s_vec))
        s = s_vec

        if alg_type == 2:
            phi_apicobasal = getApicobasalCoordinate(
                mesh, boundary_markers, BASE, degree=1)
            k_vec = np.empty((ndof_local, 3))
            k_vec_tot = df.project(
                df.grad(phi_apicobasal), W).vector().get_local()
#                k_vec[:,0] = k_vec_tot[self.W.sub(0).dofmap().dofs()]
#                k_vec[:,1] = k_vec_tot[self.W.sub(1).dofmap().dofs()]
#                k_vec[:,2] = k_vec_tot[self.W.sub(2).dofmap().dofs()]
            k_vec[:, 0] = k_vec_tot[0::3]
            k_vec[:, 1] = k_vec_tot[1::3]
            k_vec[:, 2] = k_vec_tot[2::3]
    #        print("differenza %s" % np.linalg.norm(s-s_vec))
            k = k_vec
        else:
            pass
# TODO           k = np.tile(centerline[None,:],[ndof_local,1])

#            if self.meshtype == "coarse":
#                v2d = df.vertex_to_dof_map(self.V)
#                for i in range(self.ndof_local):
#                    if self.mesh.coordinates()[i,2] < 18.2e-3:
#                        s_vec[v2d[i],:] = np.array([self.mesh.coordinates()[i,0],self.mesh.coordinates()[i,1],0.])

        phi_vec = phi_transmural.vector().get_local()

        s_norm = np.empty(ndof_local)
        theta = np.empty(ndof_local)
        kp_tilde = np.empty((ndof_local, 3))
        kp = np.empty((ndof_local, 3))
        f_tilde = np.empty((ndof_local, 3))
        f = np.empty((ndof_local, 3))
        n = np.empty((ndof_local, 3))
#            k = self.centerline
        for i in range(ndof_local):
            s_norm[i] = np.sqrt(np.inner(s[i, :], s[i, :]))
            s[i, :] = s[i, :] / s_norm[i]
            kp_tilde[i, :] = k[i, :] - np.inner(k[i, :], s[i, :]) * s[i, :]
            kp[i, :] = kp_tilde[i, :] / \
                np.sqrt(np.inner(kp_tilde[i, :], kp_tilde[i, :]))
            f_tilde[i, :] = np.cross(s[i, :], kp[i, :])
            f_tilde[i, :] = f_tilde[i, :] / \
                np.sqrt(np.inner(f_tilde[i, :], f_tilde[i, :]))
            theta[i] = (theta_endo + (theta_epi - theta_endo)
                        * phi_vec[i]) * np.pi / 180.0
            f[i, :] = f_tilde[i, :] + np.sin(theta[i]) * np.cross(s[i, :], f_tilde[i, :]) + \
                2.2 * (np.sin(theta[i] * .5))**2 * \
                np.cross(s[i, :], np.cross(s[i, :], f_tilde[i, :]))
            f[i, :] = - f[i, :] / np.sqrt(np.inner(f[i, :], f[i, :]))
            n[i, :] = np.cross(f[i, :], s[i, :])

#                print('f.s = %e \t f.n = %e \t s.n = % e' % (
#                        np.inner(f[i,:],s[i,:]),
#                        np.inner(f[i,:],n[i,:]),
#                        np.inner(s[i,:],n[i,:])))

        f_vec = np.empty(ndof_local * 3)
#            f_vec[self.W.sub(0).dofmap().dofs()] = f[:,0]
#            f_vec[self.W.sub(1).dofmap().dofs()] = f[:,1]
#            f_vec[self.W.sub(2).dofmap().dofs()] = f[:,2]
        for i in range(3):
            f_vec[i::3] = f[:, i]

        s_vec = np.empty(ndof_local * 3)
#            s_vec[self.W.sub(0).dofmap().dofs()] = s[:,0]
#            s_vec[self.W.sub(1).dofmap().dofs()] = s[:,1]
#            s_vec[self.W.sub(2).dofmap().dofs()] = s[:,2]
        for i in range(3):
            s_vec[i::3] = s[:, i]

        n_vec = np.empty(ndof_local * 3)
#            n_vec[self.W.sub(0).dofmap().dofs()] = n[:,0]
#            n_vec[self.W.sub(1).dofmap().dofs()] = n[:,1]
#            n_vec[self.W.sub(2).dofmap().dofs()] = n[:,2]
        for i in range(3):
            n_vec[i::3] = n[:, i]

        f = df.Function(W)
        s = df.Function(W)
        n = df.Function(W)
        f.vector().set_local(f_vec)
        f.vector().apply("insert")
        s.vector().set_local(s_vec)
        s.vector().apply("insert")
        n.vector().set_local(n_vec)
        n.vector().apply("insert")


#        f.vector().set_local()
    if df.MPI.rank(df.MPI.comm_world) == 0:
        print("fibers generated!", flush=True)

    f.rename("f", "f")
    s.rename("s", "s")
    n.rename("n", "n")
#        out_file = df.File("geometry_prolate/fibers_and_sheets.pvd")
#        out_file << (phi,0)
#        out_file << (f,0)
#        out_file << (s,0)
    if output_dir:
        xdmf = df.XDMFFile("{}/fibers.xdmf".format(output_dir))
        xdmf.parameters["functions_share_mesh"] = True
        xdmf.parameters["flush_output"] = True
        xdmf.write(phi_transmural, 0)
        if alg_type == 2:
            xdmf.write(phi_apicobasal, 0)
        xdmf.write(f, 0)
        xdmf.write(s, 0)
        xdmf.write(n, 0)
        xdmf.close()

    return f, s, n
