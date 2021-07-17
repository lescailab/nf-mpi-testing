#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 12:51:41 2019

@author: barnafi
"""


class PhysicalParameters:
    def __init__(self, **kwargs):
        keys = kwargs.keys()
        def setter(_x, default): return kwargs[_x] if _x in keys else default

        # Simluation specifications
        self.dim = setter("dim", 2)
        self.dt = kwargs["dt"]
        self.t0 = kwargs['t0']  # Initial time
        self.sim_time = kwargs["sim_time"]
        self.Niter = (self.sim_time - self.t0)/self.dt
        self.ys_degree = setter("ys_degree", 2)
        self.uf_degree = setter("uf_degree", 2)
        self.us_degree = setter("us_degree", 1)
        self.p_degree = setter("p_degree", 1)
        self.approach = setter("approach", "Default")  # Default/Burtschell
        self.active_stress = setter("active_stress", False)

        # Physical parameters
        self.rhos0 = self.rhos = setter("rhos", 2e3)  # Solid density
        self.rhof = setter("rhof", 2e3)  # Fluid density
        self.phi0 = self.phi = setter("phi", lambda t: 0.1)  # porosity
        self.mu_f = setter("mu_f", 0.035)  # Fluid viscosity
        self.ks = setter("ks", 2e8)  # bulk modulus
        if "E" in keys and "nu" in keys:
            self.E = setter("E", 1e6)  # Young's modulus
            self.nu = setter("nu", 0.4)  # Poisson ratio
            self.lmbda = self.E*self.nu/((1+self.nu)*(1-2*self.nu))  # Lamé params
            self.mu_s = self.E/(2*(1+self.nu))  # Lamé params
        elif "lmbda" in keys and "mu_s" in keys:
            self.lmbda = setter("lmbda", 10)  # Lamé params
            self.mu_s = setter("mu_s", 5)  # Lamé params
            self.E = self.mu_s*(3*self.lmbda+2*self.mu_s)/(self.lmbda + self.mu_s)
            self.nu = self.lmbda/(2*(self.lmbda + self.mu_s))
        self.gamma = setter("gamma", 20)  # Nitsche parameter for no slip
        self.D = setter("D", 1e7)  # Permeability, i.e k_f^{-1}, constant for now

        self.AS = setter("AS", 1e3)
        self.beta = setter("beta", 10)


class PhysicalParametersSolid:
    def __init__(self, **kwargs):
        keys = kwargs.keys()
        def setter(_x, default): return kwargs[_x] if _x in keys else default

        # Simluation specifications
        self.dt = kwargs["dt"]
        self.t0 = kwargs['t0']  # Initial time
        self.sim_time = kwargs["sim_time"]
        self.Niter = int((self.sim_time - self.t0)/self.dt)
        self.ys_degree = setter("ys_degree", 2)

        # Physical parameters
        self.ks = setter("ks", 2e8)  # bulk modulus
        self.rhos = self.rhos = setter("rhos", 1e3)  # Solid density
        if "E" in keys and "nu" in keys:
            self.E = setter("E", 1e6)  # Young's modulus
            self.nu = setter("nu", 0.4)  # Poisson ratio
            self.lmbda = self.E*self.nu/((1+self.nu)*(1-2*self.nu))  # Lamé params
            self.mu_s = self.E/(2*(1+self.nu))  # Lamé params
        else:
            self.lmbda = setter("lmbda", 10)  # Lamé params
            self.mu_s = setter("mu_s", 5)  # Lamé params
            self.E = self.mu_s*(3*self.lmbda+2*self.mu_s)/(self.lmbda + self.mu_s)
            self.nu = self.lmbda/(2*(self.lmbda + self.mu_s))

        # Activation parameter
        self.AS = setter("AS", 1e3)


class OutputParameters:
    def __init__(self, **kwargs):
        keys = kwargs.keys()
        def setter(_x, default): return kwargs[_x] if _x in keys else default

        # Solution storage
        self.name = setter("name", "result")
        # Keep lists of solutions for convergence
        self.storeSolutions = setter("store_solutions", False)
        # Export solutions to xdmf (only serial)
        self.exportSolutions = setter("export_solutions", True)
        # Amount of iterations before exporting, it%saveEvery == 0
        self.saveEvery = setter("save_every", 1)
        # Verbose output
        self.verbose = setter("verbose", True)


class ExternalForces:
    def __init__(self, **kwargs):
        self.ff = kwargs['ff']  # Load on fluid equation
        self.fs = kwargs['fs']  # Load on solid equation
        self.ts = kwargs['ts']  # Neumann solid traction
        self.tf = kwargs['tf']  # Neumann fluid traction
        self.theta = kwargs['theta']  # Mass source/sink
        self.activeStress = kwargs['active_stress'] if 'active_stress' in kwargs.keys() else None
