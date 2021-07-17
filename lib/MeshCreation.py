#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 14:06:24 2019

@author: barnafi
"""


def generateSquare(Nelements, length):
    """
    Creates a square mesh of given elements and length with markers on
    the sides: left, bottom, right and top
    """
    from dolfin import UnitSquareMesh, SubDomain, MeshFunction, Measure, near
    mesh = UnitSquareMesh(Nelements, Nelements)
    # Rescale for Chapelle-Moireau comparison
    mesh.coordinates()[:] *= length

    # Subdomains: Solid
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], 0.0) and on_boundary

    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], length) and on_boundary

    class Top(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], length) and on_boundary

    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], 0.0) and on_boundary
    left, right, top, bottom = Left(), Right(), Top(), Bottom()
    LEFT, RIGHT, TOP, BOTTOM = 1, 2, 3, 4  # Set numbering
    NONE = 99  # Marker for empty boundary

    markers = MeshFunction("size_t", mesh, 1)
    markers.set_all(0)

    boundaries = (left, right, top, bottom)
    def_names = (LEFT, RIGHT, TOP, BOTTOM)
    for side, num in zip(boundaries, def_names):
        side.mark(markers, num)

    return mesh, markers, LEFT, RIGHT, TOP, BOTTOM, NONE


def generate_boundary_measure(mesh, markers, tags_list, none_tag=42):
    from dolfin import Measure

    ds = Measure('ds', domain=mesh, subdomain_data=markers,
                 metadata={'optimize': True})
    return sum([ds(i) for i in tags_list], ds(none_tag))


def prolateGeometry(filename="prolate_4mm"):

    from dolfin import XDMFFile, Mesh, MeshValueCollection, MeshTransformation
    xdmf_meshfile = "mesh/" + filename + ".xdmf"
    xdmf_meshfile_bm = "mesh/" + filename + "_bm.xdmf"
    mesh = Mesh()
    with XDMFFile(xdmf_meshfile) as infile:
        infile.read(mesh)
    mvc = MeshValueCollection("size_t", mesh, 2)
    with XDMFFile(xdmf_meshfile_bm) as infile:
        infile.read(mvc, "name_to_read")
    from dolfin import cpp
    markers = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    ENDOCARD = 20
    EPICARD = 10
    BASE = 50
    NONE = 99

    MeshTransformation.scale(mesh, 1e-3)
    return mesh, markers, ENDOCARD, EPICARD, BASE, NONE


def setMeshAndMarkers3D(Nelements, length):
    pass
