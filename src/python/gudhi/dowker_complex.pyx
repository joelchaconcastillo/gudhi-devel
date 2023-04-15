# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Joel Chacon Castillo
#
# Copyright (C) 2023 UTD
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libc.stdint cimport intptr_t

from gudhi.simplex_tree cimport *
from gudhi.simplex_tree import SimplexTree

__author__ = "Joel Chacon Castillo"
__copyright__ = "Copyright (C) 2023 UTD"
__license__ = "MIT"

cdef extern from "Dowker_complex_interface.h" namespace "Gudhi":
    cdef cppclass Dowker_complex_interface "Gudhi::dowker_complex::Dowker_complex_interface":
        Dowker_complex_interface(vector[vector[pair[size_t, double]]] nearest_landmark_table)
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, double epsilon) except +
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, double epsilon,
            unsigned limit_dimension) except +

# DowkerComplex python interface
cdef class DowkerComplex:
    """Constructs (weak) witness complex for a given table of nearest landmarks
    with respect to witnesses.
    """

    cdef Dowker_complex_interface * thisptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, nearest_landmark_table=None):
        """DowkerComplex constructor.

        :param nearest_landmark_table: A list of lists of nearest landmarks and their distances.
            `nearest_landmark_table[w][k]==(l,d)` means that l is the k-th nearest landmark to
            witness w, and d is the (squared) distance between l and w.
        :type nearest_landmark_table: list of list of pair of int and float
        """

    # The real cython constructor
    def __cinit__(self, nearest_landmark_table=None):
        if nearest_landmark_table is not None:
            self.thisptr = new Dowker_complex_interface(nearest_landmark_table)

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def __is_defined(self):
        """Returns true if DowkerComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def create_simplex_tree(self, epsilon = float('inf'), limit_dimension = -1):
        """
        :param epsilon: The maximum relaxation parameter.
            Default is set to infinity.
        :type epsilon: float
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        cdef intptr_t stree_int_ptr=stree.thisptr
        if limit_dimension != -1:
            self.thisptr.create_simplex_tree(<Simplex_tree_interface_full_featured*>stree_int_ptr,
                epsilon, limit_dimension)
        else:
            self.thisptr.create_simplex_tree(<Simplex_tree_interface_full_featured*>stree_int_ptr,
                epsilon)
        return stree
