import math

import numpy as np
from src.Force import Force
from src.Constraint import Constraint


class Node(object):
    """Node object"""

    def __init__(self, x, y, z, index):
        self.coordiantes = np.array([x, y, z], dtype=float)
        self.force = Force(0, 0, 0)
        self.constraint = Constraint(True, True, True)
        self.dof_numbers = np.array([math.inf , math.inf, math.inf]) # random values
        self.displacement = np.array([0,0,0])
        self.index = index

    def get_node_index(self):
        """
        :return: int
            Index of the node object
        """
        return self.index

    def set_constraint(self, constraint):
        """
        :param constraint: constraint object
            Set the constraint of the node
        """
        self.constraint = constraint

    def get_constraint(self):
        """
        :return: constraint object
            Constraint of the node
        """
        return self.constraint

    def set_force(self, f):
        """
        :param f: force object
            Set the force of the node
        """
        self.force = f

    def get_force(self):
        """
        :return: force object
            Get force of the node
        """
        return self.force

    def get_coordinates(self):
        """
        :return: np.array
            Coordinates of the node as 1D-Array
        """
        return self.coordiantes

    def get_dof_numbers(self):
        """
        :return: np.array
            1D-Numpy array with indices of the enumerated DOF's
        """
        return self.dof_numbers

    def set_displacement(self, array):
        """
        :param array: np.array
            Set the displacement of the node after solving
        """
        self.displacement = array

    def get_displacement(self):
        """
        :return: np.array
            1D-Np-array with the displaced coordinates
        """
        return self.displacement

    def get_displacement_at(self, x):
        """
        :param x: int
            Index of x,y or z e.g. 0, 1 or 2
        :return: float
            Displacement at index
        """
        return self.displacement[x]

    def info(self):
        """ Print some information of the node"""
        print("Coordinates are:")
        print("[%s, %s, %s] " % (self.coordiantes[0], self.coordiantes[1], self.coordiantes[2]))

    def enumerate_dofs_node(self, eqn_count):
        """
        :param eqn_count: int
            Integer index that is passed to all nodes to enumerate their DOF's
        :return: int
            Enumerated index
        """
        for i in range(0, 3):
            if self.constraint.is_free(i) == True:
                self.dof_numbers[i] = eqn_count
                eqn_count += 1
            else:
                self.dof_numbers[i] = -1
        return eqn_count