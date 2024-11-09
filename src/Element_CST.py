import math

import numpy as np
from src.Constraint import Constraint

class Element_CST(object):

    def __init__(self, e_mod, poisson, height, density, node_1, node_2, node_3):
        # Initialise passed variables
        self.e_mod = e_mod
        self.poisson = poisson
        self.height = height
        self.density = density
        self.node_1 = node_1
        self.node_2 = node_2
        self.node_3 = node_3
        # Initialise standard variables
        self.area = self.get_area()
        self.node_1.set_constraint(Constraint(True, True, False))
        self.node_2.set_constraint(Constraint(True, True, False))
        self.node_3.set_constraint(Constraint(True, True, False))
        self.mode = None
        self.dof_numbers = np.empty(9, dtype=int)

    def get_density(self):
        return self.density

    def get_e_mod(self):
        return self.e_mod

    def get_area(self): # Linear FEM Notes p.150
        width = (self.node_3.get_coordinates()) - (self.node_1.get_coordinates())
        height = (self.node_2.get_coordinates()) - (self.node_1.get_coordinates())
        area = (np.linalg.norm(np.cross(width, height))) / 2
        return area

    def get_node_1(self):
        return self.node_1

    def get_node_2(self):
        return self.node_2

    def get_node_3(self):
        return self.node_3

    def get_all_nodes(self):
        return [self.node_1, self.node_2, self.node_3]
    
    def compute_min_length(self):
        """Used for explicit time incrementation"""

        d12 = self.distance(self.node_1.get_coordinates(), self.node_2.get_coordinates())
        d13 = self.distance(self.node_1.get_coordinates(), self.node_3.get_coordinates())
        d23 = self.distance(self.node_2.get_coordinates(), self.node_3.get_coordinates())

        return min(d12, d13, d23)
    
    def distance(self, node1, node2):
        """Calculates eucledian distance between two nodes"""
        return math.sqrt((node2[0] - node1[0])**2 + 
                        (node2[1] - node1[1])**2 + 
                        (node2[2] - node1[2])**2)

    def get_number_of_nodes(self):
        return 3

    def set_mode(self, cst_mode_as_string):
        self.mode = cst_mode_as_string

    def get_mode(self):
        return self.mode

    def enumerate_dofs(self):
        for i in range(0, 3):
            self.dof_numbers[i] = self.node_1.get_dof_numbers()[i]
            self.dof_numbers[i+3] = self.node_2.get_dof_numbers()[i]
            self.dof_numbers[i+6] = self.node_3.get_dof_numbers()[i]

    def get_dof_numbers(self):
        return self.dof_numbers

    def b_operator(self): # Linear FEM Notes p.150
        factor_1 = (self.node_2.get_coordinates()[1]) - (self.node_3.get_coordinates()[1])
        factor_2 = (self.node_3.get_coordinates()[1]) - (self.node_1.get_coordinates()[1])
        factor_3 = (self.node_1.get_coordinates()[1]) - (self.node_2.get_coordinates()[1])
        factor_4 = (self.node_3.get_coordinates()[0]) - (self.node_2.get_coordinates()[0])
        factor_5 = (self.node_1.get_coordinates()[0]) - (self.node_3.get_coordinates()[0])
        factor_6 = (self.node_2.get_coordinates()[0]) - (self.node_1.get_coordinates()[0])

        row_1 = np.array([factor_1, 0, 0, factor_2, 0, 0, factor_3, 0, 0])
        row_2 = np.array([0, factor_4, 0, 0, factor_5, 0, 0, factor_6, 0])
        row_3 = np.array([factor_4, factor_1, 0, factor_5, factor_2, 0, factor_6, factor_3, 0])

        b_operator = np.array([row_1, row_2, row_3])*((1.0)/(2.0*self.area))
        # Return 3x9 Matrix
        return b_operator

    def compute_material_matrix(self): # Linear FEM Notes p.123
        if self.mode == "Plane Stress":
            row_1 = np.array([1, self.poisson, 0])
            row_2 = np.array([self.poisson, 1, 0])
            row_3 = np.array([0, 0, (1 - self.poisson) / 2])
            c_matrix = np.array([row_1, row_2, row_3])
            factor_1 = (self.e_mod) / (1 - np.power(self.poisson, 2))
            c_matrix_return = factor_1 * c_matrix
            return c_matrix_return
        elif self.mode == "Plane Strain":
            row_1 = np.array([1 - self.poisson, self.poisson, 0])
            row_2 = np.array([self.poisson, 1 - self.poisson, 0])
            row_3 = np.array([0, 0, (1 - 2*self.poisson) / 2])
            c_matrix = np.array([row_1, row_2, row_3])
            factor_1 = (self.e_mod) / ((1 + self.poisson) * (1 - 2*self.poisson))
            c_matrix = c_matrix * factor_1
            return c_matrix
        else:
            raise Exception("CST Element Mode must be Plane Stress or Plane Strain")

    def compute_element_k(self):
        b_operator = self.b_operator()
        b_operator_transposed = np.transpose(b_operator)
        c_matrix = self.compute_material_matrix()
        element_k = (b_operator_transposed @ c_matrix @ b_operator)*self.area*self.height
        return element_k

    def compute_element_stress(self):
        # Get displacements
        displacements = np.empty([9])
        for i in range(3):
            displacements[i] = self.node_1.get_displacement()[i]
            displacements[i+3] = self.node_2.get_displacement()[i]
            displacements[i+6] = self.node_3.get_displacement()[i]
        # Get b operator to calculate strains
        b_operator = self.b_operator()
        # np.abs() or otherwise visualisation makes no sense. (Maybe not correct?)
        strains = b_operator @ displacements
        stresses = self.compute_material_matrix() @ strains
        return stresses

    def compute_element_mean_stress(self):
        stress = self.compute_element_stress()
        xx = stress[0]
        yy = stress[1]
        xy = stress[2]
        # Calculate mean stress
        if (self.mode == "Plane Stress"):
            mean_stress = math.sqrt(math.pow(xx,2) + math.pow(yy,2) - xx * yy + 3*math.pow(xy,2))
        elif (self.mode == "Plane Strain"):
            mean_stress = math.sqrt((math.pow(xx, 2) + math.pow(yy, 2)) * (math.pow(self.poisson, 2) - self.poisson + 1)
					+ xx*yy * (2 * math.pow(self.poisson, 2) - 2 * self.poisson - 1)+3*math.pow(xy, 2))
        else:
            raise Exception("Invalid argument for CST Mode. Must be Plane Stress or Plane Strain!")
        return mean_stress

    def compute_displacement_magnitude(self):
        displacements = np.empty([9])
        for i in range(3):
            displacements[i] = self.node_1.get_displacement()[i]
            displacements[i + 3] = self.node_2.get_displacement()[i]
            displacements[i + 6] = self.node_3.get_displacement()[i]
        return np.linalg.norm(displacements)


