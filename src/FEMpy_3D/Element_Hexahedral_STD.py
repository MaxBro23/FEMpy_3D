"""
Linear 8-Node Hexahedron Element using the standard integration scheme as in the text books.
"""
import numpy as np
from .Element_Hexahedral_Base import Element_Hexahedral_Base

class Element_Linear_Hexahedral_STD(Element_Hexahedral_Base):

    def __init__(self, e_mod, poisson, density, node_1, node_2, node_3, node_4,
                 node_5, node_6, node_7, node_8):
        self.e_mod = e_mod
        self.poisson = poisson
        self.density = density
        self.node_1 = node_1
        self.node_2 = node_2
        self.node_3 = node_3
        self.node_4 = node_4
        self.node_5 = node_5
        self.node_6 = node_6
        self.node_7 = node_7
        self.node_8 = node_8
        self.dof_numbers = np.empty([24], dtype=int)
        # Array of coordiantes of all nodes
        pos_node_1 = self.node_1.get_coordinates()
        pos_node_2 = self.node_2.get_coordinates()
        pos_node_3 = self.node_3.get_coordinates()
        pos_node_4 = self.node_4.get_coordinates()
        pos_node_5 = self.node_5.get_coordinates()
        pos_node_6 = self.node_6.get_coordinates()
        pos_node_7 = self.node_7.get_coordinates()
        pos_node_8 = self.node_8.get_coordinates()
        self.coordinates = np.vstack((pos_node_1, pos_node_2, pos_node_3, pos_node_4,
                                   pos_node_5, pos_node_6, pos_node_7, pos_node_8))


    def compute_b_operator(self, gauss_point_index):
        temp_matrix = self.compute_inverse_jacobian(gauss_point_index) @ self.shape_func_derivates(gauss_point_index)
        # Create B-Operator
        b_operator = np.zeros((6,24))
        for i in range(8):
            # Set 1st column from top to bottom
            b_operator[0][3*i] = temp_matrix[0][i]
            b_operator[3][3*i] = temp_matrix[1][i]
            b_operator[5][3*i] = temp_matrix[2][i]
            # Set 2nd Column
            b_operator[3][3*i+1] = temp_matrix[0][i]
            b_operator[1][3*i+1] = temp_matrix[1][i]
            b_operator[4][3*i+1] = temp_matrix[2][i]
            # Set 3rd Column
            b_operator[5][3*i+2] = temp_matrix[0][i]
            b_operator[4][3*i+2] = temp_matrix[1][i]
            b_operator[2][3*i+2] = temp_matrix[2][i]

        return np.array(b_operator)

    def compute_element_stiffness_matrix_at_one_gauss_point(self, gauss_point_index):
        transposed_B = np.transpose(self.compute_b_operator(gauss_point_index))
        b_operator = self.compute_b_operator(gauss_point_index)
        stiffness_matrix = transposed_B @ self.compute_material_matrix() @ b_operator * self.compute_determinant(gauss_point_index)

        return stiffness_matrix

    def compute_element_k(self):
        k1 = self.compute_element_stiffness_matrix_at_one_gauss_point(0)
        k2 = self.compute_element_stiffness_matrix_at_one_gauss_point(1)
        k3 = self.compute_element_stiffness_matrix_at_one_gauss_point(2)
        k4 = self.compute_element_stiffness_matrix_at_one_gauss_point(3)
        k5 = self.compute_element_stiffness_matrix_at_one_gauss_point(4)
        k6 = self.compute_element_stiffness_matrix_at_one_gauss_point(5)
        k7 = self.compute_element_stiffness_matrix_at_one_gauss_point(6)
        k8 = self.compute_element_stiffness_matrix_at_one_gauss_point(7)

        global_stiffness_matrix = k1 + k2 + k3 + k4 + k5 + k6 + k7 + k8
        return global_stiffness_matrix

    def compute_element_stress(self):

        b1 = self.compute_b_operator(0)
        b2 = self.compute_b_operator(1)
        b3 = self.compute_b_operator(2)
        b4 = self.compute_b_operator(3)
        b5 = self.compute_b_operator(4)
        b6 = self.compute_b_operator(5)
        b7 = self.compute_b_operator(6)
        b8 = self.compute_b_operator(7)

        # 6 x 24 Array
        b_operator = (b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8)* 1 / 8
        # 24x1 List
        displacements = np.empty([24])

        for i in range(3):
            displacements[i] = self.node_1.get_displacement()[i]
            displacements[i + 3] = self.node_2.get_displacement()[i]
            displacements[i + 6] = self.node_3.get_displacement()[i]
            displacements[i + 9] = self.node_4.get_displacement()[i]
            displacements[i + 12] = self.node_5.get_displacement()[i]
            displacements[i + 15] = self.node_6.get_displacement()[i]
            displacements[i + 18] = self.node_7.get_displacement()[i]
            displacements[i + 21] = self.node_8.get_displacement()[i]

        strains = b_operator @ displacements
        # C Matrix 6x6
        material_matrix = self.compute_material_matrix()
        # element_stress is a 6x1 array
        element_stress = material_matrix @ strains

        return element_stress
