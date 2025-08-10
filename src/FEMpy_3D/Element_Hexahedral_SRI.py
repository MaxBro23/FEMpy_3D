"""
Linear 8-Node Hexahedron Element using the selective reduced integration scheme.
Thanks to [1] for providing a small test case on which the presented element here is based on.

[1]:
Dinh, Tien. (2022). Re: Stiffness matrix of abaqus element c3d8 and by my 8 noded isoparametric element is not matching - can anyone help?. 
Retrieved from: https://www.researchgate.net/post/Stiffness-matrix-of-abaqus-element-c3d8-and-by-my-8-noded-isoparametric-element-is-not-matching-can-anyone-help/6238a2c02aa16869d76f2637/citation/download. 

"""
import numpy as np
from .Element_Hexahedral_Base import Element_Hexahedral_Base

class Element_Linear_Hexahedral_SRI(Element_Hexahedral_Base):
    # Selective reduced integration method
    def __init__(self, e_mod, poisson, density, node_1, node_2, node_3, node_4,
                 node_5, node_6, node_7, node_8):
        self.e_mod = e_mod
        self.poisson = poisson
        self.density =density
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

    
    def compute_b_operator_mechanic(self, X, gauss_point_index):
        integration_coords = self.get_gauss_point(gauss_point_index)
        xi, eta, zeta = integration_coords[0], integration_coords[1], integration_coords[2]
        # B-Operator for mechanical behaviour of the element
        
        GN = np.array([[-0.125 * (1 + zeta) * (1 + eta), -0.125 * (1 + zeta) * (1 - eta), -0.125 * (1 - zeta) * (1 - eta), -0.125 * (1 - zeta) * (1 + eta),
                        0.125 * (1 + zeta) * (1 + eta), 0.125 * (1 + zeta) * (1 - eta), 0.125 * (1 - zeta) * (1 - eta), 0.125 * (1 - zeta) * (1 + eta)],
                        [0.125 * (1 + zeta) * (1 - xi), -0.125 * (1 + zeta) * (1 - xi), -0.125 * (1 - zeta) * (1 - xi), 0.125 * (1 - zeta) * (1 - xi), 
                        0.125 * (1 + zeta) * (1 + xi), -0.125 * (1 + zeta) * (1 + xi), -0.125 * (1 - zeta) * (1 + xi), 0.125 * (1 - zeta) * (1 + xi)],
                        [0.125 * (1 + eta) * (1 - xi), 0.125 * (1 - eta) * (1 - xi), -0.125 * (1 - eta) * (1 - xi), -0.125 * (1 + eta) * (1 - xi),
                        0.125 * (1 + eta) * (1 + xi), 0.125 * (1 - eta) * (1 + xi), -0.125 * (1 - eta) * (1 + xi), -0.125 * (1 + eta) * (1 + xi)]])

        Jac = np.dot(GN, X)

        invJ = np.linalg.inv(Jac)
        
        temp_matrix = np.dot(invJ, GN)
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


        return b_operator

    def compute_b_operator_volumetric(self, X, gauss_point_index):
        # Volumetric B matrix of a hexahedral element for the selective reduced integration method
        integration_coords = self.get_gauss_point(gauss_point_index)
        xi, eta, zeta = integration_coords[0], integration_coords[1], integration_coords[2]


        GN = np.array([[-0.125 * (1 + zeta) * (1 + eta), -0.125 * (1 + zeta) * (1 - eta), -0.125 * (1 - zeta) * (1 - eta), -0.125 * (1 - zeta) * (1 + eta),
                        0.125 * (1 + zeta) * (1 + eta), 0.125 * (1 + zeta) * (1 - eta), 0.125 * (1 - zeta) * (1 - eta), 0.125 * (1 - zeta) * (1 + eta)],
                        [0.125 * (1 + zeta) * (1 - xi), -0.125 * (1 + zeta) * (1 - xi), -0.125 * (1 - zeta) * (1 - xi), 0.125 * (1 - zeta) * (1 - xi), 
                        0.125 * (1 + zeta) * (1 + xi), -0.125 * (1 + zeta) * (1 + xi), -0.125 * (1 - zeta) * (1 + xi), 0.125 * (1 - zeta) * (1 + xi)],
                        [0.125 * (1 + eta) * (1 - xi), 0.125 * (1 - eta) * (1 - xi), -0.125 * (1 - eta) * (1 - xi), -0.125 * (1 + eta) * (1 - xi),
                        0.125 * (1 + eta) * (1 + xi), 0.125 * (1 - eta) * (1 + xi), -0.125 * (1 - eta) * (1 + xi), -0.125 * (1 + eta) * (1 + xi)]])
        
        
        Jac = np.dot(GN, X)
        invJ = np.linalg.inv(Jac)
        dNdx = np.dot(invJ, GN)

        BMat = []
        for i in range(8):
            BMat.append(np.array([
                                [dNdx[0, i], dNdx[1, i], dNdx[2, i]],
                                [dNdx[0, i], dNdx[1, i], dNdx[2, i]],
                                [dNdx[0, i], dNdx[1, i], dNdx[2, i]],
                                [0.0, 0.0, 0.0],
                                [0.0, 0.0, 0.0],
                                [0.0, 0.0, 0.0]]))

        # Stack b-operator
        BMat_combined = np.hstack(BMat)

        return 1.0 / 3.0 * BMat_combined

    
    def compute_b_operator_deviatoric(self, gauss_point_index):
        # Deviatoric B matrix of a hexahedral element for selective reduced integration method

        # Subtract volumetric part from mechanical part to get deviatoric part
        b_operator_mechanic = self.compute_b_operator_mechanic(self.coordinates, gauss_point_index)
        b_operator_volumetric = self.compute_b_operator_volumetric(self.coordinates, gauss_point_index)
        b_operator_deviatoric = (b_operator_mechanic - b_operator_volumetric)

        return b_operator_deviatoric

    def compute_element_stiffness_matrix_deviatoric(self, gauss_point_index):

        b_operator_dev = self.compute_b_operator_deviatoric(gauss_point_index)
        stiffness_matrix = np.dot(np.dot(b_operator_dev.T, self.compute_material_matrix()), b_operator_dev) * self.compute_determinant(gauss_point_index)

        return stiffness_matrix
    
    def compute_element_stiffness_matrix_volumetric(self, gauss_point_index):
        # Gauss point coordiantes
        weight = 8.0
        # Get volumetric B-Operator
        b_vol= self.compute_b_operator_volumetric(self.coordinates, gauss_point_index)
        # Calc volumetric part of the stiffness matrix
        stiffness_matrix = weight * np.dot(np.dot(b_vol.T, self.compute_material_matrix()), b_vol) * self.compute_determinant(gauss_point_index)

        return stiffness_matrix

    def compute_element_k(self):
        # Compute deviatoric part of the stiffness matrix
        k1 = self.compute_element_stiffness_matrix_deviatoric(0)
        k2 = self.compute_element_stiffness_matrix_deviatoric(1)
        k3 = self.compute_element_stiffness_matrix_deviatoric(2)
        k4 = self.compute_element_stiffness_matrix_deviatoric(3)
        k5 = self.compute_element_stiffness_matrix_deviatoric(4)
        k6 = self.compute_element_stiffness_matrix_deviatoric(5)
        k7 = self.compute_element_stiffness_matrix_deviatoric(6)
        k8 = self.compute_element_stiffness_matrix_deviatoric(7)
        # Calculate volumetric part of the stiffness matrix
        k9 = self.compute_element_stiffness_matrix_volumetric(-1)
        # Add all up
        global_stiffness_matrix = k1 + k2 + k3 + k4 + k5 + k6 + k7 + k8 + k9
        return global_stiffness_matrix

    def compute_element_stress(self):
        # Deviatoric part
        b1 = self.compute_b_operator_deviatoric(0)
        b2 = self.compute_b_operator_deviatoric(1)
        b3 = self.compute_b_operator_deviatoric(2)
        b4 = self.compute_b_operator_deviatoric(3)
        b5 = self.compute_b_operator_deviatoric(4)
        b6 = self.compute_b_operator_deviatoric(5)
        b7 = self.compute_b_operator_deviatoric(6)
        b8 = self.compute_b_operator_deviatoric(7)
        # Volumetric part
        b9 = self.compute_b_operator_volumetric(self.coordinates, -1)

        # 6 x 24 Array
        b_operator = (b1 + b2 + b3 + b4+ b5 + b6 + b7 + b8)*1/8 + b9
        
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
        element_stress = np.dot(material_matrix,strains)

        return element_stress
