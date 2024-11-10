"""
Linear 8-Node Hexahedron Element using the selective reduced integration scheme.
Thanks to [1] for providing a small test case on which the presented element here is based on.

[1]:
Dinh, Tien. (2022). Re: Stiffness matrix of abaqus element c3d8 and by my 8 noded isoparametric element is not matching - can anyone help?. 
Retrieved from: https://www.researchgate.net/post/Stiffness-matrix-of-abaqus-element-c3d8-and-by-my-8-noded-isoparametric-element-is-not-matching-can-anyone-help/6238a2c02aa16869d76f2637/citation/download. 

"""
import numpy as np
import math


class Element_Linear_Hexahedral_SRI(object):
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
        self.coordinates = np.row_stack((pos_node_1, pos_node_2, pos_node_3, pos_node_4,
                                   pos_node_5, pos_node_6, pos_node_7, pos_node_8))

    def get_e_mod(self):
        return self.e_mod

    def get_poisson(self):
        return self.poisson

    def get_density(self):
        return self.density

    def get_node_1(self):
        return self.node_1

    def get_node_2(self):
        return self.node_2

    def get_node_3(self):
        return self.node_3

    def get_node_4(self):
        return self.node_4

    def get_node_5(self):
        return self.node_5

    def get_node_6(self):
        return self.node_6

    def get_node_7(self):
        return self.node_7

    def get_node_8(self):
        return self.node_8

    def get_all_nodes(self):
        return [self.node_1, self.node_2, self.node_3, self.node_4, self.node_5, self.node_6, self.node_7, self.node_8]
    
    def get_number_of_nodes(self):
        return 8

    def get_dof_numbers(self):
        return self.dof_numbers
    
    def compute_min_length(self):
        """Used for explicit time incrementation"""

        d12 = self.distance(self.node_1.get_coordinates(), self.node_2.get_coordinates())
        d14 = self.distance(self.node_1.get_coordinates(), self.node_4.get_coordinates())
        d15 = self.distance(self.node_1.get_coordinates(), self.node_5.get_coordinates())

        return min(d12, d14, d15)
    
    def distance(self, node1, node2):
        """Calculates eucledian distance between two nodes"""
        return math.sqrt((node2[0] - node1[0])**2 + 
                        (node2[1] - node1[1])**2 + 
                        (node2[2] - node1[2])**2)

    def enumerate_dofs(self):
        for i in range(0, 3):
            self.dof_numbers[i] = self.node_1.get_dof_numbers()[i]
            self.dof_numbers[i+3] = self.node_2.get_dof_numbers()[i]
            self.dof_numbers[i+6] = self.node_3.get_dof_numbers()[i]
            self.dof_numbers[i+9] = self.node_4.get_dof_numbers()[i]
            self.dof_numbers[i+12] = self.node_5.get_dof_numbers()[i]
            self.dof_numbers[i+15] = self.node_6.get_dof_numbers()[i]
            self.dof_numbers[i+18] = self.node_7.get_dof_numbers()[i]
            self.dof_numbers[i+21] = self.node_8.get_dof_numbers()[i]

    def get_gauss_point(self, gp_index):

        if gp_index == 0:
            natural_coordinate_1 = -1 / math.sqrt(3)
            natural_coordinate_2 = -1 / math.sqrt(3)
            natural_coordinate_3 = -1 / math.sqrt(3)
        elif gp_index == 1:
            natural_coordinate_1 = 1 / math.sqrt(3)
            natural_coordinate_2 = -1 / math.sqrt(3)
            natural_coordinate_3 = -1 / math.sqrt(3)
        elif gp_index == 2:
            natural_coordinate_1 = 1 / math.sqrt(3)
            natural_coordinate_2 = 1 / math.sqrt(3)
            natural_coordinate_3 = -1 / math.sqrt(3)
        elif gp_index == 3:
            natural_coordinate_1 = -1 / math.sqrt(3)
            natural_coordinate_2 = 1 / math.sqrt(3)
            natural_coordinate_3 = -1 / math.sqrt(3)
        elif gp_index == 4:
            natural_coordinate_1 = -1 / math.sqrt(3)
            natural_coordinate_2 = -1 / math.sqrt(3)
            natural_coordinate_3 = 1 / math.sqrt(3)
        elif gp_index == 5:
            natural_coordinate_1 = 1 / math.sqrt(3)
            natural_coordinate_2 = -1 / math.sqrt(3)
            natural_coordinate_3 = 1 / math.sqrt(3)
        elif gp_index == 6:
            natural_coordinate_1 = 1 / math.sqrt(3)
            natural_coordinate_2 = 1 / math.sqrt(3)
            natural_coordinate_3 = 1 / math.sqrt(3)
        elif gp_index == 7:
            natural_coordinate_1 = -1 / math.sqrt(3)
            natural_coordinate_2 = 1 / math.sqrt(3)
            natural_coordinate_3 = 1 / math.sqrt(3)
        elif gp_index == -1:
            natural_coordinate_1 = 0.0
            natural_coordinate_2 = 0.0
            natural_coordinate_3 = 0.0
        else:
            raise Exception("Invalid Index for Gauss Point in Hexhedral Element")

        return np.array([natural_coordinate_1, natural_coordinate_2, natural_coordinate_3])

    def shape_func_derivates(self, gauss_point_index):
        natural_coordinates = self.get_gauss_point(gauss_point_index)
        # xi_1 as abbreviation for the first natural coordiante
        xi_1 = natural_coordinates[0]
        xi_2 = natural_coordinates[1]
        xi_3 = natural_coordinates[2]

        v1 = [-(0.125)*(1.0-xi_2)*(1.0-xi_3), -(0.125)*(1.0-xi_1)*(1.0-xi_3), -(0.125)*(1.0-xi_1)*(1.0-xi_2)]
        v2 = [(0.125)*(1.0-xi_2)*(1.0-xi_3), -(0.125)*(1.0+xi_1)*(1.0-xi_3), -(0.125)*(1.0+xi_1)*(1.0-xi_2)]
        v3 = [(0.125)*(1.0+xi_2)*(1.0-xi_3), (0.125)*(1.0+xi_1)*(1.0-xi_3), -(0.125)*(1.0+xi_1)*(1.0+xi_2)]
        v4 = [-(0.125)*(1.0+xi_2)*(1.0-xi_3), (0.125)*(1.0-xi_1)*(1.0-xi_3), -(0.125)*(1.0-xi_1)*(1.0+xi_2)]
        v5 = [-(0.125)*(1.0-xi_2)*(1.0+xi_3), -(0.125)*(1.0-xi_1)*(1.0+xi_3), (0.125)*(1.0-xi_1)*(1.0-xi_2)]
        v6 = [(0.125)*(1.0-xi_2)*(1.0+xi_3), -(0.125)*(1.0+xi_1)*(1.0+xi_3), (0.125)*(1.0+xi_1)*(1.0-xi_2)]
        v7 = [(0.125)*(1.0+xi_2)*(1.0+xi_3), (0.125)*(1.0+xi_1)*(1.0+xi_3), (0.125)*(1.0+xi_1)*(1.0+xi_2)]
        v8 = [-(0.125)*(1.0+xi_2)*(1.0+xi_3), (0.125)*(1.0-xi_1)*(1.0+xi_3), (0.125)*(1.0-xi_1)*(1.0+xi_2)]

        shape_function_derivatives = np.column_stack([v1, v2, v3, v4, v5, v6, v7, v8])

        return shape_function_derivatives

    def compute_jacobian(self, gauss_point_index):

        jacobian = self.shape_func_derivates(gauss_point_index) @ self.coordinates

        return jacobian

    def compute_determinant(self, gauss_point_index):
        jacobian = self.compute_jacobian(gauss_point_index)
        determinant = np.linalg.det(jacobian)
        return determinant

    def compute_inverse_jacobian(self, gauss_point_index):
        jacobian = self.compute_jacobian(gauss_point_index)
        try:
            inverse_jacobian = np.linalg.inv(jacobian)
        except:
            raise Exception("Could not calculate matrix in hexahedral element in compute_b_operator")

        return np.array(inverse_jacobian)

    def compute_material_matrix(self):

        E = self.e_mod
        nu = self.poisson

        C = ((E) / ((1 + nu) * (1 - 2 * nu))) * np.array([
            [1 - nu, nu, nu, 0, 0, 0],
            [nu, 1 - nu, nu, 0, 0, 0],
            [nu, nu, 1 - nu, 0, 0, 0],
            [0, 0, 0, (1 - 2 * nu) / 2, 0, 0],
            [0, 0, 0, 0, (1 - 2 * nu) / 2, 0],
            [0, 0, 0, 0, 0, (1 - 2 * nu) / 2]
        ])

        return C
    
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

    def compute_element_mean_stress(self):
        """Von Mises mean stress"""
        stress = self.compute_element_stress()

        xx = np.abs(stress[0])
        yy = np.abs(stress[1])
        zz = np.abs(stress[2])
        xy = np.abs(stress[3])
        yz = np.abs(stress[4])
        zx = np.abs(stress[5])

        mean_stress = np.sqrt(0.5 * ((xx - yy)**2 + (yy - zz)**2 + (zz - xx)**2 + 6 * (xy**2 + yz**2 + zx**2)))

        return mean_stress

    def compute_displacement_magnitude(self):
        # Computes the magnitude of the displacement
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

        magnitude = np.linalg.norm(displacements)

        return magnitude
