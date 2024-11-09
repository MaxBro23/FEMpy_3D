import math
import numpy as np


class Element_Truss(object):

    def __init__(self, e_mod, area, density, node_1, node_2, node_1_idx, node_2_idx):
        self.e_mod = e_mod
        self.area = area
        self.density = density
        self.node_1 = node_1
        self.node_2 = node_2
        self.node_1_idx = node_1_idx
        self.node_2_idx = node_2_idx
        self.dof_numbers = np.empty([6], dtype=int)  # random values

    def get_length(self):
        """
        :return: float
            Eucledian distance between nodes
        """
        n1 = self.node_1.get_coordinates()
        n2 = self.node_2.get_coordinates()
        return math.sqrt(math.pow(n1[0] - n2[0],2) + math.pow(n1[1] - n2[1],2) + math.pow(n1[2] - n2[2],2))
    
    def compute_min_length(self):
        """Used for explicit time incrementation. Specific ethod name was necessary"""
        return self.get_length()
    
    def get_e_mod(self):
        """
        :return: float
            E-Modulus/Young's Modulus of the truss element
        """
        return self.e_mod

    def get_area(self):
        """
        :return: float
            Area of the truss element
        """
        return self.area
    
    def get_density(self):
        """
        :return: float
            Density of the truss
        """
        return self.density

    def get_node_1(self):
        """
        :return: node object
            Return 1st node object
        """
        return self.node_1

    def get_node_2(self):
        """
        :return: node object
            Return 2nd node
        """
        return self.node_2

    def get_all_nodes(self):
        """
        :return: list
            Return list with both node objects
        """
        return [self.node_1, self.node_2]

    def enumerate_dofs(self):
        """Enumerate the free DOF's
        :return: None
        """
        for i in range(0, 3):
            self.dof_numbers[i] = self.node_1.get_dof_numbers()[i]
            self.dof_numbers[i+3] = self.node_2.get_dof_numbers()[i]

    def get_dof_numbers(self):
        """
        :return: np.array
            1D-np.array with DOF numbers
        """
        return self.dof_numbers

    def get_number_of_nodes(self):
        """
        :return: float
            Number of nodes of the element
        """
        return 2.0


    def compute_element_k(self):
        """
        :return: 2D-np.array
            Return global element stiffness tensor
        """
        k_local = np.array([[1, -1], [-1, 1]])
        EAoverL = self.e_mod * self.area / self.get_length()
        k_local = k_local*EAoverL
        cosX = (self.node_2.get_coordinates()[0] - self.node_1.get_coordinates()[0]) / self.get_length()
        cosY = (self.node_2.get_coordinates()[1] - self.node_1.get_coordinates()[1]) / self.get_length()
        cosZ = (self.node_2.get_coordinates()[2] - self.node_1.get_coordinates()[2]) / self.get_length()
        T_transposed = np.array([[cosX, 0], [cosY, 0], [cosZ, 0], [0, cosX], [0, cosY], [0, cosZ]])
        k_global = np.matmul(np.matmul(T_transposed, k_local), T_transposed.transpose())
        return k_global

    def compute_element_m(self):
        """
        :return: 2D-np.array
            Return global element mass tensor
        """
        m_local = np.array([[2, 1], [1, 2]])
        temp_value = (self.density * self.area * self.get_length()) / (6)
        m_local = m_local * temp_value
        cosX = (self.node_2.get_coordinates()[0] - self.node_1.get_coordinates()[0]) / self.get_length()
        cosY = (self.node_2.get_coordinates()[1] - self.node_1.get_coordinates()[1]) / self.get_length()
        cosZ = (self.node_2.get_coordinates()[2] - self.node_1.get_coordinates()[2]) / self.get_length()
        T_transposed = np.array([[cosX, 0], [cosY, 0], [cosZ, 0], [0, cosX], [0, cosY], [0, cosZ]])
        m_global = np.matmul(np.matmul(T_transposed, m_local), T_transposed.transpose())
        return m_global

    def compute_element_mean_stress(self):
        """
        :return: float
            Return absolute Von Mises mean stress of the truss element
        """
        # Get Coords and displacement of nodes
        x1 = self.node_1.get_coordinates()
        x2 = self.node_2.get_coordinates()
        u1 = self.node_1.get_displacement()
        u2 = self.node_2.get_displacement()

        # Get length
        L0 = self.get_length()

        # Unit vector along element
        e = (x2 - x1) / L0

        # Length difference
        delta_L = np.dot((u2 - u1), e)

        # Strain
        strain = delta_L / L0

        # Stress
        stress = self.e_mod * strain

        return np.abs(stress)

    def compute_displacement_magnitude(self):
        """
        :return: float
            Return magnitude of displacement in all three dimensions
        """
        return np.linalg.norm([self.node_1.get_displacement(), self.node_2.get_displacement()])




