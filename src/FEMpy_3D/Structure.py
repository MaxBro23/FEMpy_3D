from datetime import datetime
import numpy as np
from .Node import Node
from .Element_Truss import Element_Truss
from .Element_CST import Element_CST
from .Element_Hexahedral_SRI import Element_Linear_Hexahedral_SRI
from .Element_Hexahedral import Element_Linear_Hexahedral


class Structure(object):
    """Structure object with all nodes and elements."""
    nodes = np.empty(0, dtype=object)
    cst_nodes = np.empty(0, dtype=object)
    cst_nodes_coordinates_list = []
    linear_hexa_nodes = np.empty(0, dtype=object)
    linear_hexa_nodes_coordinates_list = []
    nodes_idx = np.empty(0, dtype=list)
    truss_elements = np.empty(0, dtype=object)
    cst_elements = np.empty(0, dtype=object)
    linear_hexa_elements = np.empty(0, dtype=object)
    all_elements = np.empty(0, dtype=object)

    def add_node(self, x, y, z):
        """
        :param x: float
            X-coordinate of the node
        :param y: float
            Y-coordinate of the node
        :param z: float
            Z-coordinate of the node
        :return:
        """
        node = Node(x, y, z, index=len(self.nodes))
        self.nodes = np.append(self.nodes, node, axis=None)
        return node

    def get_node(self, i):
        """
        :param i: int
            Index of the node
        :return: object
            Node object
        """
        node = self.nodes[i]
        return node

    def get_all_nodes(self):
        """
        :return: np.array
            Np-array with all nodes
        """
        return self.nodes

    def get_number_of_nodes(self):
        """
        :return: int
            Number of all nodes of the structure
        """
        return len(self.nodes)

    def get_truss_element(self, i):
        """
        :param i: int
            Index of the truss
        :return: object
            Truss object
        """
        return self.truss_elements[i]

    def get_cst_element(self, i):
        """
        :param i:
            Index of the CST Element
        :return: object
            CSTElement
        """
        return self.cst_elements[i]

    def get_all_elements(self):
        """
        :return: np.array
             Array with all elements as objects
        """
        return self.all_elements

    def get_linear_hexa_element(self, i):
        """
        :param i:
            Index of the hexa element
        :return: object
            Hexahedral-Element object
        """
        return self.linear_hexa_elements[i]
    
    def get_number_of_elements(self):
        """
        :return: int
            Number of elements in the structure
        """
        return len(self.all_elements)

    def get_coordinates_of_nodes_as_list(self):
        """
        :return: list
            List with coordinates of all nodes of the structure
        """
        coordinates_list = []
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            coords = node.get_coordinates()
            coordinates_list.append(coords)
        return coordinates_list

    def get_displaced_coordiantes_of_nodes_as_list(self):
        """
        :return: list
            List with displaced coordiantes of all nodes
        """
        displaced_coordinates_list = []
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            displacements = node.get_displacement()
            coordinates = node.get_coordinates()
            displaced_coordinates = coordinates + displacements
            displaced_coordinates_list.append(displaced_coordinates)
        return displaced_coordinates_list

    def get_displaced_coordinates_of_CST_nodes_as_list(self):
        """
        :return: list
            List with displaced coordinates of all CST-Elements nodes
        """
        displaced_coordinates_list = []
        for i in range(len(self.cst_nodes)):
            node = self.cst_nodes[i]
            displacements = node.get_displacement()
            coordinates = node.get_coordinates()
            displaced_coordinates = coordinates + displacements
            displaced_coordinates_list.append(displaced_coordinates)
        return displaced_coordinates_list

    def get_displaced_coordinates_of_linear_hexahedral_nodes_as_list(self):
        """
        :return: list
            List with displaced coordinates of all Hexahedral Elements
        """
        displaced_coordinates_list = []
        for i in range(len(self.linear_hexa_nodes)):
            node = self.linear_hexa_nodes[i]
            displacements = node.get_displacement()
            coordinates = node.get_coordinates()
            displaced_coordinates = coordinates + displacements
            displaced_coordinates_list.append(displaced_coordinates)
        return displaced_coordinates_list

    def add_truss(self, e_mod, area, density, node1, node2):
        """
        :param e_mod: float
            E-Modulus/Young's Modulus of the truss
        :param area: float
            Area of the truss
        :param density: float
            Density of the truss
        :param node1: int
            Index of the first node of the truss
        :param node2: int
            Index of the second node of the truss
        :return: object
            Truss Element
        """
        truss = Element_Truss(e_mod, area, density, self.nodes[node1], self.nodes[node2], node_1_idx=node1, node_2_idx=node2)
        self.nodes_idx = np.append(self.nodes_idx, [node1, node2], axis=None)
        self.truss_elements = np.append(self.truss_elements, truss, axis=None)
        self.all_elements = np.append(self.all_elements, truss, axis=None)
        return truss

    def add_CST(self, e_mod, poisson_ratio, height, density, node_1, node_2, node_3):
        """
        :param e_mod: float
            E-Modulus/Young's Modulus of the CST
        :param poisson_ratio: float
            Poisson ratio of the CST
        :param height: float
            Height of the CST
        :param density: float
            Density of the CST
        :param node_1: int
            Node index of the first node
        :param node_2: int
            Node index of the second node
        :param node_3: int
            Node index of the third node
        :return: object
            CST Element
        """
        cst = Element_CST(e_mod, poisson_ratio, height, density, self.nodes[node_1], self.nodes[node_2], self.nodes[node_3])
        self.nodes_idx = np.append(self.nodes_idx, [node_1, node_2, node_3], axis=None)
        self.cst_nodes_coordinates_list.append(self.nodes[node_1].get_coordinates())
        self.cst_nodes_coordinates_list.append(self.nodes[node_2].get_coordinates())
        self.cst_nodes_coordinates_list.append(self.nodes[node_3].get_coordinates())
        self.cst_nodes = np.append(self.cst_nodes, cst.get_node_1(), axis=None)
        self.cst_nodes = np.append(self.cst_nodes, cst.get_node_2(), axis=None)
        self.cst_nodes = np.append(self.cst_nodes, cst.get_node_3(), axis=None)
        self.all_elements = np.append(self.all_elements, cst, axis=None)
        self.cst_elements = np.append(self.cst_elements, cst, axis=None)
        return cst

    def add_linear_hexahedral_SRI(self, e_mod, poisson_ratio, density, node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8):
        """ Adds linear hexahedral element with selective reduced integration scheme
        :param e_mod: float
            E-Modulus/Young's Modulus of the Hexahedron
        :param poisson_ratio: float
            Poisson ratio of the Hexahedron
        :param density: float
            Density of the Hexahedron
        :param node_1: int
            Index of the first node
        :param node_2: int
            Index of the second node
        :param node_3: int
            Index of the third node
        :param node_4: int
            Index of the fourth node
        :param node_5: int
            Index of the fifth node
        :param node_6: int
            Index of the sixth node
        :param node_7: int
            Index of the seventh node
        :param node_8: int
            Index of the eigth node
        :return: object
            Hexahedron Element with SRI algorithm
        """
        hexa = Element_Linear_Hexahedral_SRI(e_mod, poisson_ratio, density, self.nodes[node_1], self.nodes[node_2], self.nodes[node_3],self.nodes[node_4], self.nodes[node_5], self.nodes[node_6],self.nodes[node_7], self.nodes[node_8])
        self.nodes_idx = np.append(self.nodes_idx, [node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8], axis=None)
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_1].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_2].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_3].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_4].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_5].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_6].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_7].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_8].get_coordinates())
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_1(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_2(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_3(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_4(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_5(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_6(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_7(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_8(), axis=None)
        self.all_elements = np.append(self.all_elements, hexa, axis=None)
        self.linear_hexa_elements = np.append(self.linear_hexa_elements, hexa, axis=None)
        return hexa
    
    def add_linear_hexahedral(self, e_mod, poisson_ratio, density, node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8):
        """ Adds linear hexahedral element with standard integration algorithm as thought in text books
        :param e_mod: float
            E-Modulus/Young's Modulus of the Hexahedron
        :param poisson_ratio: float
            Poisson ratio of the Hexahedron
        :param density: float
            Density of the Hexahedron
        :param node_1: int
            Index of the first node
        :param node_2: int
            Index of the second node
        :param node_3: int
            Index of the third node
        :param node_4: int
            Index of the fourth node
        :param node_5: int
            Index of the fifth node
        :param node_6: int
            Index of the sixth node
        :param node_7: int
            Index of the seventh node
        :param node_8: int
            Index of the eigth node
        :return: object
            Hexahedron Element with standard algorithm
        """
        hexa = Element_Linear_Hexahedral(e_mod, poisson_ratio, density, self.nodes[node_1], self.nodes[node_2], self.nodes[node_3],self.nodes[node_4], self.nodes[node_5], self.nodes[node_6],self.nodes[node_7], self.nodes[node_8])
        self.nodes_idx = np.append(self.nodes_idx, [node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8], axis=None)
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_1].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_2].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_3].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_4].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_5].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_6].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_7].get_coordinates())
        self.linear_hexa_nodes_coordinates_list.append(self.nodes[node_8].get_coordinates())
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_1(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_2(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_3(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_4(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_5(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_6(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_7(), axis=None)
        self.linear_hexa_nodes = np.append(self.linear_hexa_nodes, hexa.get_node_8(), axis=None)
        self.all_elements = np.append(self.all_elements, hexa, axis=None)
        self.linear_hexa_elements = np.append(self.linear_hexa_elements, hexa, axis=None)


    def get_max_stress(self):
        """
        :return: float
            Return max von mises stress present in the structure
        """
        mean_stress_list = []
        for i in range(len(self.all_elements)):
            mean_stress = self.all_elements[i].compute_element_mean_stress()
            mean_stress_list.append(mean_stress)
        max_mean_stress = max(mean_stress_list)
        return max_mean_stress


    def set_CST_mode(self, cst_mode_as_string):
        """ Set the mode for all CST elements
        :param cst_mode_as_string: str
            Either "Plane Strain" or "Plane Stress"
        :return: None
        """
        for i in range(len(self.cst_elements)):
            if cst_mode_as_string == "Plane Strain":
                self.cst_elements[i].set_mode(cst_mode_as_string)
            elif cst_mode_as_string == "Plane Stress":
                self.cst_elements[i].set_mode(cst_mode_as_string)
            else:
                raise Exception("CST Mode must be Plane Strain or Plane Stress!")

    def get_number_of_truss_elements(self):
        """
        :return: int
            Return number of truss elements
        """
        return len(self.truss_elements)

    def get_number_of_CST_elements(self):
        """
        :return: int
            Return number of cst elements
        """
        return len(self.cst_elements)

    def get_number_of_linear_hexahedral_elements(self):
        """
        :return: int
            Return number of hexahedral elements
        """
        return len(self.linear_hexa_elements)

    def get_number_of_equations(self):
        """
        :return: int
            Return number of equations (number of free DOF'S) of the structure
        """
        eqn = 0
        for i in range(0, len(self.nodes)):
            eqn = self.nodes[i].enumerate_dofs_node(eqn)
        for i in range(0, len(self.all_elements)):
            self.all_elements[i].enumerate_dofs()
        return eqn

    def assemble_global_load_vector(self, number_of_equations):
        """
        :param number_of_equations: int
            Number of equations (free DOF'S)
        :return: np.array
            Global load vector of the system
        """
        load_vector = np.zeros(number_of_equations)
        for i in range(0, len(self.nodes)):
            for j in range(0, 3):
                if self.nodes[i].get_constraint().is_free(j) == True:
                    dof = self.nodes[i].get_dof_numbers()[j]
                    load_vector[int(dof)] = self.nodes[i].get_force().get_component(j)
        return load_vector

    def assemble_global_stiffness_tensor(self, number_of_equations):
        """
        :param number_of_equations: int
            Number of equations (free DOF'S)
        :return: np.array
            Return the stiffess matrix of the structure
        """
        # Assign dtype=int otherwise inverting the matrix will fail later due to long floating point numbers
        global_stiffness_tensor = np.zeros(shape=[number_of_equations,number_of_equations])
        for i in range(len(self.all_elements)):
            element_dof_number = self.all_elements[i].get_dof_numbers()
            element_stiffness_matrix_glob = self.all_elements[i].compute_element_k()
            for j in range(len(element_stiffness_matrix_glob[0])): # rows
                for k in range(len(element_stiffness_matrix_glob[0])): # columns
                    if ((element_dof_number[j] != -1) and (element_dof_number[k] != -1)):
                        global_stiffness_tensor[element_dof_number[j]][element_dof_number[k]] = global_stiffness_tensor[element_dof_number[j]][element_dof_number[k]] + element_stiffness_matrix_glob[j][k]
        return global_stiffness_tensor

    def assemble_global_mass_tensor(self, number_of_equations):
        """
        :param number_of_equations: int
            Number of equations (free DOF'S)
        :return: np.array
            Returns the mass matrix of the structure
        """
        global_mass_tensor = np.zeros(shape=(number_of_equations, number_of_equations))
        for i in range(len(self.all_elements)):
            element_dof_number = self.all_elements[i].get_dof_numbers()
            element_mass_matrix_glob = self.all_elements[i].compute_element_m()
            for j in range(len(global_mass_tensor[0])):  # rows
                for k in range(len(global_mass_tensor[0])):  # columns
                    if ((element_dof_number[j] != -1) and (element_dof_number[k] != -1)):
                        global_mass_tensor[element_dof_number[j]][element_dof_number[k]] += element_mass_matrix_glob[j][k]
        return global_mass_tensor

    def set_displacements_to_nodes(self, displacements):
        """ Applies the calculated displacement to all nodes
        :param displacements: np.array
            Array containg the displacements of the nodes
        :return: None
        """
        for i in range(len(self.nodes)):
            displ = np.empty(3, dtype=float)
            for j in range(len(displ)):
                if (self.nodes[i].get_dof_numbers()[j] != -1):
                    displ[j] = displacements[int(self.nodes[i].get_dof_numbers()[j])]
                else:
                    displ[j] = 0
            self.nodes[i].set_displacement(displ)

    def get_max_distance(self):
        """ Compute the max distance within the structure
        :return: float
            Max distance
        """
        distance = 0
        # Compares all nodes to the first node, which might not result in the largest distance !!!
        first_node_coords = self.nodes[0].get_coordinates()
        for i in range(len(self.nodes)):
            next_node_coords = self.nodes[i].get_coordinates()
            difference = first_node_coords - next_node_coords
            x = np.power(difference[0], 2)
            y = np.power(difference[1], 2)
            z = np.power(difference[2], 2)
            temp_value = np.sqrt(x + y +z)
            if temp_value > distance:
                distance = temp_value
        return distance


    def solve_direct_stiffness_method(self, verbose=False):
        """ Solve the system of equations using the direct stiffness method in which
        Ku=F is solved by inverting K and solving for u.
        :param verbose: bool
            Enable/Disable printing option
        :return: None
        """
        np.set_printoptions(threshold=np.inf)
        np.set_printoptions(linewidth=np.inf)
        # Get number of Equations of the system
        number_of_equations = self.get_number_of_equations()
        # Assemble global load vector
        load_vector = self.assemble_global_load_vector(number_of_equations)
        print("Global Load Vector is:") if verbose else None
        print(load_vector,"\n") if verbose else None
        # Assemble global stiffness tensor
        global_stiffness_tensor = self.assemble_global_stiffness_tensor(number_of_equations)
        print("Global Stiffness Tensor:") if verbose else None
        print(global_stiffness_tensor,"\n") if verbose else None
        # Invert stiffness tensor
        inverse_stiffness_tensor = np.linalg.inv(global_stiffness_tensor)
        # Calculate displacements
        displacements = inverse_stiffness_tensor @ load_vector
        print("\n"+"Displacements:") if verbose else None
        print(displacements,"\n") if verbose else None
        # Set dsplacement values to nodes
        self.set_displacements_to_nodes(displacements)

    def solve_explicit_euler(self, dt=None, number_timesteps=None, scale_dt=1, verbose=True):
        """ Solve the system of equations by using the explicit euler scheme with increments over the time.
        Usage of mass matrix is not implemented yet, therefore, the algorithm is equal to using the direct stiffness algorithm.
        Even though the smallest ciritcal timestep is calculated, it is still possible that the solution
        diverges. In this case, the direct stiffness method should be applied.
        :param dt: float
            time step size dt
        :param number_timesteps: int
            Number of timesteps for time discretisation
        :return:
            None
        """
        # Get number of Equations of the system
        number_of_equations = self.get_number_of_equations()
        # Get load vector and stiffness tensor
        load_vector = self.assemble_global_load_vector(number_of_equations)
        global_stiffness_tensor = self.assemble_global_stiffness_tensor(number_of_equations)
        # Create displacement vector with same shape of load_vector and fill with 0's
        displacements = np.zeros_like(load_vector)

        if dt is None:
            # Calculate min length and wave propagation velocity 
            l_min, c_speed = np.inf, np.inf
            for i in range(len(self.all_elements)):
                length = self.all_elements[i].compute_min_length()
                e_mod = self.all_elements[i].get_e_mod()
                density = self.all_elements[i].get_density()
                # Save new variables
                if length < l_min:
                    l_min = length
                if np.sqrt(e_mod/density) < c_speed:
                    c_speed = np.sqrt(e_mod/density)

            # Calculate critical time step size and scale it to ensure solution will converge
            dt = (l_min / c_speed) * scale_dt
            print("Calculated min. timestep size: ",dt)
        
        if number_timesteps is None:
            number_timesteps = 1 / dt
        
        # Solver iteratively with explicit euler
        last_printed_percentage = -1 
        for step in range(int(number_timesteps)):

            # Given that F=Ku, calculate difference between currently present internal forces and externally applied forces.
            # du/dt represents the ‘change’ of u in a time step, based on the remaining forces that are not yet balanced by the current configuration u
            du_dt = load_vector - np.dot(global_stiffness_tensor, displacements)
            # Refresh displacements with explicit euler time discretization
            # u_new = u_old + dt*(du/dt)
            displacements = displacements + dt * du_dt

            percentage_done = (step / number_timesteps) * 100
            rounded_percentage = int(percentage_done)  # Ganzzahliger Prozentwert

            # Print for each % done
            if rounded_percentage > last_printed_percentage:
                print(f"Percentage done: {rounded_percentage}%")
                last_printed_percentage = rounded_percentage 
            
        print("Final Displacement Vector:", displacements) if verbose else None
        # Set dsplacement values to nodes
        self.set_displacements_to_nodes(displacements) 

    def print_info(self, filepath="Structure_information.txt"):
        """ Print all available information of the structure
        :param filepath: str
            Path of the text file with .txt ending
        :return: None
        """
        now = datetime.now()
        with open(filepath, "w") as file:
            file.write("----------" + "\n" +
                       "File contains structural information, analysis results and information of nodes and elements"+"\n" +
                       "File was created: " + str(now) + "\n"+
                       "----------" + "\n" + "\n")
            file.write("Listing structural information" + "\n")
            file.write("Number of nodes: " + str(self.get_number_of_nodes()) + "\n")
            file.write("Number of elements: " + str(self.get_number_of_elements()) + "\n")
            file.write("   Number of truss elements: " + str(self.get_number_of_truss_elements()) + "\n")
            file.write("   Number of CST elements: " + str(self.get_number_of_CST_elements()) + "\n")
            file.write("   Number of linear hexahedral elements: " + str(self.get_number_of_linear_hexahedral_elements()) + "\n")
            file.write("\n")
            file.write("Listing analysis results \n")
            # Truss Elements
            if len(self.truss_elements) != 0.0:
                displacements_x, displacements_y, displacements_z , stress_list= [], [], [], []
                for i in range(len(self.truss_elements)):
                    truss = self.truss_elements[i]
                    displacements_x.append(truss.get_node_1().get_displacement()[0])
                    displacements_x.append(truss.get_node_2().get_displacement()[0])
                    displacements_y.append(truss.get_node_1().get_displacement()[1])
                    displacements_y.append(truss.get_node_2().get_displacement()[1])
                    displacements_z.append(truss.get_node_1().get_displacement()[2])
                    displacements_z.append(truss.get_node_2().get_displacement()[2])
                    stress_list.append(truss.compute_element_mean_stress())

                file.write("Truss elements: " + "\n") 
                file.write("   Max. displacement in x: " + str(np.max(displacements_x)) + "\n")
                file.write("   Min. displacement in x: " + str(np.min(displacements_x)) + "\n")
                file.write("   Max. displacement in y: " + str(np.max(displacements_y)) + "\n")
                file.write("   Min. displacement in y: " + str(np.min(displacements_y)) + "\n")
                file.write("   Max. displacement in z: " + str(np.max(displacements_z)) + "\n")
                file.write("   Min. displacement in z: " + str(np.min(displacements_z)) + "\n")
                file.write("   Max. Von-Mises stress: " + str(np.max(stress_list)) + "\n")
                file.write("   Min. Von-Mises stress: " + str(np.min(stress_list)) + "\n")
                file.write("   Mean. Von-Mises stress: " + str(np.mean(stress_list)) + "\n")
                stress_list.sort()
                file.write("   Median. Von-Mises stress: " + str(np.median(stress_list)) + "\n")

            # CST Elements
            if len(self.cst_elements) != 0.0:
                displacements_x, displacements_y, displacements_z , stress_list= [], [], [], []
                for i in range(len(self.cst_elements)):
                    cst = self.cst_elements[i]
                    displacements_x.append(cst.get_node_1().get_displacement()[0])
                    displacements_x.append(cst.get_node_2().get_displacement()[0])
                    displacements_x.append(cst.get_node_3().get_displacement()[0])
                    displacements_y.append(cst.get_node_1().get_displacement()[1])
                    displacements_y.append(cst.get_node_2().get_displacement()[1])
                    displacements_y.append(cst.get_node_3().get_displacement()[1])
                    displacements_z.append(cst.get_node_1().get_displacement()[2])
                    displacements_z.append(cst.get_node_2().get_displacement()[2])
                    displacements_z.append(cst.get_node_3().get_displacement()[2])
                    stress_list.append(cst.compute_element_mean_stress())

                file.write("CST elements: " + "\n") 
                file.write("   Max. displacement in x: " + str(np.max(displacements_x)) + "\n")
                file.write("   Min. displacement in x: " + str(np.min(displacements_x)) + "\n")
                file.write("   Max. displacement in y: " + str(np.max(displacements_y)) + "\n")
                file.write("   Min. displacement in y: " + str(np.min(displacements_y)) + "\n")
                file.write("   Max. displacement in z: " + str(np.max(displacements_z)) + "\n")
                file.write("   Min. displacement in z: " + str(np.min(displacements_z)) + "\n")
                file.write("   Max. Von-Mises stress: " + str(np.max(stress_list)) + "\n")
                file.write("   Min. Von-Mises stress: " + str(np.min(stress_list)) + "\n")
                file.write("   Mean. Von-Mises stress: " + str(np.mean(stress_list)) + "\n")
                stress_list.sort()
                file.write("   Median. Von-Mises stress: " + str(np.median(stress_list)) + "\n")
            # Linear Hexahedral Elements
            if len(self.linear_hexa_elements) != 0.0:
                displacements_x, displacements_y, displacements_z , stress_list= [], [], [], []
                for i in range(len(self.linear_hexa_elements)):
                    hexa = self.linear_hexa_elements[i]
                    displacements_x.append(hexa.get_node_1().get_displacement()[0])
                    displacements_x.append(hexa.get_node_2().get_displacement()[0])
                    displacements_x.append(hexa.get_node_3().get_displacement()[0])
                    displacements_x.append(hexa.get_node_4().get_displacement()[0])
                    displacements_x.append(hexa.get_node_5().get_displacement()[0])
                    displacements_x.append(hexa.get_node_6().get_displacement()[0])
                    displacements_x.append(hexa.get_node_7().get_displacement()[0])
                    displacements_x.append(hexa.get_node_8().get_displacement()[0])
                    
                    displacements_y.append(hexa.get_node_1().get_displacement()[1])
                    displacements_y.append(hexa.get_node_2().get_displacement()[1])
                    displacements_y.append(hexa.get_node_3().get_displacement()[1])
                    displacements_y.append(hexa.get_node_4().get_displacement()[1])
                    displacements_y.append(hexa.get_node_5().get_displacement()[1])
                    displacements_y.append(hexa.get_node_6().get_displacement()[1])
                    displacements_y.append(hexa.get_node_7().get_displacement()[1])
                    displacements_y.append(hexa.get_node_8().get_displacement()[1])

                    displacements_z.append(hexa.get_node_1().get_displacement()[2])
                    displacements_z.append(hexa.get_node_2().get_displacement()[2])
                    displacements_z.append(hexa.get_node_3().get_displacement()[2])
                    displacements_z.append(hexa.get_node_4().get_displacement()[2])
                    displacements_z.append(hexa.get_node_5().get_displacement()[2])
                    displacements_z.append(hexa.get_node_6().get_displacement()[2])
                    displacements_z.append(hexa.get_node_7().get_displacement()[2])
                    displacements_z.append(hexa.get_node_8().get_displacement()[2])

                    stress_list.append(hexa.compute_element_mean_stress())

                file.write("Linear hexahedral elements: " + "\n") 
                file.write("   Max. displacement in x: " + str(np.max(displacements_x)) + "\n")
                file.write("   Min. displacement in x: " + str(np.min(displacements_x)) + "\n")
                file.write("   Max. displacement in y: " + str(np.max(displacements_y)) + "\n")
                file.write("   Min. displacement in y: " + str(np.min(displacements_y)) + "\n")
                file.write("   Max. displacement in z: " + str(np.max(displacements_z)) + "\n")
                file.write("   Min. displacement in z: " + str(np.min(displacements_z)) + "\n")
                file.write("   Max. Von-Mises stress: " + str(np.max(stress_list)) + "\n")
                file.write("   Min. Von-Mises stress: " + str(np.min(stress_list)) + "\n")
                file.write("   Mean. Von-Mises stress: " + str(np.mean(stress_list)) + "\n")
                stress_list.sort()
                file.write("   Median. Von-Mises stress: " + str(np.median(stress_list)) + "\n")

            file.write("\n"+"Listing information of nodes and elements" + "\n" + "Nodes:" + "\n")
            # Node information
            for i in range(len(self.nodes)):
                node = self.nodes[i]
                coordinates = node.get_coordinates()
                file.write("   Node_ID " + str(i) + ":" + "\n" + 
                           "      Coordinates: (%.3f %.3f %.3f) " % (coordinates[0], coordinates[1], coordinates[2]) + "\n")
                displacement = node.get_displacement()
                file.write("      Displacement: (%.3e %.3e %.3e) " % (displacement[0], displacement[1], displacement[2]) + "\n")
                constraint = node.get_constraint()
                file.write(f"      Constraints: {constraint.is_free(0), constraint.is_free(1), constraint.is_free(2)}" + "\n")
                force = node.get_force()
                file.write(f"      Force: {force.get_component(0), force.get_component(1), force.get_component(2)}" + "\n")

            # Element information
            file.write("\n")
            file.write("Elements:" + "\n")
            for i in range(len(self.all_elements)):
                element = self.all_elements[i]
                nodes_list = element.get_all_nodes()

                file.write("   Element_ID " + str(i) + ":" + "\n" +
                           "      Node_ID's: (" )
                # Print node ID's
                for j, _node in enumerate(nodes_list):
                    if j < len(nodes_list)-1:
                        file.write(str(_node.get_node_index())+ ", ")
                    else:
                        file.write(str(_node.get_node_index()))

                file.write(")"+ "\n")
                file.write("      Von-Mises stress: " + str(element.compute_element_mean_stress()) + "\n")
                # Write stress tensor for CST and Hexa
                if element.get_number_of_nodes() == 8 or element.get_number_of_nodes() == 3:
                    file.write("      Element stress tensor: " + str(element.compute_element_stress()) + "\n")

                
    def create_cube_mesh(self, structure, width_x, length_y, height_z, elements_per_edge, e_mod, poisson, density, SRI=False):
        """Create a cube-shaped mesh with hexahedral elements
        :param structure: object
            Structure object to which the hexahedral elements are added
        :param width_x: float
            Width of the geometry
        :param length_y: float
            Depth of the geometry
        :param height_z: float
            Heigth of the geometry
        :param elements_per_edge: int
            Elements per edge
        :param e_mod: float
            E-Modulus/Young's Modulus of the elements
        :param poisson: float
            Poisson ratio of the elements
        :param density: float
            Density of the elements
        :param SRI: bool
            Flag to either create standard hexahedral elements or hexa's using selective reduced integration
        :return: object
            Structure object wtth the hexahedral elements
        """

        # Create nodes
        nodes = []
        for k in range (elements_per_edge+1): # z-axis
            for j in range(elements_per_edge+1): # y-axis
                for i in range(elements_per_edge +1): # x-axis
                    x = i * (width_x/elements_per_edge)
                    y = j * (length_y/elements_per_edge)
                    z = k * (height_z/(elements_per_edge))
                    nodes.append([x,y,z])
                    structure.add_node(x, y, z)
        nodes = np.array(nodes)

        # Create hexahedral elements
        elements = []
        for k in range(elements_per_edge): # z-axis
            for j in range(elements_per_edge): # y-axis
                for i in range(elements_per_edge): # x-axis
                    # Indices of the 8 corners of a hexahedral element
                    # Nummerrated clockwise, bottom side first then top
                    n0 = i + j*(elements_per_edge+1)+k*(elements_per_edge+1)**2
                    n1 = n0 +1
                    n2 = n1 + (elements_per_edge+1)
                    n3 = n0 + (elements_per_edge + 1)
                    n4 = n0 + (elements_per_edge + 1)**2
                    n5 = n1 + (elements_per_edge + 1)**2
                    n6 = n2 + (elements_per_edge + 1)**2
                    n7 = n3 + (elements_per_edge + 1)**2

                    elements.append([n0, n1, n2, n3, n4, n5, n6, n7])
                    if SRI == False:
                        structure.add_linear_hexahedral(e_mod, poisson, density, n0, n1, n2, n3, n4, n5, n6, n7)
                    else:
                        structure.add_linear_hexahedral_SRI(e_mod, poisson, density, n0, n1, n2, n3, n4, n5, n6, n7)

        elements = np.array(elements)

        return structure


    def create_tri_mesh(self, structure, width_x, length_y, elements_per_edge, e_mod, poisson, density, height_cst, mode="Plane Sress", verbose=False):
        """Create a 2D-Mesh with CST-Elements

        :param structure: object
            Structure object from Structure class to which the elements should be added
        :param width_x: float
            Width of the 2D-Mesh in x
        :param length_y: float
            Length of the 2D-Mesh in y
        :param elements_per_edge: int
            Elements per edge in x and y
        :param e_mod: float
            E-Modulus of the CST's
        :param poisson: float
            Poisson ratio of the CST'S
        :param density: float
            Density of the CST's
        :param heigth_cst: float
            Height of the CST in z-direction
        :param mode: str
            Plane Stress or Plane Strain
        """
        
        # Nodes per edge
        n_points_per_side = elements_per_edge + 1

        # X und Y coordinates of mesh
        x = np.linspace(0, width_x, n_points_per_side)
        y = np.linspace(0, length_y, n_points_per_side)

        # Meshgrid
        X, Y = np.meshgrid(x, y)

        # Change all nodes to 1D-Array
        nodes = np.vstack([X.flatten(), Y.flatten()]).T

        for node in nodes:
            structure.add_node(node[0], node[1], 0)

        # Create triangles by indices of nodes
        elements = []

        # Create 2 triangles for each quader in the mesh
        for i in range(elements_per_edge):
            for j in range(elements_per_edge):
                # Indices of the four edge nodes of a quader
                p1 = i * n_points_per_side + j  # lower left corner
                p2 = p1 + 1  # lower right corner
                p3 = p1 + n_points_per_side  # upper left corner
                p4 = p3 + 1  # upper right corner

                # Add triangles
                elements.append([p1, p2, p4])  # First triangle (lower one)
                elements.append([p1, p4, p3])  # Second triangle (upper one)

        # Numpy Array
        elements = np.array(elements)

        for element in elements:
            structure.add_CST(e_mod, poisson, height_cst, density, element[0], element[1], element[2])

        # Print number of created elements
        n_elements = elements.shape[0]
        print(f"Number of CST-Elements: {n_elements}") if verbose else None

        # Set Mode of CST
        structure.set_CST_mode(str(mode))
        
        return structure