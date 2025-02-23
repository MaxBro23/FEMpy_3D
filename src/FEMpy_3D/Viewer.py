import numpy as np
import pyvista
from pyvista import CellType


class Viewer(object):

    def __init__(self, structure, plotter):
        """
        :param structure: object
            Structure object from Structure.py
        :param plotter: object
            Plotter object from PyVista
        """
        self.structure = structure
        self.plotter = plotter
        self.cmap = "viridis" #  tab20b  viridis  cividis

    def visualise_constraints_by_cones(self, scale=0.02):
        """Visualise all constraints by cones.
        Constraints in x in red colour, in y in green and in z in blue.
        :return: None
        """
        height = self.structure.get_max_distance()*scale
        cone_x_list = np.empty(0, dtype=object)
        cone_y_list = np.empty(0, dtype=object)
        cone_z_list = np.empty(0, dtype=object)
        for i in range(self.structure.get_number_of_nodes()):
            node = self.structure.get_node(i)
            constraint = node.get_constraint()
            coords = node.get_coordinates()
            if constraint.is_free(0) == False:
                cone_x = pyvista.Cone(center=[coords[0] - (height / 2), coords[1], coords[2]], direction=(1.0, 0, 0), height=height, radius=height/2)
                cone_x_list = np.append(cone_x_list, cone_x, axis=None)
            if constraint.is_free(1) == False:
                cone_y = pyvista.Cone(center=[coords[0], coords[1] - (height / 2), coords[2]], direction=(0.0, 1.0, 0), height=height, radius=height/2)
                cone_y_list = np.append(cone_y_list, cone_y, axis=None)
            if constraint.is_free(2) == False:
                cone_z = pyvista.Cone(center=[coords[0], coords[1], coords[2] - (height / 2)], direction=(0.0, 0.0, 1.0), height=height, radius=height/2)
                cone_z_list = np.append(cone_z_list, cone_z, axis=None)
        # Add all constraints to the 3D Plotter Window
        if len(cone_x_list) != 0:
            for cone in cone_x_list:
                self.plotter.add_mesh(cone, color='red', specular=1.0, specular_power=10)
        if len(cone_y_list) != 0:
            for cone in cone_y_list:
                self.plotter.add_mesh(cone, color='green', specular=1.0, specular_power=10)
        if len(cone_z_list) != 0:
            for cone in cone_z_list:
                self.plotter.add_mesh(cone, color='blue', specular=1.0, specular_power=10)

    def visualise_forces_by_arrows(self, scale=0.1):
        """Visualise all forces by arrows in red
        :return: None
        """
        height = self.structure.get_max_distance()*scale
        arrow_list = np.empty(0, dtype=object)
        for i in range(self.structure.get_number_of_nodes()):
            force = self.structure.get_node(i).get_force().get_force_vector()
            for j in range(len(force)):
                if force[j] != 0:
                    coords = self.structure.get_node(i).get_coordinates()
                    if j == 0:
                        if force[j] > 0:
                            arrow = pyvista.Arrow(start=[coords[0]-height, coords[1], coords[2]], direction=(1.0, 0.0, 0.0), scale=height)
                            arrow_list = np.append(arrow_list, arrow, axis=None)
                        else:
                            arrow = pyvista.Arrow(start=[coords[0]+height, coords[1], coords[2]], direction=(-1.0, 0.0, 0.0), scale=height)
                            arrow_list = np.append(arrow_list, arrow, axis=None)
                    elif j == 1:
                        if force[j] > 0:
                            arrow = pyvista.Arrow(start=[coords[0], coords[1]-height, coords[2]], direction=(0.0, 1.0, 0.0), scale=height)
                            arrow_list = np.append(arrow_list, arrow, axis=None)
                        else:
                            arrow = pyvista.Arrow(start=[coords[0], coords[1]+height, coords[2]], direction=(0.0, -1.0, 0.0), scale=height)
                            arrow_list = np.append(arrow_list, arrow, axis=None)
                    else:
                        if force[j] > 0:
                            arrow = pyvista.Arrow(start=[coords[0], coords[1], coords[2]-height], direction=(0.0, 0.0, 1.0),scale=height)
                            arrow_list = np.append(arrow_list, arrow, axis=None)
                        else:
                            arrow = pyvista.Arrow(start=[coords[0], coords[1], coords[2]+height], direction=(0.0, 0.0, -1.0), scale=height)
                            arrow_list = np.append(arrow_list, arrow, axis=None)
        # Add all forces to the 3D Plotter Window
        if len(arrow_list) != 0:
            for arrow in arrow_list:
                self.plotter.add_mesh(arrow, color="red")
        else:
            print("No Forces applied to the structure")

    def visualise_truss_elements(self):
        """Visualise all truss elements
        :return: None
        """
        # Get number of truss elements
        edges = self.structure.get_number_of_truss_elements()
        # Initialise holder for node indices
        edge_node_list = np.empty(0, dtype=list)
        # Get the indices for all nodes that connect to truss elements
        for i in range(edges):
            truss_node_idx_1 = self.structure.get_truss_element(i).node_1_idx
            truss_node_idx_2 = self.structure.get_truss_element(i).node_2_idx
            edge_node_list = np.append(edge_node_list, [truss_node_idx_1, truss_node_idx_2], axis=None)
        # Reshape it from a 1D array to a 2D array where every inner array contains the nodes of one truss
        edge_node_list_reshaped = np.reshape(edge_node_list, (-1, 2))
        # Get a list of the coordinates of all nodes
        nodes = self.structure.get_coordinates_of_nodes_as_list()
        # We must "pad" the edges to indicate to vtk how many points per edge
        padding = np.empty(edges, dtype=int)
        padding[:] = 2
        # Build an array such that [[Num of Nodes per Element, Idx Node 1, Idx Node 2],...]
        edges_w_padding = np.vstack((padding, edge_node_list_reshaped.T)).T
        # Has to be of an integer type otherwise a bug will occur
        # See also: https://github.com/pyvista/pyvista/issues/810
        edges_w_padding = edges_w_padding.astype('int32')
        # Create the Poly Data
        mesh = pyvista.PolyData(nodes, edges_w_padding)
        # Add the truss element mesh to the Plotter
        self.plotter.add_mesh(mesh, render_lines_as_tubes=True, style='wireframe', line_width=10, color='grey', opacity=0.4)

    def visualise_truss_deformation(self, field):
        """Visualise all truss elements in deformed configuration
        :param field: str
            Field name for colouring the elements e.g. "Displacement" colours all elements based on their displacement
        :return: None
        """
        # Get number of truss elements
        edges = self.structure.get_number_of_truss_elements()
        # Initialise holder for node indices
        edge_node_list = np.empty(0, dtype=list)
        # Get the indices for all nodes that connect to truss elements
        for i in range(edges):
            truss_node_idx_1 = self.structure.get_truss_element(i).node_1_idx
            truss_node_idx_2 = self.structure.get_truss_element(i).node_2_idx
            edge_node_list = np.append(edge_node_list, [truss_node_idx_1, truss_node_idx_2], axis=None)
        # Reshape it from a 1D array to a 2D array where every inner array contains the nodes of one truss
        edge_node_list_reshaped = np.reshape(edge_node_list, (-1, 2))
        # Get a list of the coordinates of all nodes
        nodes = self.structure.get_displaced_coordiantes_of_nodes_as_list()
        # We must "pad" the edges to indicate to vtk how many points per edge
        padding = np.empty(edges, dtype=int)
        padding[:] = 2
        # Build an array such that [[Num of Nodes per Element, Idx Node 1, Idx Node 2],...]
        edges_w_padding = np.vstack((padding, edge_node_list_reshaped.T)).T
        # Has to be of an integer type otherwise a bug will occur
        # See also: https://github.com/pyvista/pyvista/issues/810
        edges_w_padding = edges_w_padding.astype('int32')
        # Create the Poly Data
        mesh = pyvista.PolyData(nodes, edges_w_padding)
        # Add the deformed truss element mesh to the Plotter
        field_list = []
        for i in range(self.structure.get_number_of_truss_elements()):
            element = self.structure.get_truss_element(i)
            if field == "Displacement":
                field_value = element.compute_displacement_magnitude()
            elif field == "Mean Stress":
                field_value = element.compute_element_mean_stress()
            field_list.append(field_value)
        #mean_stress_list = np.array(mean_stress_list) / self.structure.get_max_stress()
        #print("MEANNNNN", mean_stress_list)
        self.plotter.add_mesh(mesh, render_lines_as_tubes=True, style='wireframe', line_width=10, show_scalar_bar=True, scalars=field_list, cmap=self.cmap)#   , color="red"

    def visualise_CST_elements(self):
        """Visualise all cst elements
        :return: None
        """
        # Get coordinates from nodes
        node_coordinates_list = self.structure.cst_nodes_coordinates_list
        # Create empty array which holds the indices that are used to pick the correct nodes from node_coordinates_list
        faces = []
        index_0 = 0
        index_1 = 1
        index_2 = 2
        for i in range(len(self.structure.cst_elements)):
            index_array = [3, index_0, index_1, index_2]
            faces.append(index_array)
            index_0 += 3
            index_1 += 3
            index_2 += 3
        mesh = pyvista.PolyData(node_coordinates_list, faces)
        #self.plotter.add_point_labels(mesh.points, range(mesh.n_points))
        self.plotter.add_mesh(mesh, show_edges=True, color='grey', opacity=0.4)

    def visualise_CST_deformation(self, field):
        """Visualise all cst elements in deformed configuration
        :param field: str
            Field name for colouring the elements e.g. "Displacement" colours all elements based on their displacement
        :return: None
        """
        # Get a list of the coordinates of all nodes
        nodes = self.structure.get_displaced_coordinates_of_CST_nodes_as_list()
        # Create empty array which holds the indices that are used to pick the correct nodes from node_coordinates_list
        faces = []
        index_0 = 0
        index_1 = 1
        index_2 = 2
        for i in range(len(self.structure.cst_elements)):
            index_array = [3, index_0, index_1, index_2]
            faces.append(index_array)
            index_0 += 3
            index_1 += 3
            index_2 += 3
        mesh = pyvista.PolyData(nodes, faces)
        # Calculate stresses for color map
        field_list = []
        for i in range(self.structure.get_number_of_CST_elements()):
            element = self.structure.get_cst_element(i)
            if field == "Displacement":
                field_value = element.compute_displacement_magnitude()
            elif field == "Mean Stress":
                field_value = element.compute_element_mean_stress()
            field_list.append(field_value)
        #mean_stress_list = np.array(mean_stress_list) / self.structure.get_max_stress()
        self.plotter.add_mesh(mesh, show_edges=True, show_scalar_bar=True, scalars=field_list, cmap=self.cmap)#

    def visualise_linear_hexahedral_elements(self):
        """Visualise all linear hexahedral elements
        :return: None
        """
        # Get coordinates from nodes
        node_coordinates_list = self.structure.linear_hexa_nodes_coordinates_list
        # Create empty array which holds the indices that are used to pick the correct nodes from node_coordinates_list
        faces = []
        index_0 = 0
        index_1 = 1
        index_2 = 2
        index_3 = 3
        index_4 = 4
        index_5 = 5
        index_6 = 6
        index_7 = 7
        for i in range(self.structure.get_number_of_linear_hexahedral_elements()):
            index_array = [8, index_0, index_1, index_2, index_3, index_4, index_5, index_6, index_7]
            faces.append(index_array)
            index_0 += 8
            index_1 += 8
            index_2 += 8
            index_3 += 8
            index_4 += 8
            index_5 += 8
            index_6 += 8
            index_7 += 8
        # create the unstructured grid directly from the numpy arrays
        cell_type = []
        for i in range(self.structure.get_number_of_linear_hexahedral_elements()):
            cell_type.append(CellType.HEXAHEDRON)
        mesh = pyvista.UnstructuredGrid(faces, cell_type, node_coordinates_list)
        self.plotter.add_mesh(mesh, show_edges=True, color='lightgrey', pickable=True, opacity=1.0) #

    def visualise_linear_hexahedral_deformation(self, field):
        """Visualise all hexahedral elements in deformed configuration
        :param field: str
            Field name for colouring the elements e.g. "Displacement" colours all elements based on their displacement
        :return: None
        """
        # Get coordinates from nodes
        node_coordinates_list = self.structure.get_displaced_coordinates_of_linear_hexahedral_nodes_as_list()
        # Create empty array which holds the indices that are used to pick the correct nodes from node_coordinates_list
        faces = []
        index_0 = 0
        index_1 = 1
        index_2 = 2
        index_3 = 3
        index_4 = 4
        index_5 = 5
        index_6 = 6
        index_7 = 7
        for i in range(self.structure.get_number_of_linear_hexahedral_elements()):
            index_array = [8, index_0, index_1, index_2, index_3, index_4, index_5, index_6, index_7]
            faces.append(index_array)
            index_0 += 8
            index_1 += 8
            index_2 += 8
            index_3 += 8
            index_4 += 8
            index_5 += 8
            index_6 += 8
            index_7 += 8
        # create the unstructured grid directly from the numpy arrays
        cell_type = []
        for i in range(self.structure.get_number_of_linear_hexahedral_elements()):
            cell_type.append(CellType.HEXAHEDRON)
        mesh = pyvista.UnstructuredGrid(faces, cell_type ,node_coordinates_list)
        field_list = []
        for i in range(self.structure.get_number_of_linear_hexahedral_elements()):
            element = self.structure.get_linear_hexa_element(i)
            if field == "Displacement":
                field_value = element.compute_displacement_magnitude()
            elif field == "Mean Stress":
                field_value = element.compute_element_mean_stress()
            field_list.append(field_value)
        self.plotter.add_mesh(mesh,  show_edges=True, show_scalar_bar=True, scalars=field_list, cmap=self.cmap) # gist_rainbow

    def visualise_point_labels(self):
        """Visualise the index of each node
        :return: None
        """
        # Also interesting:
        # https://tutorial.pyvista.org/tutorial/03_figures/bonus/e_labels.html
        element_list = self.structure.get_all_elements()
        nodes_coordinates = []
        nodes_index_list = []
        for element in element_list:
            for node in element.get_all_nodes():
                nodes_coordinates.append(node.get_coordinates())
                node_index = node.get_node_index()
                if node_index not in nodes_index_list:
                    nodes_index_list.append(node_index)
                else:
                    nodes_index_list.append('')
        self.plotter.add_point_labels(nodes_coordinates, nodes_index_list)


