import sys
import os
sys.path.extend([str(os.getcwd())])
from src.FEMpy_3D.Structure import Structure
from src.FEMpy_3D.Viewer import Viewer
from src.FEMpy_3D.Constraint import Constraint
from src.FEMpy_3D.Force import Force
import pyvista

class struct_class(object):

    def struct_function(self):
        # Set dimensions of the cube
        x, y, z = 20, 20, 20
        # Define elements per edge
        elements_per_edge = 10
        # Define material constants
        e_mod, poisson, density =  200000, 0.35, 7.85E-09

        # Get structure object
        s = Structure()
        # Create Mesh
        structured_mesh = s.create_cube_mesh(structure=s, width_x=x, length_y=y, height_z=z, elements_per_edge=elements_per_edge, e_mod=e_mod, poisson=poisson, density=density)

        # Apply constraints and force by identifying nodes by their coordinates
        fixed_fixed_fixed = Constraint(False, False, False)
        force = Force(0.0, 0.0, -1000.0)
        for i in range(len(structured_mesh.get_all_nodes())):
            node = structured_mesh.get_node(i)
            # Set constraint
            if node.get_coordinates()[2] == 0.0:
                structured_mesh.get_node(i).set_constraint(fixed_fixed_fixed)
            # Apply force
            if node.get_coordinates()[2] == z:
                structured_mesh.get_node(i).set_force(force)

        return structured_mesh

if __name__ == "__main__":
    # Create instance of the Class which contains the structures
    structure_class = struct_class()
    # Create Structure
    structured_mesh = structure_class.struct_function()
    # Solve the structure
    structured_mesh.solve_direct_stiffness_method(verbose=False)
    #structured_mesh.solve_explicit_euler()
    # Print all information to txt file
    structured_mesh.print_info()

    # Visualise
    # Initialise the 3D Plotter Window from pyVista
    plotter = pyvista.Plotter()
    # Initialise the viewer object
    viewer = Viewer(structured_mesh, plotter)

    # Create the mesh containg the hexa elements
    viewer.visualise_linear_hexahedral_elements()

    # # Define result field of interest for visualisation
    field = "Displacement"
    field = "Mean Stress"
    # # Draw Hexahedral Displacements
    viewer.visualise_linear_hexahedral_deformation(field)


    # # Add all constraints to the 3D Plotter Window
    viewer.visualise_constraints_by_cones()
    # # Add all forces to the 3D Plotter Window
    viewer.visualise_forces_by_arrows()

    # # Open the GUI application
    plotter.show_axes()
    plotter.title = "Cube Mesh FEA"
    plotter.camera_position = [-1, -1, 0.4]
    #viewer.visualise_point_labels()
    plotter.show()




