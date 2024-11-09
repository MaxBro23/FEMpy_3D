import sys
import os
sys.path.extend([str(os.getcwd())])
from src.Constraint import Constraint
from src.Force import Force
from src.Structure import Structure
from src.Viewer import Viewer
import pyvista

class FE_Modell(object):

    def structure(self):
        # Create structure Object
        s = Structure()

        # CST Properties
        e_mod = 200000.0
        height_cst = 0.2
        poisson_ratio = 0.35
        density = 7.85E-09

        # Define dimensions of rectangle
        width = 10
        height = 10
        n_elements_per_side = 10

        # Create Force
        # 10000 / nodes_per_edge = Force per Node
        force = Force(5000.0, -5000.0, 0.0)

        # Define constraints
        fixed_fixed_fixed = Constraint(False, False, False)

        s = s.create_tri_mesh(s, width, height, n_elements_per_side, e_mod, poisson_ratio, density, height_cst, "Plane Stress")

        nodes = s.get_all_nodes()
        # Set Constraints based on coordinates
        for i, node in enumerate(nodes):
            # Force
            if node.get_coordinates()[0] == 0 and node.get_coordinates()[1] == height:
                s.get_node(i).set_force(force)
            # Constraint
            if node.get_coordinates()[1] == 0:
                s.get_node(i).set_constraint(fixed_fixed_fixed)
            if node.get_coordinates()[0] == width:
                s.get_node(i).set_constraint(fixed_fixed_fixed)

        return s


if __name__ == "__main__":
    # Create instance of the Class which contains the structures
    fem = FE_Modell()
    # Create Structure
    s = fem.structure()
    
    # Solve
    s.solve_direct_stiffness_method()
    #s.solve_explicit_euler(verbose=False)

    # Print all possible information of the structure
    s.print_info()

    # Initialise the 3D Plotter Window from pyVista
    plotter = pyvista.Plotter()
    # Initialise the viewer object
    viewer = Viewer(s, plotter)

    # Draw CST Elements
    viewer.visualise_CST_elements()

    field = "Mean Stress"
    #field = "Displacement"
    # Draw CST Displacements
    viewer.visualise_CST_deformation(field)


    # Add all constraints to the 3D Plotter Window
    viewer.visualise_constraints_by_cones()
    # Add all forces to the 3D Plotter Window
    viewer.visualise_forces_by_arrows()
    # Show all Node ID's
    #viewer.visualise_point_labels()

    plotter.show_axes()
    plotter.title = "CST Analysis"
    plotter.camera_position = [0.0, -0.5, 10]

    plotter.show()


