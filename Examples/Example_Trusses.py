import sys
import os
sys.path.extend([str(os.getcwd())])
import math
from src.FEMpy_3D.Constraint import Constraint
from src.FEMpy_3D.Force import Force
from src.FEMpy_3D.Structure import Structure
from src.FEMpy_3D.Viewer import Viewer
import pyvista

class FE_Modell():

    def create_structure(self):
        # Initialise variables
        area = math.pi * 10**2
        e_mod = 2e5
        density = 7850
        # Create structure object
        struct = Structure()

        # Create force objects
        value = -2000000.0
        force_1 = Force(0.0, value, value)

        # Create and add nodes to structure object
        n0 = struct.add_node(2.5, 0.0, 10.0)
        n1 = struct.add_node(-5.0, -5.0, 0)
        n2 = struct.add_node(10.0, -5.0, 0)
        n3 = struct.add_node(5.0, 5.0, 0)

        # Apply forces to nodes
        n0.set_force(force_1)

        # Apply BC's
        n1.set_constraint(Constraint(False, True, False))
        n2.set_constraint(Constraint(True, False, False))
        n3.set_constraint(Constraint(False, False, False))

        # Create 1-D truss-elements
        struct.add_truss(e_mod, area, density, 0, 1)
        struct.add_truss(e_mod, area, density, 0, 2)
        struct.add_truss(e_mod, area, density, 0, 3)

        struct.add_truss(e_mod, area, density, 1, 2)
        struct.add_truss(e_mod, area, density, 2, 3)
        struct.add_truss(e_mod, area, density, 3, 1)

        return struct

if __name__ == "__main__":
    # Create instance of the Class which contains the structures
    fem = FE_Modell()
    # Select one method of the class to initialise the structure
    struct = fem.create_structure()
    # Solve the system by direct stiffness method
    struct.solve_direct_stiffness_method(verbose=True)
    # Print all possible information of the structure
    struct.print_info()

    # Initialise the 3D Plotter Window from pyVista
    plotter = pyvista.Plotter()
    # Initialise the viewer object
    viewer = Viewer(struct, plotter)
    # Create the mesh containg the truss elements
    viewer.visualise_truss_elements()
    field = "Mean Stress"
    # Calculate and generate the deformed truss mesh
    viewer.visualise_truss_deformation(field)
    # Add all constraints to the 3D Plotter Window
    viewer.visualise_constraints_by_cones()
    # Add all forces to the 3D Plotter Window
    viewer.visualise_forces_by_arrows()
    # Vis point labels
    #viewer.visualise_point_labels()

    # Open the GUI application
    plotter.show_axes()
    plotter.title = "Linear Static FEA of Truss Elements using Python"
    plotter.camera_position = [-1.5, -1, 0.7]
    plotter.enable_parallel_projection()

    plotter.show()
