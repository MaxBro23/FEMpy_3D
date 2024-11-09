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

        # Define constraints
        fixed_fixed_fixed = Constraint(False, False, False)

        # Create Force
        force = Force(5000.0, -5000.0, 0.0)

        # CST Properties
        e_mod = 200000
        height = 0.15
        poisson_ratio = 0.35
        density = 7.85E-09

        #Add 4 Nodes
        s.add_node(0.0, 0.0, 0.0)
        s.add_node(50.0, 0.0, 0.0)
        s.add_node(50.0, 50.0, 0.0)
        s.add_node(0.0, 50.0, 0.0)

        # Add CST Elements
        s.add_CST(e_mod, poisson_ratio, height, density, 0, 1, 2)
        s.add_CST(e_mod, poisson_ratio, height, density, 0, 2, 3)

        # Set Mode of CST
        #s.set_CST_mode("Plane Strain")
        s.set_CST_mode("Plane Stress")


        # Set Constraints
        s.get_node(0).set_constraint(fixed_fixed_fixed)
        s.get_node(1).set_constraint(fixed_fixed_fixed)
        s.get_node(2).set_constraint(fixed_fixed_fixed)

        # Set Force
        s.get_node(3).set_force(force)

        return s


if __name__ == "__main__":
    # Create instance of the Class which contains the structures
    fem = FE_Modell()
    # Create Structure
    s = fem.structure()
    # Solve
    s.solve_direct_stiffness_method(verbose=True)

    # Print all possible information of the structure
    s.print_info()

    # Initialise the 3D Plotter Window from pyVista
    plotter = pyvista.Plotter()
    # Initialise the viewer object
    viewer = Viewer(s, plotter)

    # Draw CST Elements
    viewer.visualise_CST_elements()

    field = "Mean Stress"
    # Draw CST Displacements
    viewer.visualise_CST_deformation(field)


    # Add all constraints to the 3D Plotter Window
    viewer.visualise_constraints_by_cones()
    # Add all forces to the 3D Plotter Window
    viewer.visualise_forces_by_arrows()
    # Show all Node ID's
    viewer.visualise_point_labels()

    plotter.show_axes()
    plotter.title = "CST Analysis"
    plotter.camera_position = [1, -1, 0.4]

    plotter.show()