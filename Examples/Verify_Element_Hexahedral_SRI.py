import sys
import os
sys.path.extend([str(os.getcwd())])
from src.Constraint import Constraint
from src.Force import Force
from src.Structure import Structure
from src.Viewer import Viewer
import pyvista


class struct_class(object):

    def struct_function(self):
        e_mod = 200000.0
        poisson_ratio = 0.35
        density = 7.85E-09

        s = Structure()

        n0 = s.add_node(0.0, 0.0, 0.0)
        n1 = s.add_node(10.0, 0.0, 0.0)
        n2 = s.add_node(10.0, 10.0, 0.0)
        n3 = s.add_node(0.0, 10.0, 0.0)
        n4 = s.add_node(0.0, 0.0, 10.0)
        n5 = s.add_node(10.0, 0.0, 10.0)
        n6 = s.add_node(10.0, 10.0, 10.0)
        n7 = s.add_node(0.0, 10.0, 10.0)

        constraint_fff = Constraint(False, False, False)

        force = Force(10000.0, 0.0, 0.0)

        n1.set_force(force)
        n2.set_force(force)
        n5.set_force(force)
        n6.set_force(force)

        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 0, 1, 2, 3, 4, 5, 6, 7)

        n0.set_constraint(constraint_fff)
        n3.set_constraint(constraint_fff)
        n4.set_constraint(constraint_fff)
        n7.set_constraint(constraint_fff)

        return s

if __name__ == "__main__":
    # Create instance of the Class which contains the structures
    structure_class = struct_class()
    # Create Structure
    s = structure_class.struct_function()
    # solve
    s.solve_direct_stiffness_method(verbose=True)
    # Print all possible information of the structure to a file
    s.print_info()

    # Initialise the 3D Plotter Window from pyVista
    plotter = pyvista.Plotter()
    # Initialise the viewer object
    viewer = Viewer(s, plotter)

    #Draw Hexahedron Element
    #viewer.visualise_linear_hexahedral_elements()
    field = "Mean Stress"
    # Draw Hexahedron Deformation
    viewer.visualise_linear_hexahedral_deformation(field)
    # Add all constraints to the 3D Plotter Window
    viewer.visualise_constraints_by_cones()
    # Add all forces to the 3D Plotter Window
    viewer.visualise_forces_by_arrows()
    # Vis point labels
    viewer.visualise_point_labels()

    plotter.show_axes()
    plotter.title = "Linear Static FEA Test for a 3D Hexahedral Element"
    plotter.show()

