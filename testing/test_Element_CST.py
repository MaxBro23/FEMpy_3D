import sys
import os
sys.path.extend([str(os.getcwd())])
import numpy as np
from src.FEMpy_3D.Constraint import Constraint
from src.FEMpy_3D.Force import Force
from src.FEMpy_3D.Structure import Structure


def test_build_structure():

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

        s.solve_direct_stiffness_method(verbose=False)

        # Check input
        assert(s.get_cst_element(0).get_height() == 0.15)
        assert(s.get_cst_element(0).get_e_mod() == 200000)
        assert(s.get_cst_element(0).get_density() == 7.85E-09)
        
        # Checksture components
        assert(s.get_number_of_nodes() == 4)
        assert(s.get_number_of_elements() == 2)
        assert(s.get_number_of_truss_elements() == 0)
        assert(s.get_number_of_CST_elements() == 2)
        assert(s.get_number_of_linear_hexahedral_elements() == 0)

        # Check results
        displacements_x, displacements_y, displacements_z , stress_list= [], [], [], []
        for i in range(len(s.cst_elements)):
            cst = s.cst_elements[i]
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

        assert(np.isclose(np.max(displacements_x), 0.14624999999999996, rtol=1e-05))
        assert(np.isclose(np.min(displacements_x), 0.0, rtol=1e-05))
        assert(np.isclose(np.max(displacements_y), 0.0, rtol=1e-05))
        assert(np.isclose(np.min(displacements_y), -0.14625, rtol=1e-05))
        assert(np.isclose(np.max(displacements_z), 0.0, rtol=1e-05))
        assert(np.isclose(np.min(displacements_z), 0.0, rtol=1e-05))
        assert(np.isclose(np.max(stress_list), 1171.8930554164629, rtol=1e-05))
        assert(np.isclose(np.min(stress_list), 0.0, rtol=1e-05))
        assert(np.isclose(np.mean(stress_list), 585.9465277082314, rtol=1e-05))
        assert(np.isclose(np.median(stress_list), 585.9465277082314, rtol=1e-05))