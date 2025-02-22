import sys
import os
sys.path.extend([str(os.getcwd())])
import math
import numpy as np
from src.FEMpy_3D.Constraint import Constraint
from src.FEMpy_3D.Force import Force
from src.FEMpy_3D.Structure import Structure

# Test structure from Example_Trusses.py
def test_create_structure():
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

    struct.solve_direct_stiffness_method(verbose=False)

    # Check input
    assert(struct.get_truss_element(0).get_area() == math.pi * 10**2)
    assert(struct.get_truss_element(0).get_e_mod() == 2e5)
    assert(struct.get_truss_element(0).get_density() == 7850)
    
    # Check structure components
    assert(struct.get_number_of_nodes() == 4)
    assert(struct.get_number_of_elements() == 6)
    assert(struct.get_number_of_truss_elements() == 6)
    assert(struct.get_number_of_CST_elements() == 0)
    assert(struct.get_number_of_linear_hexahedral_elements() == 0)

    # Check results from .txt file
    displacements_x, displacements_y, displacements_z , stress_list= [], [], [], []
    for i in range(len(struct.truss_elements)):
        truss = struct.truss_elements[i]
        displacements_x.append(truss.get_node_1().get_displacement()[0])
        displacements_x.append(truss.get_node_2().get_displacement()[0])
        displacements_y.append(truss.get_node_1().get_displacement()[1])
        displacements_y.append(truss.get_node_2().get_displacement()[1])
        displacements_z.append(truss.get_node_1().get_displacement()[2])
        displacements_z.append(truss.get_node_2().get_displacement()[2])
        stress_list.append(truss.compute_element_mean_stress())

    assert(np.isclose(np.max(displacements_x), 0.23528257776171488, rtol=1e-05))
    assert(np.isclose(np.min(displacements_x), 0.0, rtol=1e-05))
    assert(np.isclose(np.max(displacements_y), 0.0, rtol=1e-05))
    assert(np.isclose(np.min(displacements_y), -1.01109266830992, rtol=1e-05))
    assert(np.isclose(np.max(displacements_z), 0.0, rtol=1e-05))
    assert(np.isclose(np.min(displacements_z), -0.24025456298651515, rtol=1e-05))
    assert(np.isclose(np.max(stress_list), 7142.296653499749, rtol=1e-05))
    assert(np.isclose(np.min(stress_list), 1882.2606220937191, rtol=1e-05))
    assert(np.isclose(np.mean(stress_list), 4087.2079820219137, rtol=1e-05))
    assert(np.isclose(np.median(stress_list), 3391.899453273931, rtol=1e-05))


    
    



