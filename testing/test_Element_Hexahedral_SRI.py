import sys
import os
sys.path.extend([str(os.getcwd())])
import numpy as np
from src.FEMpy_3D.Constraint import Constraint
from src.FEMpy_3D.Force import Force
from src.FEMpy_3D.Structure import Structure



def test_build_structure():
        
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

    n8 = s.add_node(20.0, 0.0, 0.0)
    n9 = s.add_node(20.0, 10.0, 0.0)
    n10 = s.add_node(20.0, 10.0, 10.0)
    n11 = s.add_node(20.0, 0.0, 10.0)

    constraint_fff = Constraint(False, False, False)

    force = Force(10000.0, 500.0, 0.0)

    n8.set_force(force)
    n9.set_force(force)
    n10.set_force(force)
    n11.set_force(force)

    s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 0, 1, 2, 3, 4, 5, 6, 7)
    s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 1, 8, 9, 2, 5, 11, 10, 6)

    n0.set_constraint(constraint_fff)
    n3.set_constraint(constraint_fff)
    n4.set_constraint(constraint_fff)
    n7.set_constraint(constraint_fff)

    # solve
    s.solve_direct_stiffness_method(verbose=False)

    # Check input
    assert(s.get_linear_hexa_element(0).get_poisson() == 0.35)
    assert(s.get_linear_hexa_element(0).get_e_mod() == 200000.0)
    assert(s.get_linear_hexa_element(0).get_density() == 7.85E-09)
    
    # Checksture components
    assert(s.get_number_of_nodes() == 12)
    assert(s.get_number_of_elements() == 2)
    assert(s.get_number_of_truss_elements() == 0)
    assert(s.get_number_of_CST_elements() == 0)
    assert(s.get_number_of_linear_hexahedral_elements() == 2)

    # Check results
    displacements_x, displacements_y, displacements_z , stress_list= [], [], [], []
    for i in range(len(s.linear_hexa_elements)):
        hexa = s.linear_hexa_elements[i]
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

    assert(np.isclose(np.max(displacements_x), 0.05349611364166054, rtol=1e-05))
    assert(np.isclose(np.min(displacements_x), 0.0, rtol=1e-05))
    assert(np.isclose(np.max(displacements_y), 0.04382577858409017, rtol=1e-05))
    assert(np.isclose(np.min(displacements_y), 0.0, rtol=1e-05))
    assert(np.isclose(np.max(displacements_z), 0.006767311382999561, rtol=1e-05))
    assert(np.isclose(np.min(displacements_z), -0.00676731138299959, rtol=1e-05))
    assert(np.isclose(np.max(stress_list), 387.8418998233544, rtol=1e-05))
    assert(np.isclose(np.min(stress_list), 352.1175750696611, rtol=1e-05))
    assert(np.isclose(np.mean(stress_list), 369.9797374465078, rtol=1e-05))
    assert(np.isclose(np.median(stress_list), 369.9797374465078, rtol=1e-05))