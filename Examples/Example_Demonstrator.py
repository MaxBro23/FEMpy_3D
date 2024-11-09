"""Technology demonstrator for using multiple element types together.
The structure is built by manually adding all nodes and elements. 
Automated structure generation algorithms are available in Structure.py. 
"""
import sys
import os
sys.path.extend([str(os.getcwd())])
from src.Constraint import Constraint
from src.Force import Force
from src.Structure import Structure
from src.Viewer import Viewer
import pyvista

class Bridge(object):

    def bridge(self):
        e_mod = 250000
        poisson_ratio = 0.35
        density = 7.85E-09

        s = Structure()

        # Define constraints
        fixed_fixed_fixed = Constraint(False, False, False)
        free_fixed_fixed = Constraint(True, False, False)
        free_fixed_free = Constraint(True, False, True)

        # Create the first 8 nodes of the first element manually
        s.add_node(0.0, 0.0, 0.0)
        s.add_node(1.0, 0.0, 0.0)
        s.add_node(1.0, 1.0, 0.0)
        s.add_node(0.0, 1.0, 0.0)
        s.add_node(0.0, 0.0, 1.0)
        s.add_node(1.0, 0.0, 1.0)
        s.add_node(1.0, 1.0, 1.0)
        s.add_node(0.0, 1.0, 1.0)
        # Creates nodes for 20 Hexa-Elements in 1st row
        for i in range(2, 21, 1):
            s.add_node(i, 0.0, 0.0)
            s.add_node(i, 1.0, 0.0)
            s.add_node(i, 0.0, 1.0)
            s.add_node(i, 1.0, 1.0)
        # Create first two elements manually
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 0, 1, 2, 3, 4, 5, 6, 7)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 1, 8, 9, 2, 5, 10, 11, 6)
        # Create Elements
        for i in range(8, 80, 4):
            s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, i, i + 4, i + 5, i + 1, i + 2, i + 6, i + 7, i + 3)

        # Create second row of hexa elements
        s.add_node(1, 2, 0)
        s.add_node(0, 2, 0)
        s.add_node(1, 2, 1)
        s.add_node(0, 2, 1)
        s.add_node(2, 2, 0)
        s.add_node(2, 2, 1)
        s.add_node(3, 2, 0)
        s.add_node(3, 2, 1)
        # Create more nodes automatically
        for i in range(3, 20, 1):
            s.add_node(i + 1, 2, 0)
            s.add_node(i + 1, 2, 1)
        # Add the first three elements manually
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 3, 2, 84, 85, 7, 6, 86, 87)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 2, 9, 88, 84, 6, 11, 89, 86)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 9, 13, 90, 88, 11, 15, 91, 89)
        # Automate element adding
        increment_2 = 2
        increment_4 = 4
        for i in range(13, 81, 4):
            s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, i, i + 4, 90 + increment_2, 88 + increment_2,
                                    11 + increment_4, 15 + increment_4, 91 + increment_2, 89 + increment_2)
            increment_2 += 2
            increment_4 += 4

        # Create third row of hexa elements
        s.add_node(1, 3, 0)
        s.add_node(0, 3, 0)
        s.add_node(1, 3, 1)
        s.add_node(0, 3, 1)
        s.add_node(2, 3, 0)
        s.add_node(2, 3, 1)
        s.add_node(3, 3, 0)
        s.add_node(3, 3, 1)
        # Create more nodes automatically
        for i in range(3, 20, 1):
            s.add_node(i + 1, 3, 0)
            s.add_node(i + 1, 3, 1)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 85, 84, 126, 127, 87, 86, 128, 129)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 84, 88, 130, 126, 86, 89, 131, 128)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 88, 90, 132, 130, 89, 91, 133, 131)
        # Automate element adding
        increment_2 = 2
        for i in range(90, 124, 2):
            s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, i, i + 2, 132 + increment_2, 130 + increment_2,
                                    89 + increment_2, 91 + increment_2, 133 + increment_2, 131 + increment_2)
            increment_2 += 2

        # Create fourth row of hexa elements
        s.add_node(1, 4, 0)
        s.add_node(0, 4, 0)
        s.add_node(1, 4, 1)
        s.add_node(0, 4, 1)
        s.add_node(2, 4, 0)
        s.add_node(2, 4, 1)
        s.add_node(3, 4, 0)
        s.add_node(3, 4, 1)
        # Create more nodes automatically
        for i in range(3, 20, 1):
            s.add_node(i + 1, 4, 0)
            s.add_node(i + 1, 4, 1)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 127, 126, 168, 169, 129, 128, 170, 171)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 126, 130, 172, 168, 128, 131, 173, 170)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 130, 132, 174, 172, 131, 133, 175, 173)
        increment_2 = 2
        for i in range(132, 166, 2):
            s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, i, i + 2, 174 + increment_2, 172 + increment_2,
                                    131 + increment_2, 133 + increment_2, 175 + increment_2, 173 + increment_2)
            increment_2 += 2

        # Create fifth row of hexa elements
        s.add_node(1, 5, 0)
        s.add_node(0, 5, 0)
        s.add_node(1, 5, 1)
        s.add_node(0, 5, 1)
        s.add_node(2, 5, 0)
        s.add_node(2, 5, 1)
        s.add_node(3, 5, 0)
        s.add_node(3, 5, 1)
        # Create more nodes automatically
        for i in range(3, 20, 1):
            s.add_node(i + 1, 5, 0)
            s.add_node(i + 1, 5, 1)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 169, 168, 210, 211, 171, 170, 212, 213)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 168, 172, 214, 210, 170, 173, 215, 212)
        s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, 172, 174, 216, 214, 173, 175, 217, 215)
        increment_2 = 2
        for i in range(174, 208, 2):
            s.add_linear_hexahedral_SRI(e_mod, poisson_ratio, density, i, i + 2, 216 + increment_2, 214 + increment_2,
                                    173 + increment_2, 175 + increment_2, 217 + increment_2, 215 + increment_2)
            increment_2 += 2

        """ Add CST Elements """
        e_mod = 250000
        height = 0.00001

        # Node Nr 252
        s.add_node(0.0, 0.0, -5.0)
        # Nr 253
        s.add_node(1.0, 0.0, -5.0)
        # Nr 254
        s.add_node(1.0, 1.0, -5.0)
        s.add_CST(e_mod, poisson_ratio, height, density, 252, 253, 254)
        # Nr 255
        s.add_node(0.0, 1.0, -5.0)
        s.add_CST(e_mod, poisson_ratio, height, density, 252, 254, 255)
        # Create more nodes automatically
        for i in range(2, 21, 1):
            s.add_node(i, 0, -5)
            s.add_node(i, 1, -5)
        s.add_CST(e_mod, poisson_ratio, height, density, 253, 256, 257)
        s.add_CST(e_mod, poisson_ratio, height, density, 253, 257, 254)
        s.add_CST(e_mod, poisson_ratio, height, density, 256, 258, 259)
        s.add_CST(e_mod, poisson_ratio, height, density, 256, 259, 257)
        increment_2 = 0
        for i in range(258, 292, 2):
            s.add_CST(e_mod, poisson_ratio, height, density, i, 260 + increment_2, 261 + increment_2)
            s.add_CST(e_mod, poisson_ratio, height, density, i, 261 + increment_2, 259 + increment_2)
            increment_2 += 2
        # Node Nr 294
        s.add_node(1.0, 2.0, -5.0)
        # Nr 295
        s.add_node(0.0, 2.0, -5.0)
        s.add_CST(e_mod, poisson_ratio, height, density, 255, 254, 294)
        s.add_CST(e_mod, poisson_ratio, height, density, 255, 294, 295)
        # Nr 296
        s.add_node(2.0, 2.0, -5.0)
        s.add_CST(e_mod, poisson_ratio, height, density, 254, 257, 296)
        s.add_CST(e_mod, poisson_ratio, height, density, 254, 296, 294)
        # Nr 297
        s.add_node(3.0, 2.0, -5.0)
        s.add_CST(e_mod, poisson_ratio, height, density, 257, 259, 297)
        s.add_CST(e_mod, poisson_ratio, height, density, 257, 297, 296)
        # Nr 298
        s.add_node(4.0, 2.0, -5.0)
        s.add_CST(e_mod, poisson_ratio, height, density, 259, 261, 298)
        s.add_CST(e_mod, poisson_ratio, height, density, 259, 298, 297)
        # Nr 299
        for i in range(5, 21, 1):
            s.add_node(i, 2, -5)
        increment_2 = 2
        increment_1 = 1
        for i in range(261, 292, 2):
            s.add_CST(e_mod, poisson_ratio, height, density, i, 261 + increment_2, 298 + increment_1)
            s.add_CST(e_mod, poisson_ratio, height, density, i, 298 + increment_1, 297 + increment_1)
            increment_2 += 2
            increment_1 += 1

        # third row
        # Nr 315
        for i in range(0, 21, 1):
            s.add_node(i, 3, -5)
        s.add_CST(e_mod, poisson_ratio, height, density, 295, 294, 316)
        s.add_CST(e_mod, poisson_ratio, height, density, 295, 315, 316)

        s.add_CST(e_mod, poisson_ratio, height, density, 294, 296, 317)
        s.add_CST(e_mod, poisson_ratio, height, density, 294, 317, 316)

        s.add_CST(e_mod, poisson_ratio, height, density, 296, 297, 318)
        s.add_CST(e_mod, poisson_ratio, height, density, 296, 318, 317)

        s.add_CST(e_mod, poisson_ratio, height, density, 297, 298, 319)
        s.add_CST(e_mod, poisson_ratio, height, density, 297, 319, 318)

        s.add_CST(e_mod, poisson_ratio, height, density, 298, 299, 320)
        s.add_CST(e_mod, poisson_ratio, height, density, 298, 320, 319)

        increment_1 = 0
        for i in range(299, 314, 1):
            s.add_CST(e_mod, poisson_ratio, height, density, i, 300 + increment_1, 321 + increment_1)
            s.add_CST(e_mod, poisson_ratio, height, density, i, 321 + increment_1, 320 + increment_1)
            increment_1 += 1

        # fourth row of tetras
        # Node Nr 335
        for i in range(0, 21, 1):
            s.add_node(i, 4, -5)
        increment_1 = 0
        for i in range(315, 335, 1):
            s.add_CST(e_mod, poisson_ratio, height, density, i, 316 + increment_1, 337 + increment_1)
            s.add_CST(e_mod, poisson_ratio, height, density, i, 337 + increment_1, 336 + increment_1)
            increment_1 += 1
        # Fifth row of tetras
        # Node Nr 357
        for i in range(0, 21, 1):
            s.add_node(i, 5, -5)
        increment_1 = 0
        for i in range(336, 356, 1):
            s.add_CST(e_mod, poisson_ratio, height, density, i, 337 + increment_1, 358 + increment_1)
            s.add_CST(e_mod, poisson_ratio, height, density, i, 358 + increment_1, 357 + increment_1)
            increment_1 += 1

        """  Add Truss Elements """
        area = 0.01
        density = 7.85E-09
        e_mod = 250000
        #
        s.add_truss(e_mod, area, density, 0, 262)
        s.add_truss(e_mod, area, density, 3, 263)
        s.add_truss(e_mod, area, density, 85, 299)
        s.add_truss(e_mod, area, density, 127, 320)
        s.add_truss(e_mod, area, density, 169, 341)
        s.add_truss(e_mod, area, density, 211, 362)

        s.add_truss(e_mod, area, density, 40, 262)
        s.add_truss(e_mod, area, density, 41, 263)
        s.add_truss(e_mod, area, density, 104, 299)
        s.add_truss(e_mod, area, density, 146, 320)
        s.add_truss(e_mod, area, density, 188, 341)
        s.add_truss(e_mod, area, density, 230, 362)

        s.add_truss(e_mod, area, density, 40, 282)
        s.add_truss(e_mod, area, density, 41, 283)
        s.add_truss(e_mod, area, density, 104, 309)
        s.add_truss(e_mod, area, density, 146, 330)
        s.add_truss(e_mod, area, density, 188, 351)
        s.add_truss(e_mod, area, density, 230, 372)

        s.add_truss(e_mod, area, density, 80, 282)
        s.add_truss(e_mod, area, density, 81, 283)
        s.add_truss(e_mod, area, density, 124, 309)
        s.add_truss(e_mod, area, density, 166, 330)
        s.add_truss(e_mod, area, density, 208, 351)
        s.add_truss(e_mod, area, density, 250, 372)

        # Left Side
        s.get_node(0).set_constraint(fixed_fixed_fixed)
        s.get_node(3).set_constraint(fixed_fixed_fixed)
        s.get_node(85).set_constraint(fixed_fixed_fixed)
        s.get_node(127).set_constraint(fixed_fixed_fixed)
        s.get_node(169).set_constraint(fixed_fixed_fixed)
        s.get_node(211).set_constraint(fixed_fixed_fixed)
        s.get_node(4).set_constraint(fixed_fixed_fixed)
        s.get_node(7).set_constraint(fixed_fixed_fixed)
        s.get_node(87).set_constraint(fixed_fixed_fixed)
        s.get_node(129).set_constraint(fixed_fixed_fixed)
        s.get_node(171).set_constraint(fixed_fixed_fixed)
        s.get_node(213).set_constraint(fixed_fixed_fixed)

        # Right Side
        s.get_node(250).set_constraint(fixed_fixed_fixed)
        s.get_node(208).set_constraint(fixed_fixed_fixed)
        s.get_node(166).set_constraint(fixed_fixed_fixed)
        s.get_node(124).set_constraint(fixed_fixed_fixed)
        s.get_node(81).set_constraint(fixed_fixed_fixed)
        s.get_node(80).set_constraint(fixed_fixed_fixed)
        s.get_node(251).set_constraint(fixed_fixed_fixed)
        s.get_node(209).set_constraint(fixed_fixed_fixed)
        s.get_node(167).set_constraint(fixed_fixed_fixed)
        s.get_node(125).set_constraint(fixed_fixed_fixed)
        s.get_node(83).set_constraint(fixed_fixed_fixed)
        s.get_node(82).set_constraint(fixed_fixed_fixed)

        # Outer CST's
        # From left to right
        s.get_node(252).set_constraint(fixed_fixed_fixed)
        s.get_node(253).set_constraint(free_fixed_fixed)
        s.get_node(256).set_constraint(free_fixed_fixed)
        s.get_node(258).set_constraint(free_fixed_fixed)
        s.get_node(260).set_constraint(free_fixed_fixed)
        s.get_node(262).set_constraint(free_fixed_fixed)
        s.get_node(264).set_constraint(free_fixed_fixed)
        s.get_node(266).set_constraint(free_fixed_fixed)
        s.get_node(268).set_constraint(free_fixed_fixed)
        s.get_node(270).set_constraint(free_fixed_fixed)
        s.get_node(272).set_constraint(free_fixed_fixed)
        s.get_node(274).set_constraint(free_fixed_fixed)
        s.get_node(276).set_constraint(free_fixed_fixed)
        s.get_node(278).set_constraint(free_fixed_fixed)
        s.get_node(280).set_constraint(free_fixed_fixed)
        s.get_node(282).set_constraint(free_fixed_fixed)
        s.get_node(284).set_constraint(free_fixed_fixed)
        s.get_node(286).set_constraint(free_fixed_fixed)
        s.get_node(288).set_constraint(free_fixed_fixed)
        s.get_node(290).set_constraint(free_fixed_fixed)
        s.get_node(292).set_constraint(free_fixed_fixed)
        # right side
        s.get_node(293).set_constraint(fixed_fixed_fixed)
        s.get_node(314).set_constraint(fixed_fixed_fixed)
        s.get_node(335).set_constraint(fixed_fixed_fixed)
        s.get_node(356).set_constraint(fixed_fixed_fixed)
        s.get_node(377).set_constraint(fixed_fixed_fixed)
        s.get_node(292).set_constraint(fixed_fixed_fixed)
        # From right to left
        s.get_node(377).set_constraint(fixed_fixed_fixed)
        s.get_node(376).set_constraint(free_fixed_fixed)
        s.get_node(375).set_constraint(free_fixed_fixed)
        s.get_node(374).set_constraint(free_fixed_fixed)
        s.get_node(373).set_constraint(free_fixed_fixed)
        s.get_node(372).set_constraint(free_fixed_fixed)
        s.get_node(371).set_constraint(free_fixed_fixed)
        s.get_node(370).set_constraint(free_fixed_fixed)
        s.get_node(369).set_constraint(free_fixed_fixed)
        s.get_node(368).set_constraint(free_fixed_fixed)
        s.get_node(367).set_constraint(free_fixed_fixed)
        s.get_node(366).set_constraint(free_fixed_fixed)
        s.get_node(365).set_constraint(free_fixed_fixed)
        s.get_node(364).set_constraint(free_fixed_fixed)
        s.get_node(363).set_constraint(free_fixed_fixed)
        s.get_node(362).set_constraint(free_fixed_fixed)
        s.get_node(361).set_constraint(free_fixed_fixed)
        s.get_node(360).set_constraint(free_fixed_fixed)
        s.get_node(359).set_constraint(free_fixed_fixed)
        s.get_node(358).set_constraint(free_fixed_fixed)
        s.get_node(357).set_constraint(fixed_fixed_fixed)
        # Left side
        s.get_node(336).set_constraint(fixed_fixed_fixed)
        s.get_node(315).set_constraint(fixed_fixed_fixed)
        s.get_node(295).set_constraint(fixed_fixed_fixed)
        s.get_node(255).set_constraint(fixed_fixed_fixed)

        s.get_node(40).set_constraint(free_fixed_free)
        s.get_node(42).set_constraint(free_fixed_free)
        s.get_node(44).set_constraint(free_fixed_free)
        s.get_node(46).set_constraint(free_fixed_free)
        s.get_node(36).set_constraint(free_fixed_free)
        s.get_node(38).set_constraint(free_fixed_free)
        s.get_node(34).set_constraint(free_fixed_free)
        s.get_node(32).set_constraint(free_fixed_free)
        s.get_node(30).set_constraint(free_fixed_free)
        s.get_node(28).set_constraint(free_fixed_free)
        s.get_node(26).set_constraint(free_fixed_free)
        s.get_node(24).set_constraint(free_fixed_free)
        s.get_node(22).set_constraint(free_fixed_free)
        s.get_node(20).set_constraint(free_fixed_free)
        s.get_node(16).set_constraint(free_fixed_free)
        s.get_node(18).set_constraint(free_fixed_free)
        s.get_node(12).set_constraint(free_fixed_free)
        s.get_node(14).set_constraint(free_fixed_free)
        s.get_node(8).set_constraint(free_fixed_free)
        s.get_node(12).set_constraint(free_fixed_free)
        s.get_node(1).set_constraint(free_fixed_free)
        s.get_node(5).set_constraint(free_fixed_free)
        s.get_node(48).set_constraint(free_fixed_free)
        s.get_node(50).set_constraint(free_fixed_free)
        s.get_node(52).set_constraint(free_fixed_free)
        s.get_node(54).set_constraint(free_fixed_free)
        s.get_node(56).set_constraint(free_fixed_free)
        s.get_node(58).set_constraint(free_fixed_free)
        s.get_node(60).set_constraint(free_fixed_free)
        s.get_node(62).set_constraint(free_fixed_free)
        s.get_node(64).set_constraint(free_fixed_free)
        s.get_node(66).set_constraint(free_fixed_free)
        s.get_node(68).set_constraint(free_fixed_free)
        s.get_node(70).set_constraint(free_fixed_free)
        s.get_node(72).set_constraint(free_fixed_free)
        s.get_node(74).set_constraint(free_fixed_free)
        s.get_node(76).set_constraint(free_fixed_free)
        s.get_node(78).set_constraint(free_fixed_free)

        # Forgotten Node
        s.get_node(10).set_constraint(free_fixed_free)

        s.get_node(230).set_constraint(free_fixed_free)
        s.get_node(231).set_constraint(free_fixed_free)
        s.get_node(229).set_constraint(free_fixed_free)
        s.get_node(228).set_constraint(free_fixed_free)
        s.get_node(232).set_constraint(free_fixed_free)
        s.get_node(233).set_constraint(free_fixed_free)
        s.get_node(227).set_constraint(free_fixed_free)
        s.get_node(226).set_constraint(free_fixed_free)
        s.get_node(225).set_constraint(free_fixed_free)
        s.get_node(224).set_constraint(free_fixed_free)
        s.get_node(223).set_constraint(free_fixed_free)
        s.get_node(222).set_constraint(free_fixed_free)
        s.get_node(221).set_constraint(free_fixed_free)
        s.get_node(220).set_constraint(free_fixed_free)
        s.get_node(219).set_constraint(free_fixed_free)
        s.get_node(218).set_constraint(free_fixed_free)
        s.get_node(217).set_constraint(free_fixed_free)
        s.get_node(216).set_constraint(free_fixed_free)
        s.get_node(215).set_constraint(free_fixed_free)
        s.get_node(214).set_constraint(free_fixed_free)
        s.get_node(212).set_constraint(free_fixed_free)
        s.get_node(210).set_constraint(free_fixed_free)
        s.get_node(235).set_constraint(free_fixed_free)
        s.get_node(234).set_constraint(free_fixed_free)
        s.get_node(236).set_constraint(free_fixed_free)
        s.get_node(237).set_constraint(free_fixed_free)
        s.get_node(239).set_constraint(free_fixed_free)
        s.get_node(238).set_constraint(free_fixed_free)
        s.get_node(241).set_constraint(free_fixed_free)
        s.get_node(240).set_constraint(free_fixed_free)
        s.get_node(242).set_constraint(free_fixed_free)
        s.get_node(243).set_constraint(free_fixed_free)
        s.get_node(244).set_constraint(free_fixed_free)
        s.get_node(245).set_constraint(free_fixed_free)
        s.get_node(246).set_constraint(free_fixed_free)
        s.get_node(247).set_constraint(free_fixed_free)
        s.get_node(248).set_constraint(free_fixed_free)
        s.get_node(249).set_constraint(free_fixed_free)

        # s.get_node(290).set_constraint(fixed_fixed_fixed)
        # s.get_node(376).set_constraint(fixed_fixed_fixed)
        # s.get_node(253).set_constraint(fixed_fixed_fixed)
        # s.get_node(358).set_constraint(fixed_fixed_fixed)

        # Apply Force to nodes
        force = Force(0.0, 0.0, -750.0)

        s.get_node(42).set_force(force)
        s.get_node(43).set_force(force)
        s.get_node(105).set_force(force)
        s.get_node(147).set_force(force)
        s.get_node(189).set_force(force)
        s.get_node(231).set_force(force)

        s.set_CST_mode("Plane Stress")
        return s


if __name__ == "__main__":
    # Create instance of the Class which contains the structures
    bridge = Bridge()
    # Create Structure
    s = bridge.bridge()
    # solve
    s.solve_direct_stiffness_method()
    #s.solve_explicit_euler()
    
    # Print infos to file
    s.print_info()

    # Initialise the 3D Plotter Window from pyVista
    plotter = pyvista.Plotter()
    # Initialise the viewer object
    viewer = Viewer(s, plotter)

    # Draw Hexahedral Elements
    viewer.visualise_linear_hexahedral_elements()
    # Draw CST Elements
    viewer.visualise_CST_elements()
    # Draw Truss Elements
    viewer.visualise_truss_elements()

    # Define result field of interest for visualisation
    field = "Displacement"
    #field = "Mean Stress"
    # Draw Hexahedral Displacements
    viewer.visualise_linear_hexahedral_deformation(field)
    # Draw CST Displacements
    viewer.visualise_CST_deformation(field)
    # Draw Truss Displacements
    viewer.visualise_truss_deformation(field)

    # Add all constraints to the 3D Plotter Window
    # viewer.visualise_constraints_by_cones()
    # Add all forces to the 3D Plotter Window
    viewer.visualise_forces_by_arrows()
    # Show all Node ID's
    # viewer.visualise_point_labels()

    plotter.show_axes()
    plotter.title = "Linear Static 3D FEA using a custom FE-Solver in Python"
    plotter.camera_position = [1, -1, 0.4]
    plotter.enable_parallel_projection()

    plotter.show()