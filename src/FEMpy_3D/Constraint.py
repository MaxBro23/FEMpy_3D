import numpy as np


class Constraint(object):
    """Constraint class. Defined wether or not, a node is free or fixed in x,y or z-direction.
    """

    def __init__(self, ux, uy, uz):
        self.constraints = np.array([ux, uy, uz])

    def is_free(self, i):
        """
        :param i: int
            Index for x, y or z-direction
        :return: bool
            Return true if direction is free or false otherwise
        """
        return self.constraints[i]

    def info(self):
        """Print information of constraint object
        :return: None
        """
        for i in range(0, 3, 1):
            if self.constraints[i] == True:
                print("Free")
            else:
                print("Fixed")
