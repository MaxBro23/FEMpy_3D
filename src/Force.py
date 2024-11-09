import numpy as np


class Force(object):

    def __init__(self, x, y, z):
        self.force_vector = np.array([x, y, z])

    def get_component(self, i):
        """
        :param i: int
            Index of the component e.g. 0, 1 or 2
        :return:
             Component of the force in x,y or z e.g. 0,1 or 2
        """
        return self.force_vector[i]

    def get_force_vector(self):
        """
        :return: np-array
            Returns this force object array
        """
        return self.force_vector

    def info(self):
        """Print information of force components"""
        print("Force in x: ", self.force_vector[0])
        print("Force in y: ", self.force_vector[1])
        print("Force in z: ", self.force_vector[2])
