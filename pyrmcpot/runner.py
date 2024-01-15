import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class Pyrmcpot_driver:
    def __init__(self):
        pass

    def run(self,*args,**kwargs):
        pass

'''
TODO: Implement callable like dictionary
'''
class Pyrmcpot_params:

    def __init__(self):
        self.nind = 9
        self.ngen = 500
        self.mintype = "global"
        self.algo = "pso"

if __name__ == '__main__':
    params = Pyrmcpot_params()
    driver = Pyrmcpot_driver(**params)

    driver.run()