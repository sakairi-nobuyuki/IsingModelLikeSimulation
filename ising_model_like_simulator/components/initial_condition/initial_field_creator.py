# coding: utf-8

from typing import Dict, Any
import cupy as cp

from ising_model_like_simulator.data_structures import ConfigDataclass

class InitialFieldCreator:
    def __init__(self, config: ConfigDataclass) -> None:
        self.config = config
        
        if self.config.dim == 1:
            if self.config.type == "binary":
                self.__create = self.__create_one_dimensional_binary_array

    def create(self) -> cp.ndarray:
        return self.__create()

    def __create_one_dimensional_binary_array(self) -> None:
        return cp.random.randint(0, 2, size=self.config.length)
        