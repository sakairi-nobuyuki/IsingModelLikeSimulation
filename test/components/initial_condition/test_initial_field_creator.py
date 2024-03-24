# coding: utf-8

import cupy as cp

from ising_model_like_simulator.components.initial_condition import InitialFieldCreator
from ising_model_like_simulator.data_structures import ConfigDataclass

class TestInitialFieldCreator:

    def test_one_dimensional_binary_array(self) -> None:
        config = ConfigDataclass()
        assert isinstance(config, ConfigDataclass)
        field_creator = InitialFieldCreator(config)

        order_parameter = field_creator.create()

        assert isinstance(order_parameter, cp.ndarray)
        assert cp.max(order_parameter) == 1
        assert cp.min(order_parameter) == 0