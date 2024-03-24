# coding: utf-8

from dataclasses import dataclass

@dataclass
class ConfigDataclass:
    dim: int = 1
    type: str = "binary"
    length: int = 100
    gap: float = 2.0   # Gap is distance between two sheets / longitudinal length of a sheet.
