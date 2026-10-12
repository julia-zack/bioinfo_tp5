"""Ejercicio 8: lo mismo que el 7 para la Ca2+-ATPasa del retículo (SERCA1a
de conejo), los tres pares posibles entre:
  1SU4  E1·2Ca2+ (dos Ca unidos)
  4XOU  Ca2-E1·Mg-AMPPCP (análogo de ATP)
  3B9B  E2·BeF3- (análogo de fosfoenzima E2-P)
Misma proteína y misma numeración, así que los CA se emparejan por número."""

from itertools import combinations
from pathlib import Path

from pairs import compare

STRUCTURES = ["1SU4", "4XOU", "3B9B"]


def run():
    for a, b in combinations(STRUCTURES, 2):
        print(f"\n######## {a} vs {b}")
        compare(a, b, figs_dir=Path(__file__).parent, prefix="ex8_")
