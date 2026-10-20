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
# Apodos que se muestran junto al código en gráficos y salidas.
LABELS = {"1SU4": "con Ca", "4XOU": "Ca + ATP", "3B9B": "sin Ca"}

# Dominios de la SERCA1a según la bibliografía (no figuran en el archivo):
# A (actuador), P (fosforilación), N (unión de nucleótido) y M (membrana, 10 hélices).
DOMAINS = {"A": list(range(1, 41)) + list(range(124, 244)),
           "P": list(range(330, 360)) + list(range(605, 738)),
           "N": list(range(360, 605)),
           "M": list(range(41, 124)) + list(range(244, 330)) + list(range(738, 995))}


def run():
    for a, b in combinations(STRUCTURES, 2):
        print(f"\n######## {a} ({LABELS[a]}) vs {b} ({LABELS[b]})")
        compare(a, b, figs_dir=Path(__file__).parent, prefix="ex8_", domains=DOMAINS, labels=LABELS)
