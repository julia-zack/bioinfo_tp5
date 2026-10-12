"""Ejercicio 7: adenilato quinasa de E. coli, cerrada (1AKE, con Ap5A) vs
abierta (4AKE), solo carbonos alfa, cadena A de cada una."""

from pathlib import Path

from pairs import compare

CLOSED, OPEN = "1AKE", "4AKE"


def run():
    return compare(CLOSED, OPEN, figs_dir=Path(__file__).parent, prefix="ex7_")
