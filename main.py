"""Entry point for the TP. Runs one exercise at a time, by name:

    venv/bin/python main.py 2b
    venv/bin/python main.py 6v
    venv/bin/python main.py          # lists the available exercises

The names match the numbering in the TP PDF. The 3D views are interactive
and live in notebooks instead: ex3_nglview, ex5_superposicion, ex6_vecinos,
ex7_visualizacion.
"""

import sys

from exercises.ex1 import ex1
from exercises.ex2 import ex2
from exercises.ex4 import ex4
from exercises.ex5 import ex5
from exercises.ex6 import ex6
from exercises.ex7 import ex7
from exercises.ex8 import ex8
from exercises.ex9 import ex9

EXERCISES = {
    "1":     ex1.run,     # cargar PDB y mmCIF, explorar el header
    "2a":    ex2.a,       # residuo hemo, id y full_id; otros hetero
    "2b":    ex2.b,       # cadenas, residuos, átomos, HETATM, agua
    "2c":    ex2.c,       # secuencia de cada cadena
    "2d":    ex2.d,       # PDB con solo los CA
    "4":     ex4.run,     # Fe-O, Fe-O-O, N-Fe-O-O en oxiMb y oxiHb
    "5":     ex5.run,     # superponer dos mioglobinas
    "6i":    ex6.i,       # SASA por residuo
    "6ii":   ex6.ii,      # contactos por CA, CB y centro de masa
    "6iii":  ex6.iii,     # correlación contactos vs SASA
    "6iv":   ex6.iv,      # vecinos polares / no polares en esferas
    "6v":    ex6.v,       # vecinos vs pKa de Asp
    "7":     ex7.run,     # adenilato quinasa abierta vs cerrada
    "7-diedros": ex7.dihedral_explainer,  # figura que explica cómo leer los diedros
    "8":     ex8.run,     # SERCA, los tres pares
    "9":     ex9.run,     # shadow contact map
}


def usage():
    print("usage: python3 main.py <exercise>")
    print("available exercises:")
    for name in EXERCISES:
        print(f"  {name}")


def main():
    if len(sys.argv) != 2:
        usage()
        return
    name = sys.argv[1]
    if name not in EXERCISES:
        print(f"'{name}' is not a known exercise.\n")
        usage()
        return
    EXERCISES[name]()


if __name__ == "__main__":
    main()
