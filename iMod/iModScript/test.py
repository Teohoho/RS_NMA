from NMAFlex import NMAFlex
import mdtraj as md

a = md.load("TrimerChainA.min.inpcrd",top="TrimerChainA.prmtop")

NMAFlex(a, 0, "test.flex", "../iMOD_v1.04_Linux/bin/imode_gcc")
