from NMAFlex import NMAFlex
import mdtraj as md

a = md.load("TrimerNoGlyc.inpcrd",top="TrimerNoGlyc.prmtop")

NMAFlex(a, 25, "Trimer", "../iMOD_v1.04_Linux/bin/imode_gcc", deleteTemp=False, modes=[1])
