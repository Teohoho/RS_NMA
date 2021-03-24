from NMAFlex import NMAtoFlex
import mdtraj as md

a = md.load("TrimerNoGlyc.inpcrd",top="TrimerNoGlyc.prmtop")

#NMAFlex(a, 25, "Trimer", "../iMOD_v1.04_Linux/bin/imode_gcc", deleteTemp=False, modes=[1])
test = NMAtoFlex(MDTrajObj=a, iModExec="../iMOD_v1.04_Linux/bin/imode_gcc", tempRoot="./ClassTest")

#test.GetFlex(Output="ClassTest",Threshold=25, modes=[1,20])
test.GetFlexScaled(Output="ClassTest", GenerateVMD=True)
