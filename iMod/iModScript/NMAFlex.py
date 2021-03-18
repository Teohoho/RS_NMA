import mdtraj as md
import numpy as np
import subprocess
import sys

class NMAtoFlex:
    def __init__(self, MDTrajObj, iModExec, tempRoot="./temp"):
        """
        Class that takes an MDTraj Trajectory object, ideally a prmtop/inpcrd pair 
        so the atom indices don't change, and runs iMod on the coordinates, then
        generates a flex file compatible with the Robosample package.

        Parameters
        ----------

        MDTrajObj:  MDTrajTrajectory Object
                    Trajectory object for which the NMA will be done.
        iModExec:   str
                    Path to iMod executable (imode_gcc)
        tempRoot:   str
                    Path where temporary iMod-generated files will be stored
        """
    
        ## Generate a PDB
        PDBName = "{}.pdb".format(tempRoot)
        MDTrajObj.save_pdb(PDBName, force_overwrite=True)
    
        ## Run the actual iMod command
        DEBUG = 1
        if (DEBUG == 0):
            iModCommand = [iModExec, PDBName, "-o {}".format(tempRoot), "--save_fixfile"]
            iMod_sub = subprocess.Popen(iModCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            iMod_sub.communicate()
        
        ## Generate Appropriate Arrays containing information
        ## Parse the fixfile, which is the same regardless of mode
        flexfields = [("ResIx", int), ("Phi", float), ("Chi", float), ("Psi", float)]
        FixIn = open("{}.fix".format(tempRoot), "r").readlines()
        FixArray = np.zeros(0, dtype=flexfields)
        for Line in FixIn:
            LineSplit = Line.strip("\n").split()
            if (len(LineSplit) == 4):
                FixArray = np.append(FixArray,np.array([(LineSplit[0], LineSplit[1], LineSplit[2], LineSplit[3])], dtype=FixArray.dtype)) 
        ## Parse the Evec file, generating a list of lists, one for each mode
        EvecIn = open("{}_ic.evec".format(tempRoot), "r").readlines()
        for lineIx in range(len(EvecIn)):
            EvecIn[lineIx] = EvecIn[lineIx].split()
        ## Get the number of computed modes
        NOfModes = int(EvecIn[1][3])
        print ("Computed a number of {} modes".format(NOfModes))
        EvecIn = EvecIn[2:]
        EvecModes = []

        ## We retrieve the absolute value of the Evec contributions
        CurrentMode=[]
        for lineIx in range(len(EvecIn)):
            if (EvecIn[lineIx][0] == "****"):
                if (len(CurrentMode) > 0):
                    EvecModes.append(CurrentMode)
                    CurrentMode=[]
                continue
            if (len(EvecIn[lineIx][0]) > 2):
                for EvecValue in EvecIn[lineIx]:
                    CurrentMode.append(abs(float(EvecValue)))
        EvecModes.append(CurrentMode)

        if (DEBUG == 1):
            print(len(EvecModes))
            print((EvecModes[0][0]))
            print((EvecModes[19][0]))

        self.FixArray = FixArray
        self.EvecModes = EvecModes
        self.MDTrajObj = MDTrajObj
        self.DEBUG = DEBUG

    def GetFlex(self, Output, Threshold=25, modes=[1]):

        """
        Function that uses the previously generated arrays to generate the
        Robosample-compatible flex files

        Parameters
        ----------

        Output:     str
                    Root name for generated flex files.

        Threshold:  int
                    After ranking angles by contribution, generate a flex file
                    using the top <Threshold> Phi/Psi angles.

        modes:      list of ints
                    Generate flex files corresponding to these modes.
        """

        for ModeIx in modes:
            ModeArray = self.FixArray.copy()
            
            ##Assign vectors to dihedral angles
            counter = 0
            for aIx in range(ModeArray.shape[0]):
                if (ModeArray[aIx]["Phi"] != 0):
                    ModeArray[aIx]["Phi"] = self.EvecModes[ModeIx-1][counter]
                    counter +=1
                if (ModeArray[aIx]["Psi"] != 0):
                    ModeArray[aIx]["Psi"] = self.EvecModes[ModeIx-1][counter]
                    counter +=1

            ##Sort Arrays by Phi/Psi, in decreasing order
            sortFlexPhi = np.sort(ModeArray, order="Phi")
            sortFlexPhi = np.flip(sortFlexPhi)
            sortFlexPsi = np.sort(ModeArray, order="Psi")
            sortFlexPsi = np.flip(sortFlexPsi)

            PhiRes = set([])
            PsiRes = set([])
            
            ## We get top Threshold Phi angles
            Residues = 0
            for i in range(ModeArray.shape[0]):
                if (Residues >= Threshold):
                    break
                if (sortFlexPhi[i]["Phi"] != 0):
                    PhiRes.add(sortFlexPhi[i]["ResIx"])
                    Residues = Residues + 1

            ## Then we get top Threshold Psi angles
            Residues = 0
            for i in range(ModeArray.shape[0]):
                if (Residues >= Threshold):
                    break
                if (sortFlexPhi[i]["Psi"] != 0):
                    PsiRes.add(sortFlexPsi[i]["ResIx"])
                    Residues = Residues + 1 

            PhiRes = list(PhiRes)
            PsiRes = list(PsiRes)
            PhiRes.sort()
            PsiRes.sort()
            if (self.DEBUG == 1):
                print ("PHI Res: {} \nPSI Res: {}".format(PhiRes,PsiRes))

            
        ## Extract the associated PHI/PSI angles from the above residues
        ## using the MDTrajObj that the user loads.
    
            FlexOut = open("{}.{}.flex".format(Output,ModeIx), "w")
            FlexOut.write("##MODE {}##\n".format(ModeIx))
            FlexOut.write("##PHI ANGLES##\n")
            for ResIx in PhiRes:
                NIndex = self.MDTrajObj.topology.select("resid {} and name N".format(ResIx))[0]
                CAIndex = self.MDTrajObj.topology.select("resid {} and name CA".format(ResIx))[0]
                if (str(self.MDTrajObj.topology.residue(ResIx))[0:3] != "PRO"):
                    FlexOut.write("{} {} Pin \n".format(NIndex,CAIndex))
            FlexOut.write("##PSI ANGLES##\n")
            for ResIx in PsiRes:
                CAIndex = self.MDTrajObj.topology.select("resid {} and name CA".format(ResIx))[0]
                CIndex = self.MDTrajObj.topology.select("resid {} and name C".format(ResIx))[0]
                FlexOut.write("{} {} Pin \n".format(CAIndex,CIndex))
