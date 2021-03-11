import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys

def NMAFlex(MDTrajObj, Threshold, Output, iModExec, tempRoot="./temp", modes=[1], deleteTemp=True):

    """
        Function that takes MDTraj Trajectory Object and a Threshold
        and returns a flex file which corresponds to the top Threshold
        highest and lowest eigenvectors.
        
        Parameters
        ----------
        
        MDTrajObj:  MDTrajTrajectory Object
                    Trajectory object for which the NMA will be done.
                
        Threshold: int
                    Threshold for dihedral selection         
    
        Output: str
                    Root of the flex file tot which to write.

        iModExed:str
                    Path to iMod executable (imode_gcc)

        tempdir: str
                    Where to store the temp files that are generated 
        modes: list of str
                    List of modes to generate flex files for
        deleteTemp: bool
                    Whether or not to remove all temp files generated during 
                    flex file generation
    """        

    ## Generate a PDB
    PDBName = "{}.pdb".format(tempRoot)
    MDTrajObj.save_pdb(PDBName, force_overwrite=True)

    ## Run the actual iMod command
    iModCommand = [iModExec, PDBName, "-o {}".format(tempRoot), "--save_fixfile"]
    iMod_sub = subprocess.Popen(iModCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    iMod_sub.communicate()

    ##Parse the fixfile
    flexfields = [("ResIx", float), ("Phi", float), ("Chi", float), ("Psi", float)]
    FixIn = open("{}.fix".format(tempRoot), "r").readlines()
    FixArray = np.zeros(4)
    for Line in FixIn:
        LineSplit = Line.strip("\n").split()
        if (len(LineSplit) == 4):
            FixArray = np.vstack((FixArray,np.array(LineSplit, dtype=int)))
    FixArray = np.delete(FixArray, 0, axis=0)

    ##Parse the Evec file
    counter = 0
    EvecIn = open("{}_ic.evec".format(tempRoot), "r").readlines()
    EvecNice = []
    for ModeIx in modes:
        for lineIx in range(len(EvecIn)):
            EvecIn[lineIx] = EvecIn[lineIx].split()
            if (EvecIn[lineIx][0] == "****"):
                counter = counter + 1
                continue
            if (counter == ModeIx):
                if (len(EvecIn[lineIx][0]) > 1):
                    for EvecValue in EvecIn[lineIx]:
                        EvecNice.append(EvecValue)
    
    ##Assign vectors to dihedral angles
    counter = 0
    for aIx in range(FixArray.shape[0]):
        if (FixArray[aIx][1] == 1):
            FixArray[aIx][1] = EvecNice[counter]
            counter +=1
        if (FixArray[aIx][3] == 1):
            FixArray[aIx][3] = EvecNice[counter]
            counter +=1

    ## We compute the Eigenvector distribution so we 
    ## can choose a proper number of dihedrals to set as flexible.
    ## In order to choose the highest contribution EVectors,
    ## we will choose those that are not in the most populated bin.
    PHIs = FixArray[::,1]
    PHIHisto = np.histogram(PHIs)
    PSIs = FixArray[::,3]
    PSIHisto = np.histogram(PSIs)

    PHIThreshold = np.array((PHIHisto[1][np.argmax(PHIHisto[0])], PHIHisto[1][np.argmax(PHIHisto[0]) + 1]))
    PSIThreshold = np.array((PSIHisto[1][np.argmax(PSIHisto[0])], PSIHisto[1][np.argmax(PSIHisto[0]) + 1]))

    ## Parse the FixArray and add ResIDs with values
    ## outside of any of the two thresholds to a list
    RelevantRes = []
    for FixLine in FixArray:
        if ((FixLine[1] > PHIThreshold[0] and FixLine[1] < PHIThreshold[1]) and (FixLine[3] > PSIThreshold[0] and FixLine[3] < PSIThreshold[1])):
            pass
        else:
            RelevantRes.append(int(FixLine[0]))

    print("Predicted hinge residues: {}".format(RelevantRes))

    ## Extract the associated PHI/PSI angles from the above residues
    ## using the MDTrajObj that the user loads.

    FlexOut = open(Output, "w")
    for ResIx in RelevantRes:
        NIndex = MDTrajObj.topology.select("resid {} and name N".format(ResIx))[0]
        CAIndex = MDTrajObj.topology.select("resid {} and name CA".format(ResIx))[0]
        CIndex = MDTrajObj.topology.select("resid {} and name C".format(ResIx))[0]
        if (str(MDTrajObj.topology.residue(ResIx))[0:3] == "PRO"):
            FlexOut.write("{} {} Pin \n".format(CAIndex, CIndex))
        else:
            FlexOut.write("{} {} Pin \n{} {} Pin\n".format(NIndex,CAIndex,CAIndex,CIndex))

    ## Now we delete all the temp files we created during the process
    if (deleteTemp is True):
        for FileExt in [".fix", ".log", ".pdb", "_ic.evec", "_model.pdb"]:
            procList   = ["rm", "{}{}".format(tempRoot,FileExt)]
            removeTemp = subprocess.Popen(procList, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    print ("Done!")
