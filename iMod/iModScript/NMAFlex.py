import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import sklearn.preprocessing as skpp
import subprocess
import sys

#class NMAFlex:
def NMAFlex(MDTrajObj, Threshold, Output, iModExec, tempRoot="./temp", modes=[0], deleteTemp=True):

    """
        Function that takes MDTraj Trajectory Object and a Threshold
        and returns a flex file which corresponds to the top Threshold
        highest and lowest eigenvectors.
        
        Parameters
        ----------
        
        MDTrajObj:  MDTrajTrajectory Object
                    Trajectory object for which the NMA will be done.
                
        Threshold: int
                    Threshold for dihedral selection. NOT USED!      
    
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
    FixArray = np.zeros(1, dtype=flexfields)
#    print ("FixArray:")
#    print (FixArray)
    for Line in FixIn:
        LineSplit = Line.strip("\n").split()
        if (len(LineSplit) == 4):
#            FixArray = np.vstack((FixArray,np.array(LineSplit, dtype=[("ResIx", float), ("Phi", float), ("Chi", float), ("Psi", float)])))
      #      print ("LineSplit:")
      #      print (np.array(LineSplit, dtype=FixArray.dtype))
            #FixArray = np.append(FixArray,np.array(LineSplit, dtype=FixArray.dtype))
            FixArray = np.append(FixArray,np.array([(LineSplit[0], LineSplit[1], LineSplit[2], LineSplit[3])], dtype=FixArray.dtype))
      #      print ("FixArray:")
     #       print (FixArray)
    FixArray = np.delete(FixArray, 0, axis=0)
    #print (FixArray)
    FixArrayEmpty = FixArray

    ##Parse the Evec file
    EvecInRaw = open("{}_ic.evec".format(tempRoot), "r").readlines()
    for ModeIx in modes:
        EvecIn = EvecInRaw
        counter = 0
        EvecNice = []
        FixArray = FixArrayEmpty
        for lineIx in range(len(EvecIn)):
            EvecIn[lineIx] = EvecIn[lineIx].split()
            if (EvecIn[lineIx][0] == "****"):
                counter = counter + 1
                continue
            if (counter == ModeIx):
                if (len(EvecIn[lineIx][0]) > 1):
                    for EvecValue in EvecIn[lineIx]:
                        EvecNice.append(EvecValue)
    ##Get absolute value of Evec, normalize values
        for EvecIx in range(len(EvecNice)):
            EvecNice[EvecIx] = abs(float(EvecNice[EvecIx]))
    ##Switch to array, normalize
        EvecNice = np.array(EvecNice)
        EvecNorm = skpp.MinMaxScaler()
        #print (EvecNice)
            
    ##Assign vectors to dihedral angles
        counter = 0
#        print (FixArray.shape)
        for aIx in range(FixArray.shape[0]):
            #print (FixArray[aIx])
         #   print (FixArray[aIx][1])
            if (FixArray[aIx]["Phi"] == 1):
                FixArray[aIx]["Phi"] = EvecNice[counter]
                counter +=1
            if (FixArray[aIx]["Psi"] == 1):
                FixArray[aIx]["Psi"] = EvecNice[counter]
                counter +=1
    
        ## We compute the Eigenvector distribution so we 
        ## can choose a proper number of dihedrals to set as flexible.
        ## In order to choose the highest contribution EVectors,
        ## we will choose those that are not in the most populated bin.

        ##We don't use histograms now
        if (0):
            PHIs = FixArray[::,1]
            PHIHisto = np.histogram(PHIs)
            PSIs = FixArray[::,3]
            PSIHisto = np.histogram(PSIs)
    
        ##We choose the first Threshold values
        sortFlexPhi = np.sort(FixArray, order="Phi")
        sortFlexPhi = np.flip(sortFlexPhi)
        #print (sortFlexPhi)
        sortFlexPsi = np.sort(FixArray, order="Psi")
        sortFlexPsi = np.flip(sortFlexPsi)
        #print (sortFlexPsi)

        #PHIThreshold = np.array((PHIHisto[1][np.argmax(PHIHisto[0])], PHIHisto[1][np.argmax(PHIHisto[0]) + 1]))
        #PSIThreshold = np.array((PSIHisto[1][np.argmax(PSIHisto[0])], PSIHisto[1][np.argmax(PSIHisto[0]) + 1]))
        PHIThreshold = Threshold
        PSIThreshold = Threshold 
    
        ## Parse the FixArray and add ResIDs with values
        ## outside of any of the two thresholds to a list
        if (0):
            RelevantRes = []
            for FixLine in FixArray:
                if ((FixLine[1] > PHIThreshold[0]) and (FixLine[3] > PSIThreshold[0])):
                    pass
                else:
                    RelevantRes.append(int(FixLine[0]))
        elif (1):
            Residues = 0
            RelevantRes = set()
            ## We get top Threshold Phi angles
            for i in range(FixArray.shape[0]):
                if (Residues >= Threshold):
                    break
                if (sortFlexPhi[i]["Phi"] != 0):
                    #print ("added!")
                    RelevantRes.add(int(sortFlexPhi[i]["ResIx"]))
                    Residues = Residues + 1
            ## Then we get top Threshold Psi angles
            Residues = 0
            for i in range(FixArray.shape[0]):
                if (Residues >= Threshold):
                    break
                if (sortFlexPhi[i]["Psi"] != 0):
                    #print ("added!")
                    RelevantRes.add(int(sortFlexPhi[i]["ResIx"]))
                    Residues = Residues + 1

            RelevantRes = list(RelevantRes)
            RelevantRes.sort()
            print (RelevantRes)


    
        print("Predicted hinge residues for mode {}: {}".format(ModeIx, RelevantRes))
    
        ## Extract the associated PHI/PSI angles from the above residues
        ## using the MDTrajObj that the user loads.
    
        FlexOut = open("{}.{}.flex".format(Output,ModeIx), "w")
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
