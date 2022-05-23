"""
    This script generates conformers for the PROTAC linkers.
    RDKit is employed and Python3 is required.
"""

import  os
import  sys

from  rdkit import Chem
from  gen_conformers.getListConformers import getListConformers
from  gen_conformers.getListLines import getListLines



def  main( inputFilePath, 
        outputFolderPath, 
        numConformers,
        maxAttempts,
        pruneRMSThresh
        ):
    listLinesSMILES  =  getListLines( inputFilePath )
    smilesLine  =  ""
    listSMILESCanonical = []
    for line in listLinesSMILES:
        if 0 != len( line.split( "\t" ) ):
            smilesLine  =  line.split( "\t" )[ 0 ]
            amgenID = line.split("\t")[1]
            if 0 == len( smilesLine ):
                print( "empty SMILES line" )
                return
            molecule = Chem.rdmolfiles.MolFromSmiles( smilesLine )
            smilesCanonical = Chem.rdmolfiles.MolToSmiles(molecule, True)
            if False == (smilesCanonical in listSMILESCanonical):
                listSMILESCanonical.append(smilesCanonical)
            else:
                continue
            if None == molecule:
                raise FileNotFoundError
            [ listIDs,  listConformers ]  =  getListConformers( molecule, numConformers, maxAttempts, pruneRMSThresh )
            subFolderPath = os.path.join(outputFolderPath, amgenID)
            if False == os.path.exists(subFolderPath):
                os.mkdir(subFolderPath)
            else:
                maxIndex = 0 # number of already present structures with the given ID (i.e., folder name)
                for folder in os.listdir(outputFolderPath):
                    if 2 == len(folder.split("_"))\
                    and folder.split("_")[0] == amgenID:
                        if maxIndex < int(folder.split("_")[1]):
                            maxIndex = int(folder.split("_")[1])
                maxIndex = maxIndex + 1
                subFolderPath = os.path.join(outputFolderPath, amgenID + "_" + str(maxIndex))
                os.mkdir(subFolderPath)
            for cID in listIDs:
                outputFileName  =  amgenID + "_conf_" + str( cID ) + ".pdb"
                outputFileNamePath  =  os.path.join( subFolderPath, outputFileName )
                if True == os.path.exists( outputFileNamePath ):
                    os.remove( outputFileNamePath )
                success  =  Chem.rdmolfiles.MolToPDBFile( listConformers, outputFileNamePath, cID )
            with open("stat_amgen_conf_all.csv", "a") as out:
                line = amgenID + "," + str(pruneRMSThresh) + "," + str(len(listIDs)) + "\n"
                out.write(line)


if "__main__" == __name__:
    if 6 != len( sys.argv ):
        print( "please provide three arguments:\n(i). input SMILES-file (full path);\n(ii). output folder (full path)\n(iii). number of conformers to generate\n(iv). maximal number of attempts to try embedding a new conformation\n(v). RMS value to retain only those conformations that are at least that far apart from each other" )
        exit()
    inputFilePath  =  sys.argv[ 1 ]
    outputFolderPath  =  sys.argv[ 2 ]
    numConformers  =  int( sys.argv[ 3 ] )
    maxAttempts  =  int( sys.argv[ 4 ] )
    pruneRMSThresh  =  float( sys.argv[ 5 ] )
    if False == os.path.exists( outputFolderPath ):
        os.mkdir( outputFolderPath )
    try:
        main( inputFilePath, outputFolderPath, numConformers, maxAttempts, pruneRMSThresh )
    except  FileNotFoundError:
        print( f'input file {inputFilePath} is not found' )
        exit()


