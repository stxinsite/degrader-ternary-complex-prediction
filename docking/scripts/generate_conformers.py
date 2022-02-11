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
    if 0 < len( listLinesSMILES ):
        smilesLine  =  listLinesSMILES[ 0 ]
        if 0 != len( smilesLine.split( "\t" ) ):
            smilesLine  =  smilesLine.split( "\t" )[ 0 ]
    if 0 == len( smilesLine ):
        print( "empty SMILES line" )
        return
    molecule = Chem.rdmolfiles.MolFromSmiles( smilesLine )
    if None == molecule:
        raise FileNotFoundError
    [ listIDs,  listConformers ]  =  getListConformers( molecule, numConformers, maxAttempts, pruneRMSThresh )
    for cID in listIDs:
        outputFileName  =  "conformer_" + str( cID ) + ".pdb"
        outputFileNamePath  =  os.path.join( outputFolderPath, outputFileName )
        if True == os.path.exists( outputFileNamePath ):
            os.remove( outputFileNamePath )
        success  =  Chem.rdmolfiles.MolToPDBFile( listConformers, outputFileNamePath, cID )


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


