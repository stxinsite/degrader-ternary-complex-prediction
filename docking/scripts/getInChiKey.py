"""
    This script generates InChi Key for a given SMILES string
    RDKit is employed and Python3 is required.
"""

import  os
import  sys

from  rdkit import Chem


def getListLines(file1):
    listLines = [line.rstrip("\t\n") for line in open(file1, "r")]
    return  listLines


def  main( inputFilePath, 
        outputFilePath, 
        ):
    listLinesSMILES  =  getListLines( inputFilePath )
    smilesLine  =  ""
    if 0 < len( listLinesSMILES ):
        smilesLine  =  listLinesSMILES[ 0 ]
    if 0 == len( smilesLine ):
        print( "empty SMILES line" )
        return
    molecule = Chem.rdmolfiles.MolFromSmiles( smilesLine )
    if None == molecule:
        raise FileNotFoundError
    inchiKey = Chem.inchi.MolToInchiKey(molecule)
    print(inchiKey)
    with open(outputFilePath, "w") as out:
        out.write( inchiKey + " " + smilesLine )


if "__main__" == __name__:
    if 3 != len( sys.argv ):
        print( "please provide three arguments:\n(i). input SMILES-file (full path);\n(ii). output folder (full path)" )
        exit()
    inputFilePath  =  sys.argv[ 1 ]
    outputFilePath  =  sys.argv[ 2 ]
    if True == os.path.exists(outputFilePath):
        os.remove(outputFilePath)
    try:
        main( inputFilePath, outputFilePath )
    except  FileNotFoundError:
        print( f'input file {inputFilePath} is not found' )
        exit()


