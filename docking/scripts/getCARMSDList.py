"""
    This script calculates CA-RMSD for a pair of input PDB files.
    RDKit is employed and Python3 is required.
"""

import  os
import  sys

from  rdkit import Chem
from rdkit.Chem import rdMolAlign


def getListLines(file1):
    listLines = [line.rstrip("\t\n") for line in open(file1, "r")]
    return listLines


def  main( inputFilePath1, 
           inputFilePath2,
         ):
    molecule1 = Chem.rdmolfiles.MolFromPDBFile(inputFilePath1,
                                               sanitize=False,
                                               removeHs=True)
    molecule2 = Chem.rdmolfiles.MolFromPDBFile(inputFilePath2,
                                               sanitize=False,
                                               removeHs=True)
    listAtomIDs1 = []
    for atom in molecule1.GetAtoms():
        if "CA" == atom.GetPDBResidueInfo().GetName().strip():
            atomID = atom.GetIdx()
            listAtomIDs1.append(atomID)
    listAtomIDs2 = []
    for atom in molecule2.GetAtoms():
        if "CA" == atom.GetPDBResidueInfo().GetName().strip():
            atomID = atom.GetIdx()
            listAtomIDs2.append(atomID)
    atom_map = list(zip(listAtomIDs1, listAtomIDs2))
    rmsd = Chem.rdMolAlign.AlignMol(molecule1,
                                    molecule2,
                                    atomMap=atom_map,
                                    maxIters=100)
    return rmsd

if "__main__" == __name__:
    if 4 != len(sys.argv):
        print("please enter:\n(i). reference PDB file;\n(ii). text file with a list of paths for input PDB files;\n(iii). output CSV file to save results;\n")
        sys.exit()
    referencePDB = sys.argv[1]
    if False == os.path.exists(referencePDB):
        print(f'file {referencePDB} does not exist')
        sys.exit()
    inputPDB = sys.argv[2]
    if False == os.path.exists(inputPDB):
        print(f'file {inputPDB} does not exist')
        sys.exit()
    listLines = getListLines(inputPDB)
    outputFile = sys.argv[3]
    with open(outputFile, "a") as out:
        for linePDB in listLines:
             if True == os.path.exists(linePDB):
                  rmsd = main(referencePDB, linePDB)
                  line = linePDB + "," + str(rmsd) + "\n"
                  out.write(line)
