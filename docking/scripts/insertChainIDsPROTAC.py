"""This changes chain IDs for a PROTAC molecule in the input PDB file
   in accordance with the provided text file that contains two columns:
   chainID      atomName
   <...> \tab   <...>
   <...> \tab   <...>
   ...   \tab   ...
   """

import os
import sys

def getListLines(file1):
    listLines = [line.rstrip("\t\n") for line in open(file1, "r")]
    return listLines

def getPairListAtomNamesDictChainID(file1):
    listLines = getListLines(file1)
    result1 = []
    result2 = {}
    for line in listLines:
        atomName = line.split("\t")[1]
        result1.append(atomName)
        chainID = line.split("\t")[0]
        result2[atomName] = chainID
    return [result1, result2]

def main(inputPDB, outputPDB, infoPROTAC, moleculeName):
    listLines = getListLines(inputPDB)
    [listAtomNames, dictAtomNameChainID] = getPairListAtomNamesDictChainID(infoPROTAC)
    with open(outputPDB, "a") as out:
        for line in listLines:
            if ("ATOM  " != line[0:6]\
            and "HETATM" != line[0:6])\
            or "CLA" in line\
            or "POT" in line\
            or "SOLV" in line:
                continue
            atomName = line[12:16].strip()
            residueName = line[17:20].strip()
            if True == (atomName in listAtomNames)\
            and moleculeName == residueName:
                chainID = dictAtomNameChainID[atomName]
                lineToSave = line[0:21] + chainID + line[22:len(line)] + "\n"
                out.write(lineToSave)
            else:
                out.write(line + "\n")

if "__main__" == __name__:
    if 5 != len(sys.argv):
        print(f'please provide:\n(i). input PDB file;\n(ii). output PDB file;\n(iii). text file with PROTAC info;\n(iv). PROTAC molecule name (residue name from PDB file);')
        sys.exit()
    inputPDB = sys.argv[1]
    outputPDB = sys.argv[2]
    infoPROTAC = sys.argv[3]
    moleculeName = sys.argv[4]
    if False == os.path.exists(inputPDB):
        print(f'file {inputPDB} does not exist')
        sys.exit()
    if False == os.path.exists(infoPROTAC):
        print(f'file {infoPROTAC} does not exist')
        sys.exit()
    if True == os.path.exists(outputPDB):
        os.remove(outputPDB)
    main(inputPDB, outputPDB, infoPROTAC, moleculeName)
