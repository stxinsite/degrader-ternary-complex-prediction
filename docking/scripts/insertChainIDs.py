"""This allows to insert/substitute Chain IDs in an input PDB file
   in the specified lines of the file
   """

import os
import sys

def getListLines(file1):
    listLines = [line.rstrip("\t\n") for line in open(file1, "r")]
    return listLines

def main(inputFileName, chainID, startLine, finishLine, outputFileName):
    with open(outputFileName, "a") as out:
        listLines = getListLines(inputFileName)
        for i in range(0, len(listLines)):
            line = listLines[i]
            if len(line) >= 21\
            and ("HOH" == line[17:20]\
            or "WAT" == line[17:20]\
            or "Cl-" == line[17:20]\
            or "K+" == line[17:20].strip()\
            or "ZN" == line[17:20].strip()\
            or "TER" == line[0:3]\
            or "H" == (line[12:16].strip())[0]):
                continue
            if startLine <= i\
            and finishLine > i\
            and ("ATOM  " == line[0:6]\
            or "HETATM" == line[0:6] ):
                lineToSave = line[0:21] + chainID + line[22:len(line)] + "\n"
                out.write(lineToSave)
            else:
                out.write(line + "\n")

if "__main__" == __name__:
    if 6 != len(sys.argv):
        print("please provide:\n(i). input PDB;\n(ii). chain ID to insert;\n(iii). start line index;\n(iv). finish line index;\n(v). output PDB;")
        sys.exit()
    inputPDB = sys.argv[1]
    chainID = sys.argv[2]
    startLine = int(sys.argv[3])
    finishLine = int(sys.argv[4])
    outputPDB = sys.argv[5]
    if False == os.path.exists(inputPDB):
        print(f'file {inputPDB} does not exist')
        sys.exit()
    if " " == chainID:
        print(f'chainID= {chainID} in not valid')
        sys.exit()
    if 0 > startLine\
    or 0 > finishLine\
    or startLine > finishLine:
        print(f'values startLine= {startLine} and finishLine= {finishLine} are not valid')
        sys.exit()
    if True == os.path.exists(outputPDB):
        os.remove(outputPDB)
    main(inputPDB, chainID, startLine, finishLine, outputPDB)
