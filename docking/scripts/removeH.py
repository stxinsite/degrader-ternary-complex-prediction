import os
import sys

def getListLines(file1):
	listLines = [line.rstrip("\t\n") for line in open(file1, "r")]
	return listLines

def main(inputFile, outputFile):
	listLines = getListLines(inputFile)
	with open(outputFile, "a") as out:
		for line in listLines:
			if 78 <= len(line)\
			and "H" == line[77]:
				continue
			out.write(line + "\n")

if "__main__" == __name__:
	if 3 != len(sys.argv):
		print("please provide:\n(i). input PDB file;\n(ii). output PDB file;")
		sys.exit()
	inputPDB = sys.argv[1]
	if False == os.path.exists(inputPDB):
		print(f'file {inputPDB} does not exist')
		sys.exit()
	outputPDB = sys.argv[2]
	if True == os.path.exists(outputPDB):
		os.remove(outputPDB)
	main(inputPDB, outputPDB)
