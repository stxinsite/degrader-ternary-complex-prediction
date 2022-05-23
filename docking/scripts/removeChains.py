import os
import sys

def getListLines(file1):
	listLines = [line.rstrip("\t\n") for line in open(file1, "r")]
	return listLines

def main(inputFolder):
	for file1 in os.listdir(inputFolder):
		if "_fixed.pdb" == file1[len(file1) - 10:len(file1)]\
		and "mini_" == file1[0:5]:
			file1Path = os.path.join(inputFolder, file1)
			fileOutPath = file1Path[ 0: len(file1Path) - 4] + ".AB.pdb"
			if True == os.path.exists(fileOutPath):
				os.remove(fileOutPath)
			with open(fileOutPath, "a") as out:
				listLines = getListLines(file1Path)
				for line in listLines:
					if 22 <= len(line):
						if "C" != line[21]\
						and "D" != line[21]:
							if "TRN" != line[17:20]:
								out.write(line + "\n")
							else:
								lineToSave = line[0:17] + "FWZ" + line[20:len(line)] + "\n"
								out.write(line + "\n")

if "__main__" == __name__:
	inputFolder = sys.argv[1]
	if False == os.path.exists(inputFolder):
		print(f'folder {inputFolder} does not exist')
		sys.exit()
	main(inputFolder)
