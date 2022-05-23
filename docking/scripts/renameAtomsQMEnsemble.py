import os

src = os.getcwd()

def getListLines(file1):
	listLines = [line.rstrip("\t\n") for line in open(file1, "r")]
	return listLines

for file1 in os.listdir(src):
	if "conformer_" in file1:
		listLines = getListLines(file1)
		outputFile = file1.split(".pd")[0] + "_fix.pdb"
		if True == os.path.join(outputFile):
			os.remove(outputFile)
		with open(file1.split(".pd")[0] + "_fix.pdb", "a") as out:
			for i in range(0, len(listLines)):
				line = listLines[i]
				if 0 == i:
					lineToSave = line[0:14] + "1" + line[15:17] + "UNL" + line[20:len(line)] + "\n"
					out.write(lineToSave)
					continue
				if 1 == i:
					lineToSave = line[0:14] + "2" + line[15:17] + "UNL" + line[20:len(line)] + "\n"
					out.write(lineToSave)
					continue
				if 14 == i:
					lineToSave = line[0:14] + "1" + line[15:17] + "UNL" + line[20:len(line)] + "\n"
					out.write(lineToSave)
					continue
				if 15 == i:
					lineToSave = line[0:14] + "3" + line[15:17] + "UNL" + line[20:len(line)] + "\n"
					out.write(lineToSave)
					continue
				if 16 == i:
					lineToSave = line[0:14] + "4" + line[15:17] + "UNL" + line[20:len(line)] + "\n"
					out.write(lineToSave)
					continue
				if 17 == i:
					lineToSave = line[0:14] + "5" + line[15:17] + "UNL" + line[20:len(line)] + "\n"
					out.write(lineToSave)
					continue
				out.write(line[0:17] + "UNL" + line[20:len(line)] + "\n")
