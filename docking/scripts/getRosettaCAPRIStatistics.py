import os
import sys

def getListLines(file1):
	listLines = [line.rstrip("\t\n") for line in open(file1, "r")]
	return listLines

PATTERN_FOLDER = "job_" # job_1 = with HDX, job_h = without HDX
PATTERN_FILE = "score.rerank.sc"

def main(inputFolder, outputFile, outputFileValues, topNArray = [10]):
	with open(outputFileValues, "a") as outValues:
		lineToSave = "topN,fNat,lRMSD,iRMSD\n"
		outValues.write(lineToSave)
	with open(outputFile, "a") as out:
		lineToSave = "topN,category,num\n"
		out.write(lineToSave)
		for i in range(0, len(topNArray)):
			topN = topNArray[i]
			count = 0
			for folder1 in os.listdir(inputFolder):
				if PATTERN_FOLDER == folder1[0:4]:
					folder1Path = os.path.join(inputFolder, folder1)
					inputFile = os.path.join(folder1Path, PATTERN_FILE)
					if False == os.path.exists(inputFile):
						continue
					listLines = getListLines(inputFile)
					topIncorrect = 0
					topAcceptable = 0
					topMedium = 0
					topHigh = 0
					for j in range(1, topN):
						if j >= len(listLines):
							break
						line = listLines[j]
						rank = int(float(line.split(",")[2]))
						if 0 == rank:
							topIncorrect = topIncorrect + 1
						if 1 == rank:
							topAcceptable = topAcceptable + 1
						if 2 == rank:
							topMedium = topMedium + 1
						if 3 == rank:
							topHigh = topHigh + 1
						fNat = float(line.split(",")[3])
						iRMSD = float(line.split(",")[5])
						lRMSD = float(line.split(",")[1])
						with open(outputFileValues, "a") as outValues:
							lineToSave = str(topN) + ",fNat," + str(fNat) +"\n"
							outValues.write(lineToSave)
							lineToSave = str(topN) + ",L-RMSD," + str(lRMSD) +"\n"
							outValues.write(lineToSave)
							lineToSave = str(topN) + ",I-RMSD," + str(iRMSD) +"\n"
							outValues.write(lineToSave)
					lineToSave = inputFile + "," + str(topN) + ",High," + str(topHigh) + "\n"
					out.write(lineToSave)
					lineToSave = inputFile + "," + str(topN) + ",Medium," + str(topMedium) + "\n"
					out.write(lineToSave)
					lineToSave = inputFile + "," + str(topN) + ",Acceptable," + str(topAcceptable) + "\n"
					out.write(lineToSave)
					lineToSave = inputFile + "," + str(topN) + ",Incorrect," + str(topIncorrect) + "\n"
					out.write(lineToSave)

if "__main__" == __name__:
	inputFolder = sys.argv[1]
	outputFileCategories = sys.argv[2]
	outputFileValues = sys.argv[3]
	if True == os.path.exists(outputFileCategories):
		os.remove(outputFileCategories)
	if True == os.path.exists(outputFileValues):
		os.remove(outputFileValues)
	if False == os.path.exists(inputFolder):
		print(f'folder {inputFolder} does not exist')
		sys.exit()
	main(inputFolder, outputFileCategories, outputFileValues, [10, 50, 100, 500, 1000])
