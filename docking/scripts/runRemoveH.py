import os

src = os.getcwd()

for file1 in os.listdir(src):
	if ".pdb" == file1[len(file1)-4:len(file1)]:
		outFile = file1[0:len(file1)-4] + ".noH.pdb"
		outFilePath = os.path.join(src, outFile)
		if True == os.path.exists(outFilePath):
			os.remove(outFilePath)
		file1Path = os.path.join(src, file1)
		cmd = "python removeH.py " + file1Path + " " + outFilePath
		os.system(cmd)
