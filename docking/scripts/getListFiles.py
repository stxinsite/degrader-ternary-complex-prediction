import os
import sys

if "__main__" == __name__:
	fol = sys.argv[1]
	output = sys.argv[2]
	if True == os.path.exists(output):
		os.remove(output)
	with open(output, "a") as out:
		for file1 in os.listdir(fol):
			if ".pdb" in file1:
				file1Path = os.path.join(fol, file1)
				out.write(file1Path + "\n")
