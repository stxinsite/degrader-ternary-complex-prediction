import  os
import  sys


def  main(inputFolder, inputFolderRef, outputFile):
    with open( outputFile, "a") as out:
        out.write("#!/bin/bash\n" )
        out.write("#SBATCH --nodes=1 --ntasks=1\n")
        out.write("#SBATCH --job-name tCAPRI.main\n")
        out.write("#SBATCH --partition=project\n")
        out.write("#SBATCH --qos=maxjobs\n")
        out.write("#SBATCH --output=tCAPRI.main_%A_%a.out\n")
        out.write("#SBATCH --error=tCAPRI.main_%A_%a.err\n")
        out.write("#SBATCH --array=1-473\n")
        out.write("#SBATCH --cpus-per-task=2\n\n")
        out.write("# >>> conda initialize >>>\n")
        out.write("# !! Contents within this block are managed by 'conda init' !!\n")
        out.write("__conda_setup=\"$('/bgfs01/insite/taras.dauzhenka/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)\"\n")
        out.write("if [ $? -eq 0 ]; then\n")
        out.write("    eval \"$__conda_setup\"\n")
        out.write("else\n")
        out.write("    if [ -f \"/bgfs01/insite/taras.dauzhenka/anaconda3/etc/profile.d/conda.sh\" ]; then\n")
        out.write("        . \"/bgfs01/insite/taras.dauzhenka/anaconda3/etc/profile.d/conda.sh\"\n")
        out.write("    else\n")
        out.write("        export PATH=\"/bgfs01/insite/taras.dauzhenka/anaconda3/bin:$PATH\"\n")
        out.write("    fi\n")
        out.write("fi\n")
        out.write("unset __conda_setup\n")
        out.write("# <<< conda initialize <<<\n\n")
        out.write("source /bgfs01/insite/taras.dauzhenka/anaconda3/etc/profile.d/conda.sh\n")
        out.write("conda init bash\n\n")
        out.write("arrayFiles=(")
        count = 0
        for fileRef in os.listdir(inputFolderRef):
            if "_ABXYZ.pdb" == fileRef[len(fileRef) - 10:len(fileRef)]\
            and True == ("contact.cluster_rep" in fileRef):
                fileRefPath = os.path.join(inputFolderRef, fileRef)
                for file1 in os.listdir(inputFolder):
                    if "mini_" == file1[0:5]\
                    and "_fixed.AB.pdb" == file1[len(file1) - 13:len(file1)]:
                        file1Path = os.path.join(inputFolder, file1)
                        out.write(fileRefPath + "+" + file1Path + " ")
                        count = count + 1
        print(count)
        out.write(")\n\n")
        out.write("i=${SLURM_ARRAY_TASK_ID}\n")
        out.write("IFS=+\n")
        line = "python /bgfs01/insite/taras.dauzhenka/repos/ternary/ternary_CAPRI/getCAPRITernary.py ${arrayFiles[$i]} 6HAX_cst1.tern0.3.mini.ETcapri 6HAX_protac_info_ETC.txt A_B_X_Y_Z $i\n"
        out.write( line )
        out.write("\n")
        out.write("wait\necho\necho \"All done.\n")


if "__main__" == __name__:
    if 4 != len(sys.argv):
        print("wrong input")
        sys.exit()
    inputFolder = sys.argv[1]
    inputFolderRef = sys.argv[2]
    outputFile = sys.argv[3]
    if False == os.path.exists(inputFolder):
        print("wrong input folder")
        sys.exit()
    if False == os.path.exists(inputFolderRef):
        print("wrong input reference folder")
        sys.exit()
    if True == os.path.exists(outputFile):
        os.remove(outputFile)
    main(inputFolder, inputFolderRef, outputFile)

