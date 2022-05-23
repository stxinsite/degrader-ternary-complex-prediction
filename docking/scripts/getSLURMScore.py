import  os
import  sys

def  main( inputFolder, inputFolderPARAMS, slurmScript ):
    if True == os.path.exists( slurmScript ):
        os.remove( slurmScript )
    with open( slurmScript, "a" ) as out:
        line  = "#!/bin/bash\n\
#SBATCH --nodes=1 --ntasks=1\n\
#SBATCH --job-name SiTX-0040002_score.0.3A__taras.dauzhenka\n\
#SBATCH --partition=project\n\
#SBATCH --qos=maxjobs\n\
#SBATCH --output=/bgfs01/insite/taras.dauzhenka/data/SM2_VHL_SiTX/006.minimize/SiTX-0040002/0.3A/job_1/score_dock_%A_%a.out\n\
#SBATCH --error=/bgfs01/insite/taras.dauzhenka/data/SM2_VHL_SiTX/006.minimize/SiTX-0040002/0.3A/job_1/score_dock_%A_%a.err\n\
#SBATCH --array=0-29\n\
#SBATCH --cpus-per-task=2\n\n"
        out.write( line )
        listTernary  =  []
        listPARAMS   =  []
        listLog      =  []
        listScore    =  []
        for file1 in os.listdir( inputFolder ):
            if "6hax_" == file1[ 0 : 5 ]\
            and ".noH.pdb" == file1[ len(file1)-8:len(file1) ]:
                ternaryFilePath  =  os.path.join( inputFolder, file1 )
                paramsFile  =  file1 + "_protac.pdb_lig.params"
                paramsFilePath  =  os.path.join( inputFolderPARAMS, paramsFile )
                logFile  =  file1.split( ".noH.pd" )[ 0 ] + ".log"
                logFilePath  =  os.path.join( inputFolder, logFile )
                scoreFile  =  file1.split( ".noH.pd" )[ 0 ] + ".score"
                scoreFilePath  =  os.path.join( inputFolder, scoreFile )
                listLog.append( logFilePath )
                listTernary.append( ternaryFilePath )
                listPARAMS.append( paramsFilePath )
                listScore.append( scoreFilePath )
        line  =  "i=${SLURM_ARRAY_TASK_ID}\n"
        out.write( line + "\n" )
        line  =  "arrayTernary=("
        for ternary in listTernary:
            line  =  line + " " + ternary
        line  =  line + ")\n"
        out.write( line )
        line  =  "arrayPARAMS=("
        for params in listPARAMS:
            line  =  line + " " + params
        line  =  line + ")\n"
        out.write( line )
        line  =  "arrayLog=("
        for log in listLog:
            line  =  line + " " + log
        line  =  line + ")\n"
        out.write( line )
        line  =  "arrayScoreFile=("
        for score in listScore:
            line  =  line + " " + score
        line  =  line + ")\n"
        out.write( line )
        out.write( "\n" )
#        out.write( "cd /bgfs01/insite/taras.dauzhenka/data/SM2_VHL_SiTX/006.minimize/SiTX-0040002/0.3A/job_1/\n" )
        cmd  =  "time /software/ROSETTA/rosetta_bin_linux_2019.35.60890_bundle/main/source/bin/score_jd2.static.linuxgccrelease -database $ROSETTA_DATABASE -s ${arrayTernary[$i]} -extra_res_fa ${arrayPARAMS[$i]} -out:file:scorefile ${arrayScoreFile[$i]} > ${arrayLog[$i]}\n"
        out.write( cmd )
        out.write( "\nwait\necho\necho \"All done.\"\n" )
        print( len( listTernary ) )



if "__main__" == __name__:
    if 4 != len( sys.argv ):
        print( "please provide:\n(i). input folder with PDB files;\n(ii). input folder with PARAMS files;\n(iii). file to save slurm script\n" )
        sys.exit()
    inputFolder  =  sys.argv[ 1 ]
    inputFolderPARAMS = sys.argv[ 2 ]
    slurmScript  =  sys.argv[ 3 ]
    if False == os.path.exists( inputFolder )\
    or False == os.path.isdir( inputFolder )\
    or False == os.path.exists( inputFolderPARAMS ):
        print( f'folder {inputFolder} does not exist' )
        sys.exit()
    if True == os.path.exists( slurmScript ):
        os.remove( slurmScript )
    main( inputFolder, inputFolderPARAMS, slurmScript )
