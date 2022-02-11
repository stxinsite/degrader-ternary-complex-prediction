"""
    This reads all the *_lig.pdb files
    and tranfsorms them into the corresponding SDF files
    that will be further pipelined into the
    ROSETTA: molecule_to_params.py
"""

import  os
import  sys
import  shutil


def  main( inputFolder ):
    for file1 in os.listdir( inputFolder ):
        if "_lig.pdb" == file1[ len( file1 ) - 8 : len( file1 ) ]:
            file1Path  =  os.path.join( inputFolder, file1 )
            outputFileName  =  file1.split( ".pd" )[ 0 ] + ".sdf"
            outputFileNamePath  =  os.path.join( inputFolder, outputFileName )
            if True == os.path.exists( outputFileNamePath ):
                os.remove( outputFileNamePath )
            cmd  =  "python /bgfs01/insite/taras.dauzhenka/data/6HAX/scripts/fixBondsInModel.py " + file1Path + " " + outputFileNamePath
            os.system( cmd )


if "__main__" == __name__:
    if 2 != len( sys.argv ):
        print( "please provide a full path to folder with *_lig.pdb files" )
        sys.exit()
    folderPath  =  sys.argv[ 1 ]
    if False == os.path.exists( folderPath )\
    or False == os.path.isdir( folderPath ):
        print( f'folder {folderPath} does not exist' )
        sys.exit()
    main( folderPath )

