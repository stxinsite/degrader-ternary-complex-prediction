"""
    takes an input folder and generates a list of PDB files in it
"""

import  os
import  sys
from gen_conformers.getListLines import getListLines


def  main( inputFolder, outputFile ):
    with open( outputFile, "a" ) as out:
        for file1 in os.listdir( inputFolder ):
            if "pdb" == file1.split( "." )[ len( file1.split( "." ) ) - 1 ]:
                filePath  =  os.path.join( inputFolder, file1 )
                line  =  filePath + "\n"
                out.write( line )


if "__main__" == __name__:
    if 3 != len( sys.argv ):
        print( "please provide:\n(i). folder with PDB files;\n(ii). output text file;" )
        exit()
    inputFolder  =  sys.argv[ 1 ]
    outputFile   =  sys.argv[ 2 ]
    if False == os.path.exists( inputFolder )\
    or False == os.path.isdir( inputFolder ):
        print( f'folder {inputFolder} does not exist' )
        exit()
    if True == os.path.exists( outputFile ):
        os.remove( outputFile )
    main( inputFolder, outputFile )
