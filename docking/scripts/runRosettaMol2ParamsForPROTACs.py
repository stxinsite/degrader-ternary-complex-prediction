"""
    call the Rosetta mol_to_params.py for every *_lig.sdf file in the speficied folder;
    This will create the necessary PDB and PARAMS files to be used in further
    Rosetta:minimize_ppi software
"""

import  os
import  sys


def  main( inputFolder ):
    for file1 in os.listdir( inputFolder ):
        if "_lig.sdf" == file1[ len( file1 ) - 8 : len( file1 ) ]:
            file1Path  =  os.path.join( inputFolder, file1 )
            outFilePrefix  =  file1.split( "_lig.sd" )[ 0 ] + "_lig"
            outFilePrefixPath  =  os.path.join( inputFolder, outFilePrefix )
            cmd  =  "$MOL2PARAMS " + file1Path + " -n TRN --clobber -p " + outFilePrefixPath
            os.system( cmd )

if "__main__" == __name__:
    if 2 != len( sys.argv ):
        print( "please provide a full path to folder with the *_lig.sdf files" )
        sys.exit()
    inputFolder  =  sys.argv[ 1 ]
    main( inputFolder )
