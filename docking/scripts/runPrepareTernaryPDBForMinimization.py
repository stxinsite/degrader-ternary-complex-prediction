import  os
import  sys


def  main( inputFolder ):
    for file1 in os.listdir( inputFolder ):
        if 13 < len( file1 )\
        and "_lig_0001.pdb" == file1[ len( file1 ) - 13 : len( file1 ) ]:
            inputPROTACPath  =  os.path.join( inputFolder, file1 )
            inputTernary  =  file1.split( "_lig_0001" )[ 0 ] + ".pdb"
            inputTernaryPath  =  os.path.join( inputFolder, inputTernary )
            outputTernary  =  file1.split( "_lig_0001" )[ 0 ] + "_fixed.pdb"
            outputTernaryPath  =  os.path.join( inputFolder, outputTernary )
            cmd  =  "python scripts/prepareTernaryPDBForMinimization.py " + inputTernaryPath + " " + inputPROTACPath + " " + outputTernaryPath
            os.system( cmd )

if "__main__" == __name__:
    if 2 != len( sys.argv ):
        print( "please provide:\n(i). path to the folder with ternary PDB files and PROTAC PDB files" )
        sys.exit()
    inputFolder  =  sys.argv[ 1 ]
    if False == os.path.exists( inputFolder )\
    or False == os.path.isdir( inputFolder ):
        print( f'folder {inputFolder} does not exist' )
        sys.exit()
    main( inputFolder )

