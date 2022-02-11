import  os
import  sys


def  getListLines( file1 ):
    listLines  =  [ line.rstrip( "\t\n" ) for line in open( file1, "r" ) ]
    return  listLines


def  main( fileTernary, fileLigand, fileTernaryOut ):
    listLinesTernary  =  getListLines( fileTernary )
    listLinesLigand   =  getListLines( fileLigand )
    listLinesTernaryFixed  =  []
    if True == os.path.exists( fileTernaryOut ):
        os.remove( fileTernaryOut )
    with open( fileTernaryOut, "a" ) as out:
        for line in listLinesTernary:
            if "LG" != line[ 17:19 ]:
                out.write( line + "\n" )
        for line in listLinesLigand:
            out.write( line + "\n" )

if "__main__" == __name__:
    if 4 != len( sys.argv ):
        print( "please provide:\n(i). input PDB with ternary structure;\n(ii). input PDB file ligand/PROTAC;\n(iii). output PDB to save ternary structure;" )
        sys.exit()
    inputTernary  =  sys.argv[ 1 ]
    inputPROTAC   =  sys.argv[ 2 ]
    outputTernary =  sys.argv[ 3 ]
    if False == os.path.exists( inputTernary ):
        print( f'file {inputTernary} does not exist' )
        sys.exit()
    if False == os.path.exists( inputPROTAC ):
        print( f'file {inputPROTAC} does not exist' )
        sys.exit()
    if True == os.path.exists( outputTernary ):
        os.remove( outputTernary )
    main( inputTernary, inputPROTAC, outputTernary )
