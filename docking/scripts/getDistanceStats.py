"""
    outputs distances between ends of the linker.
    make sure that the input folder contains only PDBs with the conformations of a linker.
"""

import  os
import  sys
import  math
from gen_conformers.getListLines import getListLines


def  main( inputPath, atomName1, atomName2, outputFile ):
    with open( outputFile, "a" ) as out:
        listDistances  =  []
        for file1 in os.listdir( inputPath ):
            if ".pdb" == file1[ len( file1 ) - 4 : len( file1 ) ]:
                file1Path  =  os.path.join( inputPath, file1 )
                listLinesPDB  =  getListLines( file1Path )
                p1  =  []
                p2  =  []
                for line in listLinesPDB:
                    if 3 == len( p1 )\
                    and 3 == len( p2 ):
                        break
                    if "HETATM" == line[ 0:6 ]:
                        atomName  =  line[ 12:17 ].strip()
                        if atomName == atomName1\
                        or atomName == atomName2:
                            x  =  float( line[ 30:39 ].strip() )
                            y  =  float( line[ 39:47 ].strip() )
                            z  =  float( line[ 47:55 ].strip() )
                            if 0 == len( p1 ):
                                p1  =  [ x, y, z ]
                            else:
                                p2  =  [ x, y, z ]
                if 3 == len( p1 )\
                and 3 == len( p2 ):
                    dist2  =  ( ( p1[ 0 ] - p2[ 0 ] ) ** 2 +
                            ( p1[ 1 ] - p2[ 1 ] ) ** 2 + 
                            ( p1[ 2 ] - p2[ 2 ] ) ** 2 )
                    lineToSave  =  file1Path + "\t" + str( math.sqrt( dist2 ) ) + "\n"
                    out.write( lineToSave )
                    listDistances.append( math.sqrt( dist2 ) )
        meanDist  =  0
        sdDist  =  0
        for dist in listDistances:
            meanDist  =  meanDist + dist
        meanDist  =  meanDist / len( listDistances )
        for dist in listDistances:
            sdDist  =  sdDist + ( dist - meanDist ) ** 2
        sdDist  =  math.sqrt( sdDist / len( listDistances ) )
        out.write( str( meanDist ) + "\t" + str( sdDist ) + "\t" + str( len( listDistances ) ) + "\n" )



if "__main__" == __name__:
    if 5 != len( sys.argv ):
        print( "please provide:\n(i). path to the folder with input PDB files;\n(ii). atom name 1;\n(iii). atom name 2;\n(iv). output file." )
        exit()
    inputPath  =  sys.argv[ 1 ]
    atomName1  =  sys.argv[ 2 ]
    atomName2  =  sys.argv[ 3 ]
    outputFile  =  sys.argv[ 4 ]
    if False == os.path.exists( inputPath )\
    or False == os.path.isdir( inputPath ):
        print( f'folder {inputPath} does not exist' )
        exit()
    if True == os.path.exists( outputFile ):
        os.remove( outputFile )
    main( inputPath, atomName1, atomName2, outputFile )
