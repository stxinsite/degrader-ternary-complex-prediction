"""
    changes all HETATM into ATOM labels
"""

import  os


def  getListLines( file1 ):
    listLines  =  [ line.rstrip( "\t\n" ) for line in open( file1, "r" ) ]
    return  listLines


def  changeHETATM( inputFile, outputFile ):
    if True == os.path.exists( outputFile ):
        os.remove( outputFile )
    with open( outputFile, "a" ) as out:
        listLines  =  getListLines( inputFile )
        for line in listLines:
            if "HETATM" == line[ 0:6 ]:
                lineToSave  =  "ATOM  " + line[ 6 : len( line ) ] + "\n"
                out.write( lineToSave )
            else:
                out.write( line + "\n" )
