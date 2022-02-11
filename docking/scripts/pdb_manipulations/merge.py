""" merge two PDB files into a single one.
    This is a step required by ROSETTADock: docking_protocol.linuxgccrelease
    """

import  os


# to avoid dependences, assume that the user supplies legitimate PDB files
# and do the job on a low level:
def  getListLines( file1 ):
    try:
        listLines  =  [ line.rstrip( "\t\n" ) for line in open( file1, "r" ) ]
        return  listLines
    except FileNotFoundError:
        print( 'file {file1} not found.' )
    return  []


def  mergePDBs( pdb1, pdb2, outPDB ):
    listLinesPDB1  =  getListLines( pdb1 )
    listLinesPDB2  =  getListLines( pdb2 )
    if True == os.path.exists( outPDB ):
        os.remove( outPDB )
    with open( outPDB, "a" ) as out:
        for line in listLinesPDB1:
            if "ATOM  " == line[ 0:6 ]\
            or "HETATM" == line[ 0:6 ]:
                lineToSave  =  line + "\n"
                out.write( lineToSave )
        for line in listLinesPDB2:
            if "ATOM  " == line[ 0:6 ]\
            or "HETATM" == line[ 0:6 ]:
                lineToSave  =  line + "\n"
                out.write( lineToSave )

