""" get chain ID from a PDB file
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


def  getChainIDPDB( file1 ):
    listLines  =  getListLines( file1 )
    for line in listLines:
        if ( "ATOM  " == line[ 0 : 6 ]\
        or "HETATM" == line[ 0 : 6 ] )\
        and 54 <= len( line ):
            return  line[ 21 ]
    return  ""

def  getChainIDPDBLigand( file1, chainIDProtein = 'A' ):
    listLines  =  getListLines( file1 )
    for line in listLines:
        if "HETATM" == line[ 0 : 6 ]\
        and 54 <= len( line ):
            if " " != line[ 21 ]:
                return  line[ 21 ]
            else:
                return  chr( ord( chainIDProtein ) + 1 )
    return  ""

def  fixChainIDs( file1, outputFilePath ):
    listLines  =  getListLines( file1 )
    chainIDcurrent  =  ""
    chainIDnext     =  ""
    chainID         =  ""
    listChainIDs  =  []
    chainIDChanged  =  False
    if True == os.path.exists( outputFilePath ):
        os.remove( outputFilePath )
    with open( outputFilePath, "a" ) as out:
        for i in range( 0, len( listLines ) - 1 ):
            line  =  listLines[ i ]
            lineNext  =  listLines[ i + 1 ]
            if ( "ATOM  " == line[ 0 : 6 ]\
            or "HETATM" == line[ 0 : 6 ] )\
            and 54 <= len( line ):
                chainIDcurrent  =  line[ 21 ]
                chainIDnext     =  lineNext[ 21 ]
                #chainID         =  chainIDnext
                if 0 == len( listChainIDs ):
                    listChainIDs.append( chainIDnext )
                if chainIDcurrent != chainIDnext:
                    if False == ( chainIDnext in listChainIDs ):
                        listChainIDs.append( chainIDnext )
                        chainIDChanged  =  False
                    else:
                        chainIDChanged  =  True
                        chainID  =  chr( ord( chainIDnext ) + 1 )
                if False == chainIDChanged:
                    lineToSave  =  line + "\n"
                    out.write( lineToSave )
                elif chainIDcurrent != chainIDnext:
                    lineToSave  =  line + "\n"
                    out.write( lineToSave )
                    lineToSave  =  lineNext[ 0:21 ] + chainID + lineNext[ 22 : len( lineNext ) ] + "\n"
                    out.write( lineToSave )
                    i  =  i + 1
                    if False == ( chainID in listChainIDs ):
                        listChainIDs.append( chainID )
                else:
                    lineToSave  =  lineNext[ 0:21 ] + chainID + lineNext[ 22 : len( lineNext ) ] + "\n"
                    out.write( lineToSave )
                    if False == ( chainID in listChainIDs ):
                        listChainIDs.append( chainID )
    return  listChainIDs
