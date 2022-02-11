def  getListLines( file1 ):
    listLines  =  []
    try:
        listLines  =  [ line.rstrip( "\t\n" ) for line in open( file1, "r" ) ]
    except FileNotFoundError:
        print( f'file {file1} not found' )
    return  listLines
