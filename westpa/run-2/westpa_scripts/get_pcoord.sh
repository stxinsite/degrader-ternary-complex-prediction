#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT/init_files

COMMAND="parm $WEST_SIM_ROOT/common_files/complex.parm7\n"
COMMAND="${COMMAND} trajin $WEST_STRUCT_DATA_REF \n"
COMMAND="${COMMAND} nativecontacts name SVI  ':34,35,36,37,39,40,41,42,43,45,46,47,48,49,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95 & !@H*'   ':116,118,119,120,121,122,123,124,125,126,128,129 & !@H*'  byresidue mindist distance 4.5   out Contacts_Smarca_VHL.dat   writecontacts ContactsInfo_SVI.dat\n" 
COMMAND="${COMMAND} nativecontacts name SPI  ':34,35,36,37,39,40,41,42,43,45,46,47,48,49,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95 & !@H*'   ':FWZ & !@H*'  byresidue mindist distance 4.5   out Contacts_Smarca_Protac.dat   writecontacts ContactsInfo_SPI.dat\n" 
COMMAND="${COMMAND} nativecontacts name VPI  ':116,118,119,120,121,122,123,124,125,126,128,129 & !@H*'   ':FWZ & !@H*'   byresidue mindist distance 4.5   out Contacts_VHL_Protac.dat   writecontacts ContactsInfo_VPI.dat\n" 
COMMAND="${COMMAND} go\n"
echo -e ${COMMAND} | $CPPTRAJ

sleep 1 

paste <(cat Contacts_Smarca_VHL.dat | tail -n 1 | gawk '{print $2+$3}')   <(paste Contacts_Smarca_Protac.dat Contacts_VHL_Protac.dat  | tail -n 1  |  gawk '{print $2+$3+$6+$7}')   >>$WEST_PCOORD_RETURN
paste <(cat Contacts_Smarca_VHL.dat | tail -n 1 | gawk '{print $2+$3}')   <(paste Contacts_Smarca_Protac.dat Contacts_VHL_Protac.dat  | tail -n 1  |  gawk '{print $2+$3+$6+$7}')   >>pcoord.dat

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
