#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/{complex.parm7,rest.in} .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF/seg.rst ./parent.rst
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF ./parent.rst
fi


export CUDA_DEVICES=(`echo $CUDA_VISIBLE_DEVICES_ALLOCATED | tr , ' '`)
export CUDA_VISIBLE_DEVICES=${CUDA_DEVICES[$WM_PROCESS_INDEX]}

echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES_ALLOCATED = " $CUDA_VISIBLE_DEVICES_ALLOCATED
echo "RUNSEG.SH: WM_PROCESS_INDEX = " $WM_PROCESS_INDEX
echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES


$PMEMD -O -i md.in   -p complex.parm7  -c parent.rst -r seg.rst -x seg.nc  -o seg.log  -inf seg.info  -e seg.mden

COMMAND="parm $WEST_SIM_ROOT/common_files/complex.parm7\n"
COMMAND="${COMMAND} trajin parent.rst\n"
COMMAND="${COMMAND} trajin seg.rst\n"
COMMAND="${COMMAND} autoimage \n"
COMMAND="${COMMAND} nativecontacts name SVI  ':34,35,36,37,39,40,41,42,43,45,46,47,48,49,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95 & !@H*'   ':116,118,119,120,121,122,123,124,125,126,128,129 & !@H*'  byresidue mindist distance 4.5   out Contacts_Smarca_VHL.dat   writecontacts ContactsInfo_SVI.dat\n"
COMMAND="${COMMAND} nativecontacts name SPI  ':34,35,36,37,39,40,41,42,43,45,46,47,48,49,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95 & !@H*'   ':FWZ & !@H*'  byresidue mindist distance 4.5   out Contacts_Smarca_Protac.dat   writecontacts ContactsInfo_SPI.dat\n"
COMMAND="${COMMAND} nativecontacts name VPI  ':116,118,119,120,121,122,123,124,125,126,128,129 & !@H*'   ':FWZ & !@H*'   byresidue mindist distance 4.5   out Contacts_VHL_Protac.dat   writecontacts ContactsInfo_VPI.dat\n"
COMMAND="${COMMAND} go\n"

# This wasn't working at one random iteration, thus we want to monitor

sleep 2
echo -e $COMMAND | $CPPTRAJ > cpptraj1.out
if [ ! -f "Contacts_Smarca_VHL.dat" ] || [ ! -f "Contacts_Smarca_Protac.dat" ] || [ ! -f "Contacts_VHL_Protac.dat" ] ; then
sleep 2
echo -e $COMMAND | $CPPTRAJ > cpptraj2.out
fi

sleep 1

paste <(cat Contacts_Smarca_VHL.dat | tail -n 2 | gawk '{print $2+$3}')   <(paste Contacts_Smarca_Protac.dat Contacts_VHL_Protac.dat  | tail -n 2  |  gawk '{print $2+$3+$6+$7}')   >>$WEST_PCOORD_RETURN
paste <(cat Contacts_Smarca_VHL.dat | tail -n 2 | gawk '{print $2+$3}')   <(paste Contacts_Smarca_Protac.dat Contacts_VHL_Protac.dat  | tail -n 2  |  gawk '{print $2+$3+$6+$7}')   >>pcoord.dat

# Clean up
rm -f md.in rest.in parent.rst complex.parm7 seg.log seg.info seg.mden 

