#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT || exit 1

cp west.h5  west_backup.h5

ITER=$(printf "%06d" $WEST_CURRENT_ITER)
ITER_PREV=$(printf "%06d" $(($WEST_CURRENT_ITER-1)))
tar -cf seg_logs/$ITER.tar seg_logs/$ITER-*.log
rm  -f  seg_logs/$ITER-*.log
rm -rf  traj_segs/$ITER_PREV/*/*rst


