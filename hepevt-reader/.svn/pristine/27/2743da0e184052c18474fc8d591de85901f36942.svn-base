#!/bin/bash
#
#

# location of data
datadir=`pwd`

# basename of the output of GENIE in hepevt format: basename.hepevt.gz
# will also be used as base for i3 files: basename.pi.i3.gz
basename=genie-output
runnum_label="000001"

i3_dir=/path/to/icesim
i3_env=${i3_dir}/build/env-shell.sh

$i3_env python GENIE-p1.py $datadir $basename $runnum_label

$i3_env python GENIE-p2.py $datadir $basename $runnum_label

$i3_env python GENIE-p3.py $datadir $basename $runnum_label
