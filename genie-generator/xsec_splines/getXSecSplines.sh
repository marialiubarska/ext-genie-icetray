#!/bin/bash

# The Following code grabs the external files for
# the cross-sections needed for the IceCube
# implementation of GENIE.

# Check to see which terminal based download
# programs are available
d1=`type -t wget`
d2=`type -t curl`

# Location of the two important files.
file1='http://w3.iihe.ac.be/~mdier/icecube/genie/xsec_splines/gxspl-water-v2.6.0.xml'
file2='http://w3.iihe.ac.be/~mdier/icecube/genie/xsec_splines/gxspl-water-v2.6.0.root'

if [ -n "$d1" ]; then
    echo "Using wget to download the GENIE cross-section files"
    wget $file1
    wget $file2
elif [ -n "$d2" ]; then
    echo "Using curl to download the GENIE cross-section files"
    curl $file1 -O
    curl $file2 -O
else
    echo 'The following terminal download programs could not be found: wget, curl'
    echo 'Download the files manually.'
fi




