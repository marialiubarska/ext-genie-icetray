# GENIE generation app for IceCube

Test (works, I did not check contents of output file yet):
```
$MYGENIE/genie/gicevgen -r 1 -o test -n 5 -e 10,100 -t 1000080160[0.95],1000010010[0.05] -z 0,90 -a 0,360 -R 1300 -L 5000 --nu-type NuMu --nu-fraction 0.7 --gamma 2 --cross-sections /home/mliubar/Software/genie_workspace/genie-generator/xsec_splines/GENIE_2_12_8_Water_splines.xml --force-singleprob-scale
```