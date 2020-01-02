from icecube import icetray, dataio

from oneweightupdator import OneWeightUpdator 

infile=sys.argv[1]
outfile=sys.argv[2]
print infile
print outfile

tray = I3Tray()

tray.AddModule("I3Reader","reader",
               Filename = infile)

tray.AddModule(OneWeightUpdator,"update")

tray.AddModule("I3Writer","writer",Filename=outfile)

tray.AddModule("TrashCan","mycan")

tray.Execute()

tray.Finish()

print "Done!"
