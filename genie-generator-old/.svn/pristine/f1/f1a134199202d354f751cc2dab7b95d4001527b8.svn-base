from icecube import icetray, dataclasses

from oneweightcalc import UpdateOneWeightFrame

class OneWeightUpdator(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox("OutBox")

    def Configure(self):
        pass
        
    def Physics(self, frame):
        UpdateOneWeightFrame(frame)
        self.PushFrame(frame)
