### For W analyses ######################

class wModel(PhysicsModel):
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        # --- Signal Strength as only POI --- 
        self.modelBuilder.doVar("r[1,0,20]");
        self.modelBuilder.doSet("POI","r")
    def getYieldScale(self,bin,process):
        "Return the name of a RooAbsReal to scale this yield by or the two special values 1 and 0 (don't scale, and set to zero)"
        # Important - the expected yields for nEWK and nSig should already satisfy nEWK = cewk*nSig
        # that way they will both scale with r
        if self.DC.isSignal[process] or process=="EWK": return "r"
        else: return 1
