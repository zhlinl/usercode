first: effTrig.root effId.root effTQ.root

%.root: %.py ../../python/tagAndProbe.py
	cmsRun $*.py >& $*.log&

all: effTrigJpsiMC.root effTrigUpsMC.root effTrigJpsiTnP.root effTrig.root

