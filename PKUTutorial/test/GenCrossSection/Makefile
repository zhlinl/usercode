CXXFLAGS = `root-config --cflags`
LDFLAGS = -L/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.25.02-cms6/lib `root-config --libs` -lRooFit

all: genXSupsilon1S.root genXSupsilon2S.root genXSupsilon3S.root genXSupsilon1Srap.root genXSupsilon2Srap.root genXSupsilon3Srap.root
genUpsilon1S.root: Upsilon1S/res/CMSSW_1.stdout
	hadd -f genUpsilon1S.root Upsilon1S/res/genUpsilon1S*.root

genUpsilon2S.root: Upsilon2S/res/CMSSW_1.stdout
	hadd -f genUpsilon2S.root Upsilon2S/res/genUpsilon2S*.root

genUpsilon3S.root: Upsilon3S/res/CMSSW_1.stdout
	hadd -f genUpsilon3S.root Upsilon3S/res/genUpsilon3S*.root

genXSupsilon1S.root: genXSupsilon1S.C genUpsilon1S.root
	root -l -b -q genXSupsilon1S.C

genXSupsilon2S.root: genXSupsilon2S.C genUpsilon2S.root
	root -l -b -q genXSupsilon2S.C

genXSupsilon3S.root: genXSupsilon3S.C genUpsilon3S.root
	root -l -b -q genXSupsilon3S.C

genXSupsilon1Srap.root: genXSupsilon1Srap.C genUpsilon1S.root
	root -l -b -q genXSupsilon1Srap.C

genXSupsilon2Srap.root: genXSupsilon2Srap.C genUpsilon2S.root
	root -l -b -q genXSupsilon2Srap.C

genXSupsilon3Srap.root: genXSupsilon3Srap.C genUpsilon3S.root
	root -l -b -q genXSupsilon3Srap.C


