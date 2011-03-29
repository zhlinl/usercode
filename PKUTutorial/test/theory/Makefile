CXXFLAGS = `root-config --cflags`
LDFLAGS = -L/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.25.02-cms6/lib `root-config --libs` -lRooFit

all: theoryPlot1S.gif theoryPlot2S.gif theoryPlot3S.gif theoryPlot1Srap.gif theoryPlot2Srap.gif theoryPlot3Srap.gif

theoryPlot1S.gif: theoryPlot1S.C ../GenCrossSection/genXSupsilon1S.root
	root -l -b -q theoryPlot1S.C

theoryPlot2S.gif: theoryPlot2S.C ../GenCrossSection/genXSupsilon2S.root
	root -l -b -q theoryPlot2S.C

theoryPlot3S.gif: theoryPlot3S.C ../GenCrossSection/genXSupsilon3S.root
	root -l -b -q theoryPlot3S.C

theoryPlot1Srap.gif: theoryPlot1Srap.C ../GenCrossSection/genXSupsilon1Srap.root
	root -l -b -q theoryPlot1Srap.C

theoryPlot2Srap.gif: theoryPlot2Srap.C ../GenCrossSection/genXSupsilon2Srap.root
	root -l -b -q theoryPlot2Srap.C

theoryPlot3Srap.gif: theoryPlot3Srap.C ../GenCrossSection/genXSupsilon3Srap.root
	root -l -b -q theoryPlot3Srap.C


