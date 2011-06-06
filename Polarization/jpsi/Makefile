CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)
ROOTCINT=rootcint

CXX+= -O3

CXX+=-I/usr/local/include/root
LIBS+=-L/usr/local/lib/root

FIT=AcceptanceMaps.o LifetimeModel.o MassModel.o PolarizationModel.o \
	CompositeModelBuilder.o 

#RooPolarizationPdf.o RooPolarizationConstraint.o

#DATA=

%.o : %.cc
	$(CXX) -c $<

all: polarizationFit oneBinTest injectAcceptanceMaps extractFSRShape \
	extractPromptLifetimeShape extractBackgroundShapes \
	extractNonPromptLifetimeShape makeAcceptancePdfs	
#	PlotMass PlotLifetime PlotPolarization \
#	PlotFinalMass PlotFinalLifetime \
#	PlotParameters	\


makeAcceptancePdfs: makeAcceptancePdfs.cc $(FIT) 
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit 

polarizationFit: polarizationFit.cc $(FIT) 
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit -rRooPolarizationFit.o -rRooPolarizationConstraint.o

oneBinTest: oneBinTest.cc $(FIT) 
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

extractFSRShape: extractFSRShape.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

injectAcceptanceMaps: injectAcceptanceMaps.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

extractPromptLifetimeShape: extractPromptLifetimeShape.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

extractNonPromptLifetimeShape: extractNonPromptLifetimeShape.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

extractBackgroundShapes: extractBackgroundShapes.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

#PlotMass: PlotMass.cc $(FIT)
#	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

#PlotLifetime: PlotLifetime.cc $(FIT)
#	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

#PlotPolarization: PlotPolarization.cc $(FIT)
#	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

#PlotFinalMass: PlotFinalMass.cc $(FIT)
#	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

#PlotFinalLifetime: PlotFinalLifetime.cc $(FIT)
#	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

#PlotParameters: PlotParameters.cc $(FIT)
#	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

clean: 
	rm polarizationFit oneBinTest injectAcceptanceMaps extractFSRShape \
	extractPromptLifetimeShape extractBackgroundShapes \
	extractNonPromptLifetimeShape *.o #PlotMass PlotLifetime \
#	PlotFinalMass PlotFinalLifetime \
#	PlotPolarization *.o