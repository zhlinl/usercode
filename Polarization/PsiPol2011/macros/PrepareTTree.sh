############ INPUTS ####################

# Define JobID
JobID=FrameworkTest_5Dec2012

Cdir=$PWD
# input arguments
for nState in 5; do		#1,2,3,Upsi(1S,2S,3S); 4,Jpsi 5,PsiPrime
for FidCuts in 11;do 		#defines the set of cuts to be used, see macros/polFit/effsAndCuts.h

rapMin=1     #takes bins, not actual values
rapMax=2     #if you only want to process 1 y bin, rapMax = rapMin
ptMin=1     #takes bins, not acutal values
ptMax=4     #if you only want to process 1 pt bin, ptMax = ptMin
Plotting=2   #plotting macro: 1 = plot all, 2 = plot mass, 3 = plot lifetime sidebands, 4 = plot lifetime singal region

rejectCowboys=true
RequestTrigger=true
MC=false
f_BG_zero=false
doCtauUncer=false

# input files
# In case of more input Files: define inputTreeX and adapt the line starting with inputTrees, at the moment up to 4 files implemented
#inputTree1=/scratch/knuenz/Polarization/RootInput/Psi/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Jpsi.root
#inputTree1=/scratch/knuenz/Polarization/RootInput/Psi/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Psi.root
if [ ${nState} -eq 4 ] 
then
inputTree1=/afs/cern.ch/user/z/zhlinl/work/polarization/TTree/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Jpsi.root
fi

if [ ${nState} -eq 5 ] 
then
inputTree1=/afs/cern.ch/user/z/zhlinl/work/polarization/TTree/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Psi.root
fi

################ EXECUTABLES #################
 
#following flags decide if the step is executed (1) or not (0):
execute_runData=0			           #independent of rapMin, rapMax, ptMin, ptMax
execute_runWorkspace=0			     #independent of rapMin, rapMax, ptMin, ptMax
execute_runMassFit=0			       #can be executed for different pt and y bins
execute_runLifetimeFit=0         #can be executed for different pt and y bins
execute_runPlotMassLifetime=0    #can be executed for different pt and y bins
execut_PlotFitPar=0              #independent of rapMin, rapMax, ptMin, ptMax
execute_runBkgHistos=0           #can be executed for different pt and y bins
execute_PlotCosThetaPhiBG=1 		 #This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)

#execute_PlotCosThetaPhiDistribution=0	#For each set of cuts, you can choose different values for FracLSB and nSigma

#################################

# Make directories
CutDir=${Cdir}/DataFiles/SetOfCuts${FidCuts}_${JobID}

WorkDir=${CutDir}/Psi$[nState-3]S                                  
mkdir -p ${CutDir}                                       
mkdir -p ${WorkDir}                                      
cp ../interface/commonVar_Psi$[nState-3]S.h ${WorkDir}/commonVar.h 

mkdir -p DataFiles
mkdir -p ${WorkDir}/tmpFiles/backupWorkSpace
mkdir -p ${WorkDir}/Figures
mkdir -p ${WorkDir}/PDF
mkdir -p ${WorkDir}/Fit

# Copy files to directory
cp Makefile ${WorkDir}/Makefile
cp ../interface/rootIncludes.inc ${WorkDir}/rootIncludes.inc

cp runData.cc ${WorkDir}/runData.cc
cp PolData.C ${WorkDir}/PolData.C
cp PolData.h ${WorkDir}/PolData.h
cp polFit/effsAndCuts.h ${WorkDir}/effsAndCuts.h

cp runWorkspace.cc ${WorkDir}/runWorkspace.cc
cp createWorkspace.C ${WorkDir}/createWorkspace.C

cp runMassFit.cc ${WorkDir}/runMassFit.cc
cp massFit.cc ${WorkDir}/massFit.cc

cp runLifetimeFit.cc ${WorkDir}/runLifetimeFit.cc
cp lifetimeFit.cc ${WorkDir}/lifetimeFit.cc
cp ../interface/calculatePar.cc ${WorkDir}/calculatePar.cc
cp ../interface/RooUtils.h ${WorkDir}/RooUtils.h

cp runPlotMassLifetime.cc ${WorkDir}/runPlotMassLifetime.cc
cp PlotMassLifetime.cc ${WorkDir}/PlotMassLifetime.cc

cp PlotFitPar.cc ${WorkDir}/PlotFitPar.cc

cp runBkgHistos.cc ${WorkDir}/runBkgHistos.cc
cp bkgHistos.C ${WorkDir}/bkgHistos.C
cp calcPol.C ${WorkDir}/calcPol.C

cp PlotCosThetaPhiBG.cc ${WorkDir}/PlotCosThetaPhiBG.cc

#cp PlotCosThetaPhiDistribution.cc ${CutDir}/PlotCosThetaPhiDistribution.cc

#cp runTrimEventContent.cc ${CutDir}/runTrimEventContent.cc
#cp TrimEventContent.C ${CutDir}/TrimEventContent.C

#cp runMeanPt.cc ${CutDir}/runMeanPt.cc
#cp calcMeanPt.C ${CutDir}/calcMeanPt.C

cp ../latex/Mass_sigma.tex ${WorkDir}/Mass_sigma.tex
cp ../latex/Lifetime_fitParameter.tex ${WorkDir}/Lifetime_fitParameter.tex
cp ../latex/myStyle.tex ${WorkDir}/myStyle.tex
cp ../latex/evaluateCtau.tex ${WorkDir}/evaluateCtau.tex
cp ../latex/NumEvents.tex ${WorkDir}/NumEvents.tex

cp ../latex/cosThetaPhi_$[nState-3]S_BG.tex ${WorkDir}/cosThetaPhi_$[nState-3]S_BG.tex 
cp ../latex/cosThetaPhi_$[nState-3]S_BG_highct.tex ${WorkDir}/cosThetaPhi_$[nState-3]S_BG_highct.tex 
cp ../latex/cosThetaPhi_$[nState-3]S_NPBG.tex ${WorkDir}/cosThetaPhi_$[nState-3]S_NPBG.tex
cp ../latex/MassLifetime_Psi$[nState-3]S.tex ${WorkDir}/MassLifetime_Psi$[nState-3]S.tex
#cp ../latex/massFits.tex ${CutDir}/massFits.tex
#cp ../latex/cosThetaPhi_BG.tex ${CutDir}/cosThetaPhi_BG.tex
#cp ../latex/cosThetaPhi.tex ${CutDir}/cosThetaPhi.tex

cd ${WorkDir}

make

inputTrees="inputTree=${inputTree1} inputTree=${inputTree2} inputTree=${inputTree3} inputTree=${inputTree4}"
if [ ${execute_runData} -eq 1 ]
then
./runData ${inputTrees} rejectCowboys=${rejectCowboys} FidCuts=${FidCuts} nState=${nState} MC=${MC} RequestTrigger=${RequestTrigger}
fi

if [ ${execute_runWorkspace} -eq 1 ]
then
./runWorkspace nState=${nState}
fi

if [ ${execute_runMassFit} -eq 1 ] 
then
rootfile=fit_Psi$[nState-3]S_rap${rapMin}_pt${ptMin}.root
cp tmpFiles/backupWorkSpace/${rootfile} tmpFiles/${rootfile}
cp runMassFit runMassFit_$[nState-3]S_rap${rapMin}_pt${ptMin}
./runMassFit_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState}
rm runMassFit_$[nState-3]S_rap${rapMin}_pt${ptMin}
fi

if [ ${execute_runLifetimeFit} -eq 1 ] 
then
cp runLifetimeFit runLifetimeFit_$[nState-3]S_rap${rapMin}_pt${ptMin}
./runLifetimeFit_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState}
rm runLifetimeFit_$[nState-3]S_rap${rapMin}_pt${ptMin}
fi

if [ ${execute_runPlotMassLifetime} -eq 1 ]
then
cp runPlotMassLifetime runPlotMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin}
#./runPlotMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} Plotting=${Plotting}
rm runPlotMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin}
pdflatex MassLifetime_Psi$[nState-3]S.tex
mv MassLifetime_Psi$[nState-3]S.pdf PDF/MassLifetime_Psi$[nState-3]S.pdf
fi

if [ ${execut_PlotFitPar} -eq 1 ]
then
./PlotFitPar nState=${nState} doCtauUncer=${doCtauUncer}
pdflatex Lifetime_fitParameter.tex
pdflatex Mass_sigma.tex
pdflatex evaluateCtau.tex
pdflatex evaluateCtau.tex
pdflatex NumEvents.tex
pdflatex NumEvents.tex
mv Lifetime_fitParameter.pdf PDF/Lifetime_fitParameter.pdf
mv Mass_sigma.pdf PDF/Mass_sigma.pdf
mv evaluateCtau.pdf PDF/evaluateCtau.pdf
mv NumEvents.pdf PDF/NumEvents.pdf
fi

if [ ${execute_runBkgHistos} -eq 1 ]
then
cp runBkgHistos runBkgHistos_$[nState-3]S_rap${rapMin}_pt${ptMin}
./runBkgHistos_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} MC=${MC} f_BG_zero=${f_BG_zero} doCtauUncer=${doCtauUncer}
rm runBkgHistos_$[nState-3]S_rap${rapMin}_pt${ptMin}
fi

if [ ${execute_PlotCosThetaPhiBG} -eq 1 ]
then
./PlotCosThetaPhiBG nState=${nState}
pdflatex cosThetaPhi_$[nState-3]S_BG.tex
pdflatex cosThetaPhi_$[nState-3]S_BG_highct.tex
pdflatex cosThetaPhi_$[nState-3]S_NPBG.tex
mv cosThetaPhi_$[nState-3]S_BG.pdf PDF/cosThetaPhi_$[nState-3]S_BG.pdf
mv cosThetaPhi_$[nState-3]S_BG_highct.pdf PDF/cosThetaPhi_$[nState-3]S_BG_highct.pdf
mv cosThetaPhi_$[nState-3]S_NPBG.pdf PDF/cosThetaPhi_$[nState-3]S_NPBG.pdf
fi

#if [ ${execute_PlotCosThetaPhiDistribution} -eq 1 ]
#then
#./PlotCosThetaPhiDistribution 1nState AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent=DataPath
#./PlotCosThetaPhiDistribution 2nState AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent=DataPath
#./PlotCosThetaPhiDistribution 3nState AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent=DataPath
#cp cosThetaPhi.tex AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent/cosThetaPhi.tex
#cd AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent
#pdflatex "\newcommand\TreeBinID{Ups1S}\input{cosThetaPhi.tex}"
#mv cosThetaPhi.pdf PDF/cosThetaPhi_1SUps.pdf
#pdflatex "\newcommand\TreeBinID{Ups2S}\input{cosThetaPhi.tex}"
#mv cosThetaPhi.pdf PDF/cosThetaPhi_2SUps.pdf
#pdflatex "\newcommand\TreeBinID{Ups3S}\input{cosThetaPhi.tex}"
#mv cosThetaPhi.pdf PDF/cosThetaPhi_3SUps.pdf
#cd ..
#fi


#rm runData
#rm runWorkspace
#rm runMassFit
#rm runLifetimeFit
#rm runPlotMassLifetime
#rm PlotFitPar
#rm runBkgHistos
#rm PlotCosThetaPhiBG
rm *.aux
rm *.log
#rm *.tex
#rm *.so
#rm *.d
#
#cd ../..
#
#mkdir -p tmp
#mv ${WorkDir}/*.cc tmp/
#mv ${WorkDir}/*.C tmp/
#mv ${WorkDir}/*.h tmp/
#mv ${WorkDir}/*.inc tmp/
#mv ${WorkDir}/Makefile tmp/

done
done

