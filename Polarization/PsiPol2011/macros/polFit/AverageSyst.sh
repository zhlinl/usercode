#!/bin/sh

homedir=`pwd`
cd ..
cd ..
basedir=`pwd`
cd macros/polFit
storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir

########## INPUTS ##########

for nState in 1 2 3;do



SystID=MassSigmaDep

nSystematics=3

JobID1=BiasCorrectionAug12_1S_3Sig
JobID2=BiasCorrectionAug12_2S_3Sig
JobID3=BiasCorrectionAug12_3S_3Sig
JobID4=
JobID5=
JobID6=
JobID7=
JobID8=
JobID9=


#SystID=TotalSyst
#
#nSystematics=7
#
#JobID1=BestSyst_Bkg
#JobID2=BestSyst_FrameworkI
#JobID3=BestSyst_Param
#JobID4=BestSyst_Sig_NoUnpol
#JobID5=BestSyst_TnP
#JobID6=ConstSyst
#JobID7=BestSyst_28p_SQRT12
#JobID8=
#JobID9=
 


ptBinMin=1
ptBinMax=10

########################################

cd ${homedir}

touch AverageSystematics.cc
make

mkdir Systematics
mkdir Systematics/${SystID}

SystDir=Systematics/${SystID}/AverageSyst

mkdir ${SystDir}


./AverageSystematics ${JobID1}=JobID1 ${JobID2}=JobID2 ${JobID3}=JobID3 ${JobID4}=JobID4 ${JobID5}=JobID5 ${JobID6}=JobID6 ${JobID7}=JobID7 ${JobID8}=JobID8 ${JobID9}=JobID9 ${SystID}=SystID ${storagedir}=storagedir ${basedir}=basedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${nState}nState ${nSystematics}nSystematics


done


