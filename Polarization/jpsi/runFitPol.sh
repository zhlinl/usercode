#!/bin/bash

source /afs/ihep.ac.cn/users/z/zhangll/workspace/rootset.sh 9

scen=CS
#sub=lambda_tilde
#sub=lambda_phi
sub=F_invar
dirstruct=/afs/ihep.ac.cn/users/z/zhangll/workspace/work/polarization/Polarization/
cd ${dirstruct}
dataType=Run2010AB
project=${dirstruct}project/${dataType}/${scen}/
mkdir -p ${project}${sub}/

mapdir=/publicfs/cms/user/zhangll/data/polarization/
#SETTINGS:::MAPS
geomAccPR=${mapdir}geomAccHistos_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded_zeroBinsCorrected.root
recoEffPR=${mapdir}recoEffHistos_ATLASPT_20March2011_phiFolded_zeroBinsCorrected.root
trigEffPR=${mapdir}trigEffHistos_ATLASPT_DoubleMu0_20March2011_phiFolded_zeroBinsCorrected.root
geomAccNP=${mapdir}geomAccHistos_NP_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded_zeroBinsCorrected.root
recoEffNP=${mapdir}recoEffHistos_NP_ATLASPT_19March2011_phiFolded_zeroBinsCorrected.root
trigEffNP=${mapdir}trigEffHistos_NP_ATLASPT_DoubleMu0_19March2011_phiFolded_zeroBinsCorrected.root

for rap_ in 1 2;do
for pT_  in 1 2 3 4 5 6 7 8;do
#for rap_ in 1;do
#for pT_  in 6;do
if [ $pT_ -eq 8 ] && [ $rap_ -eq 1 ]; then
continue
fi

python polarizationFit.py 	--workspaceName=pola${dataType}_rap${rap_}_pt${pT_} ${dirstruct}tree${dataType}.root --treeName=data --fitFrame=${scen} --testBin=${rap_},${pT_} --acceptanceMap=${geomAccPR},${geomAccNP} --efficiencyMap=${recoEffPR},${recoEffNP},${trigEffPR},${trigEffNP} #>& ${project}log_Fit_${dataType}_rap${rap_}_pt${pT_} &

#cp ${project}pola${dataType}_rap${rap_}_pt${pT_}.root ${dirstruct}pola${dataType}_rap${rap_}_pt${pT_}.root
python polarizationFitSimple.py 	--workspaceName=pola${dataType}_rap${rap_}_pt${pT_} --fitFrame=${scen} --testBin=${rap_},${pT_} --acceptanceMap=${geomAccPR},${geomAccNP} --recoEfficiencyMap=${recoEffPR},${recoEffNP} --trigEfficiencyMap=${trigEffPR},${trigEffNP} --lambdaPhiSub=${sub}
#--noNonPrompt
### old verison fit ###
#python polarizationFitSimple.py 	--workspaceName=pola${dataType}_rap${rap_}_pt${pT_} --fitFrame=${scen} --testBin=${rap_},${pT_} --acceptanceMap=${geomAccPR} --efficiencyMap=${recoEffPR},${recoEffNP},${trigEffPR},${trigEffNP} --lambdaPhiSub=lambda_tilde #--noNonPrompt

cp ${dirstruct}pola${dataType}_rap${rap_}_pt${pT_}-${scen}.root ${dirstruct}pola${dataType}_rap${rap_}_pt${pT_}-${scen}_test.root
python makePolPlots.py pola${dataType}_rap${rap_}_pt${pT_}-${scen}_test.root --plotPol --testBin=${rap_},${pT_} 
#--pedagogical

mv ${dirstruct}pola${dataType}_rap${rap_}_pt${pT_}.root ${project}
mv ${dirstruct}pola${dataType}_rap${rap_}_pt${pT_}-${scen}.root ${project}${sub}/
mv ${dirstruct}pola${dataType}_rap${rap_}_pt${pT_}-${scen}_test.root ${project}${sub}/
#rm pola${dataType}_rap${rap_}_pt${pT_}.root

done
done

#./extractPromptLifetimeShape /home/zhlinl/cluster1/data/polarization/Spring2010/TTree_pol_Mu0Track0Jpsi_MCprompt.root /home/zhlinl/cluster1/data/polarization/Spring2010/TTree_pol_Mu0Track0Jpsi_B0ToJPsiMuMu.root /home/zhlinl/cluster1/data/polarization/Spring2010/TTree_pol_Mu0Track0Jpsi_BpToJPsiMuMu.root /home/zhlinl/cluster1/data/polarization/Spring2010/TTree_pol_Mu0Track0Jpsi_BsToJPsiMuMu.root /home/zhlinl/cluster1/data/polarization/Spring2010/TTree_pol_Mu0Track0Jpsi_LambdaBToJPsiMuMu.root 
