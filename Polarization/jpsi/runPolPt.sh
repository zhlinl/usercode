#!/bin/bash
source /afs/ihep.ac.cn/users/z/zhangll/workspace/rootset.sh 9

dirstruct=/afs/ihep.ac.cn/users/z/zhangll/workspace/work/polarization/Polarization/
cd ${dirstruct}

#root -l -b -q 'plotPol.cc+("Run2010AB","CS","lambda_phi")'
#root -l -b -q 'plotPol.cc+("Run2010AB","HX","lambda_phi")'
#root -l -b -q 'plotFinal.cc+("Run2010AB","CS","lambda_phi")'
#root -l -b -q 'plotFinal.cc+("Run2010AB","HX","lambda_phi")'
#
#root -l -b -q 'plotPol.cc+("Run2010AB","CS","lambda_tilde")'
#root -l -b -q 'plotPol.cc+("Run2010AB","HX","lambda_tilde")'
root -l -b -q 'plotFinal.cc+("Run2010AB","CS","lambda_tilde",1)'
#root -l -b -q 'plotFinal.cc+("Run2010AB","HX","lambda_tilde")'
#
#root -l -b -q 'plotPol.cc+("Run2010AB","CS","F_invar")'
#root -l -b -q 'plotPol.cc+("Run2010AB","HX","F_invar")'
root -l -b -q 'plotFinal.cc+("Run2010AB","CS","F_invar",2)'
#root -l -b -q 'plotFinal.cc+("Run2010AB","HX","F_invar")'

###project=project/backup/
###
###scen=CS
###dataType=Run2010A
###
###for rap_ in 1 2;do
####for pT_  in 1 3 4 8;do
###for pT_  in 1 2 3 4 5 6 7 8;do
###if [ $pT_ -eq 8 ] && [ $rap_ -eq 1 ]; then
###continue
###fi
###cp ${project}polaRun2010A_rap${rap_}_pt${pT_}-${scen}_test.root polaRun2010A_${rap_}_${pT_}-${scen}.root
###done
###done
###
####python  makePolPt.py --fitFrame=${scen} --rapBins=1 polaRun2010A_1_1-${scen}.root polaRun2010A_1_2-${scen}.root polaRun2010A_1_3-${scen}.root polaRun2010A_1_4-${scen}.root polaRun2010A_1_5-${scen}.root polaRun2010A_1_6-${scen}.root polaRun2010A_1_7-${scen}.root
###python  makePolPt.py  --fitFrame=${scen} --rapBins=1  polaRun2010A_1_2-${scen}.root polaRun2010A_1_5-${scen}.root 
###python  makePolPt.py --fitFrame=${scen} --rapBins=2 polaRun2010A_2_2-${scen}.root  polaRun2010A_2_5-${scen}.root polaRun2010A_2_6-${scen}.root polaRun2010A_2_7-${scen}.root
####python  makePolPt.py --fitFrame=${scen} --rapBins=2 polaRun2010A_2_1-${scen}.root polaRun2010A_2_2-${scen}.root polaRun2010A_2_3-${scen}.root polaRun2010A_2_4-${scen}.root polaRun2010A_2_5-${scen}.root polaRun2010A_2_6-${scen}.root polaRun2010A_2_7-${scen}.root polaRun2010A_2_8-${scen}.root 
###
###for rap_ in 1 2;do
###for pT_  in 1 2 3 4 5 6 7 8;do
###if [ $pT_ -eq 8 ] && [ $rap_ -eq 1 ]; then
###continue
###fi
###rm  polaRun2010A_${rap_}_${pT_}-${scen}.root
###done
###done
###
