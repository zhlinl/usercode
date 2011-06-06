#!/bin/bash

source /afs/ihep.ac.cn/users/z/zhangll/workspace/rootset.sh 9

dirstruct=/afs/ihep.ac.cn/users/z/zhangll/workspace/work/polarization/Polarization/
cd ${dirstruct}

dataType=Run2010AB

mkdir -p pic/${dataType}/log

for rap_ in 1 2;do
for pT_ in 1 2 3 4 5 6 7 8;do
#for pT_ in 4;do

if [ $pT_ -eq 8 ] && [ $rap_ -eq 1 ]; then
continue
fi

root -l -q -b plotLM.cc\(${rap_},${pT_}\) >& pic/${dataType}/log/log_rap${rap_}_pt${pT_} &
#root -l -q -b plotLM.cc+\(${rap_},${pT_}\) 

done
done
