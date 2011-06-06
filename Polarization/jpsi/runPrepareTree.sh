#!/bin/bash

source /afs/ihep.ac.cn/users/z/zhangll/workspace/rootset.sh 9

mapdir=/publicfs/cms/user/zhangll/data/polarization/
dirstruct=/afs/ihep.ac.cn/users/z/zhangll/workspace/work/polarization/Polarization/
project=${dirstruct}project/

#python PrepareFitTree.py --output=treeFall10_test.root	--mc  ${mapdir}TTree_pol_noTriggerFilter_MCPromptFall10_23Mar2011.root

#python PrepareFitTree.py --output=treeRun2010A.root  ${mapdir}TTree_pol_noTriggerFilter_Run2010A-Nov4ReReco_v1_06Dec2010.root

python PrepareFitTree.py --output=treeRun2010B.root  ${mapdir}TTree_pol_noTriggerFilter_Run2010B-Nov4ReReco_v1_06Dec2010.root


