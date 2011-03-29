lumiTotal: lumiMinimumBiasJun14 lumiMuJun14 lumiMu lumiMuOnia
	grep -A 5 "=== Total :" lumiMinimumBiasJun14 lumiMuJun14 lumiMu lumiMuOnia | awk 'BEGIN{tot=0}{tot+=\$6;print \$0}END{print "Total: "tot"/mub"}'
lumiMinimumBiasJun14: MinimumBiasJun14/res/lumiSummary.json
	lumiCalc.py -c frontier://LumiProd/CMS_LUMI_PROD -i MinimumBiasJun14/res/lumiSummary.json -hltpath HLT_L1DoubleMuOpen recorded > lumiMinimumBiasJun14
lumiMuJun14: MuJun14/res/lumiSummary.json
	lumiCalc.py -c frontier://LumiProd/CMS_LUMI_PROD -i MuJun14/res/lumiSummary.json -hltpath HLT_L1DoubleMuOpen recorded > lumiMuJun14
lumiMu: Mu/res/lumiSummary.json 
	lumiCalc.py -c frontier://LumiProd/CMS_LUMI_PROD -i Mu/res/lumiSummary.json -hltpath HLT_L1DoubleMuOpen recorded > lumiMu
lumiMuOnia: MuOnia/res/lumiSummary.json
	lumiCalc.py -c frontier://LumiProd/CMS_LUMI_PROD -i MuOnia/res/lumiSummary.json -hltpath HLT_L1DoubleMuOpen_Tight recorded > lumiMuOnia

data.root: MinimumBiasJun14.root MuJun14.root Mu.root MuOnia.root
	hadd -f data.root MinimumBiasJun14.root MuJun14.root Mu.root MuOnia.root
MinimumBiasJun14.root:
	hadd -f MinimumBiasJun14.root MinimumBiasJun14/res/selection*.root
MuJun14.root:
	hadd -f MuJun14.root MuJun14/res/selection*.root
Mu.root:
	hadd -f Mu.root Mu/res/selection*.root
MuOnia.root:
	hadd -f MuOnia.root MuOnia/res/selection*.root

MC: JPsi.root Upsilon1S.root Upsilon2S.root Upsilon3S.root QCD.root
JPsi.root:
	hadd -f JPsi.root JPsi/res/selection*.root
Upsilon1S.root:
	hadd -f Upsilon1S.root Upsilon1S/res/selection*.root
Upsilon2S.root:
	hadd -f Upsilon2S.root Upsilon2S/res/selection*.root
Upsilon3S.root:
	hadd -f Upsilon3S.root Upsilon3S/res/selection*.root
QCD.root:
	hadd -f QCD.root QCD/res/selection*.root

copy: data.root
	lcg-cp --verbose -b -D srmv2 "file:data.root" "srm://srm-cms.cern.ch:8443//srm/managerv2?SFN=/castor/cern.ch/user/z/zgecse/Onia2MuMu-v5/Upsilon-v1/data.root"
