import FWCore.ParameterSet.Config as cms

process = cms.Process("SELECTION")

from zhlinl.UpsilonAna.selection_cff import *
selection(process, GlobalTag="GR_R_38X_V8::All", MC=False, TagTrigger="Mu3")

process.source.fileNames = cms.untracked.vstring(
    '/store/user/zgecse/Mu/Run2010A-PromptReco-v4-Onia2MuMu-v5/af687901d6eb0367c9b0c5fdf9b70ada/onia2MuMuPAT_185_1_Arw.root',
)   

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

