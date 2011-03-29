import FWCore.ParameterSet.Config as cms

process = cms.Process('tp')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/EventContent/EventContent_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/mc/Summer09/Upsilon1S/GEN-SIM-RECO/MC_31X_V3_7TeV-v1/0016/F8C2372E-ABA3-DE11-BB4A-001D0967DBE8.root'
))

process.DT = cms.EDAnalyzer("DimuonTree",
  OutputFile = cms.string("tree_tracking.root")
)

process.p = cms.Path(process.DT)

# Schedule definition
process.schedule = cms.Schedule(process.p)

