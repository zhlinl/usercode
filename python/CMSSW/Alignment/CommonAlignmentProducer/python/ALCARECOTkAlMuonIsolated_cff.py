# AlCaReco for track based alignment using isolated muon tracks
import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
ALCARECOTkAlMuonIsolatedHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    andOr = True, ## choose logical OR between Triggerbits
    eventSetupPathsKey = 'TkAlMuonIsolated',
    throw = False # tolerate triggers stated above, but not available
    )

# DCS partitions
# "EBp","EBm","EEp","EEm","HBHEa","HBHEb","HBHEc","HF","HO","RPC"
# "DT0","DTp","DTm","CSCp","CSCm","CASTOR","TIBTID","TOB","TECp","TECm"
# "BPIX","FPIX","ESp","ESm"
import DPGAnalysis.Skims.skim_detstatus_cfi
ALCARECOTkAlMuonIsolatedDCSFilter = DPGAnalysis.Skims.skim_detstatus_cfi.dcsstatus.clone(
    DetectorType = cms.vstring('TIBTID','TOB','TECp','TECm','BPIX','FPIX'),
    ApplyFilter  = cms.bool(True),
    AndOr        = cms.bool(True),
    DebugOn      = cms.untracked.bool(False)
)

import Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi
ALCARECOTkAlMuonIsolated = Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi.AlignmentTrackSelector.clone(
    filter = True, ##do not store empty events
    applyBasicCuts = True,
    ptMin = 2.0, ##GeV 
    etaMin = -3.5,
    etaMax = 3.5,
    nHitMin = 0
    )
# These unfortunately cannot be put into the clone(..): 
ALCARECOTkAlMuonIsolated.GlobalSelector.applyIsolationtest = True
ALCARECOTkAlMuonIsolated.GlobalSelector.minJetDeltaR = 0.1
ALCARECOTkAlMuonIsolated.GlobalSelector.applyGlobalMuonFilter = True
ALCARECOTkAlMuonIsolated.TwoBodyDecaySelector.applyMassrangeFilter = False
ALCARECOTkAlMuonIsolated.TwoBodyDecaySelector.applyChargeFilter = False
ALCARECOTkAlMuonIsolated.TwoBodyDecaySelector.applyAcoplanarityFilter = False

seqALCARECOTkAlMuonIsolated = cms.Sequence(ALCARECOTkAlMuonIsolatedHLT+ALCARECOTkAlMuonIsolatedDCSFilter+ALCARECOTkAlMuonIsolated)
