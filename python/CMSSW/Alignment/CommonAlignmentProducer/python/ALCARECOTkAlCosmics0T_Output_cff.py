# last update: $Date: 2010/03/17 18:17:17 $ by $Author: mussgill $

import FWCore.ParameterSet.Config as cms

# AlCaReco for track based alignment using Cosmic muon events
OutALCARECOTkAlCosmics0T_noDrop = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlCosmicsCTF0T', 
            'pathALCARECOTkAlCosmicsCosmicTF0T')
    ),
    outputCommands = cms.untracked.vstring(
        'keep *_ALCARECOTkAlCosmics*0T_*_*',
        'keep siStripDigis_DetIdCollection_*_*',
        'keep L1AcceptBunchCrossings_*_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep *_TriggerResults_*_*',
        'keep DcsStatuss_scalersRawToDigi_*_*',
        'keep Si*Cluster*_si*Clusters_*_*', # for cosmics keep original clusters
        'keep recoMuons_muons1Leg_*_*', # save muons as timing info is needed for BP corrections in deconvolution
        'keep *_MEtoEDMConverter_*_*')
)

import copy
OutALCARECOTkAlCosmics0T = copy.deepcopy(OutALCARECOTkAlCosmics0T_noDrop)
OutALCARECOTkAlCosmics0T.outputCommands.insert(0, "drop *")
