import FWCore.ParameterSet.Config as cms

process = cms.Process("read")
process.load("CondCore.DBCommon.CondDBSetup_cfi")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000000)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200000)
)
process.source = cms.Source("EmptySource",
                            numberEventsInRun = cms.untracked.uint32(1), # do not change!
                            firstRun = cms.untracked.uint32(1)
                            )

process.GlobalPositionSource = cms.ESSource("PoolDBESSource",
    process.CondDBSetup,
    # Reading from oracle (instead of Frontier) needs the following shell variable setting (in zsh):
    # export CORAL_AUTH_PATH=/afs/cern.ch/cms/DB/conddb
    # string connect = "oracle://cms_orcoff_int2r/CMS_COND_ALIGNMENT"
    # untracked uint32 authenticationMethod = 1
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('GlobalPositionRcd'),
##         tag = cms.string('GlobalAlignment_2009_v1_express')
        tag = cms.string('GlobalAlignment_v1_offline_100430V0')
    )),
    connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT')
    # connect = cms.string('sqlite_file:output.db')
)

process.GlobalPositionRcdScan = cms.EDAnalyzer("GlobalPositionRcdScan")

process.p = cms.Path(process.GlobalPositionRcdScan)


