# Auto generated configuration file
# using: 
# Revision: 1.123 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/PYTHIA6_Upsilon3S_7TeV_cff.py -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,HLT --conditions FrontierConditions_GlobalTag,STARTUP31X_V2::All --fileout GenHLT_8E29.root --number 100 --mc --no_exec --datatier GEN-SIM-RAW --eventcontent RAWSIM --processName HLT8E29
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT8E29')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('Summer09: Pythia6 generation of Upsilon(3S), 10TeV, D6T tune'),
    name = cms.untracked.string('$Source: /cvs/CMSSW/UserCode/YZheng/UpsilonAna/test/acceptance/gen_7tev_ups3s_cfg.py,v $')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000000)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("EmptySource")

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('/tmp/yzheng/filter3s.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'STARTUP3X_V15::All'
process.mugenfilter = cms.EDFilter("MCSmartSingleParticleFilter",
    Status = cms.untracked.vint32(1, 1),
    MaxDecayRadius = cms.untracked.vdouble(1500.0, 1500.0),
    MinPt = cms.untracked.vdouble(0.0, 0.0),
    ParticleID = cms.untracked.vint32(13, -13),
    MaxEta = cms.untracked.vdouble(1000, 1000),
    MinEta = cms.untracked.vdouble(-1000, -1000),
    MaxDecayZ = cms.untracked.vdouble(3000.0, 3000.0),
    MinDecayZ = cms.untracked.vdouble(-3000.0, -3000.0)
)
process.oniafilter = cms.EDFilter("PythiaFilter",
    MaxEta = cms.untracked.double(1000.0),
    Status = cms.untracked.int32(2),
    MinEta = cms.untracked.double(-1000.0),
    MinPt = cms.untracked.double(0.0),
    ParticleID = cms.untracked.int32(100553)
)
process.mumugenfilter = cms.EDFilter("MCParticlePairFilter",
    Status = cms.untracked.vint32(1, 1),
    MinPt = cms.untracked.vdouble(0.0, 0.0),
    MaxEta = cms.untracked.vdouble(1000, 1000),
    MinEta = cms.untracked.vdouble(-1000, -1000),
    ParticleCharge = cms.untracked.int32(-1),
    ParticleID1 = cms.untracked.vint32(13),
    ParticleID2 = cms.untracked.vint32(13)
)
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(0.567),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(7000.0),
    crossSection = cms.untracked.double(2850.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTJ(11)=3     ! Choice of the fragmentation function', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(51)=10042     ! CTEQ6L1 structure function chosen', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model', 
            'MSTU(21)=1     ! Check on possible errors during program execution', 
            'PARP(82)=1.8387   ! pt cutoff for multiparton interactions', 
            'PARP(89)=1960. ! sqrts for which PARP82 is set', 
            'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter', 
            'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter', 
            'PARP(90)=0.16  ! Multiple interactions: rescaling power', 
            'PARP(67)=2.5    ! amount of initial-state radiation', 
            'PARP(85)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(86)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(62)=1.25   ! ', 
            'PARP(64)=0.2    ! ', 
            'MSTP(91)=1     !', 
            'PARP(91)=2.1   ! kt distribution', 
            'PARP(93)=15.0  ! '),
        processParameters = cms.vstring('MSEL=62          ! Quarkonia NRQCD ', 
            'PMAS(296,1)  = 10.3552    ! change Upsilon(2S) mass to Upsoilon(3S) PDG2006', 
            'KFPR(461,1)  = 100553     ! change 461 to Upsilon(3S) + g', 
            'PMAS(365,1)  = 10.3600   ! change bb~ mass larger than Upsilon(3S)', 
            'PMAS(366,1)  = 10.3600   ! change bb~ mass larger than Upsilon(3S)', 
            'PMAS(367,1)  = 10.3600   ! change bb~ mass larger than Upsilon(3S)', 
            'KFDP(4214,1) = 100553     ! bb~ -> Upsilon(3S)', 
            'KFDP(4215,1) = 100553     ! bb~ -> Upsilon(3S)', 
            'KFDP(4216,1) = 100553     ! bb~ -> Upsilon(3S)', 
            'PARP(146)=3.54   ! New values for COM matrix elements', 
            'PARP(147)=0.075  ! New values for COM matrix elements', 
            'PARP(148)=0.01   ! New values for COM matrix elements', 
            'PARP(149)=0.01   ! New values for COM matrix elements', 
            'PARP(150)=0.0    ! New values for COM matrix elements', 
            'MDME(1578,1) = 0 ! 0.014000    e-              e+', 
            'MDME(1579,1) = 1 ! 0.014000    mu-             mu+', 
            'MDME(1580,1) = 0 ! 0.014000    tau-            tau+', 
            'MDME(1581,1) = 0 ! 0.008000    d               dbar', 
            'MDME(1582,1) = 0 ! 0.024000    u               ubar', 
            'MDME(1583,1) = 0 ! 0.008000    s               sbar', 
            'MDME(1584,1) = 0 ! 0.024000    c               cbar', 
            'MDME(1585,1) = 0 ! 0.425000    g               g            g', 
            'MDME(1586,1) = 0 ! 0.020000    gamma           g            g', 
            'MDME(1587,1) = 0 ! 0.185000    Upsilon         pi+          pi-', 
            'MDME(1588,1) = 0 ! 0.088000    Upsilon         pi0          pi0', 
            'MDME(1589,1) = 0 ! 0.043000    chi_0b          gamma', 
            'MDME(1590,1) = 0 ! 0.067000    chi_1b          gamma', 
            'MDME(1591,1) = 0 ! 0.066000    chi_2b          gamma', 
            'MSTP(142)=2      ! turns on the PYEVWT Pt re-weighting routine', 
            'PARJ(13)=0.750   ! probability that a c or b meson has S=1', 
            'PARJ(14)=0.162   ! probability that a meson with S=0 is produced with L=1, J=1', 
            'PARJ(15)=0.018   ! probability that a meson with S=1 is produced with L=1, J=0', 
            'PARJ(16)=0.054   ! probability that a meson with S=1 is produced with L=1, J=1', 
            'MSTP(145)=0      ! choice of polarization', 
            'MSTP(146)=0      ! choice of polarization frame ONLY when mstp(145)=1', 
            'MSTP(147)=0      ! particular helicity or density matrix component when mstp(145)=1', 
            'MSTP(148)=1      ! possibility to allow for final-state shower evolution, extreme case !', 
            'MSTP(149)=1      ! if mstp(148)=1, it determines the kinematics of the QQ~3S1(8)->QQ~3S1(8)+g branching'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters', 
            'CSAParameters'),
        CSAParameters = cms.vstring('CSAMODE = 6     ! cross-section reweighted quarkonia')
    )
)
process.ProductionFilterSequence = cms.Sequence(process.generator*process.oniafilter*process.mumugenfilter*process.mugenfilter)
process.DT = cms.EDAnalyzer("DimuonTree",
  OutputFile = cms.string("gen_7tev_ups3s_tree.root")
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen * process.DT)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step)
# special treatment in case of production filter sequence  
for path in process.paths: 
    getattr(process,path)._seq = process.ProductionFilterSequence*getattr(process,path)._seq


