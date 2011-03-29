# Auto generated configuration file
# using:
# Revision: 1.123
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v
# with command line options: Configuration/GenProduction/python/PYTHIA6_Upsilon1SWithFSR_7TeV_cff.py -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,HLT --conditions FrontierConditions_GlobalTag,STARTUP31X_V2::All --fileout GenHLT_8E29.root --number 100 --mc --no_exec --datatier GEN-SIM-RAW --eventcontent RAWSIM --processName HLT8E29 --customise=Configuration/GenProduction/HERWIGPP_custom.py
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

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
    fileName = cms.untracked.string('/tmp/yzheng/filter1s.root'),
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
#process.GlobalTag.globaltag = 'STARTUP31X_V2::All'
process.GlobalTag.globaltag = "STARTUP3X_V15::All"
process.mugenfilter = cms.EDFilter("MCSmartSingleParticleFilter",
    Status = cms.untracked.vint32(1, 1),
    MaxDecayRadius = cms.untracked.vdouble(1500.0, 1500.0),
    MinPt = cms.untracked.vdouble(0.0,0.0),
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
    ParticleID = cms.untracked.int32(553)
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
    ExternalDecays = cms.PSet(
        EvtGen = cms.untracked.PSet(
            use_default_decay = cms.untracked.bool(False),
            decay_table = cms.FileInPath('GeneratorInterface/ExternalDecays/data/DECAY_NOLONGLIFE.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/evt.pdl'),
            user_decay_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/Onia_mumu.dec'),
            list_forced_decays = cms.vstring('MyUpsilon'),
            operates_on_particles = cms.vint32(0)
        ),
        parameterSets = cms.vstring('EvtGen')
    ),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(0.139),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(7000.0),
    crossSection = cms.untracked.double(99900.0),
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
            'MDME(1034,1)=0   ! 0.025200    e- e+',
            'MDME(1035,1)=1   ! 0.024800    mu- mu+',
            'MDME(1036,1)=0   ! 0.026700    tau- tau+',
            'MDME(1037,1)=0   ! 0.015000    d dbar',
            'MDME(1038,1)=0   ! 0.045000    u ubar',
            'MDME(1039,1)=0   ! 0.015000    s sbar',
            'MDME(1040,1)=0   ! 0.045000    c cbar',
            'MDME(1041,1)=0   ! 0.774300    g g g',
            'MDME(1042,1)=0   ! 0.029000    gamma g',
            'MSTP(142)=2      ! turns on the PYEVWT Pt re-weighting routine',
            'PARJ(13)=0.750   ! probability that a c or b meson has S=1',
            'PARJ(14)=0.162   ! probability that a meson with S=0 is produced with L=1, J=1',
            'PARJ(15)=0.018   ! probability that a meson with S=1 is produced with L=1, J=0',
            'PARJ(16)=0.054   ! probability that a meson with S=1 is produced with L=1, J=1',
            'MSTP(145)=0      ! choice of polarization',
            'MSTP(146)=0      ! choice of polarization frame ONLY when mstp(145)=1',
            'MSTP(147)=0      ! particular helicity or density matrix component when mstp(145)=1',
            'MSTP(148)=1      ! possibility to allow for final-state shower evolution, extreme case !',
            'MSTP(149)=1      ! if mstp(148)=1, it determines the kinematics of the QQ~3S1(8)->QQ~3S1(8)+g branching',
            'PARP(141)=1.16   ! New values for COM matrix elements',
            'PARP(142)=0.0119 ! New values for COM matrix elements',
            'PARP(143)=0.01   ! New values for COM matrix elements',
            'PARP(144)=0.01   ! New values for COM matrix elements',
            'PARP(145)=0.05   ! New values for COM matrix elements',
            'PARP(146)=9.28   ! New values for COM matrix elements',
            'PARP(147)=0.15   ! New values for COM matrix elements',
            'PARP(148)=0.02   ! New values for COM matrix elements',
            'PARP(149)=0.02   ! New values for COM matrix elements',
            'PARP(150)=0.085  ! New values for COM matrix elements'),
        parameterSets = cms.vstring('pythiaUESettings',
            'processParameters',
            'CSAParameters'),
        CSAParameters = cms.vstring('CSAMODE = 6     ! cross-section reweighted quarkonia')
    )
)
process.ProductionFilterSequence = cms.Sequence(process.generator*process.oniafilter*process.mumugenfilter*process.mugenfilter)

process.DT = cms.EDAnalyzer("DimuonTree",
  OutputFile = cms.string("gen_7tev_ups1s_tree.root")
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen * process.DT)
#process.endjob_step = cms.Path(process.endOfProcess)
#process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step)
#process.schedule.extend([process.out_step])
#process.schedule.extend([process.endjob_step,process.out_step])
# special treatment in case of production filter sequence
for path in process.paths:
    getattr(process,path)._seq = process.ProductionFilterSequence*getattr(process,path)._seq


# Automatic addition of the customisation function

def customise(process):
        process.genParticles.abortOnUnknownPDGCode = False

        return(process)


# End of customisation function definition

process = customise(process)


