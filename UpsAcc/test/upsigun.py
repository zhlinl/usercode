# Force Upsilon->mumu decay 

import FWCore.ParameterSet.Config as cms

process = cms.Process("UPSGUN")

process.load("Configuration.Generator.PythiaUESettings_cfi")
process.load("Configuration.StandardSequences.Generator_cff")
process.load("Configuration.StandardSequences.VtxSmearedEarly7TeVCollision_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("PhysicsTools.HepMCCandAlgos.genParticleCandidatesFast_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.RandomNumberGeneratorService.evtgenproducer = cms.PSet(
    initialSeed = cms.untracked.uint32(93278151),
    engineName = cms.untracked.string('HepJamesRandom')
    )

process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.EventContent.EventContent_cff')

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    )
)

process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
     ignoreTotal=cms.untracked.int32(0),
     oncePerEventMode = cms.untracked.bool(False)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_38Y_V14::All'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)




process.source = cms.Source("EmptySource")

process.generator = cms.EDProducer("FlatRapidFlatPtGunProducer",
                              
                                   PGunParameters = cms.PSet(
    # you can request more than 1 particle
    PartID = cms.vint32(553), 
    MinRapidity = cms.untracked.double(-2.0),
    MaxRapidity = cms.untracked.double(2.0),
    # for some reason these are tracked now...this gun doesnt use eta though
    MinEta=cms.double(0),
    MaxEta=cms.double(0),
    # phi must be given in radians
    MinPhi = cms.double(-3.14159265358979323846),
    MaxPhi = cms.double( 3.14159265358979323846),
    MinPt  = cms.double(0),
    MaxPt  = cms.double(20),
    ),
                                   AddAntiParticle = cms.bool(False),
                                   Verbosity = cms.untracked.int32(0)
                                   )




process.evtgenproducer = cms.EDProducer("EvtGenProducer",
    src = cms.string('generator'),
    use_default_decay = cms.untracked.bool(False),
    decay_table = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/DECAY.DEC'),
    particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt.pdl'),                                 
    user_decay_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/mydecays.dec'),
    list_forced_decays = cms.vstring('MyUpsilon','MyUpsilon2','MyUpsilon3'),
    processParameters = cms.vstring(
        'MDME(858,1)=0                 ! Upsilon -> ee turned OFF', 
        'MDME(859,1)=0                 ! Upsilon -> mumu turned OFF', 
        'MDME(860,1)=0                 ! Upsilon -> random turned OFF',

        'MDCY(134,1) = 0', 
        'MDCY(137,1) = 0', 
        'MDCY(138,1) = 0', 
        'MDCY(135,1) = 0', 
        'MDCY(141,1) = 0', 
        'MDCY(140,1) = 0', 
        'MDCY(15,1) = 0', 
        'MDCY(123,1) = 0', 
        'MDCY(126,1) = 0', 
        'MDCY(129,1) = 0', 
        'MDCY(122,1) = 0', 
        'MDCY(125,1) = 0', 
        'MDCY(128,1) = 0', 
        'MDCY(262,1) = 0', 
        'MDCY(264,1) = 0', 
        'MDCY(263,1) = 0', 
        'MDCY(265,1) = 0', 
        'MDCY(286,1) = 0', 
        'MDCY(287,1) = 0', 
        'MDCY(124,1) = 0', 
        'MDCY(127,1) = 0', 
        'MDCY(266,1) = 0', 
        'MDCY(288,1) = 0', 
        'MDCY(267,1) = 0', 
        'MDCY(130,1) = 0', 
        'MDCY(112,1) = 0', 
        'MDCY(113,1) = 0', 
        'MDCY(114,1) = 0', 
        'MDCY(117,1) = 0', 
        'MDCY(258,1) = 0', 
        'MDCY(256,1) = 0', 
        'MDCY(257,1) = 0', 
        'MDCY(259,1) = 0', 
        'MDCY(284,1) = 0', 
        'MDCY(283,1) = 0', 
        'MDCY(118,1) = 0', 
        'MDCY(115,1) = 0', 
        'MDCY(102,1) = 0', 
        'MDCY(109,1) = 0', 
        'MDCY(103,1) = 0', 
        'MDCY(107,1) = 0', 
        'MDCY(110,1) = 0', 
        'MDCY(119,1) = 0', 
        'MDCY(120,1) = 0', 
        'MDCY(281,1) = 0', 
        'MDCY(280,1) = 0', 
        'MDCY(281,1) = 0', 
        'MDCY(108,1) = 0', 
        'MDCY(104,1) = 0', 
        'MDCY(253,1) = 0', 
        'MDCY(251,1) = 0', 
        'MDCY(250,1) = 0', 
        'MDCY(252,1) = 0', 
        'MDCY(254,1) = 0', 
        'MDCY(282,1) = 0', 
        'MDCY(285,1) = 0', 
        'MDCY(111,1) = 0', 
        'MDCY(121,1) = 0', 
        'MDCY(255,1) = 0', 
        'MDCY(261,1) = 0', 
        'MDCY(131,1) = 0', 
        'MDCY(132,1) = 0', 
        'MDCY(295,1) = 0', 
        'MDCY(268,1) = 0', 
        'MDCY(289,1) = 0', 
        'MDCY(133,1) = 0', 
        'MDCY(146,1) = 0', 
        'MDCY(147,1) = 0', 
        'MDCY(296,1) = 0', 
        'MDCY(278,1) = 0', 
        'MDCY(294,1) = 0', 
        'MDCY(148,1) = 0', 
        'MDCY(279,1) = 0', 
        'MDCY(181,1) = 0', 
        'MDCY(182,1) = 0', 
        'MDCY(84,1) = 0', 
        'MDCY(179,1) = 0', 
        'MDCY(185,1) = 0', 
        'MDCY(189,1) = 0', 
        'MDCY(187,1) = 0', 
        'MDCY(194,1) = 0', 
        'MDCY(192,1) = 0', 
        'MDCY(164,1) = 0', 
        'MDCY(169,1) = 0', 
        'MDCY(158,1) = 0', 
        'MDCY(159,1) = 0', 
        'MDCY(175,1) = 0', 
        'MDCY(155,1) = 0', 
        'MDCY(151,1) = 0', 
        'MDCY(162,1) = 0', 
        'MDCY(167,1) = 0', 
        'MDCY(163,1) = 0', 
        'MDCY(170,1) = 0', 
        'MDCY(168,1) = 0', 
        'MDCY(174,1) = 0', 
        'MDCY(172,1) = 0', 
        'MDCY(173,1) = 0', 
        'MDCY(176,1) = 0', 
        'MDCY(180,1) = 0', 
        'MDCY(186,1) = 0', 
        'MDCY(188,1) = 0', 
        'MDCY(193,1) = 0', 
        'MDCY(195,1) = 0', 
        'MDCY(196,1) = 0', 
        'MDCY(197,1) = 0', 
        'MDCY(43,1) = 0', 
        'MDCY(44,1) = 0', 
        'MDCY(269,1) = 0', 
        'MDCY(210,1) = 0', 
        'MDCY(211,1) = 0', 
        'MDCY(219,1) = 0', 
        'MDCY(227,1) = 0', 
        'MDCY(217,1) = 0', 
        'MDCY(208,1) = 0', 
        'MDCY(215,1) = 0', 
        'MDCY(143,1) = 0', 
        'MDCY(223,1) = 0', 
        'MDCY(225,1) = 0', 
        'MDCY(272,1) = 0', 
        'MDCY(291,1) = 0', 
        'MDCY(273,1) = 0', 
        'MDCY(139,1) = 0', 
        'MDCY(270,1) = 0', 
        'MDCY(290,1) = 0', 
        'MDCY(271,1) = 0', 
        'MDCY(136,1) = 0', 
        'MDCY(274,1) = 0', 
        'MDCY(292,1) = 0', 
        'MDCY(275,1) = 0', 
        'MDCY(142,1) = 0', 
        'MDCY(144,1) = 0', 
        'MDCY(145,1) = 0', 
        'MDCY(209,1) = 0', 
        'MDCY(218,1) = 0', 
        'MDCY(216,1) = 0', 
        'MDCY(224,1) = 0', 
        'MDCY(226,1) = 0', 
        'MDCY(228,1) = 0', 
        'MDCY(276,1) = 0', 
        'MDCY(277,1) = 0', 
        'MDCY(293,1) = 0', 
        'MDCY(105,1) = 0' )
)

EvtGenEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_evtgenproducer_*_*')
)

process.load("Configuration.EventContent.EventContent_cff")


process.RECOSIMEventContent.outputCommands.extend(EvtGenEventContent.outputCommands)
process.RECOSIM = cms.OutputModule("PoolOutputModule",
                                process.RECOSIMEventContent,
                                fileName = cms.untracked.string('ups1s_gun.root')
                                )

process.evtgen = cms.Path(process.evtgenproducer)



process.options = cms.untracked.PSet(
   Rethrow = cms.untracked.vstring('ProductNotFound')
)


process.pstart = cms.Path(process.generator)
process.p0 = cms.Path(process.pgen)
process.p1 = cms.Path(process.psim)
process.p2 = cms.Path(process.pdigi)
process.p3 = cms.Path(process.L1Emulator)
process.p4 = cms.Path(process.DigiToRaw)
process.p5 = cms.Path(process.RawToDigi)
process.p6 = cms.Path(process.reconstruction)

process.outpath = cms.EndPath(process.RECOSIM)

process.schedule = cms.Schedule(process.pstart,
                                process.evtgen,
                                process.p0,
                                process.p1,
                                process.p2,
                                process.p3,
                                process.p4,
                                process.p5,
                                process.p6,

process.outpath
)

process.g4SimHits.Generator.HepMCProductLabel = 'evtgenproducer'
process.genParticleCandidates.src = 'evtgenproducer'
process.genParticles.src = 'evtgenproducer'
process.VtxSmeared.src = 'evtgenproducer' 
process.genParticles.abortOnUnknownPDGCode = False
