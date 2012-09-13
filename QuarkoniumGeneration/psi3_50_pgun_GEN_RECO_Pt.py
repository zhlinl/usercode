# Auto generated configuration file
# using: 
# Revision: 1.303.2.7 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: test_hlt -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT:quarkonium_1E33_3E33,RAW2DIGI,L1Reco,RECO --conditions START42_V14A::All --datatier GEN-SIM-RECO --eventcontent RECOSIM --no_exec --filein file:step0GEN.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_quarkonium_1E33_3E33_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

    )


process.RandomNumberGeneratorService.pregenerator = cms.PSet(
        initialSeed = cms.untracked.uint32(932751),
            engineName = cms.untracked.string('HepJamesRandom')
            )

process.pregenerator = cms.EDProducer("FlatRapidMeasuredPtGunProducer",

        PGunParameters = cms.PSet(
        # you can request more than 1 particle
        PartID = cms.vint32(100443),
            MinRapidity = cms.untracked.double(-1.3),
            MaxRapidity = cms.untracked.double(1.3),
            #
            # phi must be given in radians
            #
            ##for some reason these are tracked now...this gun doesnt use eta though
            MinEta=cms.double(0),
            MaxEta=cms.double(0),

            MinPhi = cms.double(-3.14159265358979323846),
            MaxPhi = cms.double( 3.14159265358979323846),
            MinPt  = cms.double(3),
            MaxPt  = cms.double(50),
						Resonance = cms.untracked.double(4),
         ),
        AddAntiParticle = cms.bool(False),
        Verbosity = cms.untracked.int32(0)
)

process.generator = cms.EDProducer("EvtGenProducer",
                                       src = cms.string('pregenerator'),
                                       use_default_decay = cms.untracked.bool(False),
                                       decay_table = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/DECAY.DEC'),
                                       particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt.pdl'),
                                       user_decay_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/mydecays.dec'),
                                       list_forced_decays = cms.vstring('Mypsi2S'),
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
        outputCommands = cms.untracked.vstring('keep *_pregenerator_*_*')
        )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.303.2.7 $'),
    annotation = cms.untracked.string('test_hlt nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('psi_pgun_RECO.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START42_V14A::All'

process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1MuGMTParametersRcd' ),
        tag     = cms.string( 'L1MuGMTParameters_synctf_10_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
    ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1MuDTTFParametersRcd' ),
        tag     = cms.string( 'L1MuDTTFParameters_dttf11_TSC_09_17_col_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1MuCSCTFConfigurationRcd' ),
        tag     = cms.string( 'L1MuCSCTFConfiguration_90511_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCBxOrConfigRcd' ),
        tag     = cms.string( 'L1RPCBxOrConfig_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCConeDefinitionRcd' ),
        tag     = cms.string( 'L1RPCConeDefinition_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCConfigRcd' ),
        tag     = cms.string( 'L1RPCConfig_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCHsbConfigRcd' ),
        tag     = cms.string( 'L1RPCHsbConfig_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        )
    )

# Path and EndPath definitions
process.generation_step = cms.Path(process.pregenerator*process.generator*process.pgen)
#process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step])
