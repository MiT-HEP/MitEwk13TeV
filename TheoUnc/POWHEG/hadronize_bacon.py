import FWCore.ParameterSet.Config as cms
import sys

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

process = cms.Process('MakingGenBacon')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
    )

# Input source
process.source = cms.Source("LHESource",
                            fileNames = cms.untracked.vstring(sys.argv[2])
                            )

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(False),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings'
                                    )
	)
				 )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.19 $'),
	annotation = cms.untracked.string('Configuration/Generator/python/jay_pythia8_powheg_cfi.py nevts:-1'),
	name = cms.untracked.string('Applications')
)

process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("genParticles"),
    ignoreParticleIDs = cms.vuint32(
         1000022,
         1000012, 1000014, 1000016,
         2000012, 2000014, 2000016,
         1000039, 5100039,
         4000012, 4000014, 4000016,
         9900012, 9900014, 9900016,
         39,12,14,16),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(True),
    excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
    tausAsJets = cms.bool(False)
)

process.partons = cms.EDProducer("PartonSelector",
withLeptons = cms.bool(False),
src = cms.InputTag("genParticles")
)

process.ntupler = cms.EDAnalyzer('GenNtuplerMod',
  skipOnHLTFail = cms.untracked.bool(False),
  outputName    = cms.untracked.string(sys.argv[3]),
  TriggerFile   = cms.untracked.string("tt"),
  edmPVName     = cms.untracked.string('offlinePrimaryVertices'),
  edmPFCandName = cms.untracked.string('particleFlow'),

  GenInfo = cms.untracked.PSet(
    isActive            = cms.untracked.bool(True),
    edmGenEventInfoName = cms.untracked.string('generator'),
    edmLHEEventName     = cms.untracked.string('source'),
    edmGenMETName = cms.untracked.string('genMetTrue'),
    edmGenParticlesName = cms.untracked.string('genParticles'),
    fillAllGen          = cms.untracked.bool(True)
  )
)


# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_FULL', '')

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.baconSequence = cms.Path(process.partons+process.ntupler)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RECOSIMoutput_step)

process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.baconSequence)

# filter all path with the production filter sequence                                                                                 
for path in process.paths:
        getattr(process,path)._seq = process.generator * getattr(process,path)._seq
