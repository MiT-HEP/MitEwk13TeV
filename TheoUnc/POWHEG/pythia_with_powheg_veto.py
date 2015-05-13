import FWCore.ParameterSet.Config as cms
import sys

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
#from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *

IN_FILE=sys.argv[2]
OUT_FILE=sys.argv[3]
NEVTS=sys.argv[4]
NSKIP=sys.argv[5]

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
    input = cms.untracked.int32(int(NEVTS))
    )

# Input source
process.source = cms.Source("LHESource",
                            fileNames = cms.untracked.vstring(IN_FILE),
			    skipEvents = cms.untracked.unit32(int(NSKIP))
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
    emissionVeto1 = cms.untracked.PSet(),
    EV1_vetoOn = cms.bool(True),
    EV1_maxVetoCount = cms.int32(0),
    EV1_pThardMode = cms.int32(0),
    EV1_pTemtMode = cms.int32(0),
    EV1_emittedMode = cms.int32(0),
    EV1_pTdefMode = cms.int32(1),
    EV1_MPIvetoOn = cms.bool(False),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        processParameters = cms.vstring(
			'SpaceShower:pTmaxMatch = 2',
			'TimeShower:pTmaxMatch = 2',
                        'TimeShower:QEDshowerByL  = on',
                        'TimeShower:QEDshowerByQ  = on',
                        'SpaceShower:QEDshowerByQ  = on',
                        'TimeShower:alphaEMorder = 0',
                        'SpaceShower:alphaEMorder = 0',
                        'TimeShower:pTminChgL=1.e-6',
                        'TimeShower:pTminChgQ=0.8944e0',
                        'SpaceShower:pTminChgQ=0.8944e0',
			),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
                                    'processParameters'
                                    )
	)
				 )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.19 $'),
	annotation = cms.untracked.string('Configuration/Generator/python/jay_pythia8_powheg_cfi.py nevts:-1'),
	name = cms.untracked.string('Applications')
)

# Output definition                                                                                      

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string(OUT_FILE),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_FULL', '')

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RECOSIMoutput_step)

# filter all path with the production filter sequence                                                                                 
for path in process.paths:
        getattr(process,path)._seq = process.generator * getattr(process,path)._seq
