import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/homes/jbkim/analysis/cmssw_analyze/CMSSW_10_2_11_patch1/src/JTools/JAnalyze/root_files/TTJets_SingleLeptFromT_Tune_run1_lumi12868_event10047461.root'
    )
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
#add them to the VID producer
setupAllVIDIdsInModule(process,'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',setupVIDElectronSelection)

process.demo = cms.EDAnalyzer("MiniAnalyzer",
    electrons = cms.InputTag("slimmedElectrons"),
)

process.p = cms.Path(process.egmGsfElectronIDSequence * process.demo)
