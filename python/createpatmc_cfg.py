#Does primary event skimming and PFBRECO
#Author: Joosep Pata joosep.pata@cern.ch

#from Configuration.StandardSequences.Geometry_cff import *
from Configuration.Geometry.GeometryIdeal_cff import *
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
import FWCore.ParameterSet.Config as cms

## import skeleton process
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

from UserCode.TTHPAT.eventCounting import *

from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *

from FWCore.ParameterSet.VarParsing import VarParsing
import sys
print "sys.argv=",sys.argv

options = VarParsing("analysis")
options.register ("isMC", True,
  VarParsing.multiplicity.singleton,
  VarParsing.varType.bool,
  "Run on MC"
)
options.register ("doDebug", False,
  VarParsing.multiplicity.singleton,
  VarParsing.varType.bool,
  "Run in debugging mode"
)
options.register ("doSkimming", True,
  VarParsing.multiplicity.singleton,
  VarParsing.varType.bool,
  "Preselect events"
)
options.register ("doSlimming", True,
  VarParsing.multiplicity.singleton,
  VarParsing.varType.bool,
  "Drop unnecessary collections"
)


#Tag from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions?redirectedfrom=CMS.SWGuideFrontierConditions#2012_MC_production
# Latest for "53Y Releases (MC)"
options.register (
    "globalTag",
    "GR_R_52_V7::All",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global tag"
)

#needs to be disabled for crab to work, otherwise get configuration errors
#validate using edmConfigHash
if hasattr(sys, "argv") and not "-cfg" in sys.argv:
    options.parseArguments()

process = cms.Process("TTH")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

if len(options.inputFiles) == 0:
    #options.inputFiles = ["/store/mc/Summer12_DR53X/TTH_HToBB_M-125_8TeV-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/00FA9388-81FC-E111-A80D-00215E2217BE.root"]
    options.inputFiles = ["/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00037C53-AAD1-E111-B1BE-003048D45F38.root"]

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),

    #required for taus: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#53X_recommendation_for_Run_I_ana
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop recoPFTaus_*_*_*'
    )
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName=cms.untracked.string(options.outputFile),
    SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring()),
    outputCommands=cms.untracked.vstring(["drop *"])
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

if options.doDebug:
    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger = cms.Service("MessageLogger",
        destinations=cms.untracked.vstring("cout", "debug"),
        debugModules=cms.untracked.vstring("*"),
        cout=cms.untracked.PSet(threshold=cms.untracked.string("INFO")),
        debug=cms.untracked.PSet(threshold=cms.untracked.string("DEBUG")),
    )
    process.MessageLogger.cerr.FwkReport.reportEvery = 1
else:
    process.load("FWCore.MessageService.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = 10

postfix = ""
jetCorr = ["L1FastJet", "L2Relative", "L3Absolute"]
if not options.isMC:
    jetCorr += ["L2L3Residual"]

usePF2PAT(process, runPF2PAT=True, jetAlgo="AK5", runOnMC=options.isMC, postfix=postfix,
    jetCorrections=("AK5PFchs", jetCorr),
    pvCollection=cms.InputTag("goodOfflinePrimaryVertices"),
    #typeIMetCorrections = True
    typeIMetCorrections = False #Type1 MET now applied later using runMETUncertainties
)

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorPFnoPU2012
process.pfPileUp.Enable = True
process.pfPileUp.checkClosestZVertex = False

#-------------------------------------------------
# selection step 2: vertex filter
#-------------------------------------------------

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorPFnoPU2012
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel#Cleaning_Filters
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = cms.PSet(
        minNdof = cms.double(4.0),
        maxZ = cms.double(24.0),
        maxRho = cms.double(2.0)
    ),
    filter = cms.bool(True),
    src = cms.InputTag("offlinePrimaryVertices")
)

#-------------------------------------------------
# OBJECT RECO
#-------------------------------------------------
# Muons
#-------------------------------------------------

process.selectedPatMuons.cut = "pt>10 && abs(eta)<3.0"
process.pfIsolatedMuons.doDeltaBetaCorrection = True
process.pfIsolatedMuons.isolationCut = 0.2
process.patMuons.pfMuonSource = cms.InputTag("pfIsolatedMuons")
process.muonMatch.src = cms.InputTag("pfIsolatedMuons")

process.muonMatchAll = process.muonMatch.clone(
  src = cms.InputTag("pfMuons")
)
process.patMuonsAll = process.patMuons.clone(
  pfMuonSource = cms.InputTag("pfMuons"),
  genParticleMatch = cms.InputTag("muonMatchAll"),
)
process.selectedPatMuonsAll = process.selectedPatMuons.clone(
  src = cms.InputTag("patMuonsAll"),
)
process.muonsWithID = cms.EDProducer(
  "MuonIDProducer",
  muonSrc = cms.InputTag("selectedPatMuons"),
  primaryVertexSource = cms.InputTag("goodOfflinePrimaryVertices")
)
process.muonsWithIDAll = process.muonsWithID.clone(
  muonSrc = cms.InputTag("selectedPatMuonsAll")
)

process.muonSequence = cms.Sequence()
process.muonSequence += process.muonMatchAll
process.muonSequence += (
  process.patMuonsAll *
  process.selectedPatMuonsAll *
  process.muonsWithID *
  process.muonsWithIDAll
)

#-------------------------------------------------
# Electrons
# Implemented as in https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=208765
#-------------------------------------------------

#if not maxLeptonIso is None:
#    process.pfIsolatedElectrons.isolationCut = maxLeptonIso
#Use both isolated and un-isolated electrons as patElectrons.
#NB: no need to change process.electronMatch.src to pfElectrons,
#    it"s already gsfElectrons, which is a superset of the pfElectrons

#From EgammaAnalysis/ElectronTools/test/patTuple_electronId_cfg.py
process.load("EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi")
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaTrigNoIPV0 + process.mvaNonTrigV0 )
process.patElectrons.electronIDSources = cms.PSet(
  mvaTrigV0 = cms.InputTag("mvaTrigV0"),
  mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"),
  mvaTrigNoIPV0 = cms.InputTag("mvaTrigNoIPV0"),
)
process.patPF2PATSequence.replace(process.patElectrons, process.mvaID * process.patElectrons)
process.selectedPatElectrons.cut = "pt>20 && abs(eta)<3.0"
process.pfIsolatedElectrons.isolationCut = 0.2

process.pfElectrons.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFId"))
process.pfElectrons.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFId")
process.pfElectrons.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"), cms.InputTag("elPFIsoValueGamma03PFId"))
process.pfElectrons.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFId")
process.pfElectrons.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"), cms.InputTag("elPFIsoValueGamma03PFId"))

process.patElectrons.isolationValues.pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFId")
process.patElectrons.isolationValues.pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFId")
process.patElectrons.isolationValues.pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFId")
process.patElectrons.isolationValues.pfPhotons = cms.InputTag("elPFIsoValueGamma03PFId")
process.patElectrons.isolationValues.pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFId")

process.pfIsolatedElectrons.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFId"))
process.pfIsolatedElectrons.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFId")
process.pfIsolatedElectrons.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"), cms.InputTag("elPFIsoValueGamma03PFId"))

process.patElectronsAll = process.patElectrons.clone(
  pfElectronSource=cms.InputTag("pfElectrons")
)
process.selectedPatElectronsAll = process.selectedPatElectrons.clone(
  src=cms.InputTag("patElectronsAll")
)
process.electronsWithID = cms.EDProducer(
    "ElectronIDProducer",
    electronSrc = cms.InputTag("selectedPatElectrons"),
    primaryVertexSource = cms.InputTag("goodOfflinePrimaryVertices")
)
process.electronsWithIDAll = process.electronsWithID.clone(
    electronSrc = cms.InputTag("selectedPatElectronsAll")
)

process.electronSequence = cms.Sequence(
  process.patElectronsAll *
  process.selectedPatElectronsAll *
  process.electronsWithID *
  process.electronsWithIDAll
)

#---------------------------------------------
# Taus https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID
# Using the recipe Tau ID 2014 (preparation for Run II) -> 53X (recommendation for Run I analyses)
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#53X_recommendation_for_Run_I_ana
#---------------------------------------------
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

process.tauSequence = cms.Sequence(
process.recoTauClassicHPSSequence
)

#---------------------------------------------
# Trigger matching
#---------------------------------------------

process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi")

process.patTriggerSequence = cms.Sequence(
  process.patTrigger *
  process.patTriggerEvent
)

#-------------------------------------------------
# Jets
# MET corrections as https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis#Type_I_0_with_PAT
#-------------------------------------------------

##Taus are NOT removed from jets (==False)
#single-top specific
#process.pfNoTau.enable = False

process.selectedPatJets.cut = cms.string("pt>30 && abs(eta)<5.0")

#enabled by git submodule add https://github.com/latinos/UserCode-CMG-CMGTools-External CMSSW/src/CMGTools/External
process.load("CMGTools.External.pujetidsequence_cff")
process.patPF2PATSequence += process.puJetIdSqeuence

#from https://github.com/vhbb/vhbb/blob/master/HbbAnalyzer/test/patMC.py
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, "TC")
# we have already PF Met as default in PF2PAT (AR)

#A few generic options
process.patJets.addTagInfos  = True

# rho2.5 calculation
process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJetsForIsolation = process.kt4PFJets.clone(
    rParam=0.6,
    doRhoFastjet=True
)
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

# This is only needed if using the obsolete jetTools.py from VHbbAnalysis/additionalFiles
if not hasattr(process,"kt6PFJets"):
    setattr(process,"kt6PFJets", process.kt4PFJets.clone(doAreaFastjet=True, doRhoFastjet=True, rParam=0.6))

process.kt6PFJets25 = process.kt4PFJets.clone( src = "pfNoElectron"+postfix,rParam = 0.6,doRhoFastjet = True,Ghost_EtaMax = 2.5, Rho_EtaMax = 2.5 )
process.kt6PFJetsCentralNeutral = process.kt6PFJets.clone( src = cms.InputTag("pfAllNeutralHadronsAndPhotons"+postfix), Ghost_EtaMax = cms.double(3.1), Rho_EtaMax = cms.double(2.5), inputEtMin = cms.double(0.5) )

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
process.ak7PFJets = ak5PFJets.clone( rParam = cms.double(0.7) )
addJetCollection(
    process, cms.InputTag("ak7PFJets"), "AK7", "PF",
    doJTA=True, doBTagging=True, jetCorrLabel=("AK7PF", jetCorr),
    doType1MET=False, doL1Cleaning = False, doL1Counters=False,
    doJetID = False
)

from PhysicsTools.PatAlgos.tools.jetTools import *

#process.load("RecoJets.JetProducers.caSubjetFilterPFJets_cfi")
from RecoJets.JetProducers.caSubjetFilterPFJets_cfi import caSubjetFilterPFJets
process.caVHPFJets = caSubjetFilterPFJets.clone(src=cms.InputTag("pfNoElectron"+postfix),useAdjacency = cms.int32(0))

#process.load("RecoJets.JetProducers.caSubjetFilterGenJets_cfi")
from RecoJets.JetProducers.caSubjetFilterGenJets_cfi import caSubjetFilterGenJets
process.caVHGenJets = caSubjetFilterGenJets.clone()

addJetCollection(process, cms.InputTag("caVHPFJets:fat"),
    "CAVHFat", "PF",
    doJTA            = True,
    doBTagging       = True,
    jetCorrLabel     = ("AK5PF", jetCorr),
    doType1MET       = False,
    doL1Cleaning     = False,
    doL1Counters     = False,
    doJetID          = False,
    )

addJetCollection(process, cms.InputTag("caVHPFJets:sub"),
    "CAVHSub", "PF",
    doJTA            = True,
    doBTagging       = True,
    jetCorrLabel     = ("AK5PF", jetCorr),
    doType1MET       = False,
    doL1Cleaning     = False,
    doL1Counters     = False,
    genJetCollection = (cms.InputTag("caVHGenJets:sub") if options.isMC else None),
    doJetID          = False,
    )

addJetCollection(process, cms.InputTag("caVHPFJets:filter"),
    "CAVHFilter","PF",
    doJTA            = True,
    doBTagging       = True,
    jetCorrLabel     = ("AK5PF", jetCorr),
    doType1MET       = False,
    doL1Cleaning     = False,
    doL1Counters     = False,
    genJetCollection = (cms.InputTag("caVHGenJets:filter") if options.isMC else None),
    doJetID          = False,
    )

# Place appropriate jet cuts (NB: no cut on number of constituents)
defaultJetCut = cms.string("pt > 15. & abs(eta) < 5.0")
defaultFatJetCut = cms.string("pt > 100. & abs(eta) < 5.0")
process.selectedPatJets.cut = defaultJetCut
process.selectedPatJetsAK7PF.cut = defaultJetCut
#process.selectedPatJetsCAVHFatCalo.cut = defaultFatJetCut
process.selectedPatJetsCAVHFatPF.cut = defaultFatJetCut
process.selectedPatJetsCAVHSubPF.cut = cms.string("pt > 15. & abs(eta) < 5.0")
process.selectedPatJetsCAVHFilterPF.cut = cms.string("pt > 5. & abs(eta) < 5.0")

# ------------------------------------------------------------------------------
# Jet Substructure (FastJet 3)
# ------------------------------------------------------------------------------
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets

from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
ak5PrunedPFlow = ak5PFJetsPruned.clone(doAreaFastjet = cms.bool(True))

from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
ak5FilteredPFlow = ak5PFJetsFiltered.clone(doAreaFastjet = cms.bool(True))

from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsMassDropFiltered
ak5MassDropFilteredPFlow = ak5PFJetsMassDropFiltered.clone(doAreaFastjet = cms.bool(True))

#process.ca12GenJetsNoNu = ca4GenJets.clone( rParam = cms.double(1.2),src = cms.InputTag("genParticlesForJetsNoNu"))
process.ca12GenJets = ca4GenJets.clone( rParam = cms.double(1.2),src = cms.InputTag("genParticlesForJets"))
process.ca12PFJetsPFlow = ca4PFJets.clone(
    rParam = cms.double(1.2),
    src = cms.InputTag("pfNoElectron"+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax = cms.double(6.0),
    Ghost_EtaMax = cms.double(7.0)
    )
## this thing produces subjets by default
process.ca12PFJetsPrunedPFlow = ak5PrunedPFlow.clone(
    src = cms.InputTag("pfNoElectron"+postfix),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(1.2),
    jetAlgorithm = cms.string("CambridgeAachen"),
    #writeCompound = cms.bool(True), # this is used by default
    #jetCollInstanceName = cms.string("SubJets"), # this is used by default
    )
## this thing produces subjets by default
process.ca12PFJetsFilteredPFlow = ak5FilteredPFlow.clone(
    src = cms.InputTag("pfNoElectron"+postfix),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(1.2),
    jetAlgorithm = cms.string("CambridgeAachen"),
    )
## this thing produces subjets by default
process.ca12PFJetsMassDropFilteredPFlow = ak5MassDropFilteredPFlow.clone(
    src = cms.InputTag("pfNoElectron"+postfix),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(1.2),
    jetAlgorithm = cms.string("CambridgeAachen"),
    )

addJetCollection(process,
    cms.InputTag("ca12PFJetsPFlow"), # Jet collection; must be already in the event when patLayer0 sequence is executed
    "CA12", "PF",
    doJTA=True, # Run Jet-Track association & JetCharge
    doBTagging=True, # Run b-tagging
    jetCorrLabel=None,
    doType1MET=True,
    doL1Cleaning=False,
    doL1Counters=False,
    genJetCollection = (cms.InputTag("ca12GenJets") if options.isMC else None),
    doJetID = False
    )

addJetCollection(process,
    cms.InputTag("ca12PFJetsMassDropFilteredPFlow"), # Jet collection; must be already in the event when patLayer0 sequence is executed
    "CA12MassDropFiltered", "PF",
    doJTA=True, # Run Jet-Track association & JetCharge
    doBTagging=False, # Run b-tagging
    jetCorrLabel=None,
    doType1MET=True,
    doL1Cleaning=False,
    doL1Counters=False,
    #genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
    doJetID = False
    )

## adding the subjet collections which are b-tagged...
addJetCollection(process,
    cms.InputTag("ca12PFJetsMassDropFilteredPFlow", "SubJets"), # Jet collection; must be already in the event when patLayer0 sequence is executed
    "CA12MassDropFilteredSubjets", "PF",
    doJTA=True, # Run Jet-Track association & JetCharge
    doBTagging=True, # Run b-tagging
    jetCorrLabel=( "AK5PF", jetCorr ),
    doType1MET=True,
    doL1Cleaning=False,
    doL1Counters=False,
    #genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
    doJetID = False
    )

addJetCollection(process,
    cms.InputTag("ca12PFJetsFilteredPFlow", "SubJets"), # Jet collection; must be already in the event when patLayer0 sequence is executed
    "CA12FilteredSubjets", "PF",
    doJTA=True, # Run Jet-Track association & JetCharge
    doBTagging=True, # Run b-tagging
    jetCorrLabel=( "AK5PF", jetCorr ),
    doType1MET=True,
    doL1Cleaning=False,
    doL1Counters=False,
    #genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
    doJetID = False
    )

addJetCollection(process,
    cms.InputTag("ca12PFJetsPrunedPFlow", "SubJets"), # Jet collection; must be already in the event when patLayer0 sequence is executed
    "CA12PrunedSubjets", "PF",
    doJTA=True, # Run Jet-Track association & JetCharge
    doBTagging=True, # Run b-tagging
    jetCorrLabel=( "AK5PF", jetCorr ),
    doType1MET=True,
    doL1Cleaning=False,
    doL1Counters=False,
    #genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
    doJetID = False
    )

jetCutCA12    = "pt > 100."
subjetCutCA12 = "pt > 5."
process.selectedPatJetsCA12PF.cut = jetCutCA12
process.selectedPatJetsCA12MassDropFilteredPF.cut = jetCutCA12
process.selectedPatJetsCA12MassDropFilteredSubjetsPF.cut = subjetCutCA12
process.selectedPatJetsCA12FilteredSubjetsPF.cut = subjetCutCA12
process.selectedPatJetsCA12PrunedSubjetsPF.cut = subjetCutCA12

#load btag SF
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1107")

#do cleaning only against isolatted and ID"ed leptons (TODO: review the electron selection)
#process.selectedPatJets.checkOverlaps.muons.requireNoOverlaps = cms.bool(True)
#process.selectedPatJets.checkOverlaps.electrons.requireNoOverlaps = cms.bool(True)

#-------------------------------------------------
# MET uncertainty step
#-------------------------------------------------
#Embed the reference to the original jet in the jets, which is constant during the propagation
process.patJetsWithOwnRef = cms.EDProducer("PatObjectOwnRefProducer<pat::Jet>",
    src=cms.InputTag("selectedPatJets")
)

# type 1 +2 MET corrected
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")

#Note: this module causes a large memory increase when crossing the file boundary
#Reason - unknown, solution: limit processing to ~1 file.
from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties
runMEtUncertainties(process,
     electronCollection=cms.InputTag("selectedPatElectrons"),
     photonCollection=None,
     muonCollection=cms.InputTag("selectedPatMuons"),
     #FIXME: presently tau variations disabled, results in segfault
     tauCollection="", # "" means emtpy, None means cleanPatTaus
     jetCollection=cms.InputTag("selectedPatJets"),
     jetCorrLabel="L3Absolute" if options.isMC else "L2L3Residual",
     doSmearJets=options.isMC, #Note: switch this to False for the sync!
     jetCorrPayloadName="AK5PFchs",
     addToPatDefaultSequence=False
)

process.patPFMetNoPU = process.patMETs.clone(
    metSource = cms.InputTag("pfMETNoPU"),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag("genMetTrue")
)

process.pfMETNoPU = process.pfMET.clone()
process.pfMETNoPU.src=cms.InputTag("pfNoPileUp"+postfix)
#rocess.pfMETNoPU.jets = cms.InputTag("pfJets"+postfix)



process.pfNoPileUpCharge  = cms.EDFilter(
   "GenericPFCandidateSelector",
   src = cms.InputTag("pfNoPileUp"+postfix),
     cut = cms.string("charge!=0" )
)
process.pfMETNoPUCharge = process.pfMET.clone()
process.pfMETNoPUCharge.src=cms.InputTag("pfNoPileUpCharge")
process.pfMETNoPUCharge.calculateSignificance = cms.bool(False)

# load the PU JetID sequence
process.load("CMGTools.External.pujetidsequence_cff")
process.puJetId.jets =  cms.InputTag("selectedPatJets")
process.puJetMva.jets =  cms.InputTag("selectedPatJets")

# load the PU JetID sequence
process.load("CMGTools.External.pujetidsequence_cff")
process.puJetId.jets =  cms.InputTag("selectedPatJets")
process.puJetMva.jets =  cms.InputTag("selectedPatJets")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.savedGenParticles = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep++ pdgId >= 23 && pdgId <= 25", #keep W,Z,H and theirs products
    "keep++ pdgId == 22 && pt > 15", #keep gamma above 15 GeV
    "drop++   status == 2 ", #drop all non stable decay products (and daughters) [quarks are going to be added below]
    "keep++ abs(pdgId) == 15", #keep tau and its decay prods
    "keep  numberOfMothers() > 0 && abs(mother(0).pdgId) == 15", #keep first generation of tau daugthers (this is redundant I think)
    "drop  numberOfMothers() > 0 && abs(mother(0).pdgId) == {pi0}", #drop pi0 daugthers photons
    "keep  (abs(pdgId) ==13 || abs(pdgId) ==11 || abs(pdgId) ==15 ) &&  pt > 5.0", #keep leptons of decent pT
    "keep  (abs(pdgId) > 400 &&  abs(pdgId) < 600)    ||     (  (abs(pdgId) > 4000 &&  abs(pdgId) < 6000)  )",  # track-back the origin of B/D
    "keep  (  (abs(pdgId) >= 4 &&  abs(pdgId) <= 6)) ", #keep heavy quarks
    "keep ( status == 3)"  #keep event summary status3 (for pythia)
    )
)
### B Hadron truth
process.bhadrons = cms.EDProducer("MCBHadronProducer",
    quarkId = cms.uint32(5)
)

process.gen = cms.Sequence(  process.genParticlesForJets* process.caVHGenJets * process.ca12GenJets * process.bhadrons * process.savedGenParticles)

#Switch off checking for overlaps between leptons and jets
#process.selectedPatJetsNotOverlappingWithLeptonsForMEtUncertainty.checkOverlaps = cms.PSet()

process.out.outputCommands = cms.untracked.vstring([
    "drop *",
#    "keep *",

    "keep edmMergeableCounter_*_*_*", # Keep the lumi-block counter information
    "keep edmTriggerResults_TriggerResults__*", #Keep the trigger results
#    "keep *_genParticles__*", #keep all the genParticles
    "keep *_savedGenParticles__*",
    #"keep recoVertexs_offlinePrimaryVertices__*", #keep the offline PV-s
    "keep recoVertexs_goodOfflinePrimaryVertices__*", #keep the offline PV-s

    # Trigger
    "keep *_patTrigger_*_*",
    "keep *_patTriggerEvent_*_*",

    # Jets
    "keep patJets_*__*",
    "keep double_*_rho_*", #For rho-corr rel iso
    "keep recoGenJets_selectedPatJets_genJets_*", #For Jet MC smearing we need to keep the genJets
    "keep *_puJetId_*_*", # input variables
    "keep *_puJetMva_*_*", # final MVAs and working point flags

    # Muons
    "keep *_muons__*", #reco muons
    "keep patMuons_*__*",

    # Taus
    "keep patTaus_*__*",
    "keep recoPFTauDiscriminator_*__*",

    # Electrons
    "keep patElectrons_*__*",
    "keep *_electronClones__*",

    # METs
    "keep patMETs_*__*",

    #ECAL laser corr filter
    "keep bool_ecalLaserCorrFilter__*",

    #For flavour analyzer
    "keep GenEventInfoProduct_generator__*",

    #PU info
    "keep PileupSummaryInfos_addPileupInfo__*",

    ##PFCandidates
    "keep recoPFCandidates_*_pfCandidates_PAT",
    "keep recoPFMETs_pfMET__*",
    "keep recoPFMETs_pfMet__*",
    "keep recoGenMETs_genMetTrue__*",
    "keep recoPFCandidates_particleFlow__*",
    "keep recoConversions_allConversions__*",
    "keep recoVertexCompositeCandidates_generalV0Candidates_*_*",
    "keep recoTracks_generalTracks__*",
    "keep recoBeamSpot_offlineBeamSpot__*",
    "keep recoMuons_muons__*",

    "keep int_*__PAT",
    "keep ints_*__PAT",
    "keep double_*__PAT",
    "keep doubles_*__PAT",
    "keep float_*__PAT",
    "keep floats_*__PAT",

    #Hbb
    "keep *_HbbAnalyzerNew_*_*",
    "keep VHbbCandidates_*_*_*",
    "keep *_bcandidates_*_*",
    "keep *_bhadrons_*_*",
    "keep *_HLTDiCentralJet20MET80_*_*",
    "keep *_HLTDiCentralJet20MET100HBHENoiseFiltered_*_*",
    "keep *_HLTPFMHT150_*_*",
    "keep *_HLTQuadJet40_*_*",
    "keep *_HLTDoubleMu7_*_*",
    "keep *_EcalDeadCellEventFilter_*_*",
    "keep *_patType1CorrectedPFMet*_*_*",
    "keep *_patType1p2CorrectedPFMet*_*_*",
    "keep *_patMETsHT*_*_*",
    "keep patTriggerAlgorithms_patTrigger_*_*",
    "keep patTriggerConditions_patTrigger_*_*",
    "keep patTriggerObjects_patTrigger_*_*",
    "keep patTriggerFilters_patTrigger_*_*",
    "keep patTriggerPaths_patTrigger_*_*",
    "keep *_patTriggerEvent_*_*",
    "keep LHEEventProduct_*__LHE",
    "keep LHEEventProduct_*_*_LHE"
])

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring(
        ["*"]
    )
)



#-------------------------------------------------
# Paths
#-------------------------------------------------

#pre-count PV-s since opening the collection is slow
process.goodOfflinePVCount = cms.EDProducer(
    "CollectionSizeProducer<reco::Vertex>",
    src = cms.InputTag("goodOfflinePrimaryVertices")
)

process.preCalcSequences = cms.Sequence(
    process.patJetsWithOwnRef *
    process.gen
)

process.skimSequence = cms.Sequence(
    process.goodOfflinePrimaryVertices
    * process.goodOfflinePVCount
    #* process.eventFiltersSequence
    * process.patPF2PATSequence
)

process.GlobalTag.globaltag = cms.string(options.globalTag)

process.skimSequence += process.preCalcSequences
process.skimSequence += process.metUncertaintySequence
process.skimSequence += process.patTriggerSequence
process.skimSequence += process.muonSequence
process.skimSequence += process.electronSequence
process.skimSequence += process.tauSequence

#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagJetProbabilityCalibration?redirectedfrom=CMS.SWGuideBTagJetProbabilityCalibration#Calibration_in_53x_Data_and_MC
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
    tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
    connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
    tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
    connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
)

process.skimPath = cms.Path(process.skimSequence)
process.outPath = cms.EndPath(process.out)
