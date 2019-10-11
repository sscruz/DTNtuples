import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.StandardSequences.Eras import eras

import subprocess
import sys

options = VarParsing.VarParsing()

options.register('globalTag',
                 '106X_upgrade2023_realistic_v3', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global Tag")

options.register('nEvents',
                 -1, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Maximum number of processed events")

options.register('inputFolder',
                 '/eos/cms/store/group/dpg_dt/comm_dt/L1T_TDR/', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "EOS folder with input files")


options.register('secondaryInputFolder',
                 '', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "EOS folder with input files for secondary files")

options.register('applySegmentAgeing',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "If True applies ageing to RECO segments")

options.register('applyTriggerAgeing',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "If True applies ageing to trigger emulators")

options.register('ageingInput',
                 '', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Input with customised ageing, used only if non ''")

options.register('ntupleName',
                 './DTDPGNtuple_10_6_0_Phase2_Simulation.root', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Folder and name ame for output ntuple")

options.register('useRPC',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Use RPC info")

options.parseArguments()

process = cms.Process("DTNTUPLES",eras.Run2_2018)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.nEvents))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = cms.string(options.globalTag)

if options.ageingInput != "" :
    process.GlobalTag.toGet = cms.VPSet()
    process.GlobalTag.toGet.append(cms.PSet(record  = cms.string("MuonSystemAgingRcd"),
                                        connect = cms.string(options.ageingInput),
                                        tag     = cms.string("MuonSystemAging_3000fbm1")
                                        )
                               )
    
    #process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2023_realistic_v3', 'MuonSystemAging_3000fbm1,MuonSystemAgingRcd,sqlite_file:./MuonSystemAging.db')

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring()

)

#if options.inputFolder != '': 
#    files = subprocess.check_output(["ls", options.inputFolder])
#    process.source.fileNames = ["file://" + options.inputFolder + "/" + f for f in files.split()]
process.source.fileNames = ['/store/mc/PhaseIITDRSpring19DR/Mu_FlatPt2to100-pythia8-gun/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3-v2/70000/F42F882F-B3A8-4346-870D-3E62C930D076.root']


if options.secondaryInputFolder != "" :
    files = subprocess.check_output(["ls", options.secondaryInputFolder])
    process.source.secondaryFileNames = ["file://" + options.secondaryInputFolder + "/" + f for f in files.split()]

process.TFileService = cms.Service('TFileService',
        fileName = cms.string(options.ntupleName)
    )

process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D41_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

# process.DTGeometryESModule.applyAlignment = False
# process.DTGeometryESModule.fromDDD = False

process.load("Phase2L1Trigger.CalibratedDigis.CalibratedDigis_cfi") 
process.load("L1Trigger.DTPhase2Trigger.dtTriggerPhase2PrimitiveDigis_cfi")

process.CalibratedDigis.dtDigiTag = "simMuonDTDigis"
process.CalibratedDigis.scenario = 0
process.dtTriggerPhase2PrimitiveDigis.scenario = 0

#Produce RPC clusters from RPCDigi
process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
process.rpcRecHits.rpcDigiLabel = cms.InputTag('simMuonRPCDigis')
# Use RPC
process.load('Configuration.Geometry.GeometryExtended2023D38Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D38_cff')
if options.useRPC:
    process.dtTriggerPhase2PrimitiveDigis.useRPC = True
process.dtTriggerPhase2PrimitiveDigis.max_quality_to_overwrite_t0 = 10 # strict inequality
process.dtTriggerPhase2PrimitiveDigis.scenario = 0 # 0 for mc, 1 for data, 2 for slice test


process.dtTriggerPhase2AmPrimitiveDigis = process.dtTriggerPhase2PrimitiveDigis.clone()

process.load('L1Trigger.DTHoughTPG.DTTPG_cfi')

process.dtTriggerPhase2HbPrimitiveDigis = process.DTTPG.clone()
process.dtTriggerPhase2HbPrimitiveDigis.FirstBX = cms.untracked.int32(20)
process.dtTriggerPhase2HbPrimitiveDigis.LastBX = cms.untracked.int32(20)

process.load('RecoLocalMuon.Configuration.RecoLocalMuon_cff')
process.dt1DRecHits.dtDigiLabel = "simMuonDTDigis"

process.load('DTDPGAnalysis.DTNtuples.dtNtupleProducer_phase2_cfi')

process.p = cms.Path(process.dt1DRecHits
                     + process.dt4DSegments
                     + process.rpcRecHits
                     + process.CalibratedDigis
                     + process.dtTriggerPhase2AmPrimitiveDigis
                     + process.dtTriggerPhase2HbPrimitiveDigis
                     + process.dtNtupleProducer)

from DTDPGAnalysis.DTNtuples.customiseDtNtuples_cff import customiseForRunningOnMC, customiseForFakePhase2Info, customiseForAgeing
customiseForAgeing(process,"p",options.applySegmentAgeing,options.applyTriggerAgeing)
customiseForRunningOnMC(process,"p")
customiseForFakePhase2Info(process)




