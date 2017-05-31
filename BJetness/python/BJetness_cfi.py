import FWCore.ParameterSet.Config as cms

demo = cms.EDProducer('BJetness',
  # Format
  is_data   = cms.bool(False),
  # Input tags
  vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
  patElectrons        = cms.InputTag("slimmedElectrons"),
  eleMVAnonTrigIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
  muons               = cms.InputTag("slimmedMuons"),
  jets                = cms.InputTag("slimmedJets"),
  jecPayloadNamesAK4PFchsMC1 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC2 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC3 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA1 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA2 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA3 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA4 = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt"),
  jerAK4PFchs   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt"),
  jerAK4PFchsSF = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt")
)
