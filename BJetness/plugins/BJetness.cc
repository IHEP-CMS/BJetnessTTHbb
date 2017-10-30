// -*- C++ -*-
//
// Package:    FlavorJetness/BJetness
// Class:      BJetness
//
/**\class BJetness BJetness.cc FlavorJetness/BJetness/plugins/BJetness.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Authors:  Francesco Romeo & Joshuha Thomas-Wilsker
//         Created:  Sun, 30 Oct 2016 11:02:22 GMT
//
//

#include "BJetnessTTHbb/BJetness/interface/BJetness.h"

#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//Track builder infos
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
#include "TMath.h"

KalmanVertexFitter vtxFitterFV(true);
using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;

//
// Class definition
//
class BJetness : public edm::stream::EDProducer<> {

  public:

    explicit BJetness(const edm::ParameterSet&);
    ~BJetness();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    // ===== Methods aligned to the TTHbb selection =====
    void JECInitialization(std::vector<boost::shared_ptr<JetCorrectionUncertainty> > &factorisedJECsMC, std::vector<boost::shared_ptr<JetCorrectionUncertainty> > &factorisedJECsDATA);
    void GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN, const edm::EventSetup& iSetup);
    bool isGoodVertex(const reco::Vertex& vtx);
    bool is_loose_muon(const pat::Muon& mu, const reco::Vertex& vtx);
    bool is_tight_muon(const pat::Muon& mu, const reco::Vertex& vtx);
    double rel_iso_dbc_mu(const pat::Muon& lepton);
    bool is_loose_electron(const pat::Electron& ele, double rhopog);//, const reco::Vertex& vtx);
    bool is_tight_electron(const pat::Electron& ele, double rhopog);//, const reco::Vertex& vtx);
    double rel_iso_dbc_ele(const pat::Electron& lepton, double rhopog);
    double get_effarea(double eta);
    bool is_good_jet(const pat::Jet &j, double rho, double rhoJER, int vtxsize, const edm::EventSetup& iSetup, double & jetpt, boost::shared_ptr<JetCorrectionUncertainty> temp_JECUncMC, boost::shared_ptr<JetCorrectionUncertainty> temp_JECUncDATA, string sysUncSource_, string is_updown);
    //===========================================================

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    // Ensure the class recognises containers we will book into event:
    typedef std::vector<float> FlavorJetnessValues;
    //JEC systematic variations of BJetness values.
    typedef std::vector< std::vector<float> > FlavorJetnessValues_JECSysts;
    //typedef std::vector< FlavorJetnessValues > FlavorJetnessValues_JECSysts;

    bool _is_data;
    edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
    edm::EDGetTokenT<edm::View<pat::Electron> > electron_pat_;
    edm::EDGetTokenT<edm::ValueMap<bool>  > elemvanontrig_;
    edm::EDGetTokenT<edm::ValueMap<bool>  > electronTightIdMapToken_;
    edm::EDGetTokenT<edm::View<pat::Muon> > muon_h_;
    edm::EDGetTokenT<pat::JetCollection> jets_;
    edm::EDGetTokenT<double> rhopogHandle_;
    edm::EDGetTokenT<double> rhoJERHandle_;
    int    _vtx_ndof_min;
    int    _vtx_rho_max;
    double _vtx_position_z_max;
    edm::FileInPath jecPayloadNamesAK4PFchsMC1_;
    edm::FileInPath jecPayloadNamesAK4PFchsMC2_;
    edm::FileInPath jecPayloadNamesAK4PFchsMC3_;
    edm::FileInPath jecPayloadNamesAK4PFchsMCUnc_;
    edm::FileInPath jecPayloadNamesAK4PFchsDATA1_;
    edm::FileInPath jecPayloadNamesAK4PFchsDATA2_;
    edm::FileInPath jecPayloadNamesAK4PFchsDATA3_;
    edm::FileInPath jecPayloadNamesAK4PFchsDATA4_;
    edm::FileInPath jecPayloadNamesAK4PFchsDATAUnc_;
    std::string jerAK4PFchs_;
    std::string jerAK4PFchsSF_;
    boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsMC_;
    boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsMCUnc_;
    //===========================================================
    boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsMCUnc_AbsoluteStat;
    //===========================================================

    boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsDATA_;
    boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsDATAUnc_;
    //===========================================================
    boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsDATAUnc_AbsoluteStat;
    //===========================================================

    std::vector<boost::shared_ptr<JetCorrectionUncertainty> > factorisedJECsMC_;
    std::vector<boost::shared_ptr<JetCorrectionUncertainty> > factorisedJECsDATA_;


    //Methods for the BJetness variables
    void get_bjetness_vars(
                           vector<pat::Jet> evtjets, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h,
                           double& bjetnessFV_num_leps, double& bjetnessFV_npvTrkOVcollTrk, double& bjetnessFV_pvTrkOVcollTrk, double& bjetnessFV_avip3d_val, double& bjetnessFV_avip3d_sig, double& bjetnessFV_avsip3d_val, double& bjetnessFV_avsip3d_sig, double& bjetnessFV_avip2d_val, double& bjetnessFV_avip2d_sig, double& bjetnessFV_avsip2d_val, double& bjetnessFV_avsip2d_sig, double& bjetnessFV_avip1d_val, double& bjetnessFV_avip1d_sig, double& bjetnessFV_avsip1d_val, double& bjetnessFV_avsip1d_sig, double& bjetnessFV_npvTrkOVpvTrk, double &bjetnessFV_npvPtOVcollPt, double &bjetnessFV_pvPtOVcollPt, double &bjetnessFV_npvPtOVpvPt
                          );
    void get_bjetness_trkinfos(vector<pat::Jet> evtjets, const reco::Vertex& vtx, vector<Track>& jetchtrks, double& bjetness_num_pvtrks, double& bjetness_num_npvtrks, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h, double& bjetness_num_eles, double& bjetness_num_mus, vector<tuple<double, double, double> >& jetsdir, vector<Track>& jetchtrkspv, vector<Track>& jetchtrksnpv);
    bool is_goodtrk(Track trk,const reco::Vertex& vtx);
    bool is_loosePOG_jetmuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h);
    bool is_softLep_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx);
    bool is_vetoPOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat);
    bool is_loosePOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx);
    void get_avip3d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip3d_val, double& jetchtrks_avip3d_sig, double& jetchtrks_avsip3d_val, double& jetchtrks_avsip3d_sig);
    void get_avip2d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip2d_val, double& jetchtrks_avip2d_sig, double& jetchtrks_avsip2d_val, double& jetchtrks_avsip2d_sig);
    void get_avip1d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip1d_val, double& jetchtrks_avip1d_sig, double& jetchtrks_avsip1d_val, double& jetchtrks_avsip1d_sig);
    vector<TransientTrack> get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder);
};


//
// Class constructor.
//
BJetness::BJetness(const edm::ParameterSet& iConfig):
  vtx_h_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  electron_pat_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectrons"))),
  elemvanontrig_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVAnonTrigIdMap"))),
  //electronMediumIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"))),
  electronTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"))),
  muon_h_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"))),
  jets_(consumes<pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jets"))),
  rhopogHandle_(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"))),
  rhoJERHandle_(consumes<double>(edm::InputTag("fixedGridRhoAll")))
{
  // Read info in from config:
  _is_data = iConfig.getParameter<bool>("is_data");
  _vtx_ndof_min       = 4;
  _vtx_rho_max        = 2;
  _vtx_position_z_max = 24;
  jecPayloadNamesAK4PFchsMC1_     = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC1");
  jecPayloadNamesAK4PFchsMC2_     = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC2");
  jecPayloadNamesAK4PFchsMC3_     = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMC3");
  jecPayloadNamesAK4PFchsMCUnc_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsMCUnc");
  jecPayloadNamesAK4PFchsDATA1_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA1");
  jecPayloadNamesAK4PFchsDATA2_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA2");
  jecPayloadNamesAK4PFchsDATA3_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA3");
  jecPayloadNamesAK4PFchsDATA4_   = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATA4");
  jecPayloadNamesAK4PFchsDATAUnc_ = iConfig.getParameter<edm::FileInPath>("jecPayloadNamesAK4PFchsDATAUnc");
  jerAK4PFchs_                    = iConfig.getParameter<edm::FileInPath>("jerAK4PFchs").fullPath();
  jerAK4PFchsSF_                  = iConfig.getParameter<edm::FileInPath>("jerAK4PFchsSF").fullPath();

  //Instantiate JEC uncertainty sources and push into vectors:
  JECInitialization(factorisedJECsMC_ , factorisedJECsDATA_);

  // Register class products with in framework.
  produces<FlavorJetnessValues>("BJetnessValue").setBranchAlias("BJetnessValues");
  produces<FlavorJetnessValues_JECSysts>("BJetnessValueJECSysts").setBranchAlias("BJetnessValueJECSysts");
  //produces<std::vector< std::vector<float> > >("BJetnessValueJECSysts").setBranchAlias("BJetnessValueJECSysts");


}

//Class destructor
BJetness::~BJetness(){
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// ------------ Create, fill and put products into events ------------
void BJetness::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::cout << "Event number = " << iEvent.id().event() << std::endl;
  //Store bjetness vars
  auto_ptr<FlavorJetnessValues> BJetnessValues( new FlavorJetnessValues );
  auto_ptr<FlavorJetnessValues_JECSysts> BJetnessValueJECSysts( new FlavorJetnessValues_JECSysts );
  //auto_ptr<std::vector< std::vector<float> > > BJetnessValueJECSysts( new std::vector< std::vector<float> > );

  std::vector<float> bJetnessVariables;

  const char* sysUncNames[] = {"Nominal"
                              /*"Total",
                              "AbsoluteStat",
                              "AbsoluteScale",
                              "AbsoluteFlavMap",
                              "AbsoluteMPFBias",
                              "Fragmentation",
                              "SinglePionECAL",
                              "SinglePionHCAL",
                              "FlavorQCD",
                              "TimePtEta",
                              "RelativeJEREC1",
                              "RelativeJEREC2",
                              "RelativeJERHF",
                              "RelativePtBB",
                              "RelativePtEC1",
                              "RelativePtEC2",
                              "RelativePtHF",
                              "RelativeBal",
                              "RelativeFSR",
                              "RelativeStatFSR",
                              "RelativeStatEC",
                              "RelativeStatHF",
                              "PileUpDataMC",
                              "PileUpPtRef",
                              "PileUpPtBB",
                              "PileUpPtEC1",
                              "PileUpPtEC2",
                              "PileUpPtHF",
                              "PileUpMuZero",
                              "PileUpEnvelope",
                              "SubTotalPileUp",
                              "SubTotalRelative",
                              "SubTotalPt",
                              "SubTotalScale",
                              "SubTotalAbsolute",
                              "SubTotalMC",
                              "TotalNoFlavor",
                              "TotalNoFlavorNoTime",
                              "FlavorZJet",
                              "FlavorPhotonJet",
                              "FlavorPureGluon",
                              "FlavorPureQuark",
                              "FlavorPureCharm",
                              "FlavorPureBottom"*/};


  const char* updownNames[] = {"_UP_","_DOWN_"};

  std::vector<std::string> sysUncSources(sysUncNames,std::end(sysUncNames));
  std::vector<std::string> updown(updownNames,std::end(updownNames));

  //Loop over sources systematic uncertainty

  if (factorisedJECsMC_.size() != factorisedJECsDATA_.size()){
    std::cout << "BJetness::WARNING: Crash imminent. Different number of factorised JEC uuncertainty sources for data and MC." << std::endl;
    std::cout << "BJetness:: WARNING: Check sources JECInitialisation method." <<std::endl;
  }

  // Declare temp pointer to be initialised with each systematic uncertainty.
  boost::shared_ptr<JetCorrectionUncertainty> temp_JECUncMC;
  boost::shared_ptr<JetCorrectionUncertainty> temp_JECUncDATA;

  for (uint h=0; h<BJetnessValueJECSysts->size(); h++){
    (BJetnessValueJECSysts->at(h)).clear();
  }

  for(uint i=0; i<sysUncSources.size(); i++){

    uint updown_count;
    if (sysUncSources[i]=="Nominal"){
      updown_count=1;
    }
    else{
      updown_count=2;
    }

    for(uint k=0; k<updown_count; k++){
      bJetnessVariables->clear();

      cout << "sysUncSource = " << sysUncSources[i] << " ; up/down count = " << updown_count << endl;

      if (sysUncSources[i]!="Nominal"){
        //std::cout << sysUncSources.at(i-1) << std::endl;
        temp_JECUncMC = factorisedJECsMC_.at(i-1);
        temp_JECUncDATA = factorisedJECsDATA_.at(i-1);
      }

      /////
      //   Recall collections
      /////
      edm::Handle<reco::VertexCollection> vtx_h;
      iEvent.getByToken(vtx_h_, vtx_h);
      edm::Handle<edm::View<pat::Electron> > electron_pat;
      iEvent.getByToken(electron_pat_, electron_pat);
      edm::Handle<edm::View<pat::Muon> > muon_h;
      iEvent.getByToken(muon_h_, muon_h);
      edm::Handle<pat::JetCollection> jets;
      iEvent.getByToken(jets_, jets);
      edm::Handle<double> rhopogHandle;
      iEvent.getByToken(rhopogHandle_,rhopogHandle);
      double rhopog = *rhopogHandle;
      edm::Handle<double> rhoJERHandle;
      iEvent.getByToken(rhoJERHandle_,rhoJERHandle);
      double rhoJER = *rhoJERHandle;
      edm::Handle<edm::ValueMap<bool>  > mvanontrig_id_decisions;
      iEvent.getByToken(elemvanontrig_, mvanontrig_id_decisions);
      edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
      iEvent.getByToken(electronTightIdMapToken_, tight_id_decisions);

      edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);

      /////
      //   First clean the jet according the TTHbb selection
      /////
      // Require a good vertex (This part has to be clarified (do we want the first PV to be a good one?))
      if(vtx_h->empty()){

        //////////////////////////////////////////
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);


        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);

        BJetnessValueJECSysts->push_back(bJetnessVariables);
        //////////////////////////////////////////
        cout << "Event failed vertex cuts." << endl;
        cout << "Vector of bJetness vars. pushed into BJetnessValueJECSysts" << endl;
        cout << "Put bJetness object into Event" << endl;
        iEvent.put(BJetnessValueJECSysts,"BJetnessValueJECSysts");
        //iEvent.put(BJetnessValues,"BJetnessValue");
        return; // skip the event if no PV found
      }
      reco::VertexCollection::const_iterator firstgoodVertex = vtx_h->end();
      for(reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstgoodVertex; it++){
        if(isGoodVertex(*it)){
          firstgoodVertex = it;
          break;
        }
      }
      if(firstgoodVertex == vtx_h->end()){

        ////////////////////////////////////////////////
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);


        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValueJECSysts->push_back(bJetnessVariables);
        ////////////////////////////////////////////////
        cout << "Event failed first good vertex cut." << endl;
        cout << "Vector of bJetness vars. pushed into BJetnessValueJECSysts" << endl;
        cout << "Put bJetness object into Event" << endl;
        iEvent.put(BJetnessValueJECSysts,"BJetnessValueJECSysts");
        //iEvent.put(BJetnessValues,"BJetnessValue");
        return;
      }
      const reco::Vertex &PV = vtx_h->front(); //It still takes the first vertex

      //Look for muons, electrons
      vector<const reco::Candidate*> looseleps;
      vector<const reco::Candidate*> tightleps;

      //================ Muons ======================
      for(const pat::Muon &mu : *muon_h){
        if(!is_loose_muon(mu,PV)) continue;
        looseleps.push_back((const reco::Candidate*)&mu);
        if(!is_tight_muon(mu,PV)) continue;
        tightleps.push_back((const reco::Candidate*)&mu);
      }
      //====================================================

      //================ Electrons ==================

      for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
        const Ptr<pat::Electron> elPtr(electron_pat, ele - electron_pat->begin() );
        bool isPassEletrig = (*tight_id_decisions)[ elPtr ];
        const pat::Electron &lele = *ele;
        if(!(isPassEletrig)) continue;
        if(!(is_loose_electron(*ele,rhopog))) continue;
        if(!(rel_iso_dbc_ele(*ele,rhopog)<0.15)) continue;
        looseleps.push_back((const reco::Candidate*)&lele);
        if(!(is_tight_electron(*ele,rhopog))) continue;
        tightleps.push_back((const reco::Candidate*)&lele);
      }
      //====================================================


      // =================== Jets =========================
      //Get the good jets of the event
      //Iterate to access jet by decreasing b-tagging value
      int jet_pos = 0; //This counter helps to order jets
      int jet_num = 0; //This counter accounts for the number of good jets in the events
      //The definition of good jet in the event must be the same of the TTHbb analysis
      //so that jet_num corresponds to the number of jets that define the categories in the TTHbb search
      int jetb_num = 0;
      vector<pair<double,int> > jet_csv_pos;
      vector<pair<double,int> > jet_cmva_pos;


      //Check memory address is changing for each systematic.

      int jetcounter = 0;
      for(const pat::Jet &j : *jets){
        int vtxsize = vtx_h->size();
        double jetpt = 0;
        ++jetcounter;

        // temp_JECUnc's now needs to be copied into is_good_jet method so pointer to
        // JetCorrectionUncertainty points at memory block of relevant JEC unc. source.
        if(!is_good_jet(j,rhopog,rhoJER,vtxsize,iSetup,jetpt, temp_JECUncMC, temp_JECUncDATA, sysUncSources[i], updown[k])){jet_pos++; continue;}
        bool jetmatchedlepts = false;
        for(uint gl=0; gl<looseleps.size(); gl++) if(deltaR(looseleps[gl]->p4(),j.p4())<0.4) jetmatchedlepts = true;
        if(jetmatchedlepts){jet_pos++; continue;}
        double csvcurrjet = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet_csv_pos.push_back(make_pair(csvcurrjet,jet_pos));
        double cmvacurrjet = j.bDiscriminator("pfCombinedMVAV2BJetTags");
        jet_cmva_pos.push_back(make_pair(cmvacurrjet,jet_pos));
        if(csvcurrjet>0.8484) jetb_num++;
        jet_pos++;
        jet_num++;
      }
      //=======================================================

      /////
      //   Select only TTHbb events (mainly lep sel)
      /////

      //This selection should be changed for dilepton analysis.
      if(!(tightleps.size()==1 && looseleps.size()==1 && jet_num>=4 && jetb_num>=2)){


        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);


        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValueJECSysts->push_back(bJetnessVariables);
        ////////////////////////////////////////////////
        cout << "Event failed event selection." << endl;
        cout << "Vector of bJetness vars. pushed into BJetnessValueJECSysts" << endl;
        cout << "Put bJetness object into Event" << endl;
        iEvent.put(BJetnessValueJECSysts,"BJetnessValueJECSysts");
        //iEvent.put(BJetnessValues,"BJetnessValue");
        return;
      }

      /////
      //   You need to provide as input the jets selected in the event (selection according to the TTHbb analysis),
      //   which have to be ordered by decreasing b-tagging value
      /////

      sort(jet_csv_pos.rbegin(), jet_csv_pos.rend());//Order by descreasing csv value
      if(jet_num!=0){
        /////
        //   Definition of BJetness variables
        //   NOTE: Highest CSV jet excluded (start from jn=1 below) and consider up to 6th jet in an event (optimised for single lepton)
        /////
        vector<pat::Jet> evtjets; evtjets.clear();
        int maxjetnum = 6;
        if(jet_num<maxjetnum) maxjetnum = jet_num;
        for(int jn=1; jn<maxjetnum; jn++) evtjets.push_back((*jets)[jet_csv_pos[jn].second]);
        //Define the variables you want to access
        double bjetnessFV_num_leps        = -1;
        double bjetnessFV_npvTrkOVcollTrk = -1;
        double bjetnessFV_pvTrkOVcollTrk = -1;
        double bjetnessFV_npvTrkOVpvTrk =-1;
        double bjetnessFV_npvPtOVcollPt =-1;
        double bjetnessFV_pvPtOVcollPt = -1;
        double bjetnessFV_npvPtOVpvPt =-1;
        double bjetnessFV_avip3d_val      = -1;
        double bjetnessFV_avip3d_sig      = -1;
        double bjetnessFV_avsip3d_val     = -1;
        double bjetnessFV_avsip3d_sig     = -1;

        double bjetnessFV_avip2d_val      = -1;
        double bjetnessFV_avip2d_sig      = -1;
        double bjetnessFV_avsip2d_val      = -1;
        double bjetnessFV_avsip2d_sig      = -1;

        double bjetnessFV_avip1d_val      = -1;
        double bjetnessFV_avip1d_sig      = -1;
        double bjetnessFV_avsip1d_val      = -1;
        double bjetnessFV_avsip1d_sig      = -1;

        //This is the method to access the BJetness variables
        get_bjetness_vars(
          //Inputs:
          evtjets,      //Jets used to build the BJetness jets
          PV,           //Prinary vertex of the event
          *ttrkbuilder, //Transient tracker builder to measure impact parameters
          electron_pat, muon_h, //Leptons collections to count the number of electrons and muons
          //BJetness variables
          bjetnessFV_num_leps, bjetnessFV_npvTrkOVcollTrk, bjetnessFV_pvTrkOVcollTrk, bjetnessFV_avip3d_val, bjetnessFV_avip3d_sig, bjetnessFV_avsip3d_val, bjetnessFV_avsip3d_sig, bjetnessFV_avip2d_val, bjetnessFV_avip2d_sig, bjetnessFV_avsip2d_val, bjetnessFV_avsip2d_sig, bjetnessFV_avip1d_val, bjetnessFV_avip1d_sig, bjetnessFV_avsip1d_val, bjetnessFV_avsip1d_sig, bjetnessFV_npvTrkOVpvTrk, bjetnessFV_npvPtOVcollPt, bjetnessFV_pvPtOVcollPt, bjetnessFV_npvPtOVpvPt
        );
        bjetnessFV_num_leps = 2;
        //Fill the quantities for the event

        ///////////////////////////////////////////////////////////////
        //Num_of_trks
        bJetnessVariables.push_back(bjetnessFV_num_leps);
        bJetnessVariables.push_back(bjetnessFV_npvTrkOVcollTrk);
        bJetnessVariables.push_back(bjetnessFV_pvTrkOVcollTrk);
        bJetnessVariables.push_back(bjetnessFV_npvTrkOVpvTrk);
        bJetnessVariables.push_back(bjetnessFV_npvPtOVcollPt);
        bJetnessVariables.push_back(bjetnessFV_pvPtOVcollPt);
        bJetnessVariables.push_back(bjetnessFV_npvPtOVpvPt);

        //ImpactParameter
        bJetnessVariables.push_back(bjetnessFV_avip3d_val);
        bJetnessVariables.push_back(bjetnessFV_avip3d_sig);
        bJetnessVariables.push_back(bjetnessFV_avsip3d_val);
        bJetnessVariables.push_back(bjetnessFV_avsip3d_sig);

        bJetnessVariables.push_back(bjetnessFV_avip2d_val);
        bJetnessVariables.push_back(bjetnessFV_avip2d_sig);
        bJetnessVariables.push_back(bjetnessFV_avsip2d_val);
        bJetnessVariables.push_back(bjetnessFV_avsip2d_sig);

        bJetnessVariables.push_back(bjetnessFV_avip1d_val);
        bJetnessVariables.push_back(bjetnessFV_avip1d_sig);
        bJetnessVariables.push_back(bjetnessFV_avsip1d_val);
        bJetnessVariables.push_back(bjetnessFV_avsip1d_sig);


        BJetnessValues->push_back(bjetnessFV_num_leps);
        BJetnessValues->push_back(bjetnessFV_npvTrkOVcollTrk);
        BJetnessValues->push_back(bjetnessFV_pvTrkOVcollTrk);
        BJetnessValues->push_back(bjetnessFV_npvTrkOVpvTrk);
        BJetnessValues->push_back(bjetnessFV_npvPtOVcollPt);
        BJetnessValues->push_back(bjetnessFV_pvPtOVcollPt);
        BJetnessValues->push_back(bjetnessFV_npvPtOVpvPt);

        //ImpactParameter
        BJetnessValues->push_back(bjetnessFV_avip3d_val);
        BJetnessValues->push_back(bjetnessFV_avip3d_sig);
        BJetnessValues->push_back(bjetnessFV_avsip3d_val);
        BJetnessValues->push_back(bjetnessFV_avsip3d_sig);

        BJetnessValues->push_back(bjetnessFV_avip2d_val);
        BJetnessValues->push_back(bjetnessFV_avip2d_sig);
        BJetnessValues->push_back(bjetnessFV_avsip2d_val);
        BJetnessValues->push_back(bjetnessFV_avsip2d_sig);

        BJetnessValues->push_back(bjetnessFV_avip1d_val);
        BJetnessValues->push_back(bjetnessFV_avip1d_sig);
        BJetnessValues->push_back(bjetnessFV_avsip1d_val);
        BJetnessValues->push_back(bjetnessFV_avsip1d_sig);

        ///////////////////////////////////////////////////////////////
        cout << "Event passed." << endl;
        cout << "Vector of bJetness vars. pushed into BJetnessValueJECSysts" << endl;
        BJetnessValueJECSysts->push_back(bJetnessVariables);
        ///////////////////////////////////////////////////////////////

      }else{//if(jet_num!=0)

        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);
        BJetnessValues->push_back(-999);

        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);
        bJetnessVariables.push_back(-999);

        cout << "Event failed due to number of jets < 4." << endl;
        cout << "Vector of bJetness vars. pushed into BJetnessValueJECSysts" << endl;
        BJetnessValueJECSysts->push_back(bJetnessVariables);

      }
    }// End of variation up/down loop.
  }//End of uncertainty source loop.
  cout << "======= Final BJetnessValueJECSysts =======" << endl;
  cout << "No. JEC Systematics = " << BJetnessValueJECSysts->size() << endl;
  for (uint o=0; o<BJetnessValueJECSysts->size(); o++){
    cout << "=========== Systematic # " << o << " ==========" << endl;
    cout << "No. BJetness Variables = " << (BJetnessValueJECSysts->at(o)).size() << endl;
    for (uint p=0; p<(BJetnessValueJECSysts->at(o)).size(); p++){
      cout << "=========== BJetness variable # " << p << " ==========" << endl;
      cout << (BJetnessValueJECSysts->at(o)).at(p) << endl;
    }
  }
  cout << "=====================================" << endl;
  cout << "Put bJetness object into Event" << endl;
  iEvent.put(BJetnessValueJECSysts,"BJetnessValueJECSysts");
  //iEvent.put(BJetnessValues,"BJetnessValue");
}

/////
//   Methods to be aligned to the TTHbb selection
/////
//Ask for good vertices
bool BJetness::isGoodVertex(const reco::Vertex& vtx){
  if(vtx.isFake())                                   return false;
  if(vtx.ndof()<_vtx_ndof_min)                       return false;
  if(vtx.position().Rho()>_vtx_rho_max)              return false;
  if(fabs(vtx.position().Z()) > _vtx_position_z_max) return false;
  return true;
}
//Look for loose muon (definition for the jet cleaning)
bool BJetness::is_loose_muon(const pat::Muon& mu, const reco::Vertex& vtx){
  bool isloosemu = false;
  if(mu.pt()>15 &&
    TMath::Abs(mu.eta()) < 2.4 &&
    mu.isTightMuon(vtx) &&
    rel_iso_dbc_mu(mu) < 0.25
    ) isloosemu = true;
  return isloosemu;
}
bool BJetness::is_tight_muon(const pat::Muon& mu, const reco::Vertex& vtx){
  bool isloosemu = false;
  if(mu.pt()>26 &&
    TMath::Abs(mu.eta()) < 2.1 &&
    mu.isTightMuon(vtx) &&
    rel_iso_dbc_mu(mu) < 0.15
    ) isloosemu = true;
  return isloosemu;
}
double BJetness::rel_iso_dbc_mu(const pat::Muon& lepton){
  return((lepton.pfIsolationR04().sumChargedHadronPt + max(lepton.pfIsolationR04().sumNeutralHadronEt + lepton.pfIsolationR04().sumPhotonEt - 0.5 * lepton.pfIsolationR04().sumPUPt,0.0))/lepton.pt()
        );
}
//Look for loose electron (definition for the jet cleaning)
bool BJetness::is_loose_electron(const pat::Electron& ele, double rhopog){
  bool isele = false;
  if(ele.pt()>15 && TMath::Abs(ele.eta())<2.4 && !(fabs(ele.superCluster()->position().eta()) > 1.4442 && fabs(ele.superCluster()->position().eta()) < 1.5660)){
      isele = true;
  }
  return isele;
}
bool BJetness::is_tight_electron(const pat::Electron& ele, double rhopog){
  bool isele = false;
  if(ele.pt()>30 && TMath::Abs(ele.eta())<2.1 && !(fabs(ele.superCluster()->position().eta()) > 1.4442 && fabs(ele.superCluster()->position().eta()) < 1.5660)){
      isele = true;
  }
  return isele;
}
double BJetness::rel_iso_dbc_ele(const pat::Electron& el, double rhopog){
  double SumChHadPt       = el.pfIsolationVariables().sumChargedHadronPt;
  double SumNeuHadEt      = el.pfIsolationVariables().sumNeutralHadronEt;
  double SumPhotonEt      = el.pfIsolationVariables().sumPhotonEt;
  double EffArea          = get_effarea(el.superCluster()->position().eta());
  double SumNeutralCorrEt = std::max( 0.0, SumNeuHadEt+SumPhotonEt - rhopog*EffArea );
  double relIsoRhoEA = (SumChHadPt + SumNeutralCorrEt)/el.pt();
  return relIsoRhoEA;
}
double BJetness::get_effarea(double eta){
  double effarea = -1;
  if(abs(eta) < 1.0)        effarea = 0.1752;
  else if(abs(eta) < 1.479) effarea = 0.1862;
  else if(abs(eta) < 2.0)   effarea = 0.1411;
  else if(abs(eta) < 2.2)   effarea = 0.1534;
  else if(abs(eta) < 2.3)   effarea = 0.1903;
  else if(abs(eta) < 2.4)   effarea = 0.2243;
  else                      effarea = 0.2687;
  return effarea;
}


// Require good jets (according to TTHbb analysis)
bool BJetness::is_good_jet(const pat::Jet &j,double rho, double rhoJER, int vtxsize, const edm::EventSetup& iSetup, double& jetpt, boost::shared_ptr<JetCorrectionUncertainty> temp_JECUncMC_, boost::shared_ptr<JetCorrectionUncertainty> temp_JECUncDATA_, string sysUncSource_, string is_updown_){

  bool isgoodjet = true;

  //Jet Energy Corrections and Uncertainties
  double corrAK4PFchs     = 1;
  double corrAK4PFchs_Unc = 1;

  //==================================Factorised JES Unc.=========================================

  reco::Candidate::LorentzVector uncorrJetAK4PFchs = j.correctedP4(0);
  if(!_is_data){
    jecAK4PFchsMC_->setJetEta( uncorrJetAK4PFchs.eta() );
    jecAK4PFchsMC_->setJetPt ( uncorrJetAK4PFchs.pt() );
    jecAK4PFchsMC_->setJetE  ( uncorrJetAK4PFchs.energy() );
    jecAK4PFchsMC_->setRho( rho  );
    jecAK4PFchsMC_->setNPV( vtxsize  );
    jecAK4PFchsMC_->setJetA( j.jetArea() );
    corrAK4PFchs = jecAK4PFchsMC_->getCorrection();

    if(sysUncSource_!="Nominal"){
      if (is_updown_.find("_UP_") != std::string::npos) {
        temp_JECUncMC_->setJetEta( uncorrJetAK4PFchs.eta() );
        temp_JECUncMC_->setJetPt( corrAK4PFchs * uncorrJetAK4PFchs.pt() );
        corrAK4PFchs_Unc = corrAK4PFchs * (1 + fabs(temp_JECUncMC_->getUncertainty(1)));
      }
      else if (is_updown_.find("_DOWN_") != std::string::npos) {
        temp_JECUncMC_->setJetEta( uncorrJetAK4PFchs.eta() );
        temp_JECUncMC_->setJetPt( corrAK4PFchs * uncorrJetAK4PFchs.pt() );
        corrAK4PFchs_Unc = corrAK4PFchs * ( 1 - fabs(temp_JECUncMC_->getUncertainty(-1)) );
      }
    }
  }
  else {
    jecAK4PFchsDATA_->setJetEta( uncorrJetAK4PFchs.eta()    );
    jecAK4PFchsDATA_->setJetPt ( uncorrJetAK4PFchs.pt()     );
    jecAK4PFchsDATA_->setJetE  ( uncorrJetAK4PFchs.energy() );
    jecAK4PFchsDATA_->setRho    ( rho  );
    jecAK4PFchsDATA_->setNPV    ( vtxsize  );
    jecAK4PFchsDATA_->setJetA  ( j.jetArea()         );
    corrAK4PFchs = jecAK4PFchsDATA_->getCorrection();

    //====================================Factorised JES===========================================

    if(sysUncSource_!="Nominal"){
      if (is_updown_.find("_UP_") != std::string::npos) {
        temp_JECUncDATA_->setJetEta( uncorrJetAK4PFchs.eta() );
        temp_JECUncDATA_->setJetPt( corrAK4PFchs * uncorrJetAK4PFchs.pt() );
        corrAK4PFchs_Unc = corrAK4PFchs * (1 + fabs(temp_JECUncDATA_->getUncertainty(1)));
      }
      else if (is_updown_.find("_DOWN_") != std::string::npos) {
        temp_JECUncDATA_->setJetEta( uncorrJetAK4PFchs.eta() );
        temp_JECUncDATA_->setJetPt( corrAK4PFchs * uncorrJetAK4PFchs.pt() );
        corrAK4PFchs_Unc = corrAK4PFchs * ( 1 - fabs(temp_JECUncDATA_->getUncertainty(-1)) );
      }
    }
  }

  //Depending on whether the nominal or one of the systematics is running, select which weight (corrAK4PFchs_XXXXX) to use.
  float JERScaleFactor     = 1;
  float JERScaleFactorUP   = 1;
  float JERScaleFactorDOWN = 1;
  if(!_is_data) GetJER(j, corrAK4PFchs, rhoJER, true, JERScaleFactor, JERScaleFactorUP, JERScaleFactorDOWN, iSetup);


  //Acceptance
  if(sysUncSource_=="Nominal"){
    jetpt = (j.correctedJet("Uncorrected").pt()*corrAK4PFchs*JERScaleFactor);
  }
  else {
    jetpt = (j.correctedJet("Uncorrected").pt()*corrAK4PFchs_Unc*JERScaleFactor);
  }
  //std::cout << "jetpt = " << jetpt << std::endl;
  if(jetpt < 30)       isgoodjet = false; //PLEASE NOTE: that this requirement is for the SL channel, while for DL channel we require pT > 20!
  if(fabs(j.eta())>2.4) isgoodjet = false;


  //ID requirements
  if(j.neutralHadronEnergyFraction() >= 0.99) isgoodjet = false;
  if(j.chargedEmEnergyFraction()     >= 0.99) isgoodjet = false;
  if(j.neutralEmEnergyFraction()     >= 0.99) isgoodjet = false;
  if(j.numberOfDaughters()           <= 1)    isgoodjet = false;
  if(j.chargedHadronEnergyFraction() <= 0.0)  isgoodjet = false;
  if(j.chargedMultiplicity()         <= 0.0)  isgoodjet = false;

  //cout<<setw(20)<<"Jet pt,eta,phi"<<setw(20)<<jetpt<<setw(20)<<j.eta()<<setw(20)<<j.phi()<<endl;

  return isgoodjet;
}
void BJetness::JECInitialization(std::vector<boost::shared_ptr<JetCorrectionUncertainty> > &factorisedJECsMC, std::vector<boost::shared_ptr<JetCorrectionUncertainty> > &factorisedJECsDATA) {
  // Procedure for running JEC and JECUnc on crab: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite

  //AK4chs - MC: Get the factorized jet corrector parameters.
  std::vector<std::string> jecPayloadNamesAK4PFchsMC_;
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC1_.fullPath());
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC2_.fullPath());
  jecPayloadNamesAK4PFchsMC_.push_back(jecPayloadNamesAK4PFchsMC3_.fullPath());
  std::vector<JetCorrectorParameters> vParAK4PFchsMC;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFchsMC_.begin(),
  payloadEnd = jecPayloadNamesAK4PFchsMC_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFchsMC.push_back(pars);
  }
  jecAK4PFchsMC_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFchsMC) );

  jecAK4PFchsMCUnc_   = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "Total"))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "Total")))));
  //=============================================================================
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "AbsoluteStat")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "AbsoluteScale")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "AbsoluteFlavMap")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "AbsoluteMPFBias")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "Fragmentation")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "SinglePionECAL")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "SinglePionHCAL")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "FlavorQCD")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "TimePtEta")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativeJEREC1")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativeJEREC2")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativeJERHF")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativePtBB")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativePtEC1")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativePtEC2")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativePtHF")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativeBal")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativeFSR")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativeStatFSR")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativeStatEC")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "RelativeStatHF")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "PileUpDataMC")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "PileUpPtRef")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "PileUpPtBB")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "PileUpPtEC1")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "PileUpPtEC2")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "PileUpPtHF")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "PileUpMuZero")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "PileUpEnvelope")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "SubTotalPileUp")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "SubTotalRelative")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "SubTotalPt")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "SubTotalScale")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "SubTotalAbsolute")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "SubTotalMC")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "TotalNoFlavor")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "TotalNoFlavorNoTime")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "FlavorZJet")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "FlavorPhotonJet")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "FlavorPureGluon")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "FlavorPureQuark")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "FlavorPureCharm")))));
  factorisedJECsMC.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsMCUnc_.fullPath(), "FlavorPureBottom")))));
  //=============================================================================

  //AK4chs - DATA: Get the factorized jet corrector parameters.
  std::vector<std::string> jecPayloadNamesAK4PFchsDATA_;
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA1_.fullPath());
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA2_.fullPath());
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA3_.fullPath());
  jecPayloadNamesAK4PFchsDATA_.push_back(jecPayloadNamesAK4PFchsDATA4_.fullPath());
  std::vector<JetCorrectorParameters> vParAK4PFchsDATA;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNamesAK4PFchsDATA_.begin(),
	  payloadEnd = jecPayloadNamesAK4PFchsDATA_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vParAK4PFchsDATA.push_back(pars);
  }

  jecAK4PFchsDATA_    = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4PFchsDATA) );
  jecAK4PFchsDATAUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "Total"))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "Total")))));
  //=============================================================================
  //jecAK4PFchsDATAUnc_AbsoluteStat   = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "AbsoluteStat"))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "AbsoluteStat")))));
  //=============================================================================
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "AbsoluteScale")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "AbsoluteFlavMap")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "AbsoluteMPFBias")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "Fragmentation")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "SinglePionECAL")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "SinglePionHCAL")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "FlavorQCD")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "TimePtEta")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativeJEREC1")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativeJEREC2")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativeJERHF")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativePtBB")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativePtEC1")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativePtEC2")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativePtHF")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativeBal")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativeFSR")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativeStatFSR")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativeStatEC")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "RelativeStatHF")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "PileUpDataMC")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "PileUpPtRef")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "PileUpPtBB")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "PileUpPtEC1")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "PileUpPtEC2")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "PileUpPtHF")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "PileUpMuZero")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "PileUpEnvelope")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "SubTotalPileUp")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "SubTotalRelative")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "SubTotalPt")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "SubTotalScale")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "SubTotalAbsolute")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "SubTotalMC")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "TotalNoFlavor")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "TotalNoFlavorNoTime")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "FlavorZJet")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "FlavorPhotonJet")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "FlavorPureGluon")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "FlavorPureQuark")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "FlavorPureCharm")))));
  factorisedJECsDATA.push_back(boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecPayloadNamesAK4PFchsDATAUnc_.fullPath(), "FlavorPureBottom")))));
  //=============================================================================

}
void BJetness::GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN, const edm::EventSetup& iSetup){
  if(!jet.genJet()) return;
  double jetEta=fabs(jet.eta());
  double cFactorJER = 1.0;
  double cFactorJERdown = 1.0;
  double cFactorJERup = 1.0;
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Unce_AN1
  if( jetEta<0.5 ){
    cFactorJER = 1.109;
    cFactorJERdown = 1.109-0.008;
    cFactorJERup   = 1.109+0.008;
  } else if( jetEta<0.8 ){
    cFactorJER = 1.138;
    cFactorJERdown = 1.138-0.013;
    cFactorJERup   = 1.138+0.013;
  } else if( jetEta<1.1 ){
    cFactorJER = 1.114;
    cFactorJERdown = 1.114-0.013;
    cFactorJERup   = 1.114+0.013;
  } else if( jetEta<1.3 ){
    cFactorJER = 1.123;
    cFactorJERdown = 1.123-0.024;
    cFactorJERup   = 1.123+0.024;
  } else if( jetEta<1.7 ){
    cFactorJER = 1.084;
    cFactorJERdown = 1.084-0.011;
    cFactorJERup   = 1.084+0.011;
  } else if( jetEta<1.9 ){
    cFactorJER = 1.082;
    cFactorJERdown = 1.082-0.035;
    cFactorJERup   = 1.082+0.035;
  } else if( jetEta<2.1 ){
    cFactorJER = 1.140;
    cFactorJERdown = 1.140-0.047;
    cFactorJERup   = 1.140+0.047;
  } else if( jetEta<2.3 ){
    cFactorJER = 1.067;
    cFactorJERdown = 1.067-0.053;
    cFactorJERup   = 1.067+0.053;
  } else if( jetEta<2.5 ){
    cFactorJER = 1.177;
    cFactorJERdown = 1.177-0.041;
    cFactorJERup   = 1.177+0.041;
  } else if( jetEta<2.8 ){
    cFactorJER = 1.364;
    cFactorJERdown = 1.364-0.039;
    cFactorJERup   = 1.364+0.039;
  } else if( jetEta<3.0 ){
    cFactorJER = 1.857;
    cFactorJERdown = 1.857-0.071;
    cFactorJERup   = 1.857+0.071;
  } else if( jetEta<3.2 ){
    cFactorJER = 0.328;
    cFactorJERdown = 0.328-0.022;
    cFactorJERup   = 0.328+0.022;
  } else if( jetEta<5.0 ){
    cFactorJER = 1.160;
    cFactorJERdown = 1.160-0.029;
    cFactorJERup   = 1.160+0.029;
  }

  double recoJetPt = (jet.correctedJet("Uncorrected").pt())*JesSF;
  double genJetPt  = jet.genJet()->pt();
  double diffPt    = recoJetPt - genJetPt;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor res_sf;
  if(AK4PFchs){
    resolution = JME::JetResolution(jerAK4PFchs_);
    res_sf = JME::JetResolutionScaleFactor(jerAK4PFchsSF_);
  } else {
    std::cout << "BJetness.cc: WARNING - unsupported jet collection used. Currently only supported for AK4PFchs" << std::endl;
  }
  JME::JetParameters parameters;
  parameters.setJetPt(jet.pt());
  parameters.setJetEta(jet.eta());
  parameters.setRho(rhoJER);
  float relpterr = resolution.getResolution(parameters);
  if(genJetPt>0. && deltaR(jet.eta(),jet.phi(),jet.genJet()->eta(),jet.genJet()->phi())<0.2
     && (abs(jet.pt()-jet.genJet()->pt())<3*relpterr*jet.pt())) {
    JERScaleFactor     = (std::max(0., genJetPt + cFactorJER*diffPt))/recoJetPt;
    JERScaleFactorUP   = (std::max(0., genJetPt + cFactorJERup*diffPt))/recoJetPt;
    JERScaleFactorDOWN = (std::max(0., genJetPt + cFactorJERdown*diffPt))/recoJetPt;
  } else {
    JERScaleFactor     = 1.;
    JERScaleFactorUP   = 1.;
    JERScaleFactorDOWN = 1.;
  }
}
/////
//   Methods for the BJetness variables
/////
void BJetness::get_bjetness_vars(vector<pat::Jet> evtjets, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h,
                                           double& bjetnessFV_num_leps, double& bjetnessFV_npvTrkOVcollTrk, double& bjetnessFV_pvTrkOVcollTrk, double& bjetnessFV_avip3d_val, double& bjetnessFV_avip3d_sig, double& bjetnessFV_avsip3d_val, double& bjetnessFV_avsip3d_sig, double& bjetnessFV_avip2d_val, double& bjetnessFV_avip2d_sig, double& bjetnessFV_avsip2d_val, double& bjetnessFV_avsip2d_sig, double& bjetnessFV_avip1d_val, double& bjetnessFV_avip1d_sig, double& bjetnessFV_avsip1d_val, double& bjetnessFV_avsip1d_sig, double &bjetnessFV_npvTrkOVpvTrk , double &bjetnessFV_npvPtOVcollPt , double &bjetnessFV_pvPtOVcollPt , double &bjetnessFV_npvPtOVpvPt
                                          ){
  vector<Track> jetschtrks; jetschtrks.clear();
  vector<Track> jetschtrkspv;
  vector<Track> jetschtrksnpv;
  double num_pvtrks  = 0;
  double num_npvtrks = 0;
  double num_eles    = 0;
  double num_mus     = 0;
  vector<tuple<double, double, double> > jetsdir; jetsdir.clear();

  // Select jet tracks (PV and non-PV associated) + count leptons.
  get_bjetness_trkinfos(evtjets, vtx, jetschtrks, num_pvtrks, num_npvtrks, electron_pat, muon_h, num_eles, num_mus, jetsdir, jetschtrkspv, jetschtrksnpv);

  bjetnessFV_num_leps = num_eles+num_mus;


  // Construct PV & NPV fraction variables.
  if(jetschtrks.size()!=0){
    bjetnessFV_npvTrkOVcollTrk       = num_npvtrks/double(jetschtrks.size());
    bjetnessFV_pvTrkOVcollTrk       = num_pvtrks/double(jetschtrks.size());
    bjetnessFV_npvTrkOVpvTrk        = num_npvtrks/num_pvtrks;

    double ptjettrks    = 0;
    double ptjettrkspv    = 0;
    double ptjettrksnpv    = 0;
    for(uint jt=0; jt<jetschtrks.size(); jt++) ptjettrks += jetschtrks[jt].pt();
    for(uint jt=0; jt<num_pvtrks; jt++) ptjettrkspv += jetschtrkspv[jt].pt();
    for(uint jt=0; jt<num_npvtrks; jt++) ptjettrksnpv += jetschtrksnpv[jt].pt();

    if(ptjettrks!=0) bjetnessFV_npvPtOVcollPt = ptjettrksnpv/ptjettrks;
    else             bjetnessFV_npvPtOVcollPt = -997;
    if(ptjettrks!=0) bjetnessFV_pvPtOVcollPt = ptjettrkspv/ptjettrks;
    else             bjetnessFV_pvPtOVcollPt = -997;
    if(ptjettrkspv!=0) bjetnessFV_npvPtOVpvPt = ptjettrksnpv/ptjettrkspv;
    else               bjetnessFV_npvPtOVpvPt = -997;


    //================= IMPACT PARAMETERS ================
    //================= 3D Variables =================
    double ip_valtemp = 0;
    double jetchtrks_avip3d_val  = 0;
    double jetchtrks_avip3d_sig  = 0;
    double jetchtrks_avsip3d_val = 0;
    double jetchtrks_avsip3d_sig = 0;
    get_avip3d(jetschtrks, ttrkbuilder, vtx, jetsdir, jetchtrks_avip3d_val,jetchtrks_avip3d_sig,jetchtrks_avsip3d_val,jetchtrks_avsip3d_sig);
    ip_valtemp = jetchtrks_avip3d_val/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip3d_val = ip_valtemp;
    else                       bjetnessFV_avip3d_val = -996;
    ip_valtemp = jetchtrks_avip3d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip3d_sig = ip_valtemp;
    else                       bjetnessFV_avip3d_sig = -996;
    ip_valtemp = jetchtrks_avsip3d_val/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avsip3d_val = ip_valtemp;
    else                       bjetnessFV_avsip3d_val = -996;
    ip_valtemp = jetchtrks_avsip3d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avsip3d_sig = ip_valtemp;
    else                       bjetnessFV_avsip3d_sig = -996;

    //================= 2D Variables =================
    double jetchtrks_avip2d_val  = 0;
    double jetchtrks_avip2d_sig  = 0;
    double jetchtrks_avsip2d_val = 0;
    double jetchtrks_avsip2d_sig = 0;
    get_avip2d(jetschtrks, ttrkbuilder, vtx, jetsdir, jetchtrks_avip2d_val,jetchtrks_avip2d_sig,jetchtrks_avsip2d_val,jetchtrks_avsip2d_sig);
    ip_valtemp = jetchtrks_avip2d_val/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip2d_val = ip_valtemp;
    else                       bjetnessFV_avip2d_val = -996;
    ip_valtemp = jetchtrks_avip2d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip2d_sig = ip_valtemp;
    else                       bjetnessFV_avip2d_sig = -996;
    ip_valtemp = jetchtrks_avsip2d_val/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avsip2d_val = ip_valtemp;
    else                       bjetnessFV_avsip2d_val = -996;
    ip_valtemp = jetchtrks_avsip2d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avsip2d_sig = ip_valtemp;
    else                       bjetnessFV_avsip2d_sig = -996;

    //================= 1D Variables =================
    double jetchtrks_avip1d_val  = 0;
    double jetchtrks_avip1d_sig  = 0;
    double jetchtrks_avsip1d_val  = 0;
    double jetchtrks_avsip1d_sig  = 0;
    get_avip1d(jetschtrks, ttrkbuilder, vtx, jetsdir, jetchtrks_avip1d_val, jetchtrks_avip1d_sig, jetchtrks_avsip1d_val, jetchtrks_avsip1d_sig);
    ip_valtemp = jetchtrks_avip1d_val/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip1d_val = ip_valtemp;
    else                       bjetnessFV_avip1d_val = -996;
    ip_valtemp = jetchtrks_avip1d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip1d_sig = ip_valtemp;
    else                       bjetnessFV_avip1d_sig = -996;
    ip_valtemp = jetchtrks_avsip1d_val/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avsip1d_val = ip_valtemp;
    else                       bjetnessFV_avsip1d_val = -996;
    ip_valtemp = jetchtrks_avsip1d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avsip1d_sig = ip_valtemp;
    else                       bjetnessFV_avsip1d_sig = -996;

  }
  else{
    bjetnessFV_npvTrkOVcollTrk       = -998;
    bjetnessFV_pvTrkOVcollTrk       = -998;
    bjetnessFV_npvTrkOVpvTrk =-998;
    bjetnessFV_npvPtOVcollPt =-998;
    bjetnessFV_pvPtOVcollPt = -998;
    bjetnessFV_npvPtOVpvPt =-998;
    bjetnessFV_avip3d_val            = -998;
    bjetnessFV_avip3d_sig            = -998;
    bjetnessFV_avsip3d_val           = -998;
    bjetnessFV_avsip3d_sig           = -998;
    bjetnessFV_avip1d_val            = -998;
    bjetnessFV_avip1d_sig            = -998;
    bjetnessFV_avsip1d_val            = -998;
    bjetnessFV_avsip1d_sig            = -998;
  }
}





void BJetness::get_bjetness_trkinfos(vector<pat::Jet> evtjets, const reco::Vertex& vtx, vector<Track>& jetchtrks, double& bjetness_num_pvtrks, double& bjetness_num_npvtrks, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h, double& bjetness_num_eles, double& bjetness_num_mus, vector<tuple<double, double, double> >& jetsdir, vector<Track>& jetchtrkspv, vector<Track>& jetchtrksnpv){
  //Loop over evt jet
  for(uint j=0; j<evtjets.size(); j++){
    pat::Jet jet = evtjets[j];

    //Access jet daughters
    vector<CandidatePtr> jdaus(jet.daughterPtrVector());
    sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});

    for(uint jd=0; jd<jdaus.size(); jd++){
      const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
      //dR requirement
      if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
      Track trk = Track(jcand.pseudoTrack());
      bool isgoodtrk = is_goodtrk(trk,vtx);


      //Minimal conditions for a BJetness jet constituent
      if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1){
        jetchtrks.push_back(trk);
        if(jcand.fromPV()==pat::PackedCandidate::PVUsedInFit){
          jetchtrkspv.push_back(trk);
        }else{
          jetchtrksnpv.push_back(trk);
        }
        //if(jcand.fromPV()==3) bjetness_num_pvtrks++;
        //if(jcand.fromPV()==2) bjetness_num_npvtrks++;
        if(jcand.fromPV()==pat::PackedCandidate::PVUsedInFit) bjetness_num_pvtrks++;
        else bjetness_num_npvtrks++;
        jetsdir.push_back(make_tuple(jet.px(),jet.py(),jet.pz()));
        if(fabs(jcand.pdgId())==13 && is_loosePOG_jetmuon(jcand,muon_h)) bjetness_num_mus++;
        if(fabs(jcand.pdgId())==11 && is_softLep_jetelectron(jcand,electron_pat,vtx)) bjetness_num_eles++;
      }//Ch trks
    }//Loop on jet daus
  }//Loop on evt jet
}



//Check that the track is a good track
bool BJetness::is_goodtrk(Track trk,const reco::Vertex& vtx){
 bool isgoodtrk = false;
 if(trk.pt()>1 &&
   trk.hitPattern().numberOfValidHits()>=8 &&
   trk.hitPattern().numberOfValidPixelHits()>=2 &&
   trk.normalizedChi2()<5 &&
   std::abs(trk.dxy(vtx.position()))<0.2 &&
   std::abs(trk.dz(vtx.position()))<17
   ) isgoodtrk = true;
 return isgoodtrk;
}
//Loose muons (look for candidates among jet daughters)
bool BJetness::is_loosePOG_jetmuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h){
  bool ismu = false;
  for(const pat::Muon &mu : *muon_h){
    if(deltaR(jcand.p4(),mu.p4())<0.1 && fabs(jcand.pt()-mu.pt())/mu.pt()<0.05){
     if(mu.isLooseMuon()) ismu = true;
     if(ismu) break;
    }
  }
  return ismu;
}
//Loose electrons (look for candidates among jet daughters)
bool BJetness::is_softLep_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele;
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05 ){
      const HitPattern &hitPattern = lele.gsfTrack().get()->hitPattern();
      uint32_t hit = hitPattern.getHitPattern(HitPattern::TRACK_HITS, 0);
      bool hitCondition = !(HitPattern::validHitFilter(hit) && ((HitPattern::pixelBarrelHitFilter(hit) && HitPattern::getLayer(hit) < 3) || HitPattern::pixelEndcapHitFilter(hit)));
      if(!hitCondition && lele.passConversionVeto()) isele = true;
      if(isele) break;
    }
  }
  return isele;
}

bool BJetness::is_vetoPOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat){
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele;
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05){

      double ooEmooP = 999;
      if(lele.ecalEnergy()==0)                   ooEmooP = 1e30;
      else if(!std::isfinite(lele.ecalEnergy())) ooEmooP = 1e30;
      else                                       ooEmooP = fabs(1.0/lele.ecalEnergy() - lele.eSuperClusterOverP()/lele.ecalEnergy());

      if(!(fabs(lele.superCluster()->position().eta()) > 1.4442 && fabs(lele.superCluster()->position().eta()) < 1.5660)){
        if(fabs(lele.superCluster()->position().eta())<=1.4442
        && lele.full5x5_sigmaIetaIeta()<0.0115
        && fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.00749
        && fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.228
        && lele.hcalOverEcal()<0.356
        && ooEmooP<0.299
        && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=2
        && lele.passConversionVeto())
        {isele = true;}
        if(fabs(lele.superCluster()->position().eta())>1.5660
        && (lele.full5x5_sigmaIetaIeta()<0.037)
        && (lele.hcalOverEcal()<0.211)
        && (fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.00895)
        && (fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.213)
        && ooEmooP<0.15 && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=3
        && lele.passConversionVeto() )
        {isele = true;}
      }
      if(isele) break;
    }
  }
  return isele;
}

bool BJetness::is_loosePOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele;
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05){
      double ooEmooP = 999;
      if(lele.ecalEnergy()==0)                   ooEmooP = 1e30;
      else if(!std::isfinite(lele.ecalEnergy())) ooEmooP = 1e30;
      else                                       ooEmooP = fabs(1.0/lele.ecalEnergy() - lele.eSuperClusterOverP()/lele.ecalEnergy());
      if(!(fabs(lele.superCluster()->position().eta()) > 1.4442 && fabs(lele.superCluster()->position().eta()) < 1.5660)){
        if(fabs(lele.superCluster()->position().eta())<=1.4442
          && (lele.full5x5_sigmaIetaIeta()<0.011)
          && (lele.hcalOverEcal()<0.298)
          && (fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.00477)
          && (fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.222)
          && ooEmooP<0.241
          && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=1
          && lele.passConversionVeto()
          ){
            isele = true;
        }
        if(fabs(lele.superCluster()->position().eta())>1.5660
          && (lele.full5x5_sigmaIetaIeta()<0.0314)
          && (lele.hcalOverEcal()<0.101)
          && (fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.00868)
          && (fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.213)
          && ooEmooP<0.14
          && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=1
          && lele.passConversionVeto()
          ){
            isele = true;
        }
      }
      if(isele) break;
    }
  }
  return isele;
}
//Methods related to IP
void BJetness::get_avip3d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip3d_val, double& jetchtrks_avip3d_sig, double& jetchtrks_avsip3d_val, double& jetchtrks_avsip3d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avip3d_val  += valtemp;
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avip3d_sig  += valtemp;
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],jetsdirgv,vtx).second.value();
    if(valtemp==valtemp){jetchtrks_avsip3d_val += valtemp;}
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],jetsdirgv,vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avsip3d_sig += valtemp;
  }
}
void BJetness::get_avip2d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip2d_val, double& jetchtrks_avip2d_sig, double& jetchtrks_avsip2d_val, double& jetchtrks_avsip2d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteTransverseImpactParameter(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avip2d_val  += valtemp;
    valtemp = IPTools::absoluteTransverseImpactParameter(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avip2d_sig  += valtemp;
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = IPTools::signedTransverseImpactParameter(ttrks[t],jetsdirgv,vtx).second.value();
    if(valtemp==valtemp){jetchtrks_avsip2d_val += valtemp;}
    valtemp = IPTools::signedTransverseImpactParameter(ttrks[t],jetsdirgv,vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avsip2d_sig += valtemp;
  }
}


/*void BJetness::get_avip1d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip1d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  SignedTransverseImpactParameter stip;
  for(uint t=0; t<ttrks.size(); t++){
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = fabs(stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.significance());
    if(valtemp==valtemp) jetchtrks_avip1d_sig  += valtemp;
  }
}*/

void BJetness::get_avip1d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip1d_val, double& jetchtrks_avip1d_sig, double& jetchtrks_avsip1d_val, double& jetchtrks_avsip1d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  SignedTransverseImpactParameter stip;
  for(uint t=0; t<ttrks.size(); t++){
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = fabs(stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.value());
    if(valtemp==valtemp) jetchtrks_avip1d_val  += valtemp;
    valtemp = fabs(stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.significance());
    if(valtemp==valtemp) jetchtrks_avip1d_sig  += valtemp;
    valtemp = stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avsip1d_val += valtemp;
    valtemp = stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avsip1d_sig += valtemp;
  }
}

vector<TransientTrack> BJetness::get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
 vector<TransientTrack> ttrks;
 for(uint tr=0; tr<trks.size(); tr++){
  TransientTrack ttrk = ttrkbuilder.build(&trks[tr]);
  ttrks.push_back(ttrk);
 }
 return ttrks;
}
// ------------ method called once each stream before processing any runs, lumis or events  ------------
void BJetness::beginStream(edm::StreamID) {
}
// ------------ method called once each stream after processing all runs, lumis and events  ------------
void BJetness::endStream() {
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BJetness::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ===== Declare producer to PluginManager =====
// User should be able to compiled analysis code into same library as plugin.
// All such plugins must have a "DEFINE_FWK_MODULE" declaration to
// the PluginManager.
DEFINE_FWK_MODULE(BJetness);
