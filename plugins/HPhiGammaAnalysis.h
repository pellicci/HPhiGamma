
//---------- class declaration----------

class HPhiGammaAnalysis : public edm::EDAnalyzer {
public:
  explicit HPhiGammaAnalysis(const edm::ParameterSet&);
  ~HPhiGammaAnalysis();

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;
 
 //test method
  //  int countPrimaryVertex(edm::Handle<std::vector<reco::Vertex > > slimmedPV);

  bool runningOnData_;
  const edm::InputTag packedPFCandidates_;
  const edm::InputTag slimmedMuons_; 
  const edm::InputTag prunedGenParticles_;
  //const edm::InputTag packedGenParticles_;
  const edm::InputTag genParticles_;
  const edm::InputTag slimmedPhotons_;
  const edm::InputTag slimmedElectrons_;
  const edm::InputTag slimmedJets_;
  const edm::InputTag slimmedMETs_;
  const edm::InputTag slimmedMETsPuppi_;
  const edm::InputTag pvCollection_;  
  const edm::InputTag bsCollection_;  
  const edm::InputTag PileupSrc_;
  const edm::InputTag GenInfo_;
  //const edm::InputTag LHEEventProduct_; //LHE reader

  edm::LumiReWeighting Lumiweights_;

  edm::Service<TFileService> fs;

  void create_trees();

  // ----------member data ---------------------------
  TH1F* h_Events;
  TH1F* h_TriggerFilters;

  TH1F* h_pileup;

  //debug

  bool debug;
  bool verbose;
  bool trig_verbose;
  bool isTriggerStudies;

  //Counters
  int nPV;

  int nMuons10;
  int nMuons20;
  int nElectrons10;
  int nElectrons20;
  int nPhotonsChosen;
  int nPhotonsWP90_pT25_endcap;
  int nPhotonsWP90_pT15_barrel;
  int nPhotonsWP90_pT20_2p5eta3p0;
  int nPhotons38WP80;
  int nJets30;
  int nJets25;
  int _Nevents_triggered;
  int _Nevents_processed;
  int _Nevents_isTwoKaons;
  int _Nevents_isPhoton;
  int _Nevents_HiggsMatched; 
  int _Nevents_HiggsNotMatched; 
  int _Nevents_HiggsMassMatched; 
  int _Nevents_HiggsMassNotMatched; 
  int _Nevents_MesonPtNotMatched;
  int _Nevents_bestCoupleFound;
  int _Nevents_candPtFilter;
  int _Nevents_coupleIsolationFilter;
  int _Nevents_VBFVeto;

  int _NeventsFilter1;
  int _NeventsFilter2;
  int _NeventsFilter3;
  int _NeventsFilter4;
  int _NeventsFilter5;
  int _NeventsFilter6;
  int _NeventsFilter7;
  int _NeventsFilter8;
  int _NeventsFilter9;

  int nKK_found;
  int nKK_notFound;
  
  //TTree and TTree variables
  TTree *mytree;

  int run_number;
  int event_number;

  float ph_eT;
  float ph_en_sigmaUP;
  float ph_en_sigmaDW;
  float ph_en_scaleUP;
  float ph_en_scaleDW;
  float ph_eta;
  float ph_etaSC;
  float ph_phi;
  float ph_iso_ChargedHadron;
  float ph_iso_NeutralHadron;
  float ph_iso_Photon;
  float ph_iso_eArho;
  bool is_photon_wp90;
  float eTphMax;
  float photonRegressionError;
  bool isPhotonTriggerMatched;

  float firstCandPx;
  float firstCandPy;
  float firstCandPz;
  float bestFirstCand_dxy;
  float bestFirstCand_dxyErr;
  float bestFirstCand_dz;
  float bestFirstCand_dzErr;
  float secondCandPx;
  float secondCandPy;
  float secondCandPz;
  float bestSecondCand_dxy;
  float bestSecondCand_dxyErr;
  float bestSecondCand_dz;
  float bestSecondCand_dzErr;
  float firstCandEnergy;
  float secondCandEnergy;
  float firstCandEnergy_K;
  float secondCandEnergy_K;
  float firstCandEnergy_Pi;
  float secondCandEnergy_Pi;
  float _firstCandPt;
  float _firstCandEta;
  float _firstCandPhi;
  float _firstCandCharge;
  float _secondCandPt;
  float _secondCandEta;
  float _secondCandPhi;
  float _secondCandCharge;
  float _bestCouplePt;
  float _bestCoupleEta;
  float _bestCouplePhi;
  float minPDFWeight;
  float maxPDFWeight;
  float minQCDWeight;
  float maxQCDWeight;

  float _Jet_Photon_invMass;
  //float _MesonMass;
  //float _Hmass_From2K_Photon;
  float _PhiMass      ;
  float _RhoMass      ;
  float _K0starMass   ;
  float _PhiGammaMass ;
  float _RhoGammaMass;
  float _K0starGammaMass;

  float K1_sum_pT_05;
  float K1_sum_pT_05_ch;
  float K2_sum_pT_05;
  float K2_sum_pT_05_ch;
  float couple_sum_pT_05;
  float couple_sum_pT_05_ch;
  float _iso_K1;
  float _iso_K1_ch;
  float _iso_K2;
  float _iso_K2_ch;
  float _iso_couple;
  float _iso_couple_ch;

  
  float met_pT;
  float metpuppi_pT;
  
  bool isTwoProngTrigger;

  //trigger studies
  bool is_hltEG35R9Id90HE10IsoMEcalIsoFilter;
  bool is_hltEG35R9Id90HE10IsoMHcalIsoFilter;
  bool is_hltEG35R9Id90HE10IsoMTrackIsoFilter;
  bool is_hltOverlapFilterPhoton35MediumChargedIsoPFTau35;
  bool is_hltEG35R9Id90HE10IsoMR9Filter;
  bool is_hltEG35R9Id90HE10IsoMHEFilter;
  bool is_hltEG35R9Id90HE10IsoMEtFilter;
  bool is_hltEGL1EGAndTauFilter;
  bool is_hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3;
  bool is_hltL1sSingleMu22;
  bool is_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23;
  bool isIsoMuTrigger;
  bool isBestMuMu_Found;
  float bestMuMuPt;
  float bestMuMuMass;

  //Jet datamember
  
  float _jet_invMass;
  float _bestJet_pT;
  float _bestJet_eta;
  float _bestJet_phi;
  int _bestJet_nDaughters;
  float _bestJet_pTMax;
  float _bestJet_chargedEmEnergy;
  float _bestJet_neutralEmEnergy;
  float _bestJet_chargedHadEnergy;
  float _bestJet_neutralHadEnergy;
  float _bestJet_chargedEmEnergyFraction;
  float _bestJet_neutralEmEnergyFraction;
  float _bestJet_chargedHadEnergyFraction;
  float _bestJet_neutralHadEnergyFraction;
  int _bestJet_chargedHadMultiplicity;
  float _bestJet_invMass;
  float _bestJet_Photon_invMass;
  float _bestJet_JECunc;
  bool isJetTriggerMatched;
  int _bestJet_PUid;

  //MC truth
  float PU_Weight;
  float MC_Weight;
  float theta_pol;
  float deltaR_Kplus;
  float deltaR_wrong;
  float deltaR_Kminus;
  float deltaR_Piplus;
  float deltaR_Piminus;
  float genPhoton_eT;
  float genPhoton_eta;
  float genPhoton_phi;
  float genMeson_pT;
  float genMeson_eta;
  float genMeson_phi;
  float genMeson_m;
  float KplusPt;
  float KminusPt;
  float Kplus_dxy;
  float Kplus_dz;
  float Kminus_dxy;
  float Kminus_dz;
  float Kminus_eta;
  float Kplus_eta;
  float Kminus_phi;
  float Kplus_phi;

  bool is_Kplus_matched;
  bool is_Kminus_matched;
  bool is_Piplus_matched;
  bool is_Piminus_matched;
  bool is_Phi_fromH;
  bool is_Rho_fromH;
  bool is_Photon_fromH;
  bool is_photon_a_photon;
  bool is_photon_matched;
  bool is_meson_matched;
  bool is_Higgs_matched;
  bool _isPhi;
  bool _isRho;
  bool _isK0star;
  bool _isFirstCandK;
  //rho for isolation
  float rho_;

  //for VBF veto
  int nJets20;


  //Tokens
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedPFCandidatesToken_; 
  edm::EDGetTokenT<std::vector<pat::Muon> > slimmedMuonsToken_;   edm::EDGetTokenT<std::vector<pat::Jet> > slimmedJetsToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsPuppiToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > offlineSlimmedPrimaryVerticesToken_; 
  edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpotToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT<GenEventInfoProduct> GenInfoToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
  edm::EDGetTokenT<LHEEventProduct> LHEEventProduct_; //LHE reader
  //edm::EDGetTokenT<std::vector<pat::PackedGenParticle> > packedGenParticlesToken_; //GENPART
  edm::EDGetTokenT<std::vector<reco::GenParticle>> prunedGenParticlesToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjectsToken_;

  //edm::EDGetTokenT<reco::JetCorrector> jetCorrectorToken_;
  //edm::EDGetTokenT<JetCorrectionUncertainty> mJetCorrectorUnc;

//  edm::EDGetTokenT<reco::JetCorrector> jetCorrectorToken_;


  //Ele ID decisions objects
  edm::EDGetToken electronsMiniAODToken_;

  //Photon ID decisions
  edm::EDGetToken photonsMiniAODToken_;

  //rho (PU energy density)
  edm::EDGetTokenT<double> rhoToken_;

  bool verboseIdFlag_;

  //Effective areas for isolation
  EffectiveAreas   effectiveAreas_el_;
  EffectiveAreas   effectiveAreas_ph_;

};
