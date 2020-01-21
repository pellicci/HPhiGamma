
//---------- class declaration----------

class HPhiGammaAnalysis : public edm::EDAnalyzer {
public:
  explicit HPhiGammaAnalysis(const edm::ParameterSet&);
  ~HPhiGammaAnalysis();

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;

  bool runningOnData_;
  const edm::InputTag packedPFCandidates_;
  const edm::InputTag slimmedMuons_; 
  const edm::InputTag prunedGenParticles_;
  const edm::InputTag slimmedPhotons_;
  const edm::InputTag slimmedElectrons_;
  const edm::InputTag slimmedJets_;
  const edm::InputTag slimmedMETs_;
  const edm::InputTag slimmedMETsPuppi_;
  const edm::InputTag pvCollection_;  
  const edm::InputTag bsCollection_;  
  const edm::InputTag PileupSrc_;
  const edm::InputTag GenInfo_;

  edm::LumiReWeighting Lumiweights_;

  edm::Service<TFileService> fs;

  void create_trees();

  // ----------member data ---------------------------
  TH1F* h_Events;

  TH1F* h_pileup;

  //Counters
  int nPV;

  int nMuons;
  int nElectrons;
  int nPhotons;
  int nJets;
  int nJets_25;

  int _Nevents_processed;
  int _Nevents_isTwoKaons;
  int _Nevents_isPhoton;
  int _Nevents_HiggsFound; 
  int _Nevents_HiggsNotMatched; 
  int _Nevents_HiggsMassMatched; 
  int _Nevents_HiggsMassNotMatched; 
  int _Nevents_JetIsBetterForMass;

  //TTree and TTree variables
  TTree *mytree;

  int run_number;

  float K1_pT;
  float K1_eta;
  float K1_phi;
  float K1_energy;
  float K1_dxy;
  float K1_dz;
  float K1_charge;
  float K1_sum_pT_03;
  float K1_sum_pT_05;
  float K1_sum_pT_05_ch;
  float K1_pTMax;

  float K2_pT;
  float K2_eta;
  float K2_phi;
  float K2_energy;
  float K2_dxy;
  float K2_dz;
  float K2_charge;
  float K2_sum_pT_03;
  float K2_sum_pT_05;
  float K2_sum_pT_05_ch;
  float K2_pTMax;

  float ph_eT;
  float ph_eta;
  float ph_etaSC;
  float ph_phi;
  float ph_energy;
  float ph_iso_ChargedHadron;
  float ph_iso_NeutralHadron;
  float ph_iso_Photon;
  float ph_iso_eArho;

  float eTphMax;

  float firstCandPx;
  float firstCandPy;
  float firstCandPz;
  float secondCandPx;
  float secondCandPy;
  float secondCandPz;
  float firstCandEnergy;
  float secondCandEnergy;
  float _Jet_Photon_invMass;
  float _Phimass;
  float _Hmass_From2K_Photon;
    
  float met_pT;
  float metpuppi_pT;
  
  bool isTwoProngTrigger;

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
  int _bestJet_chargedHadMultiplicity;
  float _bestJet_invMass;
  float _bestJet_Photon_invMass;
  
  //MC truth
  float PU_Weight;
  float MC_Weight;

  bool is_Kplus_fromPhi;
  bool is_Kminus_fromPhi;
  bool is_Phi_fromH;
  bool is_Photon_fromH;

  bool is_photon_a_photon;
  bool is_photon_matched;
  bool _isHiggsFound;

  //rho for isolation
  float rho_;

  //Tokens
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedPFCandidatesToken_; 
  edm::EDGetTokenT<std::vector<pat::Muon> > slimmedMuonsToken_; 
  edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenParticlesToken_; 
  edm::EDGetTokenT<std::vector<pat::Jet> > slimmedJetsToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsPuppiToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > offlineSlimmedPrimaryVerticesToken_; 
  edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpotToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT<GenEventInfoProduct> GenInfoToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;

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
