
//---------- class declaration----------

class HPhiGammaTwoProngsTriggerAnalysis : public edm::EDAnalyzer {
public:
  explicit HPhiGammaTwoProngsTriggerAnalysis(const edm::ParameterSet&);
  ~HPhiGammaTwoProngsTriggerAnalysis();

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

  //debug

  bool debug;
  bool verbose;

  //Counters
  int nPV;
  int nJets;
  int nJets_25;
  int nPhotonsChosen;
  int nPhotonsOverSelection;
  int _Nevents_processed;
  int _nEvents_ZmumuFound;
  int _Nevents_isPhoton;
  
  //TTree and TTree variables
  TTree *mytree;

  int run_number;

  float currentMuMuPt;
  float bestMuMuPt;
  float bestMuMuMass;
  float firstMuEta;
  float secondMuEta;
  float firstMuPhi;
  float secondMuPhi;

  float firstCandPx;
  float firstCandPy;
  float firstCandPz;
  float secondCandPx;
  float secondCandPy;
  float secondCandPz;
  float firstCandEnergy;
  float secondCandEnergy;
  float firstCandEnergy_K;
  float secondCandEnergy_K;
  float firstCandEnergy_Pi;
  float secondCandEnergy_Pi;
  float _firstCandPt;
  float _firstCandEta;
  float _firstCandPhi;
  float _secondCandPt;
  float _secondCandEta;
  float _secondCandPhi;
  float _bestCouplePt;
  float _bestCoupleEta;
  float _bestCouplePhi;

  float _Jet_Photon_invMass;
  float _MesonMass;

  
  bool isIsoMuTrigger;
  bool isTwoProngsTrigger;
  bool isBestMuMu_Found;
  bool isBestCoupleOfTheEvent_Found;
  bool _isPhi;
  bool _isRho;

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

  //MC truth
  float PU_Weight;

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
