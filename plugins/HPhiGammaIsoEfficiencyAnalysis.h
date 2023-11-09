
//---------- class declaration----------

class HPhiGammaIsoEfficiencyAnalysis : public edm::EDAnalyzer {
public:
  explicit HPhiGammaIsoEfficiencyAnalysis(const edm::ParameterSet&);
  ~HPhiGammaIsoEfficiencyAnalysis();

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
  float tagMuEta;
  float probeMuEta;
  float tagMuPhi;
  float probeMuPhi;
  float tagMuPt;
  float probeMuPt;
  
  float sum_pT_05;
  float sum_pT_05_ch;

  bool isIsoMuTrigger;
  bool isBestMuMuFound;
  bool _hasProbeMuFiredTrigger;
  bool _hasTagMuFiredTrigger;

  float _iso;
  float _iso_ch;


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
