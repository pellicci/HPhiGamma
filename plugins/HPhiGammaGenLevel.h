//---------- class declaration----------

class HPhiGammaGenLevel : public edm::EDAnalyzer {

 public:
  explicit HPhiGammaGenLevel(const edm::ParameterSet&);
  ~HPhiGammaGenLevel();

 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;

  const edm::InputTag genParticles_;

  edm::Service<TFileService> fs;

  void create_trees();

  // ---------- member data ----------- //

  TTree *mytree;

  int event_number;
  int genH_ID_tree;
  float genH_pT_tree;
  float genH_eta_tree;
  float genH_phi_tree;
  float genH_E_tree;

  int genGamma_ID_tree;
  float genGamma_pT_tree;
  float genGamma_eta_tree;
  float genGamma_phi_tree;
  float genGamma_E_tree;

  int genPhi_ID_tree;
  float genPhi_pT_tree;
  float genPhi_eta_tree;
  float genPhi_phi_tree;
  float genPhi_E_tree;

  int genKminus_ID_tree;
  float genKminus_pT_tree;
  float genKminus_eta_tree;
  float genKminus_phi_tree;
  float genKminus_E_tree;

  int genKplus_ID_tree;
  float genKplus_pT_tree;
  float genKplus_eta_tree;
  float genKplus_phi_tree;
  float genKplus_E_tree;

  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesToken_; 
  //edm::EDGetTokenT<GenEventInfoProduct> GenInfoToken_;
};
