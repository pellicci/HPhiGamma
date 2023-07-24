//---------- class declaration----------

class HRhoGammaAOD : public edm::EDAnalyzer {

 public:
  explicit HRhoGammaAOD(const edm::ParameterSet&);
  ~HRhoGammaAOD();

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
float genH_m;       
float genMeson_m;   
float genMeson_eta; 
float genMeson_phi; 
float genMeson_pT;  
float genPhoton_eT; 
float genPhoton_eta;
float genPhoton_phi;   
float PiplusPt;     
float Piplus_eta;   
float Piplus_phi;     
float PiminusPt;    
float Piminus_eta;  
float Piminus_phi;
float theta_pol;  

  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesToken_; 
  //edm::EDGetTokenT<GenEventInfoProduct> GenInfoToken_;
};
