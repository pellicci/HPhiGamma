//ROOT includes
#include <TH1F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"
#include <stdlib.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h" 
 
#include "ZPhiGammaGenLevel.h"
 
// constructors and destructor
ZPhiGammaGenLevel::ZPhiGammaGenLevel(const edm::ParameterSet& iConfig)
{
  genParticlesToken_            = consumes<std::vector<reco::GenParticle> >(edm::InputTag("genParticles"));

  create_trees();
}

ZPhiGammaGenLevel::~ZPhiGammaGenLevel()
{
}

// ------------ method called for each event  ------------
void ZPhiGammaGenLevel::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<reco::GenParticle>  > genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  int genZ_ID    = -999;
  float genZ_pT  = -999.;
  float genZ_eta = -999.;
  float genZ_phi = -999.;
  float genZ_E   = -999.;

  int genGamma_ID    = -999;
  float genGamma_pT  = -999.;
  float genGamma_eta = -999.;
  float genGamma_phi = -999.;
  float genGamma_E   = -999.;

  int genPhi_ID    = -999;
  float genPhi_pT  = -999.;
  float genPhi_eta = -999.;
  float genPhi_phi = -999.;
  float genPhi_E   = -999.;

  int genKminus_ID    = -999;
  float genKminus_pT  = -999.;
  float genKminus_eta = -999.;
  float genKminus_phi = -999.;
  float genKminus_E   = -999.;

  int genKplus_ID    = -999;
  float genKplus_pT  = -999.;
  float genKplus_eta = -999.;
  float genKplus_phi = -999.;
  float genKplus_E   = -999.;

  genZ_ID_tree  = 0;
  genZ_pT_tree  = 0.;
  genZ_eta_tree = 0.;
  genZ_phi_tree = 0.;
  genZ_E_tree   = 0.;

  genGamma_ID_tree  = 0;
  genGamma_pT_tree  = 0.;
  genGamma_eta_tree = 0.;
  genGamma_phi_tree = 0.;
  genGamma_E_tree   = 0.;

  genPhi_ID_tree  = 0;
  genPhi_pT_tree  = 0.;
  genPhi_eta_tree = 0.;
  genPhi_phi_tree = 0.;
  genPhi_E_tree   = 0.;

  genKminus_ID_tree  = 0;
  genKminus_pT_tree  = 0.;
  genKminus_eta_tree = 0.;
  genKminus_phi_tree = 0.;
  genKminus_E_tree   = 0.;

  genKplus_ID_tree  = 0;
  genKplus_pT_tree  = 0.;
  genKplus_eta_tree = 0.;
  genKplus_phi_tree = 0.;
  genKplus_E_tree   = 0.;


  for(auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
    
    //if it is a Higgs with 2 daughters (sometimes Pythia sends a H in itself, therefore only 1 daughter)
    if(gen->pdgId() == 23 && gen->numberOfDaughters() == 2){

      //for each daughter
      for(int i = 0; i < 2; i++){
	
	//if daughters are not Phi and gamma, continue
	if(!(gen->daughter(i)->pdgId() == 22 || gen->daughter(i)->pdgId() == 333)) continue;

	//save Higgs variables
	genZ_ID  = gen->pdgId();
	genZ_pT  = gen->pt();
	genZ_eta = gen->eta();
	genZ_phi = gen->phi();
	genZ_E   = gen->energy();
	
	//if daughter(i) is a photon
	if(fabs(gen->daughter(i)->pdgId()) == 22){
	 
	  //save Gamma variables
	  genGamma_ID  = gen->daughter(i)->pdgId();
	  genGamma_pT  = gen->daughter(i)->pt();
	  genGamma_eta = gen->daughter(i)->eta();
	  genGamma_phi = gen->daughter(i)->phi();
	  genGamma_E   = gen->daughter(i)->energy();
	} 

	//if daughter(i) is a Phi
	if(gen->daughter(i)->pdgId() == 333){

	  //save Phi variables
	  genPhi_ID  = gen->daughter(i)->pdgId();
	  genPhi_pT  = gen->daughter(i)->pt();
	  genPhi_eta = gen->daughter(i)->eta();
	  genPhi_phi = gen->daughter(i)->phi();
	  genPhi_E   = gen->daughter(i)->energy();  

	  //if Phi has two daughters
	  if(gen->daughter(i)->numberOfDaughters() == 2){
	    
	    //for each Phi daughter
	    for(int j = 0; j < 2; j++){

	      //if daughter(j) is a K-
	      if(gen->daughter(i)->daughter(j)->pdgId() == -321){

		//save K- variables
		genKminus_ID = gen->daughter(i)->daughter(j)->pdgId();
		genKminus_pT = gen->daughter(i)->daughter(j)->pt();
		genKminus_eta = gen->daughter(i)->daughter(j)->eta();
		genKminus_phi = gen->daughter(i)->daughter(j)->phi();
		genKminus_E = gen->daughter(i)->daughter(j)->energy();
	      }
	      
	      //if daughter(j) is a K+
	      if(gen->daughter(i)->daughter(j)->pdgId() == 321){

		//save K+ variables
		genKplus_ID = gen->daughter(i)->daughter(j)->pdgId();
		genKplus_pT = gen->daughter(i)->daughter(j)->pt();
		genKplus_eta = gen->daughter(i)->daughter(j)->eta();
		genKplus_phi = gen->daughter(i)->daughter(j)->phi();
		genKplus_E = gen->daughter(i)->daughter(j)->energy();
	      }

	    } //j-forloop end
	  }//"if Phi has two daughters" end
	}//"if it is a Phi" end
      }//i-forloop end
    }//"if it is a Higgs" end
  }//gen-forloop end
  
  genZ_ID_tree  = genZ_ID;
  genZ_pT_tree  = genZ_pT;
  genZ_eta_tree = genZ_eta;
  genZ_phi_tree = genZ_phi;
  genZ_E_tree   = genZ_E;
  
  genGamma_ID_tree  = genGamma_ID;
  genGamma_pT_tree  = genGamma_pT;
  genGamma_eta_tree = genGamma_eta;
  genGamma_phi_tree = genGamma_phi;
  genGamma_E_tree   = genGamma_E;
  
  genPhi_ID_tree  = genPhi_ID;
  genPhi_pT_tree  = genPhi_pT;
  genPhi_eta_tree = genPhi_eta;
  genPhi_phi_tree = genPhi_phi;
  genPhi_E_tree   = genPhi_E;

  genKminus_ID_tree = genKminus_ID;
  genKminus_pT_tree = genKminus_pT;
  genKminus_eta_tree = genKminus_eta;
  genKminus_phi_tree = genKminus_phi;
  genKminus_E_tree = genKminus_E;

  genKplus_ID_tree = genKplus_ID;
  genKplus_pT_tree = genKplus_pT;
  genKplus_eta_tree = genKplus_eta;
  genKplus_phi_tree = genKplus_phi;
  genKplus_E_tree = genKplus_E;
  
  mytree->Fill();
}

//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//

void ZPhiGammaGenLevel::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen info");

  mytree->Branch("genZ_ID",&genZ_ID_tree);
  mytree->Branch("genZ_pT",&genZ_pT_tree);
  mytree->Branch("genZ_eta",&genZ_eta_tree);
  mytree->Branch("genZ_phi",&genZ_phi_tree);
  mytree->Branch("genZ_E",&genZ_E_tree);
  
  mytree->Branch("genGamma_ID",&genGamma_ID_tree);
  mytree->Branch("genGamma_pT",&genGamma_pT_tree);
  mytree->Branch("genGamma_eta",&genGamma_eta_tree);
  mytree->Branch("genGamma_phi",&genGamma_phi_tree);
  mytree->Branch("genGamma_E",&genGamma_E_tree);
  
  mytree->Branch("genPhi_ID",&genPhi_ID_tree);
  mytree->Branch("genPhi_pT",&genPhi_pT_tree);
  mytree->Branch("genPhi_eta",&genPhi_eta_tree);
  mytree->Branch("genPhi_phi",&genPhi_phi_tree);
  mytree->Branch("genPhi_E",&genPhi_E_tree);

  mytree->Branch("genKminus_ID",&genKminus_ID_tree);
  mytree->Branch("genKminus_pT",&genKminus_pT_tree);
  mytree->Branch("genKminus_eta",&genKminus_eta_tree);
  mytree->Branch("genKminus_phi",&genKminus_phi_tree);
  mytree->Branch("genKminus_E",&genKminus_E_tree);

  mytree->Branch("genKplus_ID",&genKplus_ID_tree);
  mytree->Branch("genKplus_pT",&genKplus_pT_tree);
  mytree->Branch("genKplus_eta",&genKplus_eta_tree);
  mytree->Branch("genKplus_phi",&genKplus_phi_tree);
  mytree->Branch("genKplus_E",&genKplus_E_tree);
    
}

void ZPhiGammaGenLevel::beginJob()
{
}

void ZPhiGammaGenLevel::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZPhiGammaGenLevel);
