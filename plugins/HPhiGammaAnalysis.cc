//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"
#include <stdlib.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
  
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
 
//Vertex inclusions
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/BeamSpot/interface/BeamSpot.h" 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//Electron ID stuff
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

//Photon ID stuff
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

typedef math::XYZTLorentzVector LorentzVector;

using namespace std;  
 
#include "HPhiGammaAnalysis.h"

 
// constructors and destructor
HPhiGammaAnalysis::HPhiGammaAnalysis(const edm::ParameterSet& iConfig) :
  runningOnData_(iConfig.getParameter<bool>("runningOnData")),
  verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose")),
  effectiveAreas_el_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_el")).fullPath() ),
  effectiveAreas_ph_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_ph")).fullPath() )
{
  packedPFCandidatesToken_            = consumes<std::vector<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates")); 
  slimmedMuonsToken_                  = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
  prunedGenParticlesToken_            = consumes<std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
  photonsMiniAODToken_                = consumes<std::vector<pat::Photon> > (edm::InputTag("slimmedPhotons"));
  electronsMiniAODToken_              = consumes<std::vector<pat::Electron> > (edm::InputTag("slimmedElectrons"));
  slimmedJetsToken_                   = consumes<std::vector<pat::Jet> >(edm::InputTag("slimmedJets"));
  slimmedMETsToken_                   = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETs"));
  slimmedMETsPuppiToken_              = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETsPuppi"));
  offlineSlimmedPrimaryVerticesToken_ = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));  
  offlineBeamSpotToken_               = consumes<reco::BeamSpot> (edm::InputTag("offlineBeamSpot"));
  pileupSummaryToken_                 = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
  GenInfoToken_                       = consumes<GenEventInfoProduct> (edm::InputTag("generator"));
  triggerBitsToken_                   = consumes<edm::TriggerResults> (edm::InputTag("TriggerResults","","HLT"));
  rhoToken_                           = consumes<double> (iConfig.getParameter <edm::InputTag>("rho"));

  h_Events = fs->make<TH1F>("h_Events", "Event counting in different steps", 8, 0., 8.);

  _Nevents_processed  = 0;
  _Nevents_triggered = 0;
  _Nevents_isTwoKaons   = 0;
  _Nevents_isPhoton   = 0;
  _Nevents_HiggsFound = 0;
  _Nevents_HiggsNotMatched = 0;
  _Nevents_bestCoupleFound = 0;
  _Nevents_candPtFilter = 0;
  _Nevents_coupleIsolationFilter = 0;

  debug=false;  //DEBUG datamember  

  h_pileup   = fs->make<TH1F>("pileup", "pileup", 75,0,75);

  create_trees();
}


HPhiGammaAnalysis::~HPhiGammaAnalysis()
{
}

//test method start 
/* ==========================================================================================
   countPrimaryVertex
   ------------------
   This function returns the number of primary vertex of the event
   INPUT: slimmedPV
   OUTPUT: number of valid PV
   ==========================================================================================
*/
/*int HPhiGammaAnalysis::countPrimaryVertex(edm::Handle<std::vector<reco::Vertex > > slimmedPV)
{
  int numberOfPV = 0;

  if(slimmedPV->size()<=0) return numberOfPV;
  for(reco::VertexCollection::const_iterator vtx=slimmedPV->begin();vtx!=slimmedPV->end();++vtx) 
    {
      // check that the primary vertex is not a fake one, that is the beamspot (it happens when no primary vertex is reconstructed)
      if(!vtx->isFake()) {
	numberOfPV++;
      }
    } 
  return numberOfPV;
}
//test method end
*/


// ------------ method called for each event  ------------
void HPhiGammaAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  if(debug) cout<<"Starting analyze method"<<endl;

  edm::Handle<std::vector<pat::PackedCandidate>  > PFCandidates;
  iEvent.getByToken(packedPFCandidatesToken_, PFCandidates);

  edm::Handle<std::vector<pat::Muon>  > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_, slimmedMuons);

  edm::Handle<std::vector<reco::GenParticle>  > genParticles;
  if(!runningOnData_)iEvent.getByToken(prunedGenParticlesToken_, genParticles);

  edm::Handle<std::vector<pat::Photon> > slimmedPhotons;
  iEvent.getByToken(photonsMiniAODToken_,slimmedPhotons);

  edm::Handle<std::vector<pat::Electron> > slimmedElectrons;
  iEvent.getByToken(electronsMiniAODToken_,slimmedElectrons);

  edm::Handle<std::vector<pat::Jet > > slimmedJets;
  iEvent.getByToken(slimmedJetsToken_, slimmedJets);

  edm::Handle<std::vector<pat::MET > > slimmedMETs;
  iEvent.getByToken(slimmedMETsToken_, slimmedMETs);

  edm::Handle<std::vector<pat::MET > > slimmedMETsPuppi;
  iEvent.getByToken(slimmedMETsPuppiToken_, slimmedMETsPuppi);

  edm::Handle<std::vector<reco::Vertex > > slimmedPV;
  iEvent.getByToken(offlineSlimmedPrimaryVerticesToken_, slimmedPV);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsToken_, triggerBits);

  _Nevents_processed++;

  //Retrieve the run number
  if(runningOnData_){
    run_number = iEvent.id().run();
  }

  //*************************************************************//
  //                                                             //
  //-------------------------- Vertices -------------------------//
  //                                                             //
  //*************************************************************//

  //Count the number of vertices
  nPV = 0;
  
  if(slimmedPV->size()<=0) return;
  for(reco::VertexCollection::const_iterator vtx=slimmedPV->begin();vtx!=slimmedPV->end();++vtx) {
    // check that the primary vertex is not a fake one, that is the beamspot (it happens when no primary vertex is reconstructed)
    if(!vtx->isFake()) {
      nPV++;
    }
  } 
  
  //nPV = countPrimaryVertex(slimmedPV);
  //if(nPV==0)return;  
  //cout<<"npv="<<nPV<<endl;
  // std::cout << "slimmedPV size: " << slimmedPV->size() << "   PV: " << &(slimmedPV->at(0))  << std::endl;

  //*************************************************************//
  //                                                             //
  //--------------------------- Pile Up -------------------------//
  //                                                             //
  //*************************************************************//

  PU_Weight = -1.;
  float npT = -1.;

  if(!runningOnData_){
    edm::Handle<std::vector< PileupSummaryInfo>>  PupInfo;
    iEvent.getByToken(pileupSummaryToken_, PupInfo);
  
    std::vector<PileupSummaryInfo>::const_iterator PVI; 
 
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      const int BX = PVI->getBunchCrossing();
      if(BX == 0) {
	npT  = PVI->getTrueNumInteractions();
      }
    }

    if(npT == -1) {
      std::cout << "!!!! npT = -1 !!!!" << std::endl;
      abort();
    }
    // Calculate weight using above code
    PU_Weight = Lumiweights_.weight(npT);

    // Fill histogram with PU distribution
    h_pileup->Fill(npT);
  }

  //*************************************************************//
  //                                                             //
  //-------------------------- MC Weight ------------------------//
  //                                                             //
  //*************************************************************//

  MC_Weight = -10000000.;

  if(!runningOnData_){
    edm::Handle<GenEventInfoProduct> GenInfo;
    iEvent.getByToken(GenInfoToken_, GenInfo);
    
    float _aMCatNLOweight = GenInfo->weight();
    MC_Weight = _aMCatNLOweight;

    if(MC_Weight == -10000000.) {
      std::cout << "!!!! MC_Weight = -10000000 !!!!" << std::endl;
      abort();
    }
  }



  //*************************************************************//
  //                                                             //
  //--------------------------- Trigger -------------------------//
  //                                                             //
  //*************************************************************//

  //Examine the trigger information
  isTwoProngTrigger = false;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for(unsigned int i = 0, n = triggerBits->size(); i < n; ++i){
    if(!triggerBits->accept(i)) continue;
    std::string tmp_triggername = names.triggerName(i);

    if( tmp_triggername.find("HLT_Photon35_TwoProngs35_v") != std::string::npos ){
      isTwoProngTrigger = true;
    }
  }
  if(!isTwoProngTrigger) return;
  _Nevents_triggered++;

  //*************************************************************//
  //                                                             //
  //------------------ Variable initialization ------------------//
  //                                                             //
  //*************************************************************//


  nMuons     = 0;
  nElectrons = 0;
  nPhotonsOverSelection = 0;
  nPhotonsChosen   = 0;
  nJets      = 0;
  nJets_25   = 0;
  
  //These variables will go in the tree
  ph_eT     = 0.;
  ph_eta    = 0.;
  ph_etaSC  = 0.;
  ph_phi    = 0.;
  ph_energy = 0.;
  LorentzVector ph_p4;

  ph_iso_ChargedHadron = 0.;
  ph_iso_NeutralHadron = 0.;
  ph_iso_Photon = 0.;
  ph_iso_eArho = 0.;

  is_photon_a_photon = false;
  is_photon_matched  = false;

  eTphMax  = -1000.;

  _Jet_Photon_invMass = -1.;
  _Phimass = -1.;
  _Hmass_From2K_Photon = -1.;
 
  met_pT = 0.;
  metpuppi_pT = 0.;

  _bestJet_pT=-1.;
  _bestJet_eta=-1.;
  _bestJet_phi=-1.;
  _bestJet_nDaughters=0;
  _bestJet_pTMax=-1.;
  _bestJet_chargedEmEnergy=0.;
  _bestJet_neutralEmEnergy=0.;
  _bestJet_chargedHadEnergy=0.;
  _bestJet_neutralHadEnergy=0.;
  _bestJet_chargedHadMultiplicity=0.;
  _bestJet_invMass=0.;
  _bestJet_Photon_invMass=0.;
  _isHiggsFound = false;
  _firstCandPt=0.;
  _firstCandEta=0.;
  _firstCandPhi=0.;
  _secondCandPt=0.;
  _secondCandEta=0.;
  _secondCandPhi=0.;
  _bestCouplePt=0.;
  _bestCoupleEta=0.;
  _bestCouplePhi=0.;  

  //*************************************************************//
  //                                                             //
  //----------------------------- MET ---------------------------//
  //                                                             //
  //*************************************************************//
  for(auto met = slimmedMETs->begin(); met != slimmedMETs->end(); ++met){
    met_pT = met->pt();
  }

  for(auto metpuppi = slimmedMETsPuppi->begin(); metpuppi != slimmedMETsPuppi->end(); ++metpuppi){
    metpuppi_pT = metpuppi->pt();
  }

  //*************************************************************//
  //                                                             //
  //---------------------------- Muons --------------------------//
  //                                                             //
  //*************************************************************//

  for(auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
    if(mu->pt() < 20. || !mu->CutBasedIdMedium || fabs(mu->eta()) > 2.4 || fabs(mu->muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(mu->muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
    if(!mu->PFIsoLoose) continue;
    nMuons++;
  }

  //*************************************************************//
  //                                                             //
  //-------------------------- Electrons ------------------------//
  //                                                             //
  //*************************************************************//

  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  float corr_pt = 0.;

  for(auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
    //Calculate electron p4, correct it with the Scale&Smearing correction and extract the pT
    LorentzVector el_p4 = el->p4(); // * el->userFloat("ecalTrkEnergyPostCorr")/el->energy();
    corr_pt = el_p4.pt();

    if(corr_pt < 20. || fabs(el->eta()) > 2.5 || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;

    float abseta = fabs(el->superCluster()->eta());
    float eA     = effectiveAreas_el_.getEffectiveArea(abseta);
    float el_iso   = (el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/corr_pt;
    if(el_iso > 0.35) continue;

    //-------------Conditions on loose/medium MVA electron ID-------------//
    if(el->electronID("mvaEleID-Fall17-iso-V1-wp90") == 0) continue;
    nElectrons++;
  }

  //std::cout << "Nelectrons " << nElectrons << " Nmuons " << nMuons << std::endl;

  //*************************************************************//
  //                                                             //
  //----------------------- Access MC Truth ---------------------//
  //                                                             //
  //*************************************************************//

  //In signal, identify if there's a real mu or ele from W
  is_Kplus_fromPhi = false;
  is_Kminus_fromPhi = false;
  is_Phi_fromH = false;
  is_Photon_fromH = false;

  if(!runningOnData_){
    for(auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
      if( gen->pdgId() == 321  && gen->mother()->pdgId() == 333)  is_Kplus_fromPhi  = true;
      if( gen->pdgId() == -321 && gen->mother()->pdgId() == 333)  is_Kminus_fromPhi = true;
      if( gen->pdgId() == 333  && gen->mother()->pdgId() == 25)   is_Phi_fromH      = true;
      if( gen->pdgId() == 22   && gen->mother()->pdgId() == 25)   is_Photon_fromH   = true;
    }
  }  


  //*************************************************************//
  //                                                             //
  //--------------------------- Photons -------------------------//
  //                                                             //
  //*************************************************************//

  bool cand_photon_found = false;
  float corr_et = -1.;

  for(auto photon = slimmedPhotons->begin(); photon != slimmedPhotons->end(); ++photon){

    corr_et = photon->et(); // * photon->userFloat("ecalEnergyPostCorr")/photon->energy();

    //std::cout << "photon et " << corr_et << std::endl;

    if(corr_et < 20. || fabs(photon->eta()) > 2.5) continue;
    if(photon->hasPixelSeed()) continue;   //electron veto

    //std::cout << "photon" << photon->photonID("mvaPhoID-RunIIFall17-v1-wp90") << " " << photon->photonID("mvaPhoID-RunIIFall17-v1p1-wp90") << std::endl;

    if(photon->photonID("mvaPhoID-RunIIFall17-v1p1-wp90") == 0) continue;

    float abseta = fabs(photon->superCluster()->eta());
    float eA = effectiveAreas_ph_.getEffectiveArea(abseta);
    //photon_iso = (pfIso.sumChargedHadronPt + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_))/photon->et();

    if(photon->chargedHadronIso()/corr_et > 0.3 || photon->photonIso() > 4.) continue; //|| photon->trackIso() > 6

    nPhotonsOverSelection++;

    if(corr_et < eTphMax) continue;
    eTphMax = corr_et;
    ph_iso_ChargedHadron = photon->chargedHadronIso();
    ph_iso_NeutralHadron = photon->neutralHadronIso();
    ph_iso_Photon        = photon->photonIso();
    ph_iso_eArho         = eA*rho_;

    ph_eT     = corr_et;
    ph_eta    = photon->eta();
    ph_etaSC  = photon->superCluster()->eta();
    ph_phi    = photon->phi();

    // Apply energy scale corrections to MC
    ph_energy = photon->energy(); //userFloat("ecalEnergyPostCorr");
    ph_p4     = photon->p4();  //* photon->userFloat("ecalEnergyPostCorr")/photon->energy();
    
    cand_photon_found = true;
    nPhotonsChosen++;

    float deltapTMax = 10000.;
    const float deltaRMax = 0.3;
    int   gen_mother = 0;
    int   gen_ID = 0;

    is_photon_a_photon = false;
    is_photon_matched = false;

    if(!runningOnData_){
      for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){

	//phi folding	
	float deltaPhi = fabs(ph_phi-gen->phi());
	if (deltaPhi > 3.14) deltaPhi = 6.28 - deltaPhi;
	  	      	
	float deltaR = sqrt((ph_eta-gen->eta())*(ph_eta-gen->eta())+deltaPhi*deltaPhi);
	float deltapT = fabs(ph_eT-gen->pt());

	if(deltaR > deltaRMax || deltapT > deltapTMax) continue;
	deltapTMax = deltapT;
	gen_ID = gen->pdgId();
	gen_mother = gen->mother()->pdgId();
      }
            
      if(gen_ID == 22) is_photon_a_photon = true;
      if(gen_ID != 22) std::cout << "ph gen ID = " << gen_ID << std::endl;
      if(gen_ID != 22 && fabs(gen_mother) == 24) std::cout << "ph gen ID when matched = " << gen_ID << std::endl;
      if(fabs(gen_mother) == 25) is_photon_matched = true;
    }

  }

  //Do not continue if there's no photons
  if(!cand_photon_found) return;
  _Nevents_isPhoton++;

  //  std::cout << "Nphotons " << _Nevents_isPhoton << std::endl;

  
  //*************************************************************//
  //                                                             //
  //--------------------------- N-jets --------------------------//
  //                                                             //
  //*************************************************************//
  
    
  //int nJet=1;
  int jetIndex=-1;
  int bestJet_Index=-1;
  int MCtruthIndex = -1;
  float deltaR = -1;   
  int nDaughters = 0;
  //bool isBestJetFound = false; 

  //daughters forloop variable
  int firstCandCharge;
  int secondCandCharge;
  float firstCandPt;
  float secondCandPt;
  LorentzVector firstCand_p4;
  LorentzVector secondCand_p4;
  LorentzVector couple_p4;
  LorentzVector best_firstCand_p4;
  LorentzVector best_secondCand_p4;
  LorentzVector best_couple_p4;
  float bestCoupleOfTheJet_pT = 0.;
  float deltaR_KChosen = 0.;
  float deltaR_K = 0.;
  firstCandEnergy = 0.;
  secondCandEnergy = 0.;
  firstCandPx=0.;
  firstCandPy=0.;
  firstCandPz=0.;
  secondCandPx=0.;
  secondCandPy=0.;
  secondCandPz=0.;
  float firstCandEta=0.;
  float firstCandPhi=0.;
  float secondCandEta=0.;
  float secondCandPhi=0.;
  float kMass = 0.4937;
  float candPtMin = 1.;
  bool isBestCoupleOfTheEvent_Found=false;
  bool isBestCoupleOfTheJet_Found=false;


  for (auto jet = slimmedJets->begin(); jet != slimmedJets->end(); ++jet)  //jet loop start    
    {
      if(debug) cout<<"Starting Jet forloop"<<endl;    
      jetIndex++;
      isBestCoupleOfTheJet_Found=false; //needed for debugging     
      _Jet_Photon_invMass=(jet->p4()+ph_p4).M(); //calculate inv mass
      nDaughters= jet->numberOfDaughters(); //calculate number of daughters
      
      if(debug) cout<<"Starting Pre-Filters"<<endl;    
      //-----------------------------Pre-Filters--------------------------------------------------------
      if(jet->pt() < 40. || abs(jet->eta()) > 2.5) continue;
      if(_Jet_Photon_invMass < 100.) continue; //reject jets with inv mass lower then 100 GeV
      if(jet->neutralHadronEnergyFraction() > 0.9) continue; //reject if neutralhadron-energy fraction is >0.9
      if(jet->neutralEmEnergyFraction() > 0.9) continue; //reject if neutralEm-energy fraction is >0.9, alias NO-PHOTON FILTER                              
      if(nDaughters < 2) continue; //reject if number of constituens is less then 1
      if(jet->muonEnergyFraction() > 0.8) continue; //reject if muon-energy fraction is >0.8                                             
      if(jet->chargedHadronEnergyFraction() <= 0.) continue; //reject if chargedHadron-energy fraction is 0                              
      if(jet->chargedHadronMultiplicity() == 0) continue; //reject if there are NOT charged hadrons                              
      if(jet->chargedEmEnergyFraction() > 0.8) continue; //reject if chargedEm-energy fraction is >0.8                              
     //---------------------------------------------------------------------------------------------      
      
      if(debug) cout<<"Starting access MC truth"<<endl;    
      //-------------------------------access to MC truth-------------------------------------------
      if(!runningOnData_)
	{
	  if(debug) cout<<"runningOnData=False"<<endl;
	  if(MCtruthIndex == -1)
	    {      
	      for(auto gen = genParticles->begin(); gen != genParticles->end(); ++gen) //gen particles loop start
		{ 
		  if(debug) cout<<"Accessing to genParticles"<<endl;

		  //phi folding	
		  float deltaPhi = fabs(jet->phi()-gen->phi());
		  if (deltaPhi > 3.14) deltaPhi = 6.28 - deltaPhi;

		  float R = sqrt((jet->eta()-gen->eta())*(jet->eta()-gen->eta())+deltaPhi * deltaPhi); 
		  if( gen->pdgId() == 333  && gen->mother()->pdgId() == 25 && R < 0.4)
		    {
		      if(debug) cout<<"Accessing to pdgId"<<endl;
		      MCtruthIndex=jetIndex;
		      if(debug) cout<<"MCtruth index: "<<MCtruthIndex<<endl; //MCtruthIndex and deltaR can be used out from the jet-loop 
		      deltaR=R; 
		    }
		}
	      
	    } //gen particle loop end
	}
      //------------------------------MC truth logout-----------------------------------------------
      

      //-------------------------------------daughters forloop----------------------------
      for(int firstCand_Index=0; firstCand_Index < nDaughters; firstCand_Index++) //1st loop starts
	{
	  if(debug) cout<<"Starting Jet-unpackaging forloop"<<endl;    
	  firstCandPt= slimmedJets->at(jetIndex).daughter(firstCand_Index)->pt(); //extrapolate firstCand pt
	  firstCandEta= slimmedJets->at(jetIndex).daughter(firstCand_Index)->eta(); //extrapolate firstCand eta
	  firstCandPhi= slimmedJets->at(jetIndex).daughter(firstCand_Index)->phi(); //extrapolate firstCand phi
	  
	  if(firstCandPt < candPtMin) continue; //firstCand filter if pt < candPtMin
	  
	  for(int secondCand_Index=firstCand_Index+1; secondCand_Index < nDaughters; secondCand_Index++) //2nd loop starts
	    {	  
	      
	      secondCandPt= slimmedJets->at(jetIndex).daughter(secondCand_Index)->pt(); //extrapolate secondCand pt
	      secondCandEta= slimmedJets->at(jetIndex).daughter(secondCand_Index)->eta(); //extrapolate secondCand eta
	      secondCandPhi= slimmedJets->at(jetIndex).daughter(secondCand_Index)->phi(); //extrapolate secondCand phi
	      
	      //secondCand filter if both pt are < 10GeV
	      if(secondCandPt < candPtMin) continue;
	      if(firstCandPt < 10. && secondCandPt < 10.) continue; 
	      
	      //third filter on deltaR
	      float deltaEta= firstCandEta - secondCandEta;

	      float deltaPhi = fabs(firstCandPhi - secondCandPhi);  //phi folding	
	      if (deltaPhi > 3.14) deltaPhi = 6.28 - deltaPhi;

	      deltaR_K= sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
	      if(deltaR_K > 0.02) continue;
	      
	      firstCandCharge = slimmedJets->at(jetIndex).daughter(firstCand_Index)->charge(); //extrapolate firstCand charge
	      secondCandCharge = slimmedJets->at(jetIndex).daughter(secondCand_Index)->charge(); //extrapolate secondCand charge
	      //------------------------------------DEBUG----------------------------------------------------
	      //cout<<"Particle n."<<firstCand_Index+1<<" pt= "<<firstCandPt<<" Q= "<<firstCandCharge<<endl; //debug       
	      //cout<<"Particle n."<<secondCand_Index+1<<" pt= "<<secondCandPt<<" Q= "<<secondCandCharge<<endl; //debug       
	      //--------------------------------------------------------------------------------------------
	      
	      //OPPOSITE CHARGE - FILTER
	      if(firstCandCharge*secondCandCharge >= 0) continue; //choose only opposite charges
	      
	      //PIONS INTO KAONS CORRECTION		  		  
	      firstCand_p4 = slimmedJets->at(jetIndex).daughter(firstCand_Index)->p4(); //extrapolate quadrimomentum
	      secondCand_p4 = slimmedJets->at(jetIndex).daughter(secondCand_Index)->p4();

	      firstCandPx = firstCand_p4.px(); //extrapolate px, py, pz of the first candidate
	      firstCandPy = firstCand_p4.py();
 	      firstCandPz = firstCand_p4.pz();
	      secondCandPx = secondCand_p4.px(); //extrapolate px, py, pz of the second candidate
	      secondCandPy = secondCand_p4.py();
	      secondCandPz = secondCand_p4.pz();
	      
	      firstCandEnergy = sqrt(firstCandPx*firstCandPx + firstCandPy*firstCandPy + firstCandPz*firstCandPz + kMass*kMass); //energy recalculation
	      secondCandEnergy = sqrt(secondCandPx*secondCandPx + secondCandPy*secondCandPy + secondCandPz*secondCandPz + kMass*kMass); //energy recalculation

	      firstCand_p4.SetE(firstCandEnergy); //quadrimomentum correction
	      secondCand_p4.SetE(secondCandEnergy); //quadrimomentum correction
	      couple_p4 = firstCand_p4 + secondCand_p4; //calculation of the couple-quadrimomentum after the correction
	     
	      //PHI INV MASS - FILTER
	      _Phimass=(couple_p4).M(); //calculate inv mass of the Phi candidate	
	      if(_Phimass < 1. || _Phimass > 1.05) continue; //filter on phi invariant mass	      
	      
	      if(couple_p4.pt() < 30.) continue;
	      //PT MAX OF THE JET - FILTER
	      if(couple_p4.pt() <= bestCoupleOfTheJet_pT) continue; //choose the couple with greatest pt
	      bestCoupleOfTheJet_pT = couple_p4.pt();	      
	      isBestCoupleOfTheJet_Found=true; //needed for debugging
	      isBestCoupleOfTheEvent_Found=true;	

	      bestJet_Index = jetIndex; //note the position of the chosen jet inside the vector	      	      
	      deltaR_KChosen = deltaR_K;
	      _bestJet_Photon_invMass = _Jet_Photon_invMass;
	      best_firstCand_p4 = firstCand_p4; 
	      best_secondCand_p4 = secondCand_p4;
	      best_couple_p4 = couple_p4;	      
	    }
	  //2nd loop ends
	  if(debug) cout<<"Ending Jet-unpackaging forloop"<<endl;    
	} //1st loop ends
      
      if(debug){
	if(!isBestCoupleOfTheJet_Found) cout<<"No best couple detected for jetIndex = "<<jetIndex<<endl;
      }
      
      if(jet->pt() < 25.) continue;
      nJets_25++;
      if(jet->pt() < 30.) continue;
      nJets++;
      
      if(debug) cout<<"Ending Jet forloop"<<endl;        
      
    } //jet loop end

  if(debug) cout<<"jet forloop logout"<<endl;  
  
  if(!isBestCoupleOfTheEvent_Found) 
    {
      if(debug) cout<<"No best couple detected for current event"<<endl;
      return;
    }
  _Nevents_bestCoupleFound++;      
  if(debug) cout<<"_Nevents_bestCoupleFound: "<<_Nevents_bestCoupleFound<<endl;        

  //DATAMEMBER SAVING
  _firstCandPt= best_firstCand_p4.pt();
  _firstCandEta= best_firstCand_p4.eta();
  _firstCandPhi= best_firstCand_p4.phi();
  _secondCandPt= best_secondCand_p4.pt();	      
  _secondCandEta= best_secondCand_p4.eta();
  _secondCandPhi= best_secondCand_p4.phi();
  _bestCouplePt = best_couple_p4.pt();
  _bestCoupleEta = best_couple_p4.eta();
  _bestCouplePhi = best_couple_p4.phi();
  
  _bestJet_invMass=slimmedJets->at(bestJet_Index).mass();
  _bestJet_pT=slimmedJets->at(bestJet_Index).pt();
  _bestJet_eta=slimmedJets->at(bestJet_Index).eta();
  _bestJet_phi=slimmedJets->at(bestJet_Index).phi();
  _bestJet_nDaughters=slimmedJets->at(bestJet_Index).numberOfDaughters();
  _bestJet_chargedEmEnergy=slimmedJets->at(bestJet_Index).chargedEmEnergy();
  _bestJet_neutralEmEnergy=slimmedJets->at(bestJet_Index).neutralEmEnergy();
  _bestJet_chargedHadEnergy=slimmedJets->at(bestJet_Index).chargedHadronEnergy();
  _bestJet_neutralHadEnergy=slimmedJets->at(bestJet_Index).neutralHadronEnergy();
  _bestJet_chargedHadMultiplicity=slimmedJets->at(bestJet_Index).chargedHadronMultiplicity();
 
  //PHI MASS CALCULATION
  _Phimass = (best_firstCand_p4 + best_secondCand_p4).M();
  
  //H INV MASS CALCULATION
  _Hmass_From2K_Photon = (best_firstCand_p4 + best_secondCand_p4 + ph_p4).M(); //calculate inv mass of the Higgs candidate
  
  //K-candidates and PHI ISOLATION
  K1_sum_pT_05    = 0.;
  K1_sum_pT_05_ch = 0.;
  
  K2_sum_pT_05    = 0.;
  K2_sum_pT_05_ch = 0.;

  couple_sum_pT_05 = 0.;
  couple_sum_pT_05_ch = 0.;

  _iso_K1 = 0.;
  _iso_K1_ch = 0.;
  _iso_K2 = 0.;
  _iso_K2_ch = 0.;
  _iso_couple = 0.;
  _iso_couple_ch = 0.;
    
  for(auto cand_iso = PFCandidates->begin(); cand_iso != PFCandidates->end(); ++cand_iso){
    
    if(cand_iso->pt() < 0.5) continue; //do not consider tracks with pt < 500MeV
    
    float deltaPhi_K1 = fabs(_firstCandPhi-cand_iso->phi());  //phi folding	
    if (deltaPhi_K1 > 3.14) deltaPhi_K1 = 6.28 - deltaPhi_K1;

    float deltaR_K1 = sqrt((_firstCandEta-cand_iso->eta())*(_firstCandEta-cand_iso->eta()) + deltaPhi_K1*deltaPhi_K1);
    if(deltaR_K1 < 0.02) continue;
    
    float deltaPhi_K2 = fabs(_secondCandPhi-cand_iso->phi());  //phi folding	
    if (deltaPhi_K2 > 3.14) deltaPhi_K2 = 6.28 - deltaPhi_K2;

    float deltaR_K2 = sqrt((_secondCandEta-cand_iso->eta())*(_secondCandEta-cand_iso->eta()) + deltaPhi_K2*deltaPhi_K2);
    if(deltaR_K2 < 0.02) continue;

    float deltaPhi_Couple = fabs(_bestCouplePhi-cand_iso->phi());  //phi folding	
    if (deltaPhi_Couple > 3.14) deltaPhi_Couple = 6.28 - deltaPhi_Couple;

    float deltaR_Couple = sqrt((_bestCoupleEta-cand_iso->eta())*(_bestCoupleEta-cand_iso->eta()) + deltaPhi_Couple*deltaPhi_Couple);

    if(deltaR_K1 <= 0.3) K1_sum_pT_05 += cand_iso->pt();
    if(deltaR_K2 <= 0.3) K2_sum_pT_05 += cand_iso->pt();
    if(deltaR_Couple <= 0.3) couple_sum_pT_05 += cand_iso->pt();

    if(cand_iso->charge() != 0 && (fabs(cand_iso->dxy()) >= 0.2 || fabs(cand_iso->dz()) >= 0.5) ) continue; // Requesting charged particles to come from PV
    if(deltaR_K1 <= 0.3) K1_sum_pT_05_ch += cand_iso->pt();
    if(deltaR_K2 <= 0.3) K2_sum_pT_05_ch += cand_iso->pt();
    if(deltaR_Couple <= 0.3) couple_sum_pT_05_ch += cand_iso->pt();
  }
  
  //CANDIDATES SORTING
  if(_firstCandPt < _secondCandPt)  //swap-values loop, in order to fill the tree with the candidate with max pt of the couple in firstCand branches  
    {                               //and one with min pt in secondCand branches
      float a,b,c,d;
      a = _firstCandPt;
      b = _firstCandEta;
      c = _firstCandPhi;
      d = firstCandEnergy;
      _firstCandPt = _secondCandPt;
      _firstCandEta = _secondCandEta;
      _firstCandPhi = _secondCandPhi;
      firstCandEnergy = secondCandEnergy;	      
      _secondCandPt = a;
      _secondCandEta = b;
      _secondCandPhi = c;
      secondCandEnergy = d;
    }

  //CUTS ON CANDIDATES PT
  if(_firstCandPt < 15. || _secondCandPt < 5.) return;
  _Nevents_candPtFilter++;

  //ISOLATION DATAMEMBER FOR TREE FILLING 
  _iso_K1 = K1_sum_pT_05/_firstCandPt;
  _iso_K2 = K2_sum_pT_05/_secondCandPt;
  _iso_couple = couple_sum_pT_05/_bestCouplePt;
  _iso_K1_ch = K1_sum_pT_05_ch/_firstCandPt;
  _iso_K2_ch = K2_sum_pT_05_ch/_secondCandPt;
  _iso_couple_ch = couple_sum_pT_05_ch/_bestCouplePt;

  //CUT ON PHI ISOLATION
  if(_iso_couple_ch > 1.) return;
  _Nevents_coupleIsolationFilter++;

  
  //MC TRUTH CHECK
  if(!runningOnData_)
    {
      _isHiggsFound=false; //bool initialization
      
      if(MCtruthIndex == bestJet_Index) //if the index of the best jet matches with one of the MC truth, it passes here
	{
	  if(debug) cout<<endl<<"****************TRUE HIGGS FOUND******************"<<endl;
	  if(debug) cout<<"Higgs deltaR = "<<deltaR<<endl;
	  _Nevents_HiggsFound++;
	  _isHiggsFound=true;
	}  
      else 
	{
	  _Nevents_HiggsNotMatched++;
	  if(debug) cout<<endl<<"THAT'S NOT A HIGGS!"<<endl;
	}
      
      //some prints
      
      if(debug){
	cout<<"Jet + photon inv. mass = "<<_bestJet_Photon_invMass<<endl;
	cout<<"n. of daughters: "<<_bestJet_nDaughters<<endl;
	cout<<"Best couple pT = "<<_firstCandPt + _secondCandPt<<endl;
	cout<<"Best couple DeltaR = "<<deltaR_KChosen<<endl;
	cout<<"Phi candidate inv. mass  = "<<_Phimass<<endl;
	cout<<"H inv. mass after jet-unpackaging = "<<_Hmass_From2K_Photon<<endl;
	cout<<"--------------------------------------------------"<<endl;
	cout<<"Higgs found = "<<_Nevents_HiggsFound<<",   Higgs NOT matched = "<<_Nevents_HiggsNotMatched<<endl;
	cout<<"--------------------------------------------------"<<endl<<endl;
      }
    }
  
  mytree->Fill();
  if(debug) cout<<"Ending analyze method"<<endl;
}

//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//

void HPhiGammaAnalysis::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");

  mytree->Branch("nPV",&nPV);
  mytree->Branch("isTwoProngTrigger",&isTwoProngTrigger);

  //Save run number info when running on data
  if(runningOnData_){
    mytree->Branch("run_number",&run_number);
  }

  mytree->Branch("nMuons",&nMuons);
  mytree->Branch("nElectrons",&nElectrons);
  mytree->Branch("nPhotonsOverSelection",&nPhotonsOverSelection);
  mytree->Branch("nPhotonsChosen",&nPhotonsChosen);
  mytree->Branch("nJets",&nJets);
  mytree->Branch("nJets_25",&nJets_25);
  mytree->Branch("met_pT",&met_pT);
  mytree->Branch("metpuppi_pT",&metpuppi_pT);

  mytree->Branch("photon_eT",&ph_eT);
  mytree->Branch("photon_eta",&ph_eta);
  mytree->Branch("photon_etaSC",&ph_etaSC);
  mytree->Branch("photon_phi",&ph_phi);
  mytree->Branch("photon_energy",&ph_energy);
  mytree->Branch("photon_iso_ChargedHadron",&ph_iso_ChargedHadron);
  mytree->Branch("photon_iso_NeutralHadron",&ph_iso_NeutralHadron);
  mytree->Branch("photon_iso_Photon",&ph_iso_Photon);
  mytree->Branch("photon_iso_eArho",&ph_iso_eArho);

  mytree->Branch("bestJet_pT",&_bestJet_pT);
  mytree->Branch("bestJet_eta",&_bestJet_eta);
  mytree->Branch("bestJet_phi",&_bestJet_phi);
  mytree->Branch("bestJet_nDaughters",&_bestJet_nDaughters);
  mytree->Branch("bestJet_chargedEmEnergy",&_bestJet_chargedEmEnergy);
  mytree->Branch("bestJet_neutralEmEnergy",&_bestJet_neutralEmEnergy);
  mytree->Branch("bestJet_chargedHadEnergy",&_bestJet_chargedHadEnergy);
  mytree->Branch("bestJet_neutralHadEnergy",&_bestJet_neutralHadEnergy);
  mytree->Branch("bestJet_invMass",&_bestJet_invMass);
  mytree->Branch("bestJet_Photon_invMass",&_bestJet_Photon_invMass);
  
  mytree->Branch("firstCandPt",&_firstCandPt);
  mytree->Branch("firstCandEta",&_firstCandEta);
  mytree->Branch("firstCandPhi",&_firstCandPhi);
  mytree->Branch("secondCandPt",&_secondCandPt);
  mytree->Branch("secondCandEta",&_secondCandEta);
  mytree->Branch("secondCandPhi",&_secondCandPhi);
  mytree->Branch("bestCouplePt",&_bestCouplePt);
  mytree->Branch("bestCoupleEta",&_bestCoupleEta);
  mytree->Branch("bestCouplePhi",&_bestCouplePhi);
  
  mytree->Branch("firstCandEnergy",&firstCandEnergy);
  mytree->Branch("secondCandEnergy",&secondCandEnergy);

  mytree->Branch("Phimass",&_Phimass);
  mytree->Branch("Hmass_From2K_Photon",&_Hmass_From2K_Photon);

  mytree->Branch("K1_sum_pT_05",&K1_sum_pT_05);
  mytree->Branch("K1_sum_pT_05_ch",&K1_sum_pT_05_ch);
  mytree->Branch("K2_sum_pT_05",&K2_sum_pT_05);
  mytree->Branch("K2_sum_pT_05_ch",&K2_sum_pT_05_ch);
  mytree->Branch("Couple_sum_pT_05",&couple_sum_pT_05);
  mytree->Branch("Couple_sum_pT_05_ch",&couple_sum_pT_05_ch);

  mytree->Branch("iso_K1",&_iso_K1);
  mytree->Branch("iso_K1_ch",&_iso_K1_ch);
  mytree->Branch("iso_K2",&_iso_K2);
  mytree->Branch("iso_K2_ch",&_iso_K2_ch);
  mytree->Branch("iso_couple",&_iso_couple);
  mytree->Branch("iso_couple_ch",&_iso_couple_ch);


  //Save MC info
  if(!runningOnData_){
    mytree->Branch("PU_Weight",&PU_Weight);
    mytree->Branch("MC_Weight",&MC_Weight);
    mytree->Branch("isHiggsFound",&_isHiggsFound);

    is_Kplus_fromPhi = false;
    is_Kminus_fromPhi = false;
    is_Phi_fromH = false;
    is_Photon_fromH = false;

    mytree->Branch("isKplusfromPhi",&is_Kplus_fromPhi);
    mytree->Branch("isKminusfromPhi",&is_Kminus_fromPhi);
    mytree->Branch("isPhiFromH",&is_Phi_fromH);
    mytree->Branch("isPhotonFromH",&is_Photon_fromH);

    mytree->Branch("isPhotonTrue",&is_photon_a_photon);
    mytree->Branch("isPhotonMatched",&is_photon_matched);

  }

}

void HPhiGammaAnalysis::beginJob()
{
  //Flag for PileUp reweighting
  if (!runningOnData_){ // PU reweighting for 2017
   Lumiweights_ = edm::LumiReWeighting("MCpileUp_2018_25ns_JuneProjectionFull18_PoissonOOTPU.root", "MyDataPileupHistogram.root", "pileup", "pileup");
  }
}


//*************************************************************//
//                                                             //
//------------------- Fill event loss histos ------------------//
//                                                             //
//*************************************************************//

void HPhiGammaAnalysis::endJob() 
{
  h_Events->Fill(0.5,_Nevents_processed);
  h_Events->Fill(1.5,_Nevents_triggered);
  h_Events->Fill(2.5,_Nevents_isPhoton);
  h_Events->Fill(3.5,_Nevents_bestCoupleFound);  
  h_Events->Fill(4.5,_Nevents_candPtFilter);  
  h_Events->Fill(5.5,_Nevents_coupleIsolationFilter);  

  h_Events->GetXaxis()->SetBinLabel(1,"Events processed");
  h_Events->GetXaxis()->SetBinLabel(2,"Events triggered");
  h_Events->GetXaxis()->SetBinLabel(3,"Photon requested");
  h_Events->GetXaxis()->SetBinLabel(4,"Best couple of the event found");
  h_Events->GetXaxis()->SetBinLabel(5,"Cand pT selection ");
  h_Events->GetXaxis()->SetBinLabel(6,"Phi isolation selection");
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPhiGammaAnalysis);
