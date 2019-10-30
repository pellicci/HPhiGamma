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
  _Nevents_isTwoKaons   = 0;
  _Nevents_isPhoton   = 0;
  _Nevents_isPhimass  = 0;
  _Nevents_isHmass  = 0;

  h_pileup   = fs->make<TH1F>("pileup", "pileup", 75,0,75);

  create_trees();
}

HPhiGammaAnalysis::~HPhiGammaAnalysis()
{
}

// ------------ method called for each event  ------------
void HPhiGammaAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
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
  nPV = -1;

  if(slimmedPV->size()<=0) return;
  for(reco::VertexCollection::const_iterator vtx=slimmedPV->begin();vtx!=slimmedPV->end();++vtx) {
    // check that the primary vertex is not a fake one, that is the beamspot (it happens when no primary vertex is reconstructed)
    if(!vtx->isFake()) {
      nPV++;
    }
  } 
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

  //*************************************************************//
  //                                                             //
  //------------------ Variable initialization ------------------//
  //                                                             //
  //*************************************************************//

  nMuons     = 0;
  nElectrons = 0;
  nPhotons   = 0;
  nJets      = 0;
  nJets_25   = 0;

  //These variables will go in the tree
  K1_pT     = 0.;
  K1_eta    = 0.;
  K1_phi    = 0.;
  K1_energy = 0.;
  K1_dxy    = 0.;
  K1_dz     = 0.;
  K1_charge = 0.;
  LorentzVector K1_p4;
  K1_sum_pT_03    = 0.;
  K1_sum_pT_05    = 0.;
  K1_sum_pT_05_ch = 0.;
  K1_pTMax        = -1000.;

  K2_pT     = 0.;
  K2_eta    = 0.;
  K2_phi    = 0.;
  K2_energy = 0.;
  K2_dxy    = 0.;
  K2_dz     = 0.;
  K2_charge = 0.;
  LorentzVector K2_p4;
  K2_sum_pT_03    = 0.;
  K2_sum_pT_05    = 0.;
  K2_sum_pT_05_ch = 0.;
  K2_pTMax        = -1000.;

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

  _Phimass = 0.;
  _Hmass = 0.;

  met_pT = 0.;
  metpuppi_pT = 0.;

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
    if(mu->pt() < 10. || !mu->CutBasedIdMedium || fabs(mu->eta()) > 2.4 || fabs(mu->muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(mu->muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
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

    if(corr_pt < 10. || fabs(el->eta()) > 2.5 || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;

    float abseta = fabs(el->superCluster()->eta());
    float eA     = effectiveAreas_el_.getEffectiveArea(abseta);
    float el_iso   = (el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/corr_pt;
    if(el_iso > 0.35) continue;

    //-------------Conditions on loose/medium MVA electron ID-------------//
    if(el->electronID("mvaEleID-Fall17-iso-V1-wp80") == 0) continue;
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
  //---------------------------- Kaons --------------------------//
  //                                                             //
  //*************************************************************//

  //The basic phylosophy is to get the two highest momentum isolated pions, and *then* check for the charge
  //This has little effect on signal efficiency, but large effect on background rejection

  for(auto cand = PFCandidates->begin(); cand != PFCandidates->end(); ++cand){
    
    if(cand->pt() < 20. || !cand->trackHighPurity() || fabs(cand->dxy()) >= 0.2 || fabs(cand->dz()) >= 0.5 ) continue;
    if(cand->pt() < K2_pTMax) continue;

    if(cand->pt() > K1_pTMax){

      //Move the other kaon to second position
      K2_pTMax  = K1_pTMax;
      K2_pT     = K1_pT;
      K2_eta    = K1_eta;
      K2_phi    = K1_phi;
      K2_energy = K1_energy;
      K2_p4     = K1_p4;
      K2_dxy    = K1_dxy;
      K2_dz     = K1_dz;
      K2_charge = K1_charge;

      K1_pTMax  = cand->pt();
      K1_pT     = cand->pt();
      K1_eta    = cand->eta();
      K1_phi    = cand->phi();
      K1_energy = cand->energy();
      K1_p4     = cand->p4();
      K1_dxy    = cand->dxy();
      K1_dz     = cand->dz();
      K1_charge = cand->charge();
    }
    else{
      K2_pTMax  = cand->pt();
      K2_pT     = cand->pt();
      K2_eta    = cand->eta();
      K2_phi    = cand->phi();
      K2_energy = cand->energy();
      K2_p4     = cand->p4();
      K2_dxy    = cand->dxy();
      K2_dz     = cand->dz();
      K2_charge = cand->charge();
    }

  }

  //Do NOT continue if you didn't find a pion
  if(K1_pTMax < 0. || K2_pTMax < 0.) return;
  _Nevents_isTwoKaons++;

  //std::cout << "Two kaons " << _Nevents_isTwoKaons << std::endl;

  //*************************************************************//
  //                                                             //
  //----------------------- Kaon isolation ----------------------//
  //                                                             //
  //*************************************************************//

  for(auto cand_iso = PFCandidates->begin(); cand_iso != PFCandidates->end(); ++cand_iso){
    float deltaR_K1 = sqrt((K1_eta-cand_iso->eta())*(K1_eta-cand_iso->eta())+(K1_phi-cand_iso->phi())*(K1_phi-cand_iso->phi()));
    float deltaR_K2 = sqrt((K2_eta-cand_iso->eta())*(K2_eta-cand_iso->eta())+(K2_phi-cand_iso->phi())*(K2_phi-cand_iso->phi()));

    if(deltaR_K1 <= 0.3 && deltaR_K1 >= 0.02) K1_sum_pT_03 += cand_iso->pt();
    if(deltaR_K1 <= 0.5 && deltaR_K1 >= 0.02) K1_sum_pT_05 += cand_iso->pt();

    if(deltaR_K2 <= 0.3 && deltaR_K2 >= 0.02) K2_sum_pT_03 += cand_iso->pt();
    if(deltaR_K2 <= 0.5 && deltaR_K2 >= 0.02) K2_sum_pT_05 += cand_iso->pt();

    if(cand_iso->charge() != 0 && (fabs(cand_iso->dxy()) >= 0.2 || fabs(cand_iso->dz()) >= 0.5) ) continue; // Requesting charged particles to come from PV
    if(deltaR_K1 <= 0.5 && deltaR_K1 >= 0.02) K1_sum_pT_05_ch += cand_iso->pt();
    if(deltaR_K2 <= 0.5 && deltaR_K2 >= 0.02) K2_sum_pT_05_ch += cand_iso->pt();
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

    if(corr_et < 20. || fabs(photon->eta()) > 2.5 || corr_et < eTphMax) continue;
    if(photon->hasPixelSeed()) continue;   //electron veto

    //std::cout << "photon" << photon->photonID("mvaPhoID-RunIIFall17-v1-wp90") << " " << photon->photonID("mvaPhoID-RunIIFall17-v1p1-wp90") << std::endl;

    if(photon->photonID("mvaPhoID-RunIIFall17-v1p1-wp90") == 0) continue;

    float abseta = fabs(photon->superCluster()->eta());
    float eA = effectiveAreas_ph_.getEffectiveArea(abseta);
    //photon_iso = (pfIso.sumChargedHadronPt + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_))/photon->et();
    if(photon->chargedHadronIso()/corr_et > 0.3 || photon->photonIso() > 4.) continue; //|| photon->trackIso() > 6
    ph_iso_ChargedHadron = photon->chargedHadronIso();
    ph_iso_NeutralHadron = photon->neutralHadronIso();
    ph_iso_Photon        = photon->photonIso();
    ph_iso_eArho         = eA*rho_;

    eTphMax = corr_et;

    ph_eT     = corr_et;
    ph_eta    = photon->eta();
    ph_etaSC  = photon->superCluster()->eta();
    ph_phi    = photon->phi();

    // Apply energy scale corrections to MC
    ph_energy = photon->energy(); //userFloat("ecalEnergyPostCorr");
    ph_p4     = photon->p4();  //* photon->userFloat("ecalEnergyPostCorr")/photon->energy();
    
    cand_photon_found = true;

    float deltapTMax = 10000.;
    const float deltaRMax = 0.3;
    int   gen_mother = 0;
    int   gen_ID = 0;

    is_photon_a_photon = false;
    is_photon_matched = false;

    if(!runningOnData_){
      for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
	float deltaR = sqrt((ph_eta-gen->eta())*(ph_eta-gen->eta())+(ph_phi-gen->phi())*(ph_phi-gen->phi()));
	float deltapT = fabs(ph_eT-gen->pt());

	if(deltaR > deltaRMax || deltapT > deltapTMax) continue;
	deltapTMax = deltapT;
	gen_ID = gen->pdgId();
	gen_mother = gen->mother()->pdgId();
      }
            
      if(gen_ID == 22) is_photon_a_photon = true;
      //if(gen_ID != 22) std::cout << "ph gen ID = " << gen_ID << std::endl;
      //if(gen_ID != 22 && fabs(gen_mother) == 24) std::cout << "ph gen ID when matched = " << gen_ID << std::endl;
      if(fabs(gen_mother) == 25) is_photon_matched = true;
    }

  }

  //Do not continue if there's no photons
  if(!cand_photon_found) return;
  _Nevents_isPhoton++;

  //std::cout << "Nphotons " << _Nevents_isPhoton << std::endl;

  //Only save events in a certain range
  _Phimass = (K1_p4 + K2_p4).M();
  //std::cout << "PhiMass " << _Phimass << std::endl;
  if(_Phimass < 0.1 || _Phimass > 5.) return;
  _Nevents_isPhimass++;

  _Hmass = (K1_p4 + K2_p4 + ph_p4).M();
  //std::cout << "HMass " << _Hmass << std::endl;
  if(_Hmass < 100. || _Hmass > 140.) return;
  _Nevents_isHmass++;

  //*************************************************************//
  //                                                             //
  //--------------------------- N-jets --------------------------//
  //                                                             //
  //*************************************************************//

  for (auto jet = slimmedJets->begin(); jet != slimmedJets->end(); ++jet){
    if(jet->pt() < 25.) continue;
    nJets_25++;
    if(jet->pt() < 30.) continue;
    nJets++;
  }

  mytree->Fill();
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
  mytree->Branch("nJets",&nJets);
  mytree->Branch("nJets_25",&nJets_25);
  mytree->Branch("met_pT",&met_pT);
  mytree->Branch("metpuppi_pT",&metpuppi_pT);

  mytree->Branch("K1_pT",&K1_pT);
  mytree->Branch("K1_eta",&K1_eta);
  mytree->Branch("K1_phi",&K1_phi);
  mytree->Branch("K1_energy",&K1_energy);
  mytree->Branch("K1_dxy",&K1_dxy);
  mytree->Branch("K1_dz",&K1_dz);
  mytree->Branch("K1_charge",&K1_charge);
  mytree->Branch("K1_sum_pT_03",&K1_sum_pT_03);
  mytree->Branch("K1_sum_pT_05",&K1_sum_pT_05);
  mytree->Branch("K1_sum_pT_05_ch",&K1_sum_pT_05_ch);

  mytree->Branch("K2_pT",&K2_pT);
  mytree->Branch("K2_eta",&K2_eta);
  mytree->Branch("K2_phi",&K2_phi);
  mytree->Branch("K2_energy",&K2_energy);
  mytree->Branch("K2_dxy",&K2_dxy);
  mytree->Branch("K2_dz",&K2_dz);
  mytree->Branch("K2_charge",&K2_charge);
  mytree->Branch("K2_sum_pT_03",&K2_sum_pT_03);
  mytree->Branch("K2_sum_pT_05",&K2_sum_pT_05);
  mytree->Branch("K2_sum_pT_05_ch",&K2_sum_pT_05_ch);

  mytree->Branch("photon_eT",&ph_eT);
  mytree->Branch("photon_eta",&ph_eta);
  mytree->Branch("photon_etaSC",&ph_etaSC);
  mytree->Branch("photon_phi",&ph_phi);
  mytree->Branch("photon_energy",&ph_energy);
  mytree->Branch("photon_iso_ChargedHadron",&ph_iso_ChargedHadron);
  mytree->Branch("photon_iso_NeutralHadron",&ph_iso_NeutralHadron);
  mytree->Branch("photon_iso_Photon",&ph_iso_Photon);
  mytree->Branch("photon_iso_eArho",&ph_iso_eArho);

  mytree->Branch("Phimass",&_Phimass);
  mytree->Branch("Hmass",&_Hmass);

  //Save MC info
  if(!runningOnData_){
    mytree->Branch("PU_Weight",&PU_Weight);
    mytree->Branch("MC_Weight",&MC_Weight);

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
  _Nevents_processed  = 0;
  _Nevents_isTwoKaons   = 0;
  _Nevents_isPhoton   = 0;
  _Nevents_isPhimass  = 0;
  _Nevents_isHmass  = 0;

  h_Events->Fill(0.5,_Nevents_processed);
  h_Events->Fill(1.5,_Nevents_isTwoKaons);
  h_Events->Fill(2.5,_Nevents_isPhoton);
  h_Events->Fill(3.5,_Nevents_isPhimass);
  h_Events->Fill(4.5,_Nevents_isHmass);

  h_Events->GetXaxis()->SetBinLabel(1,"Events triggered");
  h_Events->GetXaxis()->SetBinLabel(2,"Two kaons requested");
  h_Events->GetXaxis()->SetBinLabel(3,"Photon requested");
  h_Events->GetXaxis()->SetBinLabel(4,"Phi mass cut");
  h_Events->GetXaxis()->SetBinLabel(5,"Higgs mass cut");
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPhiGammaAnalysis);
