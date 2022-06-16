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
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
  
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
 
#include "HPhiGammaTriggerAnalysis.h"



// constructors and destructor
HPhiGammaTriggerAnalysis::HPhiGammaTriggerAnalysis(const edm::ParameterSet& iConfig) :
  runningOnData_(iConfig.getParameter<bool>("runningOnData")),
  verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose")),
  effectiveAreas_el_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_el")).fullPath() ),
  effectiveAreas_ph_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_ph")).fullPath() )
{
  packedPFCandidatesToken_            = consumes<std::vector<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates")); 
  slimmedMuonsToken_                  = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
  photonsMiniAODToken_                = consumes<std::vector<pat::Photon> > (edm::InputTag("slimmedPhotons"));
  offlineSlimmedPrimaryVerticesToken_ = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));  
  offlineBeamSpotToken_               = consumes<reco::BeamSpot> (edm::InputTag("offlineBeamSpot"));
  pileupSummaryToken_                 = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
  triggerBitsToken_                   = consumes<edm::TriggerResults> (edm::InputTag("TriggerResults","","HLT"));
  rhoToken_                           = consumes<double> (iConfig.getParameter <edm::InputTag>("rho"));


  h_Events = fs->make<TH1F>("h_Events", "Event counting in different steps", 8, 0., 8.);

  _Nevents_processed  = 0;
  _Nevents_isPhoton   = 0;

  debug=false;  //DEBUG datamember 
  verbose=false; 

  h_pileup   = fs->make<TH1F>("pileup", "pileup", 75,0,75);

  create_trees();
}


HPhiGammaTriggerAnalysis::~HPhiGammaTriggerAnalysis()
{
}

// ------------ method called for each event  ------------
void HPhiGammaTriggerAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  edm::Handle<std::vector<pat::PackedCandidate>  > PFCandidates;
  iEvent.getByToken(packedPFCandidatesToken_, PFCandidates);

  edm::Handle<std::vector<pat::Muon>  > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_, slimmedMuons);

  edm::Handle<std::vector<pat::Photon> > slimmedPhotons;
  iEvent.getByToken(photonsMiniAODToken_,slimmedPhotons);

  edm::Handle<std::vector<reco::Vertex > > slimmedPV;
  iEvent.getByToken(offlineSlimmedPrimaryVerticesToken_, slimmedPV);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsToken_, triggerBits);

  _Nevents_processed++;

  //Retrieve the run number
  if(runningOnData_){
    run_number = iEvent.id().run();
  }

  //some prints
  if (verbose){
    cout<<"############ Event n."<< _Nevents_processed<<" ############"<<endl;
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

    // Fill histogram with PU distribution
    h_pileup->Fill(npT);
  }


  //*************************************************************//
  //                                                             //
  //--------------------------- Trigger -------------------------//
  //                                                             //
  //*************************************************************//
  
  if (verbose) cout<<"Checking triggers..."<<endl;

  //Examine the trigger information
  isIsoMuTrigger = false;
  isPhotonTrigger = false;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for(unsigned int i = 0, n = triggerBits->size(); i < n; ++i){ //trigger forloop start
    if(!triggerBits->accept(i)) continue;
    std::string tmp_triggername = names.triggerName(i);

    if( tmp_triggername.find("HLT_IsoMu24") != std::string::npos ){ //Muon trigger
      isIsoMuTrigger = true;
      if (verbose) cout<<"IsoMu24 triggered"<<endl;
    }

    if( tmp_triggername.find("HLT_Mu17_Photon30_IsoCaloId") != std::string::npos ){ //Photon trigger
      isPhotonTrigger = true; 
      if (verbose) cout<<"Mu17_Photon30 triggered"<<endl;
    }

  } //trigger forloop end
  
  //RETURN if muon trigger does not switch on
  if(!isIsoMuTrigger) {
  if (verbose) cout<<"RETURN: IsoMu24 not triggered."<<endl<<endl;
  return; //Return only if there're not any muons
  }

  //*************************************************************//
  //                                                             //
  //------------------ Variable initialization ------------------//
  //                                                             //
  //*************************************************************//

  isBestMuMu_Found = false;
  float currentMuMuPt = -1.;
  bestMuMuPt = -1.;
  float currentMuMuMass = -1.;
  bestMuMuMass = -1.;

  nPhotonsOverSelection = 0;
  nPhotonsChosen   = 0;
  
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

  eTphMax  = -1000.;


   //*************************************************************//
  //                                                             //
  //---------------------------- Muons --------------------------//
  //                                                             //
  //*************************************************************//
  if (verbose) cout<<"Muons forloop start"<<endl;

  //for(auto firstMu= slimmedMuons->begin(); firstMu!= slimmedMuons->end(); ++firstMu){ //Muon first forloop start
  for(std::vector<pat::Muon>::size_type firstMuIndex = 0; firstMuIndex < slimmedMuons->size();firstMuIndex ++){ //Muon first forloop start
      
      if(slimmedMuons->at(firstMuIndex).pt() < 5. || !slimmedMuons->at(firstMuIndex).CutBasedIdMedium || fabs(slimmedMuons->at(firstMuIndex).eta()) > 2.4 || fabs(slimmedMuons->at(firstMuIndex).muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(slimmedMuons->at(firstMuIndex).muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
      if(!slimmedMuons->at(firstMuIndex).PFIsoLoose) continue;
    
      //for(auto secondMu= slimmedMuons->begin(); secondMu!= slimmedMuons->end(); ++secondMu){//Muon second forloop start
      for(std::vector<pat::Muon> ::size_type secondMuIndex = firstMuIndex + 1; secondMuIndex < slimmedMuons->size();secondMuIndex ++){ //Muon second forloop start
          
          if(slimmedMuons->at(secondMuIndex).pt() < 5. || !slimmedMuons->at(secondMuIndex).CutBasedIdMedium || fabs(slimmedMuons->at(secondMuIndex).eta()) > 2.4 || fabs(slimmedMuons->at(secondMuIndex).muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(slimmedMuons->at(secondMuIndex).muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;

          //at least one of the two muons must have pT > 25 GeV
          if(slimmedMuons->at(firstMuIndex).pt() < 25. && slimmedMuons->at(secondMuIndex).pt() < 25.) continue; 

          if(slimmedMuons->at(firstMuIndex).charge()*slimmedMuons->at(secondMuIndex).charge() >= 0.) continue; //take only muons with opposite charges
          currentMuMuMass = (slimmedMuons->at(firstMuIndex).p4() + slimmedMuons->at(secondMuIndex).p4()).M();
          if(currentMuMuMass < 20. || currentMuMuMass > 120.) continue; //MuMu inv mass for Z

          currentMuMuPt = (slimmedMuons->at(firstMuIndex).p4() + slimmedMuons->at(secondMuIndex).p4()).pt(); 
          if(currentMuMuPt <= bestMuMuPt) continue; //choose the pair with largest pT
          bestMuMuPt = currentMuMuPt;
          bestMuMuMass = currentMuMuMass;
          isBestMuMu_Found = true;

      } //Muon second forloop end
    }//Muon first forloop end

  if(!isBestMuMu_Found) { 
    if (verbose) cout<<"RETURN: No Z->mumu found."<<endl<<endl;
    return;
  }
  
  if(isBestMuMu_Found && verbose){
    cout<<"Muon pair found, with pT = "<<bestMuMuPt<<" and inv mass = "<<bestMuMuMass<<endl;
  } 

  //*************************************************************//
  //                                                             //
  //--------------------------- Photons -------------------------//
  //                                                             //
  //*************************************************************//
  if (verbose) cout<<"Photons forloop start"<<endl;

  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  bool cand_photon_found = false;
  float corr_et = -1.;

  for(auto photon = slimmedPhotons->begin(); photon != slimmedPhotons->end(); ++photon){ //Photons forloop start

    corr_et = photon->et(); // * photon->userFloat("ecalEnergyPostCorr")/photon->energy();

    //std::cout << "photon et " << corr_et << std::endl;

    if(corr_et < 35. || fabs(photon->eta()) > 2.5) continue;
    if(verbose) cout<<"corr_et = "<<corr_et<<endl; //FIXME
    if(photon->hasPixelSeed()) continue;   //electron veto

    //std::cout << "photon" << photon->photonID("mvaPhoID-RunIIFall17-v1-wp90") << " " << photon->photonID("mvaPhoID-RunIIFall17-v1p1-wp90") << std::endl;
    if(photon->photonID("mvaPhoID-RunIIFall17-v1p1-wp90") == 0) continue;
    if(verbose) cout<<"photonID = "<<corr_et<<endl; //FIXME

    float abseta = fabs(photon->superCluster()->eta());
    float eA = effectiveAreas_ph_.getEffectiveArea(abseta);
    //photon_iso = (pfIso.sumChargedHadronPt + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_))/photon->et();

    //if(photon->chargedHadronIso()/corr_et > 0.3 || photon->photonIso() > 4.) continue; //|| photon->trackIso() > 6

    nPhotonsOverSelection++;
    cout<<"nPhotonsOverSelection = "<<nPhotonsOverSelection<<endl; //FIXME
    if(corr_et < eTphMax) continue; //choose the best photon

    
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
  } //Photon forloops end

  //Do not continue if there are no photons
  if(!cand_photon_found) {
    if (verbose) cout<<"RETURN: No best photon found."<<endl<<endl;
    return;
  }
  if(cand_photon_found){
  cout<<"Photon found, with eT = "<<ph_eT<<endl;
  _Nevents_isPhoton++;
  }
  if(verbose) cout<<_Nevents_isPhoton<<" best photons chosen after offline selection"<<endl;

  //  std::cout << "Nphotons " << _Nevents_isPhoton << std::endl;

  if (verbose) cout<<"endl";
  if (verbose) cout<<"cand_photon_found = "<<cand_photon_found<<endl;
  mytree->Fill();
}

//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//

void HPhiGammaTriggerAnalysis::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");

  mytree->Branch("nPV",&nPV);
  mytree->Branch("isIsoMuTrigger",&isIsoMuTrigger);
  mytree->Branch("isPhotonTrigger",&isPhotonTrigger);

  //Save run number info when running on data
  if(runningOnData_){
    mytree->Branch("run_number",&run_number);
  }

  mytree->Branch("nPhotonsOverSelection",&nPhotonsOverSelection);
  mytree->Branch("nPhotonsChosen",&nPhotonsChosen);
  mytree->Branch("nEvents_isPhoton",&_Nevents_isPhoton);
  mytree->Branch("cand_photon_found",&cand_photon_found);

  mytree->Branch("MuMuPt",&bestMuMuPt);
  mytree->Branch("MuMuMass",&bestMuMuMass);
  mytree->Branch("photon_eT",&ph_eT);
  mytree->Branch("photon_eta",&ph_eta);
  mytree->Branch("photon_etaSC",&ph_etaSC);
  mytree->Branch("photon_phi",&ph_phi);
  mytree->Branch("photon_energy",&ph_energy);
  mytree->Branch("photon_iso_ChargedHadron",&ph_iso_ChargedHadron);
  mytree->Branch("photon_iso_NeutralHadron",&ph_iso_NeutralHadron);
  mytree->Branch("photon_iso_Photon",&ph_iso_Photon);
  mytree->Branch("photon_iso_eArho",&ph_iso_eArho);
}
  
void HPhiGammaTriggerAnalysis::beginJob()
{
  //Flag for PileUp reweighting
  if (!runningOnData_){ // PU reweighting for 2017
   Lumiweights_ = edm::LumiReWeighting("MCpileUp_2018_25ns_UltraLegacy_PoissonOOTPU.root", "MyDataPileupHistogram.root", "pileup", "pileup");
  }
}


//*************************************************************//
//                                                             //
//------------------- Fill event loss histos ------------------//
//                                                             //
//*************************************************************//

void HPhiGammaTriggerAnalysis::endJob() 
{
  h_Events->Fill(0.5,_Nevents_processed);
  h_Events->Fill(1.5,_Nevents_isPhoton);
  
  h_Events->GetXaxis()->SetBinLabel(1,"Events processed");
  h_Events->GetXaxis()->SetBinLabel(2,"Best photon selection");
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPhiGammaTriggerAnalysis);

