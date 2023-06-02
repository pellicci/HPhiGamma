//ROOT includes
#include <TH1F.h>
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
#include "DataFormats/Math/interface/deltaR.h"
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
 
#include "HPhiGammaIsoEfficiencyAnalysis.h"



// constructors and destructor
HPhiGammaIsoEfficiencyAnalysis::HPhiGammaIsoEfficiencyAnalysis(const edm::ParameterSet& iConfig) :
  runningOnData_(iConfig.getParameter<bool>("runningOnData")),
  verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose")),
  effectiveAreas_el_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_el")).fullPath() ),
  effectiveAreas_ph_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_ph")).fullPath() )
{
  packedPFCandidatesToken_            = consumes<std::vector<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates")); 
  slimmedMuonsToken_                  = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
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
  _nEvents_ZmumuFound = 0;

  debug=false;  //DEBUG datamember 
  verbose=false; 

  h_pileup   = fs->make<TH1F>("pileup", "pileup", 75,0,75);

  create_trees();
}


HPhiGammaIsoEfficiencyAnalysis::~HPhiGammaIsoEfficiencyAnalysis()
{
}

// ------------ method called for each event  ------------
void HPhiGammaIsoEfficiencyAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  

  edm::Handle<std::vector<pat::PackedCandidate>  > PFCandidates;
  iEvent.getByToken(packedPFCandidatesToken_, PFCandidates);

  edm::Handle<std::vector<pat::Muon>  > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_, slimmedMuons);

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
    if (debug) cout << "PU weight = "<<PU_Weight<<endl;

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

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for(unsigned int i = 0, n = triggerBits->size(); i < n; ++i){ //trigger forloop start
    if(!triggerBits->accept(i)) continue;
    std::string tmp_triggername = names.triggerName(i);

    if( tmp_triggername.find("HLT_IsoMu24") != std::string::npos ){ //Muon trigger
      isIsoMuTrigger = true;
      if (verbose) cout<<"IsoMu24 triggered"<<endl;
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
  LorentzVector tagMuP4;
  LorentzVector probeMuP4;


  //*************************************************************//
  //                                                             //
  //---------------------------- Muons --------------------------//
  //                                                             //
  //*************************************************************//
  if (verbose) cout<<"Muons forloop start"<<endl;

  //First loop over muons: TAG muon
  for(std::vector<reco::Muon>::size_type tagMuIndex = 0; tagMuIndex < slimmedMuons->size();tagMuIndex ++){ //Muon first forloop start
      
      //refuse muons not passing over basic requirements
      if(slimmedMuons->at(tagMuIndex).pt() < 25. || !slimmedMuons->at(tagMuIndex).CutBasedIdMedium || fabs(slimmedMuons->at(tagMuIndex).eta()) > 2.4 || fabs(slimmedMuons->at(tagMuIndex).muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(slimmedMuons->at(tagMuIndex).muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
      if(!slimmedMuons->at(tagMuIndex).PFIsoLoose) continue;
    
      //Second loop over muons: PROBE muon
      for(std::vector<pat::Muon> ::size_type probeMuIndex = tagMuIndex + 1; probeMuIndex < slimmedMuons->size();probeMuIndex ++){ //Muon second forloop start
          
          //refuse muons not passing over basic requirements
        if(slimmedMuons->at(probeMuIndex).pt() < 10. || !slimmedMuons->at(probeMuIndex).CutBasedIdMedium || fabs(slimmedMuons->at(probeMuIndex).eta()) > 2.4 || fabs(slimmedMuons->at(probeMuIndex).muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(slimmedMuons->at(probeMuIndex).muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue; 

          //take only muons with opposite charges
          if(slimmedMuons->at(tagMuIndex).charge() * slimmedMuons->at(probeMuIndex).charge() >= 0.) continue; 

          //take only muons pairs falling in an invariant mass range
          currentMuMuMass = (slimmedMuons->at(tagMuIndex).p4() + slimmedMuons->at(probeMuIndex).p4()).M();
          if(currentMuMuMass < 60. || currentMuMuMass > 120.) continue; //MuMu inv mass for Z

          float tagPhi = slimmedMuons->at(tagMuIndex).phi();
          float tagEta = slimmedMuons->at(tagMuIndex).eta();
          float proPhi = slimmedMuons->at(probeMuIndex).phi();
          float proEta = slimmedMuons->at(probeMuIndex).eta();

          float deltaPhi = fabs(tagPhi - proPhi);
          if (deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;

          float deltaEta = fabs(tagEta - proEta);

          float deltaR = sqrt(deltaPhi * deltaPhi + deltaEta * deltaEta);

          if (deltaR < 0.4) continue; //consider only muons with high deltaR to compute the isolation

          //choose the pair with largest pT 
          currentMuMuPt = (slimmedMuons->at(tagMuIndex).p4() + slimmedMuons->at(probeMuIndex).p4()).pt();
          if(currentMuMuPt <= bestMuMuPt) continue; 

          //Save the variables of the current pair until a better one doesn't replace it
          bestMuMuPt   = currentMuMuPt;
          bestMuMuMass = currentMuMuMass;
          tagMuP4      = slimmedMuons->at(tagMuIndex).p4();
          probeMuP4    = slimmedMuons->at(probeMuIndex).p4();
          tagMuEta     = tagEta;
          tagMuPhi     = tagPhi;
          probeMuEta   = proEta;
          probeMuPhi   = proPhi;
          isBestMuMu_Found = true;

      } //Muon second forloop end
    }//Muon first forloop end

  if(!isBestMuMu_Found) { 
    if (verbose) cout<<"RETURN: No Z->mumu found."<<endl<<endl;
    return;
  }

  _nEvents_ZmumuFound++;

  if(isBestMuMu_Found && verbose){
    cout<<"Muon pair found, with pT = "<<bestMuMuPt<<" and inv mass = "<<bestMuMuMass<<endl;
  } 

  //return if probe pT < 35 GeV (which is the HLT threshold on TwoProngs pT)
  if (probeMuP4.pt() < 35.) return;
 
  //*************************************************************//
  //                                                             //
  //--------------------------- Isolation -----------------------//
  //                                                             //
  //*************************************************************//


//PROBE MUON ISOLATION
sum_pT_05    = 0.;
sum_pT_05_ch = 0.;
_iso    = 0.;
_iso_ch = 0.;


//------------- ISOLATION -------------------------------------------------------------------------  
for(auto cand_iso = PFCandidates->begin(); cand_iso != PFCandidates->end(); ++cand_iso){ //ISOLATION FORLOOP START

  if(cand_iso->pt() < 0.5) continue; //do not consider tracks with pT < 500MeV
  //cout << "particle charge = "<<cand_iso->charge()<<endl;

  float deltaPhi = fabs(probeMuPhi - cand_iso->phi());  //phi folding 
  if (deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;

  float deltaR = sqrt((probeMuEta - cand_iso->eta())*(probeMuEta - cand_iso->eta()) + deltaPhi * deltaPhi);
  if(deltaR < 0.001) continue;

  if(deltaR <= 0.3) sum_pT_05 += cand_iso->pt(); //iso from all particles (neutral + charged)

  if(cand_iso->charge() == 0) continue;

  if(fabs(cand_iso->dxy()) >= 0.2 || fabs(cand_iso->dz()) >= 0.5) continue; // Requesting charged particles to come from PV  if(deltaR <= 0.3) sum_pT_05_ch += cand_iso->pt();
  if(deltaR <= 0.3) sum_pT_05_ch += cand_iso->pt();

} //ISOLATION FORLOOP END
  

//DATAMEMBER FOR TREE FILLING 

tagMuPt   = tagMuP4.pt();
probeMuPt = probeMuP4.pt();

_iso    = probeMuPt / (sum_pT_05 + probeMuPt);
_iso_ch = probeMuPt / (sum_pT_05_ch + probeMuPt);
 
  mytree->Fill();

}

//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//

void HPhiGammaIsoEfficiencyAnalysis::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");

  mytree->Branch("nPV",&nPV);
  mytree->Branch("isIsoMuTrigger",&isIsoMuTrigger);

  //Save run number info when running on data
  if(runningOnData_){
    mytree->Branch("run_number",&run_number);
    }

  mytree->Branch("tagMuPhi",&tagMuPhi);
  mytree->Branch("probeMuPhi",&probeMuPhi);

  mytree->Branch("tagMuEta",&tagMuEta);
  mytree->Branch("probeMuEta",&probeMuEta);
  
  mytree->Branch("tagMuPt",&tagMuPt);
  mytree->Branch("probeMuPt",&probeMuPt);

  mytree->Branch("MuMuPt",&bestMuMuPt);
  mytree->Branch("MuMuMass",&bestMuMuMass);

  mytree->Branch("sum_pT_05",&sum_pT_05);
  mytree->Branch("sum_pT_05_ch",&sum_pT_05_ch);
  
  mytree->Branch("iso",&_iso);
  mytree->Branch("iso_ch",&_iso_ch);

  //Save MC info
  if(!runningOnData_){ //NO INFO FOR DATA
    mytree->Branch("PU_Weight",&PU_Weight);
  }
}
  
void HPhiGammaIsoEfficiencyAnalysis::beginJob()
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

void HPhiGammaIsoEfficiencyAnalysis::endJob() 
{
  h_Events->Fill(0.5,_Nevents_processed);
  
  h_Events->GetXaxis()->SetBinLabel(1,"Events processed");
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPhiGammaIsoEfficiencyAnalysis);

