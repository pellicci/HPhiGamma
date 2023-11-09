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

#include <Math/PxPyPzE4D.h>
#include <TRandom3.h>


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
//---------------------------- Muons --------------------------//
//                                                             //
//*************************************************************//
if (verbose)
  cout << "Muons forloop start" << endl;

// Variable declaration
bestMuMuPt   = -1.;
bestMuMuMass = -1.;
LorentzVector tagMuP4;
LorentzVector probeMuP4;
isBestMuMuFound = false;
bool hasProbeMuFiredTrigger = false;
_hasProbeMuFiredTrigger = false;
bool hasTagMuFiredTrigger = false;
_hasTagMuFiredTrigger = false;

// First loop over muons: first muon
for (std::vector<reco::Muon>::size_type firstMuIndex = 0; firstMuIndex < slimmedMuons->size(); firstMuIndex++) {
  const reco::Muon& firstMuon = slimmedMuons->at(firstMuIndex);

  // Select the first muon with pT > 10 GeV
  if (firstMuon.pt() <= 10.) {
    if (verbose) {
      cout << "First muon rejected due to low pT: " << firstMuon.pt() << " GeV." << endl;
    }
    continue;
  }

  // Second loop over muons: second muon
  for (std::vector<reco::Muon>::size_type secondMuIndex = firstMuIndex + 1; secondMuIndex < slimmedMuons->size(); secondMuIndex++) {
    const reco::Muon& secondMuon = slimmedMuons->at(secondMuIndex);

    // Select the second muon with pT > 10 GeV
    if (secondMuon.pt() <= 10.) {
      if (verbose) {
        cout << "Second muon rejected due to low pT: " << secondMuon.pt() << " GeV." << endl;
      }
      continue;
    }

    // Select only pairs with opposite charges
    if (firstMuon.charge() * secondMuon.charge() >= 0.) {
      if (verbose) {
        cout << "Muons with the same charge rejected." << endl;
      }
      continue;
    }

    // Calculate the invariant mass of the pair
    float currentMuMuMass = (firstMuon.p4() + secondMuon.p4()).M();

    // Select only pairs with invariant mass between 60 and 120 GeV
    if (currentMuMuMass < 60. || currentMuMuMass > 120.) {
      if (verbose) {
        cout << "Muon pair rejected due to invalid invariant mass: " << currentMuMuMass << " GeV." << endl;
      }
      continue;
    }

    // Calculate the deltaR between the tracks
    float deltaPhi = fabs(firstMuon.phi() - secondMuon.phi());
    if (deltaPhi > M_PI) deltaPhi = 2 * M_PI - deltaPhi;

    float deltaEta = fabs(firstMuon.eta() - secondMuon.eta());

    float deltaR   = sqrt(deltaPhi * deltaPhi + deltaEta * deltaEta);

    // Select only pairs with deltaR > 0.4
    if (deltaR <= 0.4) {
      if (verbose) {
        cout << "Muon pair rejected due to small deltaR: " << deltaR << endl;
      }
      continue;
    }

    // Select the pair with the highest pT
    float currentMuMuPt = (firstMuon.p4() + secondMuon.p4()).Pt();
    if (currentMuMuPt <= bestMuMuPt) {
      if (verbose) {
        cout << "Muon pair rejected due to lower pT: " << currentMuMuPt << " GeV." << endl;
      }
      continue;
    }

    // Randomly assign one muon as tag and the other as probe
    bool isFirstTag = (gRandom->Rndm() < 0.5);
    const reco::Muon* tagMu;
    const reco::Muon* probeMu;
    int tagIndex;
    int probeIndex;

    if (isFirstTag) {
      tagMuP4    = firstMuon.p4();
      probeMuP4  = secondMuon.p4();
      tagMu      = &firstMuon;
      probeMu    = &secondMuon;
      tagIndex   = firstMuIndex;
      probeIndex = secondMuIndex;

    } else {
      tagMuP4    = secondMuon.p4();
      probeMuP4  = firstMuon.p4();
      tagMu      = &secondMuon;
      probeMu    = &firstMuon;
      tagIndex   = secondMuIndex;
      probeIndex = firstMuIndex;
    }

    // Check the requirements on the tag muon
    if (!slimmedMuons->at(tagIndex).CutBasedIdMedium || fabs(tagMu->eta()) > 2.4 || fabs(tagMu->muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(tagMu->muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5 || !slimmedMuons->at(tagIndex).PFIsoLoose) {
      if (verbose) {
        cout << "Tag muon failed the requirements." << endl;
      }
      continue;
    }

    // Check only the pT requirement on the probe muon (pT > 35)
    if (probeMuP4.Pt() <= 35. || !slimmedMuons->at(probeIndex).CutBasedIdMedium|| fabs(probeMu->eta()) > 2.1 || fabs(probeMu->muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(probeMu->muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) {
      if (verbose) {
        cout << "Probe muon failed the requirements." << endl;
      }
      continue;
    }

    // Verify if the probe muon fired HLT_IsoMu24
    hasProbeMuFiredTrigger = slimmedMuons->at(probeIndex).triggered("HLT_IsoMu24_v*");
    hasTagMuFiredTrigger   = slimmedMuons->at(tagIndex).triggered("HLT_IsoMu24_v*");

    // Save the variables of the current pair
    bestMuMuPt   = currentMuMuPt;
    bestMuMuMass = currentMuMuMass;
    tagMuPt      = tagMuP4.Pt();
    tagMuEta     = tagMuP4.Eta();
    tagMuPhi     = tagMuP4.Phi();
    probeMuPt    = probeMuP4.Pt();
    probeMuEta   = probeMuP4.Eta();
    probeMuPhi   = probeMuP4.Phi();

    isBestMuMuFound = true;

    // Exit the loops since a valid pair has been found
    break;
  }

  // Exit the first loop if a valid pair has been found
  if (isBestMuMuFound)
    break;
}

if (!isBestMuMuFound) {
  if (verbose) cout << "RETURN: No valid muon pair found." << endl << endl;
  return;
}

_nEvents_ZmumuFound++;

if (hasProbeMuFiredTrigger){
  _hasProbeMuFiredTrigger = hasProbeMuFiredTrigger;
}

if (hasTagMuFiredTrigger){
  _hasTagMuFiredTrigger = hasTagMuFiredTrigger;
}

if (isBestMuMuFound && verbose) {
  cout << "Muon pair found with pT = " << bestMuMuPt << " and inv mass = " << bestMuMuMass << endl;
  cout << "Tag muon: pT = " << tagMuPt << ", eta = " << tagMuEta << ", phi = " << tagMuPhi << endl;
  cout << "Probe muon: pT = " << probeMuPt << ", eta = " << probeMuEta << ", phi = " << probeMuPhi << endl;
}


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
  mytree->Branch("hasProbeMuFiredTrigger",&_hasProbeMuFiredTrigger);
  mytree->Branch("hasTagMuFiredTrigger",&_hasTagMuFiredTrigger);

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

