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
 
#include "HPhiGammaTwoProngsTriggerAnalysis.h"



// constructors and destructor
HPhiGammaTwoProngsTriggerAnalysis::HPhiGammaTwoProngsTriggerAnalysis(const edm::ParameterSet& iConfig) :
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


HPhiGammaTwoProngsTriggerAnalysis::~HPhiGammaTwoProngsTriggerAnalysis()
{
}

// ------------ method called for each event  ------------
void HPhiGammaTwoProngsTriggerAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  isIsoMuTrigger     = false;
  isTwoProngsTrigger = false;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for(unsigned int i = 0, n = triggerBits->size(); i < n; ++i){ //trigger forloop start
    if(!triggerBits->accept(i)) continue;
    std::string tmp_triggername = names.triggerName(i);

    if( tmp_triggername.find("HLT_IsoMu24") != std::string::npos ){ //Muon trigger
      isIsoMuTrigger = true;
      if (verbose) cout<<"IsoMu24 triggered"<<endl;
    }

    if( tmp_triggername.find("HLT_IsoMu24_TwoProngs35") != std::string::npos ){ //Photon trigger
      isTwoProngsTrigger = true; 
      if (verbose) cout<<"Mu24_TwoProngs35 triggered"<<endl;
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


  //*************************************************************//
  //                                                             //
  //---------------------------- Muons --------------------------//
  //                                                             //
  //*************************************************************//
  if (verbose) cout<<"Muons forloop start"<<endl;

  //First loop over muons
  for(std::vector<reco::Muon>::size_type firstMuIndex = 0; firstMuIndex < slimmedMuons->size();firstMuIndex ++){ //Muon first forloop start
      
      //refuse muons not passing over basic requirements
      if(slimmedMuons->at(firstMuIndex).pt() < 5. || !slimmedMuons->at(firstMuIndex).CutBasedIdTight || fabs(slimmedMuons->at(firstMuIndex).eta()) > 2.4 || fabs(slimmedMuons->at(firstMuIndex).muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(slimmedMuons->at(firstMuIndex).muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
      if(!slimmedMuons->at(firstMuIndex).PFIsoLoose) continue;
    
      //Second loop over muons
      for(std::vector<pat::Muon> ::size_type secondMuIndex = firstMuIndex + 1; secondMuIndex < slimmedMuons->size();secondMuIndex ++){ //Muon second forloop start
          
          //refuse muons not passing over basic requirements
          if(slimmedMuons->at(secondMuIndex).pt() < 5. || !slimmedMuons->at(secondMuIndex).CutBasedIdTight || fabs(slimmedMuons->at(secondMuIndex).eta()) > 2.4 || fabs(slimmedMuons->at(secondMuIndex).muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(slimmedMuons->at(secondMuIndex).muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
          if(!slimmedMuons->at(secondMuIndex).PFIsoLoose) continue;

          //at least one of the two muons must have pT > 25 GeV
          if(slimmedMuons->at(firstMuIndex).pt() < 25. && slimmedMuons->at(secondMuIndex).pt() < 25.) continue; 

          //take only muons with opposite charges
          if(slimmedMuons->at(firstMuIndex).charge() * slimmedMuons->at(secondMuIndex).charge() >= 0.) continue; 

          //take only muons pairs falling in an invariant mass range
          currentMuMuMass = (slimmedMuons->at(firstMuIndex).p4() + slimmedMuons->at(secondMuIndex).p4()).M();
          if(currentMuMuMass < 65. || currentMuMuMass > 115.) continue; //MuMu inv mass for Z

          //choose the pair with largest pT 
          currentMuMuPt = (slimmedMuons->at(firstMuIndex).p4() + slimmedMuons->at(secondMuIndex).p4()).pt();
          if(currentMuMuPt <= bestMuMuPt) continue; 

          //Save the variables of the current pair until a better one doesn't replace it
          bestMuMuPt   = currentMuMuPt;
          bestMuMuMass = currentMuMuMass;
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
 
 //*************************************************************//
  //                                                             //
  //--------------------------- N-tracks --------------------------//
  //                                                             //
  //*************************************************************//


  //daughters forloop variable
int firstCandCharge;
int secondCandCharge;
float firstCandPt;
float secondCandPt;
LorentzVector firstCand_p4;
LorentzVector secondCand_p4;
LorentzVector firstCand_p4_K;
LorentzVector secondCand_p4_K;
LorentzVector firstCand_p4_Pi;
LorentzVector secondCand_p4_Pi;
LorentzVector couple_p4;
LorentzVector couple_p4_K;
LorentzVector couple_p4_Pi;
LorentzVector best_firstCand_p4;
LorentzVector best_secondCand_p4;
LorentzVector best_couple_p4;
float bestCoupleOfTheEvent_pT = 0.;
float deltaR_K = 0.;
firstCandEnergy_K = 0.;
secondCandEnergy_K = 0.;
firstCandEnergy_Pi = 0.;
secondCandEnergy_Pi = 0.;
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
float PhiMass = 0.;
float RhoMass = 0.;
float kMass = 0.4937;
float PiMass = 0.13957;
float candPtMin = 1.;
isBestCoupleOfTheEvent_Found=false;
bool isPhi = false;
bool isRho = false;
_firstCandPt   = 0.;
_firstCandEta  = 0.;
_firstCandPhi  = 0.;
_secondCandPt  = 0.;
_secondCandEta = 0.;
_secondCandPhi = 0.;
_bestCouplePt  = 0.;
_bestCoupleEta = 0.;
_bestCouplePhi = 0.;
_MesonMass     = -1.;


      //-------------------------------------tracks forloop----------------------------
      if(verbose)cout<<"TRACKs:"<<endl;


  for(std::vector<pat::PackedCandidate>::size_type firstTrk_Index=0; firstTrk_Index < PFCandidates->size(); firstTrk_Index++){ //1ST LOOP STARTS

        //minimum apporach distance
        if ((PFCandidates->at(firstTrk_Index).dxy((&slimmedPV->at(0))->position())) >= 0.2 || PFCandidates->at(firstTrk_Index).dz((&slimmedPV->at(0))->position()) >= 0.5 ) continue;

        firstCandPt  = PFCandidates->at(firstTrk_Index).pt();  //extrapolate firstCand pt
        firstCandEta = PFCandidates->at(firstTrk_Index).eta(); //extrapolate firstCand eta
        firstCandPhi = PFCandidates->at(firstTrk_Index).phi(); //extrapolate firstCand phi

        if(firstCandPt < candPtMin) continue; //firstCand filter if pT < candPtMin

        for(std::vector<pat::PackedCandidate>::size_type secondTrk_Index=firstTrk_Index+1; secondTrk_Index < PFCandidates->size(); secondTrk_Index++){ //2ND LOOP STARTS

          //minimum apporach distance
          if ((PFCandidates->at(secondTrk_Index).dxy((&slimmedPV->at(0))->position())) >= 0.2 || PFCandidates->at(secondTrk_Index).dz((&slimmedPV->at(0))->position()) >= 0.5 ) continue;

          secondCandPt  = PFCandidates->at(secondTrk_Index).pt();  //extrapolate secondCand pt
          secondCandEta = PFCandidates->at(secondTrk_Index).eta(); //extrapolate secondCand eta
          secondCandPhi = PFCandidates->at(secondTrk_Index).phi(); //extrapolate secondCand phi

          //secondCand filter if both pT are < 10GeV
          if(secondCandPt < candPtMin) continue;
          if(firstCandPt < 10. && secondCandPt < 10.) continue; 
          //if(verbose) cout<<"tracks pT cut passed"<<endl; //fixme

          //third filter on deltaR ------------------------------------------------------------
          float deltaEta= firstCandEta - secondCandEta;

          float deltaPhi = fabs(firstCandPhi - secondCandPhi);  //phi folding 
          if (deltaPhi > 3.14) deltaPhi = 6.28 - deltaPhi;

          deltaR_K= sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
          //if(verbose) cout<<"deltaR = "<<deltaR_K<<endl; //fixme
          if(deltaR_K > 0.07) continue; //FIXME, it was 0.02
          //if(verbose) cout<<"deltaR cut passed"<<endl; //fixme


          //OPPOSITE CHARGE - FILTER ------------------------------------------------------------
          firstCandCharge  = PFCandidates->at(firstTrk_Index).charge(); //extrapolate firstCand charge
          secondCandCharge = PFCandidates->at(secondTrk_Index).charge(); //extrapolate secondCand charge
          if(firstCandCharge*secondCandCharge >= 0) continue; //choose only opposite charges
          if(verbose) cout<<"opposite charge cut passed"<<endl; //fixme

          //QUADRIMOMENTUM CALCULATION ------------------------------------------------------------
          firstCand_p4  = PFCandidates->at(firstTrk_Index).p4(); //extrapolate quadrimomentum
          secondCand_p4 = PFCandidates->at(secondTrk_Index).p4();

          firstCandPx   = firstCand_p4.px(); //extrapolate px, py, pz of the first candidate
          firstCandPy   = firstCand_p4.py();
          firstCandPz   = firstCand_p4.pz();
          secondCandPx  = secondCand_p4.px(); //extrapolate px, py, pz of the second candidate
          secondCandPy  = secondCand_p4.py();
          secondCandPz  = secondCand_p4.pz();

          //PIONS OR KAONS HYPOTHESIS -----------------------------------------------------------------------------------------------------------------------------------------------            
          firstCandEnergy_K   = sqrt(firstCandPx  * firstCandPx  + firstCandPy  * firstCandPy  + firstCandPz  * firstCandPz  + kMass  * kMass ); //Kaon hypothesis energy recalculation
          secondCandEnergy_K  = sqrt(secondCandPx * secondCandPx + secondCandPy * secondCandPy + secondCandPz * secondCandPz + kMass  * kMass ); //Kaon hypothesis energy recalculation
          firstCandEnergy_Pi  = sqrt(firstCandPx  * firstCandPx  + firstCandPy  * firstCandPy  + firstCandPz  * firstCandPz  + PiMass * PiMass); //Pion hypothesis energy recalculation
          secondCandEnergy_Pi = sqrt(secondCandPx * secondCandPx + secondCandPy * secondCandPy + secondCandPz * secondCandPz + PiMass * PiMass); //Pion hypothesis energy recalculation
          
          if (verbose) {
            cout<<"firstCandEnergy_K = "<<firstCandEnergy_K<<endl;
            cout<<"firstCandEnergy_Pi = "<<firstCandEnergy_Pi<<endl;
            cout<<"secondCandEnergy_K = "<<secondCandEnergy_K<<endl;
            cout<<"secondCandEnergy_Pi = "<<secondCandEnergy_Pi<<endl;
          }

          firstCand_p4_K   = firstCand_p4.SetE(firstCandEnergy_K); //Kaon hypothesis quadrimomentum correction
          secondCand_p4_K  = secondCand_p4.SetE(secondCandEnergy_K); //Kaon hypothesis quadrimomentum correction
          firstCand_p4_Pi  = firstCand_p4.SetE(firstCandEnergy_Pi); //Kaon hypothesis quadrimomentum correction
          secondCand_p4_Pi = secondCand_p4.SetE(secondCandEnergy_Pi); //Kaon hypothesis quadrimomentum correction

          couple_p4_K  = firstCand_p4_K  + secondCand_p4_K; //calculation of the couple-quadrimomentum after the correction
          couple_p4_Pi = firstCand_p4_Pi + secondCand_p4_Pi; //calculation of the couple-quadrimomentum after the correction
          
          if (verbose) {
            cout<<"KK pT = "<<couple_p4_K.pt()<<endl;
            cout<<"PiPi pT = "<<couple_p4_Pi.pt()<<endl;
          }

          //PT MAX OF THE EVENT - FILTER
          if(couple_p4_K.pt() <= bestCoupleOfTheEvent_pT) continue; //choose the couple with greatest pt
          bestCoupleOfTheEvent_pT = couple_p4_K.pt();       

         //PT CUT
          if(couple_p4_K.pt() < 38.) continue;
          if(verbose) cout<<"couplepT cut passed"<<endl; //fixme

          //PHI INV MASS - FILTER
          isPhi = false;
          isRho = false;

          PhiMass = (couple_p4_K).M(); //calculate inv mass of the Phi candidate  
          if (verbose) cout<<"mKK (before the meson mass selection) =  "<<PhiMass<<endl;
          if(PhiMass > 1. && PhiMass < 1.05) isPhi = true; //filter on Phi invariant mass  

          RhoMass = (couple_p4_Pi).M(); //calculate inv mass of the Rho candidate  
          if (verbose) cout<<"mPiPi (before the meson mass selection) =  "<<RhoMass<<endl;
          if(RhoMass > 0.5 && RhoMass < 1.) isRho = true; //filter on Rho invariant mass   

          if (!isPhi && !isRho) continue; //continue if the pair mass doesn't match any of the two mass hypothesis
          if(verbose) cout<<"mass range if the couple cut passed"<<endl; //fixme

          isBestCoupleOfTheEvent_Found = true;

          //Save if best pair has been found
          _isPhi                  = isPhi;
          _isRho                  = isRho;
          
          if(isPhi){
            best_firstCand_p4  = firstCand_p4_K; 
            best_secondCand_p4 = secondCand_p4_K;
            best_couple_p4     = couple_p4_K;
            _MesonMass         = PhiMass;
          }  
          if(isRho){
            best_firstCand_p4  = firstCand_p4_Pi; 
            best_secondCand_p4 = secondCand_p4_Pi;
            best_couple_p4     = couple_p4_Pi;
            _MesonMass         = RhoMass;

          }      
        
      } //2ND LOOP ENDS
  } //1ST LOOP ENDS


//DATAMEMBER SAVING
_firstCandPt   = best_firstCand_p4.pt();
_firstCandEta  = best_firstCand_p4.eta();
_firstCandPhi  = best_firstCand_p4.phi();
_secondCandPt  = best_secondCand_p4.pt();       
_secondCandEta = best_secondCand_p4.eta();
_secondCandPhi = best_secondCand_p4.phi();
_bestCouplePt  = best_couple_p4.pt();
_bestCoupleEta = best_couple_p4.eta();
_bestCouplePhi = best_couple_p4.phi();

if (verbose && !isBestCoupleOfTheEvent_Found) cout<<"best couple of the event not found!"<<endl;

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

//------------- ISOLATION -------------------------------------------------------------------------  
for(auto cand_iso = PFCandidates->begin(); cand_iso != PFCandidates->end(); ++cand_iso){ //ISOLATION FORLOOP START

  if(cand_iso->pt() < 0.5) continue; //do not consider tracks with pT < 500MeV

  float deltaPhi_K1 = fabs(_firstCandPhi-cand_iso->phi());  //phi folding 
  if (deltaPhi_K1 > 3.14) deltaPhi_K1 = 6.28 - deltaPhi_K1;

  float deltaR_K1 = sqrt((_firstCandEta-cand_iso->eta())*(_firstCandEta-cand_iso->eta()) + deltaPhi_K1*deltaPhi_K1);
  if(deltaR_K1 < 0.002) continue;

  float deltaPhi_K2 = fabs(_secondCandPhi-cand_iso->phi());  //phi folding  
  if (deltaPhi_K2 > 3.14) deltaPhi_K2 = 6.28 - deltaPhi_K2;

  float deltaR_K2 = sqrt((_secondCandEta-cand_iso->eta())*(_secondCandEta-cand_iso->eta()) + deltaPhi_K2*deltaPhi_K2);
  if(deltaR_K2 < 0.002) continue;

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
} //ISOLATION FORLOOP END
  
  //CANDIDATES SORTING
if(_firstCandPt < _secondCandPt)  //swap-values loop, in order to fill the tree with the candidate with max pt of the couple in firstCand branches  
  {                               //and one with min pt in secondCand branches
    float a,b,c,d;
    a = _firstCandPt;
    b = _firstCandEta;
    c = _firstCandPhi;
    d = firstCandEnergy;
    _firstCandPt    = _secondCandPt;
    _firstCandEta   = _secondCandEta;
    _firstCandPhi   = _secondCandPhi;
    firstCandEnergy = secondCandEnergy;       
    _secondCandPt    = a;
    _secondCandEta   = b;
    _secondCandPhi   = c;
    secondCandEnergy = d;
}

if (verbose) {
  cout <<"First candidate pT = "<<_firstCandPt<<endl;
  cout <<"Second candidate pT = "<<_secondCandPt<<endl;
}

//CUTS ON CANDIDATES PT
if(_firstCandPt < 20. || _secondCandPt < 5.) {
    cout<<"Final cut on candidates pT not passed, RETURN."<<endl;
    return;
}

//ISOLATION DATAMEMBER FOR TREE FILLING 
_iso_K1        = K1_sum_pT_05/_firstCandPt;
_iso_K2        = K2_sum_pT_05/_secondCandPt;
_iso_couple    = couple_sum_pT_05/_bestCouplePt;
_iso_K1_ch     = K1_sum_pT_05_ch/_firstCandPt;
_iso_K2_ch     = K2_sum_pT_05_ch/_secondCandPt;
_iso_couple_ch = couple_sum_pT_05_ch/_bestCouplePt;

//CUT ON PHI ISOLATION
if(_iso_couple_ch > 1.) {
  if(verbose) cout<<"No isolation cut passed, RETURN."<<endl;
  return;
}

  if(!isBestCoupleOfTheEvent_Found) {
    if (verbose) cout<<"RETURN: No best couple found."<<endl<<endl;
    return;
  }
 
  mytree->Fill();

}

//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//

void HPhiGammaTwoProngsTriggerAnalysis::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");

  mytree->Branch("nPV",&nPV);
  mytree->Branch("isIsoMuTrigger",&isIsoMuTrigger);
  mytree->Branch("isTwoProngsTrigger",&isTwoProngsTrigger);

  //Save run number info when running on data
  if(runningOnData_){
    mytree->Branch("run_number",&run_number);
    }

  mytree->Branch("isBestCoupleOfTheEvent_Found",&isBestCoupleOfTheEvent_Found);

  mytree->Branch("MuMuPt",&bestMuMuPt);
  mytree->Branch("MuMuMass",&bestMuMuMass);

  mytree->Branch("firstCandPt",&_firstCandPt);
  mytree->Branch("firstCandEta",&_firstCandEta);
  mytree->Branch("firstCandPhi",&_firstCandPhi);
  mytree->Branch("secondCandPt",&_secondCandPt);
  mytree->Branch("secondCandEta",&_secondCandEta);
  mytree->Branch("secondCandPhi",&_secondCandPhi);
  mytree->Branch("bestCouplePt",&_bestCouplePt);
  mytree->Branch("bestCoupleEta",&_bestCoupleEta);
  mytree->Branch("bestCouplePhi",&_bestCouplePhi);
  mytree->Branch("isPhi",&_isPhi);
  mytree->Branch("isRho",&_isRho);

  mytree->Branch("firstCandEnergy",&firstCandEnergy);
  mytree->Branch("secondCandEnergy",&secondCandEnergy);

  mytree->Branch("MesonMass",&_MesonMass);
}
  
void HPhiGammaTwoProngsTriggerAnalysis::beginJob()
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

void HPhiGammaTwoProngsTriggerAnalysis::endJob() 
{
  h_Events->Fill(0.5,_Nevents_processed);
  
  h_Events->GetXaxis()->SetBinLabel(1,"Events processed");
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPhiGammaTwoProngsTriggerAnalysis);

