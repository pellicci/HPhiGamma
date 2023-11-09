#define _USE_MATH_DEFINES
#include <cmath> 
#include <iostream>

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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

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
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//Parton distribution and QCD scale variations stuff 
#include "FWCore/Framework/interface/Run.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h" //LHE reader
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h" //LHE reader

//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h" //JEC uncertainties
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h" //JEC uncertainties
//#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

//JEC uncertainties
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
//#include "CondFormats/JetMETObjects/interface/JetResolution.h"
//#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//#include "CondFormats/DataRecord/interface/JetCorrectionsRecord.h"

//#include "CondFormats/JetMETObjects/interface/JetCorrectionsRecord.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

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
  LHEEventProduct_                    = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  triggerObjectsToken_                = consumes<std::vector<pat::TriggerObjectStandAlone> > (edm::InputTag("slimmedPatTrigger","","PAT")); //triggObj
  //packedGenParticlesToken_            = consumes<std::vector<pat::GenParticle>>(edm::InputTag("packedGenParticles", "", "PAT"));
  

  h_Events         = fs->make<TH1F>("h_Events", "Event counting in different steps", 8, 0., 8.);
  h_TriggerFilters = fs->make<TH1F>("h_TriggerFilters", "Trigger filters counting in different steps", 11, 0., 11.);

  _Nevents_processed             = 0;
  _Nevents_triggered             = 0;
  _Nevents_isTwoKaons            = 0;
  _Nevents_isPhoton              = 0;
  _Nevents_HiggsMatched          = 0;
  _Nevents_HiggsNotMatched       = 0;
  _Nevents_MesonPtNotMatched     = 0;
  _Nevents_bestCoupleFound       = 0;
  _Nevents_candPtFilter          = 0;
  _Nevents_coupleIsolationFilter = 0;
  _Nevents_VBFVeto               = 0;

  _NeventsFilter1                = 0;
  _NeventsFilter2                = 0;
  _NeventsFilter3                = 0;
  _NeventsFilter4                = 0;
  _NeventsFilter5                = 0;
  _NeventsFilter6                = 0;
  _NeventsFilter7                = 0;
  _NeventsFilter8                = 0;
  _NeventsFilter9                = 0;

  nKK_found                      = 0;
  nKK_notFound                   = 0;


  debug=false;  //DEBUG datamember 
  verbose=false;
  trig_verbose=true; 
  isTriggerStudies=true;

  h_pileup   = fs->make<TH1F>("pileup", "pileup", 75,0,75);

  create_trees();
}


HPhiGammaAnalysis::~HPhiGammaAnalysis()
{
}

//--------- Function to match reco object with trigger object ---------//
namespace{
  std::vector<const pat::TriggerObjectStandAlone*> getMatchedObjs(const float photon_eta,const float photon_phi,const float meson_eta,const float meson_phi,const std::vector<pat::TriggerObjectStandAlone>& trigObjs,const float maxDeltaR=0.1)
  {
    std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
    const float maxDR2 = maxDeltaR*maxDeltaR;
    int nTrigObj = 0;
    for(auto& trigObj : trigObjs){
      nTrigObj++;
      const float dR2_photon = reco::deltaR2(photon_eta,photon_phi,trigObj.eta(),trigObj.phi());
      if(dR2_photon<maxDR2) matchedObjs.push_back(&trigObj);

      const float dR2_meson = reco::deltaR2(meson_eta,meson_phi,trigObj.eta(),trigObj.phi());
      if(dR2_meson<maxDR2) matchedObjs.push_back(&trigObj);
    }
    return matchedObjs;
  }
}
//------------------------------------------------------------------------------//

// ------------ method called for each event  ------------
void HPhiGammaAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<pat::PackedCandidate>  > PFCandidates;
  iEvent.getByToken(packedPFCandidatesToken_, PFCandidates);

  edm::Handle<std::vector<pat::Muon>  > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_, slimmedMuons);

  edm::Handle<std::vector<reco::GenParticle>  > prunedGenParticles;
  if(!runningOnData_)iEvent.getByToken(prunedGenParticlesToken_, prunedGenParticles);

  //edm::Handle<std::vector<pat::PackedGenParticle>> packedGenParticles;  //GENPART
  //if(!runningOnData_) iEvent.getByToken(packedGenParticlesToken_, packedGenParticles);

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

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects; //triggObj
  iEvent.getByToken(triggerObjectsToken_, triggerObjects); //triggObj
  
  std::vector<pat::TriggerObjectStandAlone> unpackedTrigObjs; //triggObj

  //edm::Handle<reco::JetCorrector> jetCorr;
  //iEvent.getByToken(jetCorrectorToken_, jetCorr);

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];

  _Nevents_processed++; //This will be saved in the output tree, giving the number of processed events

  //Retrieve the run number
  if(runningOnData_){
    run_number = iEvent.id().run();
  }

  event_number = iEvent.id().event();

  //FIXMEEEEE
  bool Kplus_found  = false;
  bool Kminus_found = false;

  //cout <<endl<<"---------- event n. "<<event_number<<" ------------"<<endl;

  if (!runningOnData_){
    for(auto gen = prunedGenParticles->begin(); gen != prunedGenParticles->end(); ++gen){
      //if (gen->mother()) cout<<"pdgID = "<<gen->pdgId()<<" (mother ID = "<<gen->mother()->pdgId()<<")"<<endl;
      if( gen->pdgId() == 321  && gen->mother()->pdgId() == 333 && gen->mother()->mother()->pdgId() == 25) Kplus_found  = true;
      if( gen->pdgId() == -321 && gen->mother()->pdgId() == 333 && gen->mother()->mother()->pdgId() == 25) Kminus_found = true;
    }
    if (Kplus_found && Kminus_found){
      //cout<<"KK pair found!"<<endl;
      nKK_found ++;
    }
    else{
      //cout<<"KK not found!"<<endl;
      nKK_notFound ++;
    }
  }
  //cout <<endl<<"nKK_found = "<<nKK_found<<", nKK_notFound = "<<nKK_notFound<<endl;
  //cout <<endl<<"-----------------------------------------------------"<<endl;

   //*************************************************************//
   //                                                             //
   //--------------------------- Trigger -------------------------//
   //                                                             //
   //*************************************************************//

   //Examine the trigger information, return if the trigger doesn't switch on and count the number of events where the trigger has switched on
   isTwoProngTrigger = false;
   isIsoMuTrigger    = false;

   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   for(unsigned int i = 0, n = triggerBits->size(); i < n; ++i){
     if(!triggerBits->accept(i)) continue;
     std::string tmp_triggername = names.triggerName(i);

     if( tmp_triggername.find("HLT_Photon35_TwoProngs35_v") != std::string::npos ){
       isTwoProngTrigger = true;
     }

     if(isTriggerStudies){ //start trigger studies
       if(tmp_triggername.find("HLT_IsoMu24") != std::string::npos ){ //Muon trigger
         isIsoMuTrigger = true;
         if (verbose) cout<<"IsoMu24 triggered"<<endl;
         }
       //else return;
     } //end trigger studies
   }

   if(!isTwoProngTrigger){
    if(verbose) cout<<"Event not triggered, RETURN."<<endl;
    //return; //FIXMEEEEEE
  }
  _Nevents_triggered++;

  //Create unpacked filter names
  for(auto& trigObj : *triggerObjects){
      unpackedTrigObjs.push_back(trigObj);
      unpackedTrigObjs.back().unpackFilterLabels(iEvent,*triggerBits);
  }

  //*************************************************************//
  //                                                             //
  //-------------------------- Vertices -------------------------//
  //                                                             //
  //*************************************************************//

  //Count the number of vertices and return if there's no vertex
  nPV = 0;
  if(slimmedPV->size()<=0){
    if(verbose) cout<<"No primary vertex found, RETURN."<<endl;
    return;
  }

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
  //------------------ Variable initialization ------------------//
  //                                                             //
  //*************************************************************//

 nMuons10              = 0;
 nElectrons10          = 0;
 nMuons20              = 0;
 nElectrons20          = 0;
 nPhotonsChosen        = 0;
 nJets30               = 0;
 nJets25               = 0;
 nJets20               = 0;
 nPhotons38WP80              = 0;
 nPhotonsWP90_pT20_2p5eta3p0 = 0;
 nPhotonsWP90_pT15_barrel    = 0;
 nPhotonsWP90_pT25_endcap    = 0;
 
  //These variables will go in the tree
 ph_eT     = 0.;
 ph_eta    = 0.;
 ph_etaSC  = 0.;
 ph_phi    = 0.;
 LorentzVector ph_p4;

 ph_iso_ChargedHadron = 0.;
 ph_iso_NeutralHadron = 0.;
 ph_iso_Photon        = 0.;
 ph_iso_eArho         = 0.;
 photonRegressionError = 0.;


 eTphMax = -1000.;

 _Jet_Photon_invMass  = -1.;
 //_MesonMass           = -1.;
 //_Hmass_From2K_Photon = -1.;

 _PhiMass         = -1.;
 _RhoMass         = -1.;
 _K0starMass      = -1.;
 _PhiGammaMass    = -1.;
 _RhoGammaMass    = -1.;
 _K0starGammaMass = -1.;

 met_pT      = 0.;
 metpuppi_pT = 0.;

 _bestJet_pT                       =-1.;
 _bestJet_eta                      =-1.;
 _bestJet_phi                      =-1.;
 _bestJet_nDaughters               = 0;
 _bestJet_pTMax                    =-1.;
 _bestJet_chargedEmEnergy          = 0.;
 _bestJet_neutralEmEnergy          = 0.;
 _bestJet_chargedHadEnergy         = 0.;
 _bestJet_neutralHadEnergy         = 0.;
 _bestJet_chargedEmEnergyFraction  = 0.;
 _bestJet_neutralEmEnergyFraction  = 0.;
 _bestJet_chargedHadEnergyFraction = 0.;
 _bestJet_neutralHadEnergyFraction = 0.;
 _bestJet_chargedHadMultiplicity   = 0.;
 _bestJet_invMass                  = 0.;
 _bestJet_Photon_invMass           = 0.;
 _bestJet_JECunc                   = 0.;
 _firstCandPt                      = 0.;
 _firstCandEta                     = 0.;
 _firstCandPhi                     = 0.;
 _firstCandCharge                  = 0.;
 _secondCandPt                     = 0.;
 _secondCandEta                    = 0.;
 _secondCandPhi                    = 0.;
 _secondCandCharge                 = 0.;
 _bestCouplePt                     = 0.;
 _bestCoupleEta                    = 0.;
 _bestCouplePhi                    = 0.;  

//K-candidates and PHI ISOLATION
K1_sum_pT_05        = 0.;
K1_sum_pT_05_ch     = 0.;

K2_sum_pT_05        = 0.;
K2_sum_pT_05_ch     = 0.;

couple_sum_pT_05    = 0.;
couple_sum_pT_05_ch = 0.;

_iso_K1        = 0.;
_iso_K1_ch     = 0.;
_iso_K2        = 0.;
_iso_K2_ch     = 0.;
_iso_couple    = 0.;
_iso_couple_ch = 0.;

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
 //Count muons for each event
for(auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
  if(mu->pt() < 10. || !mu->CutBasedIdMedium || fabs(mu->eta()) > 2.4 || fabs(mu->muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(mu->muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
  if(!mu->PFIsoLoose) continue;
  nMuons10++;
  if(mu->pt() < 20.) continue;
  nMuons20++;
}

if(isTriggerStudies){ //start trigger studies
  
  isBestMuMu_Found = false;
  float currentMuMuPt = -1.;
  bestMuMuPt = -1.;
  float currentMuMuMass = -1.;
  bestMuMuMass = -1.;

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
          bestMuMuPt   = currentMuMuPt;
          bestMuMuMass = currentMuMuMass;
          isBestMuMu_Found = true;

      } //Muon second forloop end
    }//Muon first forloop end

  if(!isBestMuMu_Found) { 
    if (verbose) cout<<"RETURN: No Z->mumu found."<<endl<<endl;
    //return;
  }
  
  if(isBestMuMu_Found && verbose){
    cout<<"Muon pair found, with pT = "<<bestMuMuPt<<" and inv mass = "<<bestMuMuMass<<endl;
  } 
} //end trigger studies

  //*************************************************************//
  //                                                             //
  //-------------------------- Electrons ------------------------//
  //                                                             //
  //*************************************************************//
  //Count the number of electrons
  // Get rho value
edm::Handle< double > rhoH;
iEvent.getByToken(rhoToken_,rhoH);
rho_ = *rhoH;

float corr_pt = 0.;
ph_en_sigmaUP = 0.;
ph_en_sigmaDW = 0.;
ph_en_scaleUP = 0.;
ph_en_scaleDW = 0.;

for(auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
    //Calculate electron p4, correct it with the Scale&Smearing correction and extract the pT
    LorentzVector el_p4 = el->p4() * el->userFloat("ecalTrkEnergyPostCorr")/el->energy();
    corr_pt = el_p4.pt();

    if(corr_pt < 10. || fabs(el->eta()) > 2.5 || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;

    //float abseta = fabs(el->superCluster()->eta());
    //float eA     = effectiveAreas_el_.getEffectiveArea(abseta);
    //float el_iso   = (el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/corr_pt;
    //if(el_iso > 0.35) continue;

    //-------------Conditions on loose/medium MVA electron ID-------------//
    if(el->electronID("mvaEleID-Fall17-iso-V2-wp80") == 0) continue;
    nElectrons10++;
    if (corr_pt < 20.) continue;
    nElectrons20++;
  }

  //std::cout << "Nelectrons " << nElectrons << " Nmuons " << nMuons << std::endl;


  //*************************************************************//
  //                                                             //
  //--------------------------- Photons -------------------------//
  //                                                             //
  //*************************************************************//
  if(verbose) cout<< "PHOTONs"<<" --------------------------------"<<endl;
  
  bool cand_photon_found = false; //initialize this bool to false, return if it doesn't turn into true
  float corr_et = -1.;

  for(auto photon = slimmedPhotons->begin(); photon != slimmedPhotons->end(); ++photon){ //PHOTON FORLOOP START --------------------------------
    
    // Apply energy scale corrections, from 18Apr2023 the correction are embedded in the config file with the postReco tool
    corr_et = photon->et();// * photon->userFloat("ecalEnergyPostCorr") / photon->energy(); 

    if(photon->photonID("mvaPhoID-RunIIFall17-v2-wp90") == 0) continue; //WP90
    if(!photon->passElectronVeto()) continue; 

    if(corr_et > 20. && fabs(photon->eta()) >= 2.5   && fabs(photon->eta()) < 3.)    nPhotonsWP90_pT20_2p5eta3p0++;
    if(corr_et > 15. && fabs(photon->eta()) >= 0.    && fabs(photon->eta()) < 1.444) nPhotonsWP90_pT15_barrel++;
    if(corr_et > 25. && fabs(photon->eta()) > 1.566  && fabs(photon->eta()) < 2.5)   nPhotonsWP90_pT25_endcap++;

    if(fabs(photon->eta()) > 2.1) continue; //eta cut corresponding to the trigger's one
    if(corr_et < 38.) continue;
    if(photon->photonID("mvaPhoID-RunIIFall17-v2-wp80") == 0) continue; //WP80
    //if(photon->photonID("cutBasedPhotonID-Fall17-94X-V2-medium") == 0) continue; //WP80 FIXMEEEEEEE

    float abseta = fabs(photon->superCluster()->eta());
    float eA = effectiveAreas_ph_.getEffectiveArea(abseta);

    //photon_iso = (pfIso.sumChargedHadronPt + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_))/photon->et();

    //if(photon->chargedHadronIso()/corr_et > 0.3 || photon->photonIso() > 4.) continue; //|| photon->trackIso() > 6

    nPhotons38WP80++;

    // Apply energy scale corrections
    ph_p4 = photon->p4();// * photon->userFloat("ecalEnergyPostCorr") / photon->energy();

    if(corr_et < eTphMax) continue; //choose as best photon the one with highest eT
    
    //reinitialize the trigger matching bool
    isPhotonTriggerMatched = false;
    
    eTphMax = corr_et;
    ph_iso_ChargedHadron = photon->chargedHadronIso();
    ph_iso_NeutralHadron = photon->neutralHadronIso();
    ph_iso_Photon        = photon->photonIso();
    ph_iso_eArho         = eA * rho_;
    ph_eT     = corr_et;
    //ph_energy = photon->energy();
    ph_en_sigmaUP = photon->userFloat("energySigmaUp")/cosh(photon->eta()); //dividing for cosh(eta) means passing from energy to eT
    ph_en_sigmaDW = photon->userFloat("energySigmaDown")/cosh(photon->eta());  //dividing for cosh(eta) means passing from energy to eT
    ph_en_scaleUP = photon->userFloat("energyScaleUp")/cosh(photon->eta()); //dividing for cosh(eta) means passing from energy to eT
    ph_en_scaleDW = photon->userFloat("energyScaleDown")/cosh(photon->eta()); //dividing for cosh(eta) means passing from energy to eT
    ph_eta    = photon->eta();
    ph_etaSC  = photon->superCluster()->eta();
    ph_phi    = photon->phi();
    ph_HoverE       = photon->hadronicOverEm();
    ph_hadTowOverEm = photon->hadTowOverEm();

    photonRegressionError = photon->getCorrectedEnergyError(reco::Photon::P4type::regression2);
    if(debug) cout << "Regression2 Energy Error: " << photonRegressionError << endl;
    
    // Verify if the photon fired the HLT
    //isPhotonTriggerMatched = photon->triggered("HLT_Photon35_TwoProngs35_v*");

    cand_photon_found = true;
    nPhotonsChosen++;

/*
    if(!runningOnData_){ //ONLY FOR MC START ----------------------------------------------------------------------
      for (auto gen = prunedGenParticles->begin(); gen != prunedGenParticles->end(); ++gen){ //loop on genParticles start
	     //gen particles phi folding	
        float deltaPhi = fabs(ph_phi-gen->phi());
        if (deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;

        float deltaR = sqrt((ph_eta-gen->eta())*(ph_eta-gen->eta())+deltaPhi*deltaPhi);
        float deltapT = fabs(ph_eT-gen->pt());

        if(deltaR > deltaRMax || deltapT > deltapTMax) continue;
        deltapTMax = deltapT;
        gen_ID = gen->pdgId();
        gen_mother = gen->mother()->pdgId();
    } //loop on genParticles end

    if(gen_ID == 22) is_photon_a_photon = true;
    if(gen_ID != 22 && verbose) cout << "Photon gen ID = " << gen_ID << endl;
    if(gen_ID != 22 && fabs(gen_mother) == 24 && verbose) cout << "Photon gen ID when matched = " << gen_ID << endl;
    if(fabs(gen_mother) == 25) is_photon_matched = true;
    if(verbose) {
      cout<<"Bool: is_photon_matched = "<<is_photon_matched<<endl;
    }
  } //ONLY FOR MC END ----------------------------------------------------------------------------------------------
*/
}//PHOTON FORLOOP END -------------------------------------------------------------------------------------------------------------------------

//Return if there are no photons chosen
if(!cand_photon_found) {
  cout<<"No best photon found, RETURN."<<endl;
  return;
}

_Nevents_isPhoton++;

  //*************************************************************//
  //                                                             //
  //--------------------------- N-jets --------------------------//
  //                                                             //
  //*************************************************************//

  //int nJet=1;
int jetIndex=-1;
int bestJet_Index=-1;
//int MCtruthIndex = -1;
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
LorentzVector firstCand_p4_K;
LorentzVector secondCand_p4_K;
LorentzVector firstCand_p4_Pi;
LorentzVector secondCand_p4_Pi;
LorentzVector couple_p4;
LorentzVector couple_p4_Phi;
LorentzVector couple_p4_Rho;
LorentzVector couple_p4_K0star_PiK;
LorentzVector couple_p4_K0star_KPi;
LorentzVector best_firstCand_p4;
LorentzVector best_secondCand_p4;
LorentzVector best_couple_p4;
float bestCoupleOfTheJet_pT = 0.;
float deltaR_KChosen = 0.;
float deltaR_K = 0.;
firstCandEnergy_K = 0.;
secondCandEnergy_K = 0.;
firstCandEnergy_Pi = 0.;
secondCandEnergy_Pi = 0.;
firstCandPx=0.;
firstCandPy=0.;
firstCandPz=0.;
float firstCand_dxy=-999.;
float firstCand_dxyErr=-999.;
float firstCand_dz=-999.;
float firstCand_dzErr=-999.;
bestFirstCand_dxy=-999.;
bestFirstCand_dz=-999.;
bestFirstCand_dxyErr=-999.;
bestFirstCand_dzErr=-999.;
secondCandPx=0.;
secondCandPy=0.;
secondCandPz=0.;
float secondCand_dxy=-999.;
float secondCand_dxyErr=-999.;
float secondCand_dz=-999.;
float secondCand_dzErr=-999.;
bestSecondCand_dxy=-999.;
bestSecondCand_dz=-999.;
bestSecondCand_dxyErr=-999.;
bestSecondCand_dzErr=-999.;
float firstCandEta=0.;
float firstCandPhi=0.;
float secondCandEta=0.;
float secondCandPhi=0.;
float PhiMass = 0.;
float RhoMass = 0.;
float K0starMass = 0.;
float K0starPiKMass = 0.;
float K0starKPiMass = 0.;
float kMass  = 0.4937;
float PiMass = 0.13957;
float K0TheoryMass = 0.892;
bool isBestCoupleOfTheEvent_Found=false;
bool isPhi = false;
bool isRho = false;
bool isK0star = false;
bool isK0starPiK = false;
bool isK0starKPi = false;
_isFirstCandK = false;

std::vector<float> pt_jets_vector;   //for VBF veto
std::vector<LorentzVector> p4_jets_vector;   //for VBF veto
//std::vector<float> m_jets_vector;   //for VBF veto

//JEC uncertainties
JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar); 


//JET LOOP
  for (auto jet = slimmedJets->begin(); jet != slimmedJets->end(); ++jet) { //JET LOOP START -------------------------------------------------------- 

    jetIndex++;

    _Jet_Photon_invMass=(jet->p4()+ph_p4).M(); //calculate inv mass
    nDaughters= jet->numberOfDaughters(); //calculate number of daughters

    float NHF  = jet->neutralHadronEnergyFraction();
    float NEMF = jet->neutralEmEnergyFraction();
    float MUF  = jet->muonEnergyFraction();
    float CHF  = jet->chargedHadronEnergyFraction();
    int   CHM  = jet->chargedHadronMultiplicity();
    float CEMF = jet->chargedEmEnergyFraction();
    int   NumNeutralParticle = jet->neutralMultiplicity();
    int   pileUPjetID        = jet->userInt("pileupJetId:fullId");

    if(verbose){
    cout<<"________________________________________";
    cout<<"jet n."<<jetIndex<<endl;
    cout<<"jet pT  = "<<jet->pt()<<endl;
    cout<<"jet eta = "<<jet->eta()<<endl;
    cout<<"jet neutralHadronEnergyFraction = "<<NHF<<endl;
    cout<<"jet neutralEmEnergyFraction     = "<<NEMF<<endl;
    cout<<"jet muonEnergyFraction          = "<<MUF<<endl;
    cout<<"jet chargedHadronEnergyFraction = "<<CHF<<endl;
    cout<<"jet chargedHadronMultiplicity   = "<<CHM<<endl;
    cout<<"jet chargedEmEnergyFraction     = "<<CEMF<<endl;
    cout<<"jet numNeutralParticle          = "<<NumNeutralParticle<<endl;
    cout<<"jet pileUP                      = "<<pileUPjetID<<endl;
    //if(verbose) cout<<"nDaughters = "<<nDaughters<<endl;
    }

    /*
    //----------------------------- Pre-Filters --------------------------------------------------------
    if(jet->neutralHadronEnergyFraction() > 0.2) continue; //reject if neutralhadron-energy fraction is > 0.2
    if(jet->neutralEmEnergyFraction() > 0.9) continue; //reject if neutralEm-energy fraction is > 0.9, alias NO-PHOTON FILTER                              
    if(nDaughters < 2) continue; //reject if number of constituens is less then 1
    if(jet->muonEnergyFraction() > 0.8) continue; //reject if muon-energy fraction is > 0.8                                             
    if(jet->chargedHadronEnergyFraction() <= 0.) continue; //reject if chargedHadron-energy fraction is 0                              
    if(jet->chargedHadronMultiplicity() == 0) continue; //reject if there are NOT charged hadrons                              
    if(jet->chargedEmEnergyFraction() > 0.8) continue; //reject if chargedEm-energy fraction is > 0.8
    */

    //----------------------------- Pre-Filters --------------------------------------------------------
    bool jetID = false;

    //if(pileUPjetID <= 4) continue; //medium jet PU ID 
    
    if(abs(jet->eta())<=2.4)                        jetID = (CEMF<0.8 && CHM>0 && CHF>0 && nDaughters>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9);
    if(abs(jet->eta())>2.4 && abs(jet->eta())<=2.7) jetID = (NEMF<0.99 && NHF < 0.9);
    if(abs(jet->eta())>2.7 && abs(jet->eta())<=3.0) jetID = (NEMF>0.0 && NEMF<0.99 && NHF<0.9 && NumNeutralParticle>1);
    if(abs(jet->eta())>3.0)                         jetID = (NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10);
    
    if(verbose) cout<<"jetID = "<<jetID<<endl;
    if(jetID == 0) continue;

    //HEM 2018 correction ----------------------------------------------------------------------------------------
    float jet_pt = jet->pt();
    
    if (!runningOnData_){ 
     
      if (jet->phi() >= -1.57 && jet->phi() <= -0.87 && jet->eta() >= -2.5 && jet->eta() <= -1.3) jet_pt = 0.8 * jet_pt; 
      if (jet->phi() >= -1.57 && jet->phi() <= -0.87 && jet->eta() >= -3.0 && jet->eta() <= -2.5) jet_pt = 0.65 * jet_pt; 

    }

    //------------------------------------------------------------------------------------------------------------

    //for VBF veto ----------------------
    if(jet_pt < 20. || abs(jet->eta()) > 4.7) continue;
    nJets20++;
    p4_jets_vector.push_back(jet->p4());
    pt_jets_vector.push_back(jet->pt());
    //-----------------------------------

    //JEC and JES uncertainties ----------------------------------------------
    jecUnc->setJetEta(jet->eta());
    jecUnc->setJetPt(jet->pt()); // here you must use the CORRECTED jet pt
    double unc = jecUnc->getUncertainty(true);
    double jetPt_shifted = jet->pt() - unc; //just for JEC syst
    if (debug) cout <<"jet pT = "<<jet->pt()<<", jet pT shifted = "<<jetPt_shifted<<endl;
    //------------------------------------------------------------------------

    if(jet_pt < 40. || abs(jet->eta()) > 2.1) continue;
    if(_Jet_Photon_invMass < 100.) continue; //reject jets with inv mass lower then 100 GeV
                           
     //-------------------------------------------------------------------------------------------------      
    
    if (verbose) cout<<"    Jet at index = "<<jetIndex<<" passed the cuts:"<<endl; 

    //------------------------------- access to MC truth -------------------------------------------
    /*
    if(!runningOnData_) { //ONLY FOR MC START  
      if(MCtruthIndex == -1)
      {      
           for(auto gen = prunedGenParticles->begin(); gen != prunedGenParticles->end(); ++gen){ //GEN PARTICLES LOOP START 
	          //phi folding	
            float deltaPhi = fabs(jet->phi()-gen->phi());
            if (deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;

            float R = sqrt((jet->eta()-gen->eta())*(jet->eta()-gen->eta())+deltaPhi * deltaPhi); 
            if( (gen->pdgId() == 333 || gen->pdgId() == 113)  && gen->mother()->pdgId() == 25 && R < 0.4)
            {
              MCtruthIndex=jetIndex;
	              //if(verbose) cout<<" MC index of the jet containing the Phi from the Higgs = "<<MCtruthIndex<<endl; //MCtruthIndex and deltaR can be used out from the jet-loop 
              deltaR=R; 
            } 
	       } //GEN PARTICLE LOOP END
      } 
    } //ONLY FOR MC END //-----------------------------------------------------------------------------
    */

      //-------------------------------------daughters forloop----------------------------

      for(int firstCand_Index=0; firstCand_Index < nDaughters; firstCand_Index++){ //1ST LOOP STARTS

        //if (verbose) cout<<"Daughter n."<<firstCand_Index+1<<" pT = "<<slimmedJets->at(jetIndex).daughter(firstCand_Index)->pt()<<endl;

         //loop only over charged daughters
        if (slimmedJets->at(jetIndex).daughter(firstCand_Index)->charge() == 0) continue;
        
        if(slimmedJets->at(jetIndex).daughter(firstCand_Index)->bestTrack() == NULL) continue;

        if (abs(slimmedJets->at(jetIndex).daughter(firstCand_Index)->bestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2) continue;
        if (abs(slimmedJets->at(jetIndex).daughter(firstCand_Index)->bestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5)  continue;
        if (!(slimmedJets->at(jetIndex).daughter(firstCand_Index)->bestTrack()->quality(reco::Track::highPurity))) continue;
        
        firstCandPt  = slimmedJets->at(jetIndex).daughter(firstCand_Index)->pt();  //take firstCand pt
        firstCandEta = slimmedJets->at(jetIndex).daughter(firstCand_Index)->eta(); //take firstCand eta
        firstCandPhi = slimmedJets->at(jetIndex).daughter(firstCand_Index)->phi(); //take firstCand phi

        firstCand_dxy    = slimmedJets->at(jetIndex).daughter(firstCand_Index)->bestTrack()->dxy((&slimmedPV->at(0))->position());
        firstCand_dxyErr = slimmedJets->at(jetIndex).daughter(firstCand_Index)->bestTrack()->dxyError();
        
        firstCand_dz    = slimmedJets->at(jetIndex).daughter(firstCand_Index)->bestTrack()->dz((&slimmedPV->at(0))->position());
        firstCand_dzErr = slimmedJets->at(jetIndex).daughter(firstCand_Index)->bestTrack()->dzError();

        if(firstCandPt < 1.) continue; //firstCand filter if pT < 1 GeV

        for(int secondCand_Index=firstCand_Index+1; secondCand_Index < nDaughters; secondCand_Index++){ //2ND LOOP STARTS

          if (slimmedJets->at(jetIndex).daughter(secondCand_Index)->charge() == 0) continue;

          //minimum apporach distance
          if (slimmedJets->at(jetIndex).daughter(secondCand_Index)->bestTrack() == NULL) continue;
          //if(verbose) cout<<"Track without bestTrack() found!"<<endl;

          //cout<<"dxy = "<<slimmedJets->at(jetIndex).daughter(secondCand_Index)->bestTrack()->dxy((&slimmedPV->at(0))->position())<<endl;
          //cout<<"dz  = "<<slimmedJets->at(jetIndex).daughter(secondCand_Index)->bestTrack()->dz((&slimmedPV->at(0))->position())<<endl;
          
          if ((abs(slimmedJets->at(jetIndex).daughter(secondCand_Index)->bestTrack()->dxy((&slimmedPV->at(0))->position()))) >= 0.2 || abs(slimmedJets->at(jetIndex).daughter(secondCand_Index)->bestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5 ) continue;
          if (!(slimmedJets->at(jetIndex).daughter(secondCand_Index)->bestTrack()->quality(reco::Track::highPurity))) continue;
          
          secondCand_dxy    = slimmedJets->at(jetIndex).daughter(secondCand_Index)->bestTrack()->dxy((&slimmedPV->at(0))->position());
          secondCand_dxyErr = slimmedJets->at(jetIndex).daughter(secondCand_Index)->bestTrack()->dxyError();
          
          secondCand_dz    = slimmedJets->at(jetIndex).daughter(secondCand_Index)->bestTrack()->dz((&slimmedPV->at(0))->position());
          secondCand_dzErr = slimmedJets->at(jetIndex).daughter(secondCand_Index)->bestTrack()->dzError();

          //TRKs PT CUT --------------------------------------------------------------------------
          secondCandPt  = slimmedJets->at(jetIndex).daughter(secondCand_Index)->pt();
          if(secondCandPt < 1.) continue; //firstCand filter if pT < 1 GeV
          if(firstCandPt < 10. && secondCandPt < 10.) continue;  //filter if both cand pT are < 10GeV

          //DITRK DELTA R CUT --------------------------------------------------------------------------
          secondCandEta = slimmedJets->at(jetIndex).daughter(secondCand_Index)->eta();
          secondCandPhi = slimmedJets->at(jetIndex).daughter(secondCand_Index)->phi();
          
          float deltaEta= firstCandEta - secondCandEta;

          float deltaPhi = fabs(firstCandPhi - secondCandPhi);  //phi folding	
          if (deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;

          deltaR_K= sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
          if(deltaR_K > 0.07) continue;

          //OPPOSITE CHARGE - FILTER ------------------------------------------------------------
          firstCandCharge  = slimmedJets->at(jetIndex).daughter(firstCand_Index)->charge(); //take firstCand charge
          secondCandCharge = slimmedJets->at(jetIndex).daughter(secondCand_Index)->charge(); //take secondCand charge
          if(firstCandCharge * secondCandCharge >= 0) continue; //choose only opposite charges

          //QUADRIMOMENTUM CALCULATION ------------------------------------------------------------
          firstCand_p4  = slimmedJets->at(jetIndex).daughter(firstCand_Index)->p4(); //take quadrimomentum
          secondCand_p4 = slimmedJets->at(jetIndex).daughter(secondCand_Index)->p4();

          firstCandPx  = firstCand_p4.px(); //take px, py, pz of the first candidate
          firstCandPy  = firstCand_p4.py();
          firstCandPz  = firstCand_p4.pz();
          secondCandPx = secondCand_p4.px(); //take px, py, pz of the second candidate
          secondCandPy = secondCand_p4.py();
          secondCandPz = secondCand_p4.pz();

          //PIONS OR KAONS HYPOTHESIS -----------------------------------------------------------------------------------------------------------------------------------------------            
          firstCandEnergy_K   = sqrt(firstCandPx  * firstCandPx  + firstCandPy  * firstCandPy  + firstCandPz  * firstCandPz  + kMass  * kMass ); //Kaon hypothesis energy recalculation
          secondCandEnergy_K  = sqrt(secondCandPx * secondCandPx + secondCandPy * secondCandPy + secondCandPz * secondCandPz + kMass  * kMass ); //Kaon hypothesis energy recalculation
          firstCandEnergy_Pi  = sqrt(firstCandPx  * firstCandPx  + firstCandPy  * firstCandPy  + firstCandPz  * firstCandPz  + PiMass * PiMass); //Pion hypothesis energy recalculation
          secondCandEnergy_Pi = sqrt(secondCandPx * secondCandPx + secondCandPy * secondCandPy + secondCandPz * secondCandPz + PiMass * PiMass); //Pion hypothesis energy recalculation
          
          /*
          if (verbose) {
            cout<<"firstCandEnergy_K   = "<<firstCandEnergy_K<<endl;
            cout<<"firstCandEnergy_Pi  = "<<firstCandEnergy_Pi<<endl;
            cout<<"secondCandEnergy_K  = "<<secondCandEnergy_K<<endl;
            cout<<"secondCandEnergy_Pi = "<<secondCandEnergy_Pi<<endl;
          }
          */

          firstCand_p4_K   = firstCand_p4.SetE(firstCandEnergy_K); //Kaon hypothesis quadrimomentum correction
          secondCand_p4_K  = secondCand_p4.SetE(secondCandEnergy_K); //Kaon hypothesis quadrimomentum correction
          firstCand_p4_Pi  = firstCand_p4.SetE(firstCandEnergy_Pi); //Kaon hypothesis quadrimomentum correction
          secondCand_p4_Pi = secondCand_p4.SetE(secondCandEnergy_Pi); //Kaon hypothesis quadrimomentum correction

          couple_p4_Phi        = firstCand_p4_K  + secondCand_p4_K; //calculation of the couple-quadrimomentum after the correction
          couple_p4_Rho        = firstCand_p4_Pi + secondCand_p4_Pi; //calculation of the couple-quadrimomentum after the correction
          couple_p4_K0star_PiK = firstCand_p4_Pi + secondCand_p4_K;
          couple_p4_K0star_KPi = firstCand_p4_K  + secondCand_p4_Pi;

          if (verbose) {
            cout<<"KK   pT = "<<couple_p4_Phi.pt()<<endl;
            cout<<"PiPi pT = "<<couple_p4_Rho.pt()<<endl;
            cout<<"PiK  pT = "<<couple_p4_K0star_PiK.pt()<<endl;
            cout<<"KPi  pT = "<<couple_p4_K0star_KPi.pt()<<endl;
          }

          //DITRK PT CUT -------------------------------------------------------------------------
          if(couple_p4_Phi.pt() < 38.) {
            if(verbose) cout<<"couplePt cut NOT passed"<<endl;
            continue;
          }  
          
          //MESON INV MASS CUT -------------------------------------------------------------------------
          isPhi       = false;
          isRho       = false;
          isK0star    = false;
          isK0starKPi = false;
          isK0starPiK = false;

          //Rho tagging
          RhoMass = (couple_p4_Rho).M(); //calculate inv mass of the Rho candidate  
          if (verbose) cout<<"mPiPi (before the meson mass selection) =  "<<RhoMass<<endl;
          if(RhoMass > 0.5 && RhoMass < 1.) isRho = true;

          //K0star tagging
          K0starPiKMass = (couple_p4_K0star_PiK).M();
          K0starKPiMass = (couple_p4_K0star_KPi).M();
          if (verbose) cout<<"mPiK (before the meson mass selection) =  "<<K0starPiKMass<<endl;
          if (verbose) cout<<"mKPi (before the meson mass selection) =  "<<K0starKPiMass<<endl;

          if(abs(K0starPiKMass - K0TheoryMass) < abs(K0starKPiMass - K0TheoryMass)){
            K0starMass  = K0starPiKMass;
            isK0starPiK = true;
            }
          else{
            K0starMass  = K0starKPiMass;        
            isK0starKPi = true;
          }

          if(K0starMass > 0.6 && K0starMass < 1.) isK0star = true;

          //Phi tagging
          PhiMass = (couple_p4_Phi).M(); //calculate inv mass of the Phi candidate  
          if (verbose) cout<<"mKK (before the meson mass selection) =  "<<PhiMass<<endl;
          if(PhiMass > 1. && PhiMass < 1.05) isPhi = true;

          //Continue if no candidates
          if (!isPhi && !isRho && !isK0star) continue; //continue if the pair mass doesn't match any of the two mass hypothesis

          //save the 4-mom for the K hypothesis because they are the same for each hypothesis except for the mass value
          firstCand_p4  = firstCand_p4_K; 
          secondCand_p4 = secondCand_p4_K;
          couple_p4     = couple_p4_Phi;

          /*
          //update values of quadrimomenta
          if(isPhi){
            firstCand_p4  = firstCand_p4_K; 
            secondCand_p4 = secondCand_p4_K;
            couple_p4     = couple_p4_Phi;
          }  
          if(isRho){
            firstCand_p4  = firstCand_p4_Pi;
            secondCand_p4 = secondCand_p4_Pi;
            couple_p4     = couple_p4_Rho;
          }
          if(isK0star && isK0starPiK){
            firstCand_p4  = firstCand_p4_Pi;
            secondCand_p4 = secondCand_p4_K;
            couple_p4     = couple_p4_K0star_PiK; 
          }
          if(isK0star && isK0starKPi){
            firstCand_p4  = firstCand_p4_K;
            secondCand_p4 = secondCand_p4_Pi;
            couple_p4     = couple_p4_K0star_KPi; 
          }
          */

          // ISOLATION CUT -------------------------------------------------------------------------  
          for(auto cand_iso = PFCandidates->begin(); cand_iso != PFCandidates->end(); ++cand_iso){ //ISOLATION FORLOOP START
            
            if(debug){
              cout <<endl<<"ISO CALC DETAILS ---------------------"<<endl;
              cout << "pt cand_iso = "<<cand_iso->pt()<<endl;
            }
            if(cand_iso->pt() < 0.5) continue; //do not consider tracks with pT < 500MeV

            //calculate the deltaR between the track and the first candidate ---------------------------------------
            float deltaPhi_K1 = fabs(firstCand_p4.phi()-cand_iso->phi());  //phi folding 
            if (deltaPhi_K1 > M_PI) deltaPhi_K1 = 2*M_PI - deltaPhi_K1;

            float deltaR_K1 = sqrt((firstCand_p4.eta()-cand_iso->eta())*(firstCand_p4.eta()-cand_iso->eta()) + deltaPhi_K1*deltaPhi_K1);
            if (debug) cout << "deltaR_K1 = "<<deltaR_K1<<endl;
            if(deltaR_K1 < 0.0005) continue; //remove first candidate from the sum

            //calculate the deltaR between the track and the second candidate ---------------------------------------
            float deltaPhi_K2 = fabs(secondCand_p4.phi()-cand_iso->phi());  //phi folding  
            if (deltaPhi_K2 > M_PI) deltaPhi_K2 = 2*M_PI - deltaPhi_K2;

            float deltaR_K2 = sqrt((secondCand_p4.eta()-cand_iso->eta())*(secondCand_p4.eta()-cand_iso->eta()) + deltaPhi_K2*deltaPhi_K2);
            if (debug) cout << "deltaR_K2 = "<<deltaR_K2<<endl;
            if(deltaR_K2 < 0.0005) continue; //remove second candidate from the sum

            //calculate the deltaR between the track and the best pair ---------------------------------------
            float deltaPhi_Couple = fabs(couple_p4.phi()-cand_iso->phi());  //phi folding  
            if (deltaPhi_Couple > M_PI) deltaPhi_Couple = 2*M_PI - deltaPhi_Couple;

            float deltaR_Couple = sqrt((couple_p4.eta()-cand_iso->eta())*(couple_p4.eta()-cand_iso->eta()) + deltaPhi_Couple*deltaPhi_Couple);

            //sum pT of the tracks inside a cone of deltaR = 0.3 ---------------------------------------
            if(deltaR_K1 <= 0.3) K1_sum_pT_05 += cand_iso->pt();
            if(deltaR_K2 <= 0.3) K2_sum_pT_05 += cand_iso->pt();
            if(deltaR_Couple <= 0.3) couple_sum_pT_05 += cand_iso->pt();
            //cout<< "charge before = "<<cand_iso->charge()<<endl;

            //sum pT of the charged tracks inside a cone of deltaR = 0.3 ---------------------------------------
            if (debug) cout << "Charge = "<< cand_iso->charge()<<endl;
            if(cand_iso->charge() == 0) continue;
          // cout << "particle charge = "<<cand_iso->charge()<<endl;
            if (debug) cout << "dxy = "<<fabs(cand_iso->dxy())<<" and dz = "<< fabs(cand_iso->dz())<<endl;
            if(fabs(cand_iso->dxy()) >= 0.2 || fabs(cand_iso->dz()) >= 0.5) continue; // Requesting charged particles to come from PV
            //cout<< "charge after = "<<cand_iso->charge()<<endl;
            if(deltaR_K1 <= 0.3) K1_sum_pT_05_ch += cand_iso->pt();
            if(deltaR_K2 <= 0.3) K2_sum_pT_05_ch += cand_iso->pt();
            if (debug) cout <<"deltaR_Couple = "<<deltaR_Couple<<endl;
            if(deltaR_Couple <= 0.3){
              couple_sum_pT_05_ch += cand_iso->pt();
              if (debug) cout<<"Particle in the cone: SumPt = "<<couple_sum_pT_05_ch<<endl;
            }
          } //ISOLATION FORLOOP END

          float isoCoupleCh = couple_p4.pt()/(couple_sum_pT_05_ch + couple_p4.pt());
          if(isoCoupleCh < 0.9) {
            cout<<"No isolation cut passed."<<endl;
            continue; //FIXMEEEEEEE
          }


          //PT MAX OF THE JET - FILTER -------------------------------------------------
          if (verbose) cout<<"Current bestCoupleOfTheEvent_Pt = "<<bestCoupleOfTheJet_pT<<endl;
          
          if(couple_p4.pt() <= bestCoupleOfTheJet_pT) {
            if(verbose) cout<<"Not passed: pT lower than the current best pair of the event"<<endl;
            continue; //choose the couple with greatest pt
          }

          //reinitialize the trigger matching bool
          isJetTriggerMatched = false;

          //If passed, this is the pair with the largest pT of the event so far
          bestCoupleOfTheJet_pT = couple_p4.pt();     
          if (verbose) cout<<"couple_p4.pt() = "<<bestCoupleOfTheJet_pT<<endl;

          if(verbose) cout<<"This is the best pair so far!"<<endl<<"-------------------------"<<endl;
          isBestCoupleOfTheEvent_Found = true;

          //Save if best pair has been found
          bestJet_Index           = jetIndex; //note the position of the chosen jet inside the vector	  
          deltaR_KChosen          = deltaR_K;
          _bestJet_Photon_invMass = _Jet_Photon_invMass;
          _isPhi                  = isPhi;
          _isRho                  = isRho;
          _isK0star               = isK0star;
          _bestJet_JECunc         = unc;
          _firstCandCharge        = firstCandCharge;
          _secondCandCharge       = secondCandCharge;
          bestFirstCand_dxy       = firstCand_dxy;
          bestFirstCand_dz        = firstCand_dz;
          bestSecondCand_dxy      = secondCand_dxy;
          bestSecondCand_dz       = secondCand_dz;
          bestFirstCand_dxyErr    = firstCand_dxyErr;
          bestFirstCand_dzErr     = firstCand_dzErr;
          bestSecondCand_dxyErr   = secondCand_dxyErr;
          bestSecondCand_dzErr    = secondCand_dzErr;

          // Verify if the best jet fired the HLT
          //isJetTriggerMatched = slimmedJets->at(jetIndex).triggered("HLT_Photon35_TwoProngs35_v*");

          best_firstCand_p4  = firstCand_p4; 
          best_secondCand_p4 = secondCand_p4;
          best_couple_p4     = couple_p4;

          //save just the masses for the different hypothesis
          if(_isPhi)                      _PhiMass = PhiMass,     _PhiGammaMass    = (couple_p4_Phi + ph_p4).M();
          if(_isRho)                      _RhoMass = RhoMass,     _RhoGammaMass    = (couple_p4_Rho + ph_p4).M();
          if(_isK0star && isK0starPiK) _K0starMass = K0starMass,  _K0starGammaMass = (couple_p4_K0star_PiK + ph_p4).M(), _isFirstCandK = false;
          if(_isK0star && isK0starKPi) _K0starMass = K0starMass,  _K0starGammaMass = (couple_p4_K0star_KPi + ph_p4).M(), _isFirstCandK = true;

          if (verbose) {
          cout<<"_PhiMass         = "<<_PhiMass<<endl;
          cout<<"_RhoMass         = "<<_RhoMass<<endl;
          cout<<"_K0starMass      = "<<_K0starMass<<endl;
          cout<<"_PhiGammaMass    = "<<_PhiGammaMass<<endl;
          cout<<"_RhoGammaMass    = "<<_RhoGammaMass<<endl;
          cout<<"_K0starGammaMass = "<<_K0starGammaMass<<endl;
          cout<<"firstCandPt      = "<<firstCand_p4.pt()<<endl;
          cout<<"secondCandPt     = "<<secondCand_p4.pt()<<endl;
          cout<<"isFirstCandK     = "<<_isFirstCandK<<endl;
          }

      } //2ND LOOP ENDS
	} //1ST LOOP ENDS

  if(jet->pt() < 25.) continue;
  nJets25++;
  if(jet->pt() < 30.) continue;
  nJets30++;

} //JET LOOP END

delete jecUnc;


if(!isBestCoupleOfTheEvent_Found) 
{
  cout<<"No best couple detected for current event, RETURN."<<endl;
  return;
}
_Nevents_bestCoupleFound++;      
if(verbose) cout<<"Bool: _Nevents_bestCoupleFound: "<<_Nevents_bestCoupleFound<<endl;        

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

_bestJet_invMass                  = slimmedJets->at(bestJet_Index).mass();
_bestJet_pT                       = slimmedJets->at(bestJet_Index).pt();
_bestJet_eta                      = slimmedJets->at(bestJet_Index).eta();
_bestJet_phi                      = slimmedJets->at(bestJet_Index).phi();
_bestJet_nDaughters               = slimmedJets->at(bestJet_Index).numberOfDaughters();
_bestJet_chargedEmEnergy          = slimmedJets->at(bestJet_Index).chargedEmEnergy();
_bestJet_neutralEmEnergy          = slimmedJets->at(bestJet_Index).neutralEmEnergy();
_bestJet_chargedHadEnergy         = slimmedJets->at(bestJet_Index).chargedHadronEnergy();
_bestJet_neutralHadEnergy         = slimmedJets->at(bestJet_Index).neutralHadronEnergy();
_bestJet_chargedEmEnergyFraction  = slimmedJets->at(bestJet_Index).chargedEmEnergyFraction();
_bestJet_neutralEmEnergyFraction  = slimmedJets->at(bestJet_Index).neutralEmEnergyFraction();
_bestJet_chargedHadEnergyFraction = slimmedJets->at(bestJet_Index).chargedHadronEnergyFraction();
_bestJet_neutralHadEnergyFraction = slimmedJets->at(bestJet_Index).neutralHadronEnergyFraction();
_bestJet_chargedHadMultiplicity   = slimmedJets->at(bestJet_Index).chargedHadronMultiplicity();
_bestJet_PUid                     = slimmedJets->at(bestJet_Index).userInt("pileupJetId:fullId");

//MESON MASS CALCULATION
//_MesonMass = (best_firstCand_p4 + best_secondCand_p4).M();

//H INV MASS CALCULATION
//_PhiGammaMass = (best_firstCand_p4 + best_secondCand_p4 + ph_p4).M(); //calculate inv mass of the Higgs candidate
  

//CANDIDATES SORTING
if(_firstCandPt < _secondCandPt)  //swap-values loop, in order to fill the tree with the candidate with max pt of the couple in firstCand branches  
  {                               //and one with min pt in secondCand branches
    float a,b,c,d,e;
    a = _firstCandPt;
    b = _firstCandEta;
    c = _firstCandPhi;
    d = firstCandEnergy;
    e = _firstCandCharge;
    _isFirstCandK   = !_isFirstCandK; 
    _firstCandPt    = _secondCandPt;
    _firstCandEta   = _secondCandEta;
    _firstCandPhi   = _secondCandPhi;
    firstCandEnergy = secondCandEnergy;
    _firstCandCharge = _secondCandCharge;
    _secondCandPt     = a;
    _secondCandEta    = b;
    _secondCandPhi    = c;
    secondCandEnergy  = d;
    _secondCandCharge = e;
}

//CUTS ON CANDIDATES PT
if(_firstCandPt < 20. || _secondCandPt < 5.) {
    cout<<"Final cut on candidates pT not passed, RETURN."<<endl;
    return;
}
  _Nevents_candPtFilter++;

//ISOLATION DATAMEMBER FOR TREE FILLING 
_iso_K1        = _firstCandPt/(K1_sum_pT_05 + _firstCandPt);
_iso_K2        = _secondCandPt/(K2_sum_pT_05 + _secondCandPt);
_iso_couple    = _bestCouplePt/(couple_sum_pT_05 + _bestCouplePt);
_iso_K1_ch     = _firstCandPt/(K1_sum_pT_05_ch + _firstCandPt);
_iso_K2_ch     = _secondCandPt/(K2_sum_pT_05_ch + _secondCandPt);
_iso_couple_ch = _bestCouplePt/(couple_sum_pT_05_ch + _bestCouplePt);

//CUT ON PHI ISOLATION
if(verbose){
  cout<<endl;
  cout<<"###### ISO           = "<<_iso_couple_ch<<endl;
  cout<<"###### isPhi         = "<<_isPhi<<endl;
  cout<<"###### isRho         = "<<_isRho<<endl;
  cout<<"###### isK0star      = "<<_isK0star<<endl;
  cout<<"###### SUM pT        = "<<couple_sum_pT_05_ch<<endl;
  cout<<"###### pT leading    = "<<_firstCandPt<<endl;
  cout<<"###### pT subleading = "<<_secondCandPt<<endl;
  //cout<<"###### MesonMass     = "<<_MesonMass<<endl;
  //cout<<"###### HMass         = "<<_Hmass_From2K_Photon<<endl;
}


/*
//VBF veto -------------------------------------------------------------
int jet_i = 0;
int jet_j = 0;

if(verbose){
  cout <<"VBF veto ------------------------------------"<<endl<<endl;
  cout <<"nJets20      = "<<nJets20<<endl;
  cout <<"nJets30      = "<<nJets30<<endl;
}

if (nJets20 > 2.){

  float deltaEtaJets = -1.;
  
  for (jet_i = 0; jet_i < nJets20; jet_i++){
    
    //cout <<"deltaEta(BestJet - jet"<<jet_i<<")  = "<<abs(_bestJet_eta-eta_jets_vector.at(jet_i))<<endl;
    
    if (abs(_bestJet_eta-eta_jets_vector.at(jet_i)) == 0) {
      if(debug) cout<<"Best jet index = "<<jet_i<<endl;
      continue;
    }

    for (jet_j = jet_i + 1; jet_j < nJets20; jet_j++){

      //cout <<"deltaEta(BestJet - jet"<<jet_j<<")  = "<<abs(_bestJet_eta-eta_jets_vector.at(jet_j))<<endl;
      if (abs(_bestJet_eta-eta_jets_vector.at(jet_j)) == 0) {
      if(debug) cout<<"Best jet index = "<<jet_j<<endl;
      continue;
      }

      deltaEtaJets = abs(eta_jets_vector.at(jet_i) - eta_jets_vector.at(jet_j));

      if(debug) cout <<"deltaEta(jet"<<jet_i<<" - jet"<<jet_j<<") = "<<deltaEtaJets<<endl;
      if (deltaEtaJets > 3.) return;

    }
  }
}
if(debug) cout <<"---------------------------------------------"<<endl<<endl;
*/

//***********************************************
//**************** VBF VETO *********************
//***********************************************

//remove bestJet from vectors and those jets with a deltaR too close to the photon
bool bestJetRemoved = false;

for (int i = 0; i < static_cast<int>(p4_jets_vector.size()); i++) {
    if(verbose) cout<<"pT jets (> 20 and eta < 4.7) = "<<pt_jets_vector[i]<<" (eta jet = "<<p4_jets_vector[i].Eta()<<endl;
    if (pt_jets_vector[i] == _bestJet_pT) {
        p4_jets_vector.erase(p4_jets_vector.begin() + i);
        pt_jets_vector.erase(pt_jets_vector.begin() + i);
        i--; // to handle the reduction in vector size
        bestJetRemoved = true;
  }
}

if(!bestJetRemoved)cout<<"SOMETHING WRONG!!! No best jet removed from vectors."<<endl;

// sort the pt_jets_vector in descending order
std::sort(pt_jets_vector.begin(), pt_jets_vector.end(), std::greater<float>());

// create a vector of indices to keep track of the original position of elements in pt_jets_vector
std::vector<int> pt_sorted_indices(pt_jets_vector.size());
std::iota(pt_sorted_indices.begin(), pt_sorted_indices.end(), 0);
std::sort(pt_sorted_indices.begin(), pt_sorted_indices.end(), [&](int i, int j) {
    return pt_jets_vector[i] > pt_jets_vector[j];
});

// sort the p4_jets_vector using the indices in pt_sorted_indices
std::vector<LorentzVector> p4_jets_sorted(p4_jets_vector.size());
for (int i = 0; i < static_cast<int>(pt_sorted_indices.size()); i++) {
    p4_jets_sorted[i] = p4_jets_vector[pt_sorted_indices[i]];
}

for (int i = 0; i < static_cast<int>(p4_jets_vector.size()); i++) {
  float deltaEtaPhotonJet = fabs(ph_eta - p4_jets_vector[i].Eta());
  float deltaPhiPhotonJet = fabs(ph_phi - p4_jets_vector[i].Phi());
  if (deltaPhiPhotonJet > M_PI) deltaPhiPhotonJet = 2*M_PI - deltaPhiPhotonJet;
  float deltaR_photonJet = sqrt(deltaEtaPhotonJet * deltaEtaPhotonJet + deltaPhiPhotonJet + deltaPhiPhotonJet);
  if(verbose) cout<<"deltaR_photonJet = "<<deltaR_photonJet<<endl;
  if(verbose) cout<<"VBF jet pT       = "<<p4_jets_vector[i].Pt()<<endl;
  if (deltaR_photonJet < 0.3){
    p4_jets_vector.erase(p4_jets_vector.begin() + i);
    pt_jets_vector.erase(pt_jets_vector.begin() + i);
  }
}

bool isVBF = false;

int nExtraJet = static_cast<int>(p4_jets_vector.size());

if(nExtraJet > 1){

float deltaEtaJets = abs(p4_jets_vector[0].Eta() - p4_jets_vector[1].Eta());
float leadingJetPt = p4_jets_vector[0].Pt(); 
float mJJ          = (p4_jets_vector[0] + p4_jets_vector[1]).M();

if (debug){
cout<<"--------- VBF veto ---------------"<<endl;
cout<<"isVBF              = "<<isVBF<<endl;
cout<<"photon_eT          = "<<ph_eT<<endl;
cout<<"mesonPt            = "<<_bestCouplePt<<endl;
cout<<"bestJetPt          = "<<_bestJet_pT<<endl;
cout<<"nExtraJet          = "<<nExtraJet<<endl;
cout<<"eta leadingJet     = "<<p4_jets_vector[0].Eta()<<endl;
cout<<"eta subLeadingJet  = "<<p4_jets_vector[1].Eta()<<endl;
cout<<"deltaEtaJets       = "<<deltaEtaJets<<endl;
cout<<"leadingJetPt       = "<<leadingJetPt<<endl;
cout<<"subLeadingJetPt    = "<<p4_jets_vector[1].Pt()<<endl;
cout<<"mJJ                = "<<mJJ<<endl;
cout<<"---------------------------------"<<endl;
}

/*
//VBF definition : nJets >= 2 (without counting the jet containing the candidate meson), pT leading jet > 30, pT subleading jet > 20, deltaEta > 3, mJJ > 400
if(leadingJetPt > 30){
  isVBF = true;
  cout<<"VBF-like event: RETURN!"<<endl;
 }
else if (p4_jets_vector[0].Eta() * p4_jets_vector[1].Eta() < 0. && deltaEtaJets > 3.){
  isVBF = true;
  cout<<"VBF-like event: RETURN!"<<endl;
 } 
}
*/
/*
float deltaEtaPhotonJet1 = fabs(ph_eta - p4_jets_vector[0].Eta());
float deltaPhiPhotonJet1 = fabs(ph_phi - p4_jets_vector[0].Phi());
if (deltaPhiPhotonJet1 > M_PI) deltaPhiPhotonJet1 = 2*M_PI - deltaPhiPhotonJet1;
float deltaR_photonJet1 = sqrt(deltaEtaPhotonJet1 * deltaEtaPhotonJet1 + deltaPhiPhotonJet1 + deltaPhiPhotonJet1);

float deltaEtaPhotonJet2 = fabs(ph_eta - p4_jets_vector[1].Eta());
float deltaPhiPhotonJet2 = fabs(ph_phi - p4_jets_vector[1].Phi());
if (deltaPhiPhotonJet2 > M_PI) deltaPhiPhotonJet2 = 2*M_PI - deltaPhiPhotonJet2;
float deltaR_photonJet2 = sqrt(deltaEtaPhotonJet2 * deltaEtaPhotonJet2 + deltaPhiPhotonJet2 + deltaPhiPhotonJet2);

if(deltaEtaPhotonJet1 < 0.3 || deltaEtaPhotonJet2 < 0.3) 
*/

if((leadingJetPt < 30.) || (p4_jets_vector[0].Eta() * p4_jets_vector[1].Eta() > 0) || (p4_jets_vector[0].Eta() * p4_jets_vector[1].Eta() < 0 && deltaEtaJets < 3) || (mJJ < 300.)){
  isVBF = false;
}
else {
  isVBF = true;
  }

if(verbose){
cout<<"--------- Event n. "<<event_number<<" ---------------"<<endl;
cout<<"mesonMass       = "<<_PhiMass<<endl;
cout<<"Hmass           = "<<_PhiGammaMass<<endl;
cout<<"jet puID        = "<<_bestJet_PUid<<endl;
cout<<"photon_eT       = "<<ph_eT<<endl;
cout<<"photon_eta      = "<<ph_eta<<endl;
cout<<"mesonPt         = "<<_bestCouplePt<<endl;
cout<<"mesonEta        = "<<_bestCoupleEta<<endl;
cout<<"bestJetPt       = "<<_bestJet_pT<<endl;
cout<<"leadingJetPt    = "<<leadingJetPt<<endl;
cout<<"subLeadingJetPt = "<<p4_jets_vector[1].Pt()<<endl;
cout<<"eta_j1          = "<<p4_jets_vector[0].Eta()<<endl;
cout<<"eta_j2          = "<<p4_jets_vector[1].Eta()<<endl;
cout<<"deltaEta        = "<<deltaEtaJets<<endl;
cout<<"mJJ             = "<<mJJ<<endl;
cout<<"isVBF           = "<<isVBF<<endl;
  }
}

if(isVBF){

  cout<<"VBF-like event: RETURN!"<<endl;
  return;
}

_Nevents_VBFVeto++;





  //*************************************************************//
  //                                                             //
  //-------------------------- LHE access -----------------------//
  //                                                             //
  //*************************************************************//

 if(!runningOnData_){

  // Get the LHEEventProduct from the miniAOD
  edm::Handle<LHEEventProduct> lheEventProduct;
  iEvent.getByLabel("externalLHEProducer", lheEventProduct);

  // Get the central weight
  double centralWeight = lheEventProduct->originalXWGTUP();

  // Initialize variables for the minimum and maximum weights
  minPDFWeight = std::numeric_limits<double>::max();
  maxPDFWeight = -std::numeric_limits<double>::max();
  minQCDWeight = std::numeric_limits<double>::max();
  maxQCDWeight = -std::numeric_limits<double>::max();

  // Loop over the LHE weights and find the minimum and maximum PDF and QCD variations
  for (unsigned int i = 0; i < lheEventProduct->weights().size(); ++i) {
      std::string weightId = lheEventProduct->weights()[i].id;
      if (weightId.find("2") == 0) {
          int weightIdInt = std::stoi(weightId);
          if (weightIdInt >= 2001 && weightIdInt <= 2100) {
              double pdfWeight = lheEventProduct->weights()[i].wgt / centralWeight;
              if (pdfWeight < minPDFWeight) {
                  minPDFWeight = pdfWeight;
              }
              if (pdfWeight > maxPDFWeight) {
                  maxPDFWeight = pdfWeight;
              }
          }
          else if (weightIdInt >= 1001 && weightIdInt <= 1009) {
              double qcdWeight = lheEventProduct->weights()[i].wgt / centralWeight;
              if (qcdWeight < minQCDWeight) {
                  minQCDWeight = qcdWeight;
              }
              if (qcdWeight > maxQCDWeight) {
                  maxQCDWeight = qcdWeight;
              }
          }
      }
  }

  if(debug){
  // Print the minimum and maximum PDF and QCD variations
  cout << "Minimum PDF variation: " << minPDFWeight << endl;
  cout << "Maximum PDF variation: " << maxPDFWeight << endl;
  cout << "Minimum QCD variation: " << minQCDWeight << endl;
  cout << "Maximum QCD variation: " << maxQCDWeight << endl;
  }

}


  //*************************************************************//
  //                                                             //
  //----------------------- Access MC Truth ---------------------//
  //                                                             //
  //*************************************************************//

  //In signal, identify if there's a real mu or ele from W
  is_Kplus_matched   = false;
  is_Kminus_matched  = false;
  is_Piplus_matched  = false;
  is_Piminus_matched = false;
  //is_Phi_matched     = false;
  //is_Rho_matched     = false;
  is_photon_matched  = false;
  is_meson_matched   = false;
  is_Higgs_matched   = false; 

  Kminus_phi   = -999.;
  Kplus_phi    = -999.;
  float Piminus_phi  = -999.;
  float Piplus_phi   = -999.;
  Kminus_eta   = -999.;
  Kplus_eta    = -999.;
  float Piminus_eta  = -999.;
  float Piplus_eta   = -999.;
  deltaR_Kplus       = -999;
  deltaR_wrong       = -999;
  deltaR_Kminus      = -999.;
  deltaR_Piplus      = -999.;
  deltaR_Piminus     = -999.;
  genPhoton_eT       = -999.;
  genPhoton_eta      = -999.;
  genPhoton_phi      = -999.;
  genMeson_pT        = -999.;
  genMeson_m         = -999.;
  genMeson_eta       = -999.;
  genMeson_phi       = -999.;
  KplusPt      = -999.;
  KminusPt     = -999.;
  Kplus_dxy    = -999.;
  Kplus_dz     = -999.;
  Kminus_dxy   = -999.;
  Kminus_dz    = -999.;

  if(!runningOnData_){
    for(auto gen = prunedGenParticles->begin(); gen != prunedGenParticles->end(); ++gen){
      //if( gen->pdgId() == 333) cout<<"pdg ID = "<<gen->pdgId()<<" (mother ID = "<<gen->mother()->pdgId()<<")"<<endl;
      if( gen->pdgId() == 321  && gen->mother()->pdgId() == 333 && gen->mother()->mother()->pdgId() == 25)  Kplus_phi   = gen->phi(), Kplus_eta   = gen->eta(), KplusPt  = gen->pt();//, Kplus_dz  = gen->dz();//(&slimmedPV->at(0))->position()
      if( gen->pdgId() == -321 && gen->mother()->pdgId() == 333 && gen->mother()->mother()->pdgId() == 25)  Kminus_phi  = gen->phi(), Kminus_eta  = gen->eta(), KminusPt = gen->pt();//, Kminus_dxy = gen->dxy(), Kminus_dz = gen->dz();
      if( gen->pdgId() == 211  && gen->mother()->pdgId() == 113 && gen->mother()->mother()->pdgId() == 25)  Piplus_phi  = gen->phi(), Piplus_eta  = gen->eta();
      if( gen->pdgId() == -211 && gen->mother()->pdgId() == 113 && gen->mother()->mother()->pdgId() == 25)  Piminus_phi = gen->phi(), Piminus_eta = gen->eta();
      if( gen->pdgId() == 333  && gen->mother()->pdgId() == 25) genMeson_pT  = gen->pt(), genMeson_phi  = gen->phi(),  genMeson_eta = gen->eta(), genMeson_m = gen->mass();
      if( gen->pdgId() == 113  && gen->mother()->pdgId() == 25) genMeson_pT  = gen->pt(), genMeson_phi  = gen->phi(),  genMeson_eta = gen->eta(), genMeson_m = gen->mass();
      if( gen->pdgId() == 313  && gen->mother()->pdgId() == 25) genMeson_pT  = gen->pt(), genMeson_phi  = gen->phi(),  genMeson_eta = gen->eta(), genMeson_m = gen->mass();
      if( gen->pdgId() == 22   && gen->mother()->pdgId() == 25) genPhoton_eT = gen->pt(), genPhoton_phi = gen->phi(), genPhoton_eta = gen->eta();
    }
  } 

//MC TRUTH CHECK
if(!runningOnData_) //ONLY FOR MC START
{
    //photon matching -----------------------------------------
    float deltaPhiPhoton = fabs(ph_phi - genPhoton_phi);
    if (deltaPhiPhoton > M_PI) deltaPhiPhoton = 2*M_PI - deltaPhiPhoton;

    float deltaR_photonGenVsReco = sqrt((ph_eta - genPhoton_eta) * (ph_eta - genPhoton_eta) + deltaPhiPhoton * deltaPhiPhoton);
    //float photonRelPt = 
    if (deltaR_photonGenVsReco < 0.2) is_photon_matched = true;

    //meson matching -----------------------------------------
    float deltaPhiMeson = fabs(_bestCouplePhi - genMeson_phi);
    if (deltaPhiMeson > M_PI) deltaPhiMeson = 2*M_PI - deltaPhiMeson;

    float deltaR_mesonGenVsReco = sqrt((_bestCoupleEta - genMeson_eta) * (_bestCoupleEta - genMeson_eta) + deltaPhiMeson * deltaPhiMeson);
    if (deltaR_mesonGenVsReco < 0.3) is_meson_matched = true;    

    //meson pT matching -----------------------------------------
    if (_bestCouplePt < genMeson_pT - 4. || _bestCouplePt > genMeson_pT + 4.) _Nevents_MesonPtNotMatched ++; 

    //Higgs matching -----------------------------------------    
    if(is_photon_matched && is_meson_matched){

      if(verbose) cout<<endl<<"**************** HIGGS FOUND ******************"<<endl;
      if(verbose) cout<<"Higgs deltaR = "<<deltaR<<endl;
      _Nevents_HiggsMatched++;
      is_Higgs_matched=true;
    }
    
    else 
    {
      _Nevents_HiggsNotMatched++;
      if(verbose) cout<<endl<<"THAT'S NOT A HIGGS!"<<endl;
    
    }

//if is PhiGamma event
  if (_isPhi){

    //First cand positive and second cand negative
    if (_firstCandCharge > 0){ 
      //phi angle folding K plus
      float deltaPhi_Kplus = fabs(_firstCandPhi - Kplus_phi);
      float deltaPhi_wrong = fabs(_secondCandPhi - Kplus_phi);
      if (deltaPhi_Kplus > M_PI) deltaPhi_Kplus = 2*M_PI - deltaPhi_Kplus;
      if (deltaPhi_wrong > M_PI) deltaPhi_wrong = 2*M_PI - deltaPhi_wrong;

      //deltaR K plus
      deltaR_Kplus = sqrt((_firstCandEta - Kplus_eta) * (_firstCandEta - Kplus_eta) + deltaPhi_Kplus * deltaPhi_Kplus);
      deltaR_wrong = sqrt((_secondCandEta - Kplus_eta) * (_secondCandEta - Kplus_eta) + deltaPhi_wrong * deltaPhi_wrong);
      
      //phi angle folding K minus
      float deltaPhi_Kminus = fabs(_secondCandPhi - Kminus_phi);
      if (deltaPhi_Kminus > M_PI) deltaPhi_Kminus = 2*M_PI - deltaPhi_Kminus;

      //deltaR K minus
      deltaR_Kminus = sqrt((_secondCandEta - Kminus_eta) * (_secondCandEta - Kminus_eta) + deltaPhi_Kminus * deltaPhi_Kminus);
      cout<<endl;
      //if (_firstCandPt < 0.95*KplusPt && _firstCandPt > 1.05*KplusPt) cout<<"firstCand pT not matched"<<endl;

    }

    else{ //Second cand positive and first cand negative
      
      //phi angle folding K plus
      float deltaPhi_Kplus = fabs(_secondCandPhi - Kplus_phi);
      if (deltaPhi_Kplus > M_PI) deltaPhi_Kplus = 2*M_PI - deltaPhi_Kplus;

      //deltaR K plus
      deltaR_Kplus = sqrt((_secondCandEta - Kplus_eta) * (_secondCandEta - Kplus_eta) + deltaPhi_Kplus * deltaPhi_Kplus);

      cout<<endl;
//      if (_firstCandPt < 0.95*KminusPt && _firstCandPt > 1.05*KminusPt) cout<<"firstCandPt pT not matched"<<endl;

      }

      //phi angle folding K minus
      float deltaPhi_Kminus = fabs(_firstCandPhi - Kminus_phi);
      if (deltaPhi_Kminus > M_PI) deltaPhi_Kminus = 2*M_PI - deltaPhi_Kminus;

      //deltaR K minus
      deltaR_Kminus = sqrt((_firstCandEta - Kminus_eta) * (_firstCandEta - Kminus_eta) + deltaPhi_Kminus * deltaPhi_Kminus);


  } //if isPhi END
  
  else{ //RhoGamma event

    //First cand positive and second cand negative
    if (_firstCandCharge > 0){ 
      //phi angle folding Pi plus
      float deltaPhi_Piplus = fabs(_firstCandPhi - Piplus_phi);
      if (deltaPhi_Piplus > M_PI) deltaPhi_Piplus = 2*M_PI - deltaPhi_Piplus;

      //deltaR Pi plus
      deltaR_Piplus = sqrt((_firstCandEta - Piplus_eta) * (_firstCandEta - Piplus_eta) + deltaPhi_Piplus * deltaPhi_Piplus);

      //phi angle folding Pi minus
      float deltaPhi_Piminus = fabs(_secondCandPhi - Piminus_phi);
      if (deltaPhi_Piminus > M_PI) deltaPhi_Piminus = 2*M_PI - deltaPhi_Piminus;

      //deltaR Pi minus
      deltaR_Piminus = sqrt((_secondCandEta - Piminus_eta) * (_secondCandEta - Piminus_eta) + deltaPhi_Piminus * deltaPhi_Piminus);
      
    }

    else{ //Second cand positive and first cand negative
      
      //phi angle folding Pi plus
      float deltaPhi_Piplus = fabs(_secondCandPhi - Piplus_phi);
      if (deltaPhi_Piplus > M_PI) deltaPhi_Piplus = 2*M_PI - deltaPhi_Piplus;

      //deltaR Pi plus
      deltaR_Piplus = sqrt((_secondCandEta - Piplus_eta) * (_secondCandEta - Piplus_eta) + deltaPhi_Piplus * deltaPhi_Piplus);
      
      
      //phi angle folding Pi minus
      float deltaPhi_Piminus = fabs(_firstCandPhi - Piminus_phi);
      if (deltaPhi_Piminus > M_PI) deltaPhi_Piminus = 2*M_PI - deltaPhi_Piminus;

      //deltaR Pi minus
      deltaR_Piminus = sqrt((_firstCandEta - Piminus_eta) * (_firstCandEta - Piminus_eta) + deltaPhi_Piminus * deltaPhi_Piminus);

    }
  } // if isRho END



    //phi angle folding
    //float deltaPhi_Kpm = fabs(Kplus_phi - Kminus_phi);
    //if (deltaPhi_Kpm > M_PI) deltaPhi_Kpm = 2*M_PI - deltaPhi_Kpm;

    //float deltaR_Kpm = sqrt((Kplus_eta - Kminus_eta) * (Kplus_eta - Kminus_eta) + deltaPhi_Kpm * deltaPhi_Kpm);

    //cout<<"deltaR_Kpm   = "<<deltaR_Kpm<<endl;


    //For the polarization reweighting
    theta_pol = 0.;

    if(!runningOnData_){

    TLorentzVector mu[2];

    for(auto genParticle = prunedGenParticles->begin(); genParticle != prunedGenParticles->end(); ++genParticle){
     // last in the collection
     if( genParticle->pdgId()==321 && genParticle->mother() && genParticle->mother()->pdgId()==333) {
       const reco::GenParticle* phi = dynamic_cast<const reco::GenParticle*>(genParticle->mother());
       if( phi->mother() && phi->mother()->pdgId()==25) {
    mu[1].SetPxPyPzE(genParticle->px(),genParticle->py(),
         genParticle->pz(),genParticle->energy());
    mu[0].SetPxPyPzE(phi->px(),phi->py(),
         phi->pz(),phi->energy());
    break;
          }
        }
      }

      TVector3 phiBoost = mu[0].BoostVector();
      mu[1].Boost(-phiBoost);
      theta_pol = mu[0].Vect().Angle(mu[1].Vect());
    }


    //some prints
    
    if(verbose){
      cout<<"Photon eT = "<<ph_eT<<endl;
      //cout<<"ph_en_sigmaUP = "<< ph_en_sigmaUP<<endl;
      //cout<<"ph_en_sigmaDW = "<< ph_en_sigmaDW<<endl;
      //cout<<"ph_en_scaleUP = "<<ph_en_scaleUP<<endl;
      //cout<<"ph_en_scaleDW = "<<ph_en_scaleDW<<endl;
      cout<<"n Jets = "<<nJets25<<endl;
      cout<<"Jet + photon inv. mass = "<<_bestJet_Photon_invMass<<endl;
      cout<<"n. of daughters: "<<_bestJet_nDaughters<<endl;
      cout<<"Best couple pT = "<<_firstCandPt + _secondCandPt<<endl;
      cout<<"_firstCandPt   = "<<_firstCandPt<<endl;
      cout<<"_secondCandPt  = "<<_secondCandPt<<endl;
      cout<<"genMeson_pT    = "<<genMeson_pT<<endl;
      cout<<"genKplus pT    = "<<KplusPt<<endl;
      cout<<"genKminus pT   = "<<KminusPt<<endl;
      cout<<"trk1 dxy       = "<<bestFirstCand_dxy<<endl;
      cout<<"trk2 dxy       = "<<bestSecondCand_dxy<<endl;
      cout<<"trk1 dxyErr    = "<<bestFirstCand_dxyErr<<endl;
      cout<<"trk2 dxyErr    = "<<bestSecondCand_dxyErr<<endl;
      cout<<"deltaDxy       = "<<abs(bestFirstCand_dxy - bestSecondCand_dxy)<<endl;
      cout<<"trk1 dz        = "<<bestFirstCand_dz<<endl;
      cout<<"trk2 dz        = "<<bestSecondCand_dz<<endl;
      cout<<"trk1 dzErr     = "<<bestFirstCand_dzErr<<endl;
      cout<<"trk2 dzErr     = "<<bestSecondCand_dzErr<<endl;
      cout<<"deltaDz        = "<<abs(bestFirstCand_dz - bestSecondCand_dz)<<endl;
      cout<<"Kplus  dxy     = "<<Kplus_dxy<<endl;
      cout<<"Kplus  dz      = "<<Kplus_dz<<endl;
      cout<<"Kminus dxy     = "<<Kminus_dxy<<endl;
      cout<<"Kminus dz      = "<<Kminus_dz<<endl;
      cout<<"Best couple DeltaR = "<<deltaR_KChosen<<endl;
      //cout<<"Meson candidate inv. mass  = "<<_MesonMass<<endl;
      cout<<"isPhi = "<<isPhi<<" and isRho = "<<isRho<<endl;
      cout<<"isFirstCandK = "<<_isFirstCandK<<endl;
      //cout<<"H inv. mass = "<<_Hmass_From2K_Photon<<endl;
      cout<<"--------------------------------------------------"<<endl;
      cout<<"MC Higgs found = "<<_Nevents_HiggsMatched<<",   Higgs NOT matched = "<<_Nevents_HiggsNotMatched<<",   mesonPt not matched = "<<_Nevents_MesonPtNotMatched<<endl;
      cout<<"--------------------------------------------------"<<endl<<endl;
      }
     
  }  //ONLY FOR MC END 
 
 else //ONLY FOR DATA
 {
  cout<<"CANDIDATE HIGGS FOUND IN DATA: EVENT RECORDED!"<<endl;
  if(debug){
cout<<"Photon eT = "<<ph_eT<<endl;
cout<<"ph_en_sigmaUP = "<< ph_en_sigmaUP<<endl;
cout<<"ph_en_sigmaDW = "<< ph_en_sigmaDW<<endl;
cout<<"ph_en_scaleUP = "<<ph_en_scaleUP<<endl;
cout<<"ph_en_scaleDW = "<<ph_en_scaleDW<<endl;
  }
 }
   is_hltL1sBigORMu18erTauXXer2p1 = false;
   is_hltEG35R9Id90HE10IsoMEcalIsoFilter  = false;
   is_hltEG35R9Id90HE10IsoMHcalIsoFilter  = false;
   is_hltEG35R9Id90HE10IsoMTrackIsoFilter = false;
   is_hltOverlapFilterPhoton35MediumChargedIsoPFTau35 = false;
   is_hltEG35R9Id90HE10IsoMR9Filter = false;
   is_hltEG35R9Id90HE10IsoMHEFilter = false;
   is_hltEG35R9Id90HE10IsoMEtFilter = false;
   is_hltEGL1EGAndTauFilter = false;
   is_hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3 = false;
   is_hltL1sSingleMu22 = false; //TwoProngs leg
   is_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23 = false; //Photon leg

    int nMatchedObjs = 0;
    std::vector<const pat::TriggerObjectStandAlone*> matchedTrigObjs;

 //Trigger matching
 for (pat::TriggerObjectStandAlone obj : *triggerObjects){ // note: not "const &" since we want to call unpackPathNames
   
   //bool isAcceptedPath = false;
   obj.unpackPathNames(names);
   
   std::vector pathNamesAll = obj.pathNames(false);
   
   for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
     // Record also if the object is associated to a 'l3' filter (always true for the definition used
     // in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which means
     // that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
     bool isSuccessfulTrigger = obj.hasPathName( pathNamesAll[h], true, true );
     //std::cout << "   " << pathNamesAll[h];
     if(!isSuccessfulTrigger) continue;
     //std::cout << pathNamesAll[h] << "(L,3)" << std::endl; 
      
     //if(pathNamesAll[h].find("HLT_Photon35_TwoProngs35_v") != std::string::npos) {
     //if(!isTwoProngTrigger){
     //now match ALL objects in a cone of DR<0.1
     //it is important to match all objects as there are different ways to reconstruct the same electron
     //eg, L1 seeded, unseeded, as a jet etc
     //and so you want to be sure you get all possible objects
     matchedTrigObjs = getMatchedObjs(ph_eta,ph_phi,_bestCoupleEta,_bestCouplePhi,unpackedTrigObjs,0.1);
      }
    }
  //}   

  //if(!isTwoProngTrigger){

     for(const auto trigObj : matchedTrigObjs){
      nMatchedObjs++;
      
      auto filterIds = trigObj->filterIds();
      for (int id : filterIds) {
        if(trig_verbose) cout<<"Matched obj n."<<nMatchedObjs<<" with ID "<<id<<endl;
      }      
      //now just check if it passes the filters

      //HLT_Photon35_TwoProngs35 filters -----------------------------------------------------      
      if(trigObj->hasFilterLabel("hltL1sBigORMu18erTauXXer2p1")) {
        is_hltL1sBigORMu18erTauXXer2p1 = true;
        if(trig_verbose) cout<<"1) hltL1sBigORMu18erTauXXer2p1 passed!"<<endl;
      }
      if(trigObj->hasFilterLabel("hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3")) {
        is_hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3 = true;
        if(trig_verbose) cout<<"1) hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3 passed!"<<endl;
      }
      if(trigObj->hasFilterLabel("hltEGL1EGAndTauFilter")) {
        is_hltEGL1EGAndTauFilter = true;
        if(trig_verbose) cout<<"2) hltEGL1EGAndTauFilter passed!"<<endl;
      }
      if(trigObj->hasFilterLabel("hltEG35R9Id90HE10IsoMEtFilter")) {
        is_hltEG35R9Id90HE10IsoMEtFilter = true;
        if(trig_verbose) cout<<"3) hltEG35R9Id90HE10IsoMEtFilter passed!"<<endl;
      }
      if(trigObj->hasFilterLabel("hltEG35R9Id90HE10IsoMHEFilter")) {
        is_hltEG35R9Id90HE10IsoMHEFilter = true;
        if(trig_verbose) cout<<"4) hltEG35R9Id90HE10IsoMHEFilter passed!"<<endl;
      }
      if(trigObj->hasFilterLabel("hltEG35R9Id90HE10IsoMR9Filter")) {
        is_hltEG35R9Id90HE10IsoMR9Filter = true;
        if(trig_verbose) cout<<"5) hltEG35R9Id90HE10IsoMR9Filter passed!"<<endl;
      }
       if(trigObj->hasFilterLabel("hltEG35R9Id90HE10IsoMEcalIsoFilter")) {
         is_hltEG35R9Id90HE10IsoMEcalIsoFilter = true;
         if(trig_verbose) cout<<"6) hltEG35R9Id90HE10IsoMEcalIsoFilter passed!"<<endl;
       }
       if(trigObj->hasFilterLabel("hltEG35R9Id90HE10IsoMHcalIsoFilter")){
         is_hltEG35R9Id90HE10IsoMHcalIsoFilter = true;
         if(trig_verbose) cout<<"7) hltEG35R9Id90HE10IsoMHcalIsoFilter passed!"<<endl;
       }
       if(trigObj->hasFilterLabel("hltEG35R9Id90HE10IsoMTrackIsoFilter")){
         is_hltEG35R9Id90HE10IsoMTrackIsoFilter = true;
         if(trig_verbose) cout<<"8) hltEG35R9Id90HE10IsoMTrackIsoFilter passed!"<<endl;
       }
       if(trigObj->hasFilterLabel("hltOverlapFilterPhoton35MediumChargedIsoPFTau35")){
         is_hltOverlapFilterPhoton35MediumChargedIsoPFTau35 = true;
         if(trig_verbose) cout<<"9) hltOverlapFilterPhoton35MediumChargedIsoPFTau35 passed!"<<endl;
       }
      
      //HLT_Mu24_TwoProngs35 filters -----------------------------------------------------      
      if(trigObj->hasFilterLabel("hltL1sSingleMu22")) {
        is_hltL1sSingleMu22 = true;
        if(trig_verbose) cout<<"1) hltL1sSingleMu22 passed!"<<endl;
      }

      //HLT_Mu17_Photon30_IsoCaloId filters -----------------------------------------------------      
      if(trigObj->hasFilterLabel("hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23")) {
        is_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23 = true;
        if(trig_verbose) cout<<"1) hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23 passed!"<<endl;
      }

      cout<<endl;

     }
   //}
  
    if(is_hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3) _NeventsFilter1++;
    if(is_hltEGL1EGAndTauFilter) _NeventsFilter2++;
    if(is_hltEG35R9Id90HE10IsoMEtFilter) _NeventsFilter3++;
    if(is_hltEG35R9Id90HE10IsoMHEFilter) _NeventsFilter4++;
    if(is_hltEG35R9Id90HE10IsoMR9Filter) _NeventsFilter5++;
    if(is_hltEG35R9Id90HE10IsoMEcalIsoFilter ) _NeventsFilter6++;
    if(is_hltEG35R9Id90HE10IsoMHcalIsoFilter ) _NeventsFilter7++;
    if(is_hltEG35R9Id90HE10IsoMTrackIsoFilter) _NeventsFilter8++;
    if(is_hltOverlapFilterPhoton35MediumChargedIsoPFTau35) _NeventsFilter9++;

    cout<<"n event filter 1 : "<<_NeventsFilter1<<endl;
    cout<<"n event filter 2 : "<<_NeventsFilter2<<endl;
    cout<<"n event filter 3 : "<<_NeventsFilter3<<endl;
    cout<<"n event filter 4 : "<<_NeventsFilter4<<endl;
    cout<<"n event filter 5 : "<<_NeventsFilter5<<endl;
    cout<<"n event filter 6 : "<<_NeventsFilter6<<endl;
    cout<<"n event filter 7 : "<<_NeventsFilter7<<endl;
    cout<<"n event filter 8 : "<<_NeventsFilter8<<endl;
    cout<<"n event filter 9 : "<<_NeventsFilter9<<endl;
    cout<<"n event filter 10: "<<_NeventsFilter10<<endl;

  //cout<<endl<<"Event n = "<<event_number<<endl;
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
    mytree->Branch("event_number",&event_number);
  }
  else{
    mytree->Branch("event_number",&event_number);
  }

  mytree->Branch("isTwoProngTrigger",&isTwoProngTrigger);
  mytree->Branch("nMuons10",&nMuons10);
  mytree->Branch("nMuons20",&nMuons20);
  mytree->Branch("nElectrons10",&nElectrons10);
  mytree->Branch("nElectrons20",&nElectrons20);
  mytree->Branch("nPhotons38WP80",&nPhotons38WP80);
  mytree->Branch("nPhotonsWP90_pT20_2p5eta3p0",&nPhotonsWP90_pT20_2p5eta3p0);
  mytree->Branch("nPhotonsWP90_pT15_barrel",&nPhotonsWP90_pT15_barrel);
  mytree->Branch("nPhotonsWP90_pT25_endcap",&nPhotonsWP90_pT25_endcap);
  mytree->Branch("nPhotonsChosen",&nPhotonsChosen);
  mytree->Branch("nJets30",&nJets30);
  mytree->Branch("nJets25",&nJets25);
  mytree->Branch("met_pT",&met_pT);
  mytree->Branch("metpuppi_pT",&metpuppi_pT);

  mytree->Branch("photon_eT",&ph_eT);
  mytree->Branch("photon_eT_sigmaDW",&ph_en_sigmaDW);
  mytree->Branch("photon_eT_sigmaUP",&ph_en_sigmaUP);
  mytree->Branch("photon_eT_scaleDW",&ph_en_scaleDW);
  mytree->Branch("photon_eT_scaleUP",&ph_en_scaleUP);
  mytree->Branch("photon_eta",&ph_eta);
  mytree->Branch("photon_etaSC",&ph_etaSC);
  mytree->Branch("photon_phi",&ph_phi);
  mytree->Branch("photon_HoverE",&ph_HoverE);
  mytree->Branch("photon_hadTowOverEm",&ph_hadTowOverEm);

  //mytree->Branch("photon_iso_ChargedHadron",&ph_iso_ChargedHadron);
  //mytree->Branch("photon_iso_NeutralHadron",&ph_iso_NeutralHadron);
  //mytree->Branch("photon_iso_Photon",&ph_iso_Photon);
  mytree->Branch("photon_iso_eArho",&ph_iso_eArho);
  mytree->Branch("photonRegressionError",&photonRegressionError);
  mytree->Branch("isPhotonTriggerMatched",&isPhotonTriggerMatched);

  mytree->Branch("bestJet_pT",&_bestJet_pT);
  mytree->Branch("bestJet_eta",&_bestJet_eta);
  mytree->Branch("bestJet_phi",&_bestJet_phi);
  mytree->Branch("bestJet_nDaughters",&_bestJet_nDaughters);
  mytree->Branch("bestJet_chargedEmEnergy",&_bestJet_chargedEmEnergy);
  mytree->Branch("bestJet_neutralEmEnergy",&_bestJet_neutralEmEnergy);
  mytree->Branch("bestJet_chargedHadEnergy",&_bestJet_chargedHadEnergy);
  mytree->Branch("bestJet_neutralHadEnergy",&_bestJet_neutralHadEnergy);
  mytree->Branch("bestJet_chargedEmEnergyFraction",&_bestJet_chargedEmEnergyFraction);
  mytree->Branch("bestJet_neutralEmEnergyFraction",&_bestJet_neutralEmEnergyFraction);
  mytree->Branch("bestJet_chargedHadEnergyFraction",&_bestJet_chargedHadEnergyFraction);
  mytree->Branch("bestJet_neutralHadEnergyFraction",&_bestJet_neutralHadEnergyFraction);
  mytree->Branch("bestJet_invMass",&_bestJet_invMass);
  mytree->Branch("bestJet_Photon_invMass",&_bestJet_Photon_invMass);
  mytree->Branch("bestJet_JECunc",&_bestJet_JECunc);
  mytree->Branch("bestJet_PUid",&_bestJet_PUid);
  mytree->Branch("isJetTriggerMatched",&isJetTriggerMatched);

  mytree->Branch("firstCandCharge",&_firstCandCharge);
  mytree->Branch("firstCandPt",&_firstCandPt);
  mytree->Branch("firstCandEta",&_firstCandEta);
  mytree->Branch("firstCandPhi",&_firstCandPhi);
  mytree->Branch("firstCand_dxy",&bestFirstCand_dxy);
  mytree->Branch("firstCand_dz",&bestFirstCand_dz);
  mytree->Branch("firstCand_dxyErr",&bestFirstCand_dxyErr);
  mytree->Branch("firstCand_dzErr",&bestFirstCand_dzErr);
  mytree->Branch("secondCandCharge",&_secondCandCharge);
  mytree->Branch("secondCandPt",&_secondCandPt);
  mytree->Branch("secondCandEta",&_secondCandEta);
  mytree->Branch("secondCandPhi",&_secondCandPhi);
  mytree->Branch("secondCand_dxy",&bestSecondCand_dxy);
  mytree->Branch("secondCand_dz",&bestSecondCand_dz);
  mytree->Branch("secondCand_dxyErr",&bestSecondCand_dxyErr);
  mytree->Branch("secondCand_dzErr",&bestSecondCand_dzErr);
  mytree->Branch("bestCouplePt",&_bestCouplePt);
  mytree->Branch("bestCoupleEta",&_bestCoupleEta);
  mytree->Branch("bestCouplePhi",&_bestCouplePhi);
  mytree->Branch("isPhi",&_isPhi);
  mytree->Branch("isRho",&_isRho);
  mytree->Branch("isK0star",&_isK0star); 
  mytree->Branch("isFirstCandK",&_isFirstCandK);   

  mytree->Branch("firstCandEnergy",&firstCandEnergy);
  mytree->Branch("secondCandEnergy",&secondCandEnergy);

  //mytree->Branch("MesonMass",&_MesonMass);
  //mytree->Branch("Hmass_From2K_Photon",&_Hmass_From2K_Photon);
  mytree->Branch("_PhiMass",&_PhiMass);
  mytree->Branch("_RhoMass",&_RhoMass);
  mytree->Branch("_K0starMass",&_K0starMass);
  mytree->Branch("_PhiGammaMass",&_PhiGammaMass);
  mytree->Branch("_RhoGammaMass",&_RhoGammaMass);
  mytree->Branch("_K0starGammaMass",&_K0starGammaMass);

  //mytree->Branch("K1_sum_pT_05",&K1_sum_pT_05);
  //mytree->Branch("K1_sum_pT_05_ch",&K1_sum_pT_05_ch);
  //mytree->Branch("K2_sum_pT_05",&K2_sum_pT_05);
  //mytree->Branch("K2_sum_pT_05_ch",&K2_sum_pT_05_ch);
  mytree->Branch("Couple_sum_pT_05",&couple_sum_pT_05);
  mytree->Branch("Couple_sum_pT_05_ch",&couple_sum_pT_05_ch);

  mytree->Branch("iso_K1",&_iso_K1);
  mytree->Branch("iso_K1_ch",&_iso_K1_ch);
  mytree->Branch("iso_K2",&_iso_K2);
  mytree->Branch("iso_K2_ch",&_iso_K2_ch);
  mytree->Branch("iso_couple",&_iso_couple);
  mytree->Branch("iso_couple_ch",&_iso_couple_ch);

  mytree->Branch("is_hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3",&is_hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3);
  mytree->Branch("is_hltEGL1EGAndTauFilter",&is_hltEGL1EGAndTauFilter);
  mytree->Branch("is_hltEG35R9Id90HE10IsoMEtFilter",&is_hltEG35R9Id90HE10IsoMEtFilter);
  mytree->Branch("is_hltEG35R9Id90HE10IsoMHEFilter",&is_hltEG35R9Id90HE10IsoMHEFilter);
  mytree->Branch("is_hltEG35R9Id90HE10IsoMR9Filter",&is_hltEG35R9Id90HE10IsoMR9Filter);
  mytree->Branch("is_hltEG35R9Id90HE10IsoMEcalIsoFilter",&is_hltEG35R9Id90HE10IsoMEcalIsoFilter);
  mytree->Branch("is_hltEG35R9Id90HE10IsoMHcalIsoFilter",&is_hltEG35R9Id90HE10IsoMHcalIsoFilter);
  mytree->Branch("is_hltEG35R9Id90HE10IsoMTrackIsoFilter",&is_hltEG35R9Id90HE10IsoMTrackIsoFilter);
  mytree->Branch("is_hltOverlapFilterPhoton35MediumChargedIsoPFTau35",&is_hltOverlapFilterPhoton35MediumChargedIsoPFTau35);
  mytree->Branch("is_hltL1sSingleMu22",&is_hltL1sSingleMu22); //TwoProngs leg
  mytree->Branch("is_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23",&is_hltL1sMu5EG23IorMu5IsoEG20IorMu7EG23IorMu7IsoEG20IorMuIso7EG23); //Photon leg
 
  if(isTriggerStudies){ //start trigger studies
    mytree->Branch("isIsoMuTrigger",&isIsoMuTrigger);
    mytree->Branch("isBestMuMu_Found",&isBestMuMu_Found);
    mytree->Branch("MuMuPt",&bestMuMuPt);
    mytree->Branch("MuMuMass",&bestMuMuMass);
  } //end trigger studies

  //Save MC info
  if(!runningOnData_){ //NO INFO FOR DATA
    mytree->Branch("PU_Weight",&PU_Weight);
    mytree->Branch("MC_Weight",&MC_Weight);
    mytree->Branch("theta_pol",&theta_pol);
    mytree->Branch("minPDFWeight",&minPDFWeight);
    mytree->Branch("maxPDFWeight",&maxPDFWeight);
    mytree->Branch("minQCDWeight",&minQCDWeight);
    mytree->Branch("maxQCDWeight",&maxQCDWeight);
    mytree->Branch("isHiggsMatched",&is_Higgs_matched);
    //mytree->Branch("isKplusMatched",&is_Kplus_matched);
    //mytree->Branch("isKminusMatched",&is_Kminus_matched);
    //mytree->Branch("isPiplusMatched",&is_Piplus_matched);
    //mytree->Branch("isPiminusMatched",&is_Piminus_matched);
    //mytree->Branch("isPhiFromH",&is_Phi_fromH);
    //mytree->Branch("isRhofromH",&is_Rho_fromH);
    //mytree->Branch("isPhotonFromH",&is_Photon_fromH);
    //mytree->Branch("isPhotonTrue",&is_photon_a_photon);
    mytree->Branch("isPhotonMatched",&is_photon_matched);
    mytree->Branch("genPhoton_eT",&genPhoton_eT);
    mytree->Branch("isMesonMatched",&is_meson_matched);
    mytree->Branch("genMeson_pT",&genMeson_pT);
    mytree->Branch("genMeson_m",&genMeson_m);
    mytree->Branch("KplusPt",&KplusPt);
    mytree->Branch("KminusPt",&KminusPt);
    mytree->Branch("Kminus_eta",&Kminus_eta);
    mytree->Branch("Kplus_eta",&Kplus_eta);
    mytree->Branch("Kminus_phi",&Kminus_phi);
    mytree->Branch("Kplus_phi",&Kplus_phi);

    mytree->Branch("deltaR_Kplus",&deltaR_Kplus);
    mytree->Branch("deltaR_wrong",&deltaR_wrong);
    mytree->Branch("deltaR_Kminus",&deltaR_Kminus);
    mytree->Branch("deltaR_Piplus",&deltaR_Piplus);
    mytree->Branch("deltaR_Piminus",&deltaR_Piminus);

}

}

void HPhiGammaAnalysis::beginJob()
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

void HPhiGammaAnalysis::endJob() 
{
h_Events->Fill(0.5,_Nevents_processed);
h_Events->Fill(1.5,_Nevents_triggered);
h_Events->Fill(2.5,_Nevents_isPhoton);
h_Events->Fill(3.5,_Nevents_bestCoupleFound);  
h_Events->Fill(4.5,_Nevents_candPtFilter);  
h_Events->Fill(5.5,_Nevents_coupleIsolationFilter);
h_Events->Fill(6.5,_Nevents_VBFVeto);  
h_Events->GetXaxis()->SetBinLabel(1,"processed");
h_Events->GetXaxis()->SetBinLabel(2,"triggered");
h_Events->GetXaxis()->SetBinLabel(3,"best photon");
h_Events->GetXaxis()->SetBinLabel(4,"best pair");
h_Events->GetXaxis()->SetBinLabel(5,"trks pT");
h_Events->GetXaxis()->SetBinLabel(6,"trks iso");
h_Events->GetXaxis()->SetBinLabel(7,"VBF veto");

float scale_factor = 1./_Nevents_VBFVeto;
cout<<"scale_factor = "<<scale_factor<<endl;
h_TriggerFilters->Fill(0.5,_Nevents_VBFVeto * scale_factor);
h_TriggerFilters->Fill(1.5,_NeventsFilter1 * scale_factor);
h_TriggerFilters->Fill(2.5,_NeventsFilter2 * scale_factor);
h_TriggerFilters->Fill(3.5,_NeventsFilter3 * scale_factor);
h_TriggerFilters->Fill(4.5,_NeventsFilter4 * scale_factor);
h_TriggerFilters->Fill(5.5,_NeventsFilter5 * scale_factor);
h_TriggerFilters->Fill(6.5,_NeventsFilter6 * scale_factor);
h_TriggerFilters->Fill(7.5,_NeventsFilter7 * scale_factor);
h_TriggerFilters->Fill(8.5,_NeventsFilter8 * scale_factor);
h_TriggerFilters->Fill(9.5,_NeventsFilter9 * scale_factor);
h_TriggerFilters->Fill(10.5,_NeventsFilter10 * scale_factor);
h_TriggerFilters->GetXaxis()->SetBinLabel(1,"offline selection");
h_TriggerFilters->GetXaxis()->SetBinLabel(2,"hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3");
h_TriggerFilters->GetXaxis()->SetBinLabel(3,"hltEGL1EGAndTauFilter");
h_TriggerFilters->GetXaxis()->SetBinLabel(4,"hltEG35R9Id90HE10IsoMEtFilter");
h_TriggerFilters->GetXaxis()->SetBinLabel(5,"hltEG35R9Id90HE10IsoMHEFilter");
h_TriggerFilters->GetXaxis()->SetBinLabel(6,"hltEG35R9Id90HE10IsoMR9Filter");
h_TriggerFilters->GetXaxis()->SetBinLabel(7,"hltEG35R9Id90HE10IsoMEcalIsoFilter");
h_TriggerFilters->GetXaxis()->SetBinLabel(8,"hltEG35R9Id90HE10IsoMHcalIsoFilter");
h_TriggerFilters->GetXaxis()->SetBinLabel(9,"hltEG35R9Id90HE10IsoMTrackIsoFilter");
h_TriggerFilters->GetXaxis()->SetBinLabel(10,"hltOverlapFilterPhoton35MediumChargedIsoPFTau35");
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPhiGammaAnalysis);
