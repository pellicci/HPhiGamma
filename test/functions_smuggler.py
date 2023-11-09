#######################################################################
#                                                                     #
# <<  Hi! I smuggle bad coded functions. But, ehy, they're free! >>   #
#                                                                     #
#######################################################################                                                                                                                                                       

import ROOT
import math
import os
from array import array
import time


#------- Arrays and reader for the BDT -------#

trk1isoCh_array            = array('f', [0.])
#pairIso_array              = array('f', [0.])
mesonPt_array              = array('f', [0.])
photonEt_array             = array('f', [0.])
#photonEta_array            = array('f', [0.])
mesonEta_array             = array('f', [0.])
#nJet_array                 = array('f', [0.])
#JetNeutralEmEnergy_array   = array('f', [0.])
#JetChargedHadEnergy_array  = array('f', [0.])
#JetNeutralHadEnergy_array  = array('f', [0.])

#metPt_array      = array('f', [0.])
#bestJetPt_array  = array('f', [0.])
#dPhiGammaTrk_array = array('f', [0.])

reader = ROOT.TMVA.Reader("!Color")

#MCfromDATA correction for photons ###############################################################################
ph_ID_scale_name_2018  = "scale_factors/egammaEffi.txt_EGM2D_Pho_wp80.root_UL18.root"
ph_ID_scale_file_2018  = ROOT.TFile(ph_ID_scale_name_2018)
ph_ID_scale_histo_2018 = ROOT.TH2F()
ph_ID_scale_histo_2018 = ph_ID_scale_file_2018.Get("EGamma_SF2D")

ph_pixVeto_scale_name_2018  = "scale_factors/HasPix_SummaryPlot_UL18.root"
ph_pixVeto_scale_file_2018  = ROOT.TFile(ph_pixVeto_scale_name_2018)
ph_pixVeto_scale_histo_2018 = ROOT.TH2F()
ph_pixVeto_scale_histo_2018 = ph_pixVeto_scale_file_2018.Get("MVAID/SF_HasPix_MVAID")

class Simplified_Workflow_Handler:

    def __init__(self,signalname,dataname,isBDT):

        # Where the files are
        self.dir_bkg_input  = "rootfiles/latest_production/MC/backgrounds/"
        self.dir_sig_input  = "rootfiles/latest_production/MC/signals/"
        self.dir_data_input = "rootfiles/latest_production/dataprocess/"

        self.norm_filename  = "rootfiles/latest_production/MC/normalizations/Normalizations_table.txt"

        
        ###################################################################################
        #                                                                                 #
        #------------------------ Add BDT variables to the reader ------------------------#
        #                                                                                 #
        ###################################################################################

        reader.AddVariable("_firstTrkIsoCh",trk1isoCh_array)
        #reader.AddVariable("_coupleIso0",pairIso_array)
        reader.AddVariable("_bestCouplePt/mesonGammaMass",mesonPt_array)
        reader.AddVariable("_photonEt/mesonGammaMass",photonEt_array)
        reader.AddVariable("_bestCoupleEta",mesonEta_array)
        #reader.AddVariable("_photonEt/_HpT",photonEt_array)
        #reader.AddVariable("_bestCouplePt/_HpT",mesonPt_array)
        #reader.AddVariable("_photonEt/_HpT",photonEt_array)
        #reader.AddVariable("_JetNeutralEmEnergy",JetNeutralEmEnergy_array)
        #reader.AddVariable("_JetChargedHadEnergy",JetChargedHadEnergy_array)   
        #reader.AddVariable("_JetNeutralHadEnergy",JetNeutralHadEnergy_array)

        #reader.AddVariable("_metPt",metPt_array)
        #reader.AddVariable("_bestJetPt/mesonGammaMass",bestJetPt_array)
        #reader.AddVariable("_dPhiGammaTrk",dPhiGammaTrk_array)

        if isBDT:
            reader.BookMVA("BDT","MVA/default/weights/TMVAClassification_BDT.weights.xml")

    #Get BDT output function ###########################################################################################################

    def get_BDT_output(self,trk1isoCh,mesonPt,photonEt,mesonEta,mesonGammaMass):#,JetNeutralEmEn,JetChargedHadEn,JetNeutralHadEn):  HiggsPt  pairIso

        trk1isoCh_array[0] = trk1isoCh
        #pairIso_array[0]   = pairIso
        mesonPt_array[0]   = mesonPt/mesonGammaMass
        photonEt_array[0]  = photonEt/mesonGammaMass
        mesonEta_array[0]  = mesonEta
        #JetNeutralEmEnergy_array[0]  = JetNeutralEmEn
        #JetChargedHadEnergy_array[0] = JetChargedHadEn
        #JetNeutralHadEnergy_array[0] = JetNeutralHadEn
        #metPt_array[0]     = metPt
        #bestJetPt_array[0] = bestJetPt/mesonGammaMass
        #dPhiGammaTrk_array[0] = dPhiGammaTrk

        return reader.EvaluateMVA("BDT")

    #Photon scaling function ###########################################################################################################
    def get_photon_scale(self,ph_pt, ph_eta):

        local_ph_pt = ph_pt
        if local_ph_pt > 499.: # This is because corrections are up to 499 GeV
            local_ph_pt = 499.
        
        local_ph_eta = ph_eta
        if local_ph_eta >= 2.5:
            local_ph_eta = 2.49
        if local_ph_eta < -2.5: # Corrections reach down to eta = -2.5, but only up to eta = 2.49
            local_ph_eta = -2.5

        scale_factor_ID = ph_ID_scale_histo_2018.GetBinContent( ph_ID_scale_histo_2018.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2018.GetYaxis().FindBin(local_ph_pt) )
        ph_ID_err       = ph_ID_scale_histo_2018.GetBinError( ph_ID_scale_histo_2018.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2018.GetYaxis().FindBin(local_ph_pt) )

        etaBin = 1 if abs(local_ph_eta) < 1.48 else 4
        scale_factor_pixVeto = ph_pixVeto_scale_histo_2018.GetBinContent(etaBin)
        ph_pixVeto_err       = ph_pixVeto_scale_histo_2018.GetBinError(etaBin)
        
        scale_factor = scale_factor_ID * scale_factor_pixVeto
        tot_err      = math.sqrt( scale_factor_pixVeto * scale_factor_pixVeto * ph_ID_err * ph_ID_err + scale_factor_ID * scale_factor_ID * ph_pixVeto_err * ph_pixVeto_err )

        return scale_factor, tot_err



