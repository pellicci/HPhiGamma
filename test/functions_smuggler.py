#######################################################################
#                                                                     #
# <<  Hi! I smuggle bad coded functions. But, ehy, they're free! >>   #
#                                                                     #
#######################################################################                                                                                                                                                       

import ROOT
import math
import os
from array import array

#------- Arrays and reader for the BDT -------#

K1pT_array     = array('f', [0.])
KKIso_array    = array('f', [0.])
photonEt_array = array('f', [0.])

reader = ROOT.TMVA.Reader("!Color")

class Simplified_Workflow_Handler:

    def __init__(self,signalname,dataname,isBDT):

        # Where the files are
        self.dir_bkg_input  = "rootfiles/latest_production/MC/backgrounds/"
        self.dir_sig_input  = "rootfiles/latest_production/MC/signals/"
        self.dir_data_input = "rootfiles/latest_production/dataprocess/"

        self.norm_filename = "rootfiles/latest_production/MC/normalizations/Normalizations_table.txt"
        
        ###################################################################################
        #                                                                                 #
        #------------------------ Add BDT variables to the reader ------------------------#
        #                                                                                 #
        ###################################################################################

        reader.AddVariable("_firstCandPt",K1pT_array)
        reader.AddVariable("_coupleIso",KKIso_array)
        reader.AddVariable("_photonEt",photonEt_array)

        if isBDT:
            reader.BookMVA("BDT","MVA/default/weights/TMVAClassification_BDT.weights.xml")

    ###############################################################################################################################################

    def get_BDT_output(self,K1pT,KKIso,photonEt):

        K1pT_array[0]     = K1pT
        KKIso_array[0]    = KKIso
        photonEt_array[0] = photonEt

        return reader.EvaluateMVA("BDT")
    
