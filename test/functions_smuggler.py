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

K1iso_array     = array('f', [0.])
KKpT_array    = array('f', [0.])
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

        reader.AddVariable("_firstCandIso",K1iso_array)
        reader.AddVariable("_bestCouplePt/mass_KKg",KKpT_array)
        reader.AddVariable("_photonEt/mass_KKg",photonEt_array)

        if isBDT:
            reader.BookMVA("BDT","MVA/default/weights/TMVAClassification_BDT.weights.xml")

    ###############################################################################################################################################

    def get_BDT_output(self,K1iso,KKpT,photonEt,mass_KKg):

        K1iso_array[0]     = K1iso
        KKpT_array[0]      = KKpT/mass_KKg
        photonEt_array[0]  = photonEt/mass_KKg

        return reader.EvaluateMVA("BDT")
    
