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

trk1iso_array  = array('f', [0.])
mesoinPt_array = array('f', [0.])
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

        reader.AddVariable("_firstTrkIso",trk1iso_array)
        reader.AddVariable("_bestCouplePt/mesonGammaMass",mesoinPt_array)
        reader.AddVariable("_photonEt/mesonGammaMass",photonEt_array)

        if isBDT:
            reader.BookMVA("BDT","MVA/default/weights/TMVAClassification_BDT.weights.xml")

    ###############################################################################################################################################

    def get_BDT_output(self,trk1iso,mesonPt,photonEt,mesonGammaMass):

        trk1iso_array[0]   = trk1iso
        mesoinPt_array[0]  = mesonPt/mesonGammaMass
        photonEt_array[0]  = photonEt/mesonGammaMass

        return reader.EvaluateMVA("BDT")
    
