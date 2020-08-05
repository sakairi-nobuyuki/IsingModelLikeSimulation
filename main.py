import os
import glob
import shutil
import load_config_2d
import config_2d
import output_utils

import pprint

if __name__ == '__main__':
    ###  Set execution environment if Linux or windows


    ###  Load config file
    slimconf = load_config_2d.SheetLikeIsingModelConfig ()
    out      = output_utils.OutputUtils (len (slimconf.parameter_dict_list))


    #pprint.pprint (slimconf.parameter_dict_list)
    #for parameter_dict in slimconf.parameter_dict_list:
        #pprint.pprint (parameter_dict)
        ###slimconf.parameter_dict_list  <-- to be used
        #for param_key, param_value in parameter_dict.items ():
        #    print (param_key, param_value)



    ###  Make iterator according to the config file and iterate
    






    
