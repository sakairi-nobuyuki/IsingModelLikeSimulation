import config_2d
import itertools
import pprint

class SheetLikeIsingModelConfig:
    def __init__(self):
        self.parameter = {}
        self.parameter['n_x'] = config_2d.n_x
        self.parameter['n_y'] = config_2d.n_y
        self.parameter['n_spin_state']  = config_2d.n_spin_state
        self.parameter['IonicStrength'] = config_2d.IonicStrength
        self.parameter['Phi']   = config_2d.Phi
        self.parameter['Delta'] = config_2d.Delta
        self.parameter['L']     = config_2d.L
        self.parameter['Cexcl'] = config_2d.Cexcl
        self.parameter['kB']    = config_2d.kB
        self.parameter['T']     = config_2d.T
        self.parameter['TsaMax'] = config_2d.TsaMax
        self.parameter['TsaMin'] = config_2d.TsaMin
        self.parameter['dTsa']   = config_2d.dTsa
        self.parameter['ksa']    = config_2d.ksa
        self.parameter['PsaRef'] = config_2d.PsaRef
        self.parameter['n_trial']  = config_2d.n_trial
        self.parameter['CvdW']     = config_2d.CvdW
        self.parameter['Alpha']    = config_2d.Alpha
        self.parameter['output_file_name'] = config_2d.output_file_name
        self.parameter['output_freq']      = config_2d.output_freq

        
        self.parameter_combinations = list (itertools.product (self.parameter['n_x'], self.parameter['n_y'], \
            self.parameter['n_spin_state'], self.parameter['IonicStrength'], self.parameter['Phi'], \
            self.parameter['Delta'], self.parameter['L'], self.parameter['Cexcl'], self.parameter['kB'], \
            self.parameter['T'], self.parameter['TsaMax'], self.parameter['TsaMin'], self.parameter['dTsa'], \
            self.parameter['ksa'], self.parameter['PsaRef'], self.parameter['n_trial'], self.parameter['CvdW'], self.parameter['Alpha']))

        pprint.pprint (self.parameter_combinations)

        self.parameter_dict_list = []
        for item_list in self.parameter_combinations:
            parameter_dict = {}
            for i_item, item_key in enumerate (self.parameter.keys ()):
                if 'output' in item_key:  continue
                parameter_dict[item_key] = item_list[i_item]
            self.parameter_dict_list.append (parameter_dict)
        pprint.pprint (self.parameter_dict_list)