import argparse
import random
import re
import numpy as np

# example 1
"""
- fieldName: C3H8
  field: !ablate::mathFunctions::Linear
  startValues: [0.028115959]
"""

# example 2
"""
- fieldName: H2
  field:  2.66E-08
"""

_spread_default=0.3

class IC_Cfg:
    def __init__(self, IC_cfg_path: str, modified_IC_cfg_path: str=None, spread: float=_spread_default):
        with open(IC_cfg_path, 'r') as f:
            self.IC_cfg=f.read()
        self.IC_cfg_path=IC_cfg_path

        # default is random new name
        if not modified_IC_cfg_path:
            id=random.randint(0, int(1e5))
            modified_IC_cfg_path = re.sub(r'(\.yaml)$', r'.new{}\1'.format(id), self.IC_cfg_path)
        self.modified_IC_cfg_path=modified_IC_cfg_path

        assert spread < 1, 'Variable perturb has undesirable properties with spread>=1'
        self.spread=spread

        # These are patterns for capturing fuel ICs of type 1 & 2 as given in above examples
        # In each case there is one format placeholder (i.e. {...}) & 2 groups (i.e. (...)): 
        # * 1st group: captures the biolerplate/signature of the pattern (for restoration in a sub)
        # * 2nd group: actually captures the IC species value itself.
        # * format placeholder: exists where the species name would be, either .format() in the name or a catchall pattern 
        self.template_IC_pattern1_start=r'(?P<prefix>- +fieldName: +{}\s+field: +!ablate::mathFunctions::Linear\s+startValues:) +\[ *(?P<value>[0-9.Ee-]+) *\]'
        self.template_IC_pattern1_end=self.template_IC_pattern1_start.replace('startValues', r'.*\s*endValues')
        self.template_IC_pattern2=r'(?P<prefix>- +fieldName: +{}\s+field:) +(?P<value>[0-9.Ee-]+)' # from deprecated config version

    @staticmethod
    def rand_perm_arg(arg, spread=_spread_default):
        """ randomly permutes *one* arg """
        if type(arg) is int or type(arg) is float:
            return type(arg)(arg ** ((random.random()*2-1) * spread+1))
            # *2 because spread is two sided (i.e. +-)
        if type(arg) is bool:
            return random.random() > 0.5
            #return (float(arg) + 0.5) ** ((random.random() + 0.5) * spread * 2)
        if type(arg) in [list, tuple]:
            return [rand_perm_arg(arg_i) for arg_i in arg]

    @property
    def width(self):
        width_pattern=r'upper: +\[ *([0-9.Ee-]+) *\]'
        match=re.search(width_pattern, self.IC_cfg)
        return float(match.group(1))
    
    @width.setter
    def width(self, new_width_value):
        if not new_width_value:
            new_width_value=self.rand_perm_arg(self.width, spread=self.spread)
            print(f'performing random perturbation of width: old value={self.width}, new_value={new_width_value}')

        width_pattern1=r'upper: +\[ *([0-9.Ee-]+) *\]' # @ line 22 (in sampleDiffusionFlame.yaml)
        width_pattern2=r'end: +([0-9.Ee-]+)' # @ lines 77-112 (in sampleDiffusionFlame.yaml)
        self.IC_cfg=re.sub(width_pattern1, f'upper: [{new_width_value}]', self.IC_cfg)
        self.IC_cfg=re.sub(width_pattern2, f'end: {new_width_value}', self.IC_cfg)

    def get_nonzero_IC_species(self, include_fuel=True, include_oxidizer=True):
        assert include_fuel or include_oxidizer, 'You need to choose at least one part of ICs to search!'
        if include_fuel and include_oxidizer: # we only handle 1 case at a time, but recursion gives use both!
            return self.get_nonzero_IC_species(include_fuel=True, include_oxidizer=False) + \
             self.get_nonzero_IC_species(include_fuel=False, include_oxidizer=True)

        # Whether we are looking at oxidizer or fuel (start or end) depends where we look for the nonzero values
        generic_IC_pattern1 = self.template_IC_pattern1_start if include_fuel else self.template_IC_pattern1_end
        generic_IC_pattern1 = generic_IC_pattern1.format(r'(?P<name>[A-Z0-9]+)') # add new generic pattern group to capture species name
        #print('get_nonzero_fuel_IC_species(), generic_IC_pattern1: ', generic_IC_pattern1)

        IC_search=re.finditer(generic_IC_pattern1, self.IC_cfg) # like findall but returns match objects
        if not IC_search: raise RuntimeError(f'Couldn\'t find any species in config file!!')
        
        # We added a nested group for capturing the species name! (at group_index=2)
        # Also species must be nonzero to be reportable...
        all_species=[match.group('name') for match in IC_search if float(match.group('value'))>0]
        if include_oxidizer: assert 'O2' in all_species, 'oxygen should be in oxidizer!'
        return all_species

    def get_one_species_IC(self, species_name: str, fuel_IC=True):
        pattern1 = self.template_IC_pattern1_start if fuel_IC else self.template_IC_pattern1_end
        IC_search = re.search(pattern1.format(species_name), self.IC_cfg, re.IGNORECASE)
        if IC_search is None: raise RuntimeError(f'Couldn\'t find species named: {species_name} (with fuel_IC={fuel_IC}) in config file !')
        IC_species=float(IC_search.group('value')) # 2nd group is the species value itself
        assert 1 >= IC_species >= 0
        return IC_species

    def modify_one_species_IC(self, species_name: str, new_value: float=None, fuel_IC=True):
        if new_value is None: # default is just to randomly perturb
            new_value=self.get_one_species_IC(species_name, fuel_IC=fuel_IC)
            new_value=self.rand_perm_arg(new_value, self.spread)
        # modify example 1 type
        pattern1=self.template_IC_pattern1_start if fuel_IC else self.template_IC_pattern1_end
        pattern1=pattern1.format(species_name)
        sub1=r'\g<prefix> [{}]'.format(new_value) # first group is the boilerplate signature for restoration

        #print('modify_one_species_fuel_IC(), pattern1:', pattern1, '\nsub1: ', sub1)
        self.IC_cfg=re.sub(pattern1, sub1, self.IC_cfg)

        # # modify example 2 type
        # pattern2=self.template_IC_pattern2.format(species_name)
        # sub2=r'(?P=prefix) {}'.format(str(new_value).upper()) # first group is the boilerplate signature for restoration
        # print('modify_one_species_fuel_IC(), pattern2:', pattern2, '\nsub2: ', sub2)
        # self.IC_cfg=re.sub(pattern2, sub2, self.IC_cfg)

    def constrain_species_L1(self, fuel_IC=True):
        species_names=self.get_nonzero_IC_species(True, True)
        IC_values = np.array([self.get_one_species_IC(name, fuel_IC=fuel_IC) for name in species_names], dtype=np.float64)
        assert np.all(0<=IC_values), 'Yi\'s can\'t be negative!'

        IC_values /= IC_values.sum() # all Yis should sum to 1
        IC_values /= IC_values.sum() # all Yis should sum to 1
        assert IC_values.sum()==1
        for name, value in zip(species_names, list(IC_values)):
            self.modify_one_species_IC(name, value, fuel_IC=fuel_IC)

    def __del__(self):
        # constraint modifications to be physically valid!
        self.constrain_species_L1(True)
        self.constrain_species_L1(False)

        print('dumping perturbed cfg to: ', self.modified_IC_cfg_path)

        # dump modifications to new file
        with open(self.modified_IC_cfg_path, 'w') as f:
            f.write(self.IC_cfg)
    
if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Used to Perturb Ablate ICs from the sampleDiffusionFlame.yaml template.')
    parser.add_argument('IC_cfg_path', type=str, help='The IC file to perturb, format should match that of sampleDiffusionFlame.yaml')
    parser.add_argument('--spread', type=float, default=_spread_default, help='Perturbation spread (if perturbation is enabled).')
    parser.add_argument('--fuel_IC_species_values', type=float, nargs='*', help='The fuel IC species values to insert, should match printed order (also found in IC file).')
    parser.add_argument('--oxidizer_IC_species_values', type=float, nargs='*', help='The fuel IC species values to insert, should match printed order (also found in IC file).')
    parser.add_argument('--width', type=float, default=None, help='Desired width of the simulation.')
    args=parser.parse_args()

    ic_cfg=IC_Cfg(args.IC_cfg_path)
    ic_cfg.width=args.width # if we assign None (b/c width is unspecified), then it will make random perturbation of width

    all_species = ic_cfg.get_nonzero_IC_species()
    def set_IC_vec(new_IC_vals: list, fuel_IC=True):
        existing_ICs = ic_cfg.get_nonzero_IC_species(include_fuel=fuel_IC, include_oxidizer=not fuel_IC)
        existing_ICs = {name: ic_cfg.get_one_species_IC(name, fuel_IC=fuel_IC) for name in existing_ICs}
        print(f'Found IC species (with fuel_IC={fuel_IC}): ', existing_ICs)
        print(f'With sum={sum(existing_ICs.values())}')

        if not new_IC_vals: # if no species values are given, then default is to perturb
            new_IC_vals = [None]*len(all_species)
        for species_name, species_value in zip(all_species, new_IC_vals):
            ic_cfg.modify_one_species_IC(species_name, species_value, fuel_IC=fuel_IC)
    set_IC_vec(args.fuel_IC_species_values, fuel_IC=True)
    set_IC_vec(args.oxidizer_IC_species_values, fuel_IC=False)
    del ic_cfg
