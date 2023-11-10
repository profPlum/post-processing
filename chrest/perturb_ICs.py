import argparse
import random
import re, os
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
        self.templte_IC_pattern1=r'(- +fieldName: +{}\s+field: +!ablate::mathFunctions::Linear\s+.*?startValues:) +\[ *([0-9.E-]+) *\]'
        self.template_IC_pattern2=r'(- +fieldName: +{}\s+field:) +([0-9.Ee-]+)'

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
        width_pattern=r'upper: +\[ *([0-9.E-]+) *\]'
        match=re.search(width_pattern, self.IC_cfg)
        return float(match.group(1))
    
    @width.setter
    def width(self, new_width_value):
        if not new_width_value:
            new_width_value=self.rand_perm_arg(self.width, spread=self.spread)
            print('performing random permutation of width: old value={self.width}, new_value={new_width_value}')

        width_pattern1=r'upper: +\[ *([0-9.E-]+) *\]' # @ line 22 (in sampleDiffusionFlame.yaml)
        width_pattern2=r'end: +([0-9.E-]+)' # @ lines 77-112 (in sampleDiffusionFlame.yaml)
        self.IC_cfg=re.sub(width_pattern1, f'upper: [{new_width_value}]', self.IC_cfg)
        self.IC_cfg=re.sub(width_pattern2, f'end: {new_width_value}', self.IC_cfg)

    def get_nonzero_IC_species(self, include_fuel=True, include_oxidizer=True):
        assert include_fuel or include_oxidizer, 'You need to choose at least one part of ICs to search!'

        # NOTE: only the type two pattern is easily seperable into oxidizier & fuel
        # add new generic pattern group to capture species name
        generic_IC_pattern2 = self.template_IC_pattern2.format(r'([A-Z0-9]+)')
        print('get_nonzero_fuel_IC_species(), generic_IC_pattern2: ', generic_IC_pattern2)

        # Here we are matching repeated patterns of generic_IC_pattern2, starting at 'massFractionsFuel:'
        # This search pattern will be interrupted once massFractionsFuel sequence is over (e.g. by oxidizer)
        search_region_template_pattern = r'\s+{}: +!ablate::finiteVolume::fieldFunctions::MassFractions\s+.+\s+values:(?:\s+{})+'
        search_region_template_pattern=search_region_template_pattern.format('{}', re.sub('[()]', '', generic_IC_pattern2))
        # NOTE: we reinsert '{}' b/c that section is dependant on fuel or oxidizer search (hence it's a template)
        
        print('search_region_template_pattern: ', search_region_template_pattern)

        search_region=''
        if include_fuel:
            search_region_pattern=search_region_template_pattern.format('massFractionsFuel')
            print('Searching for fuel with pattern: ', search_region_pattern)
            search_region+=re.search(search_region_pattern, self.IC_cfg).group(0) # (group=0 is whole match)
        if include_oxidizer:
            search_region_pattern=search_region_template_pattern.format('massFractionsOxidizer')
            print('Searching for oxidizer with pattern: ', search_region_pattern)
            search_region+=re.search(search_region_pattern, self.IC_cfg, re.DOTALL).group(0)

        print('Extracted Search Region: ', search_region)

        IC_search=re.finditer(generic_IC_pattern2, search_region) # like findall but returns match objects
        if not IC_search:
            raise RuntimeError(f'Couldn\'t find any species in config file!!')
        
        # We added a nested group for capturing the species name! (at group_index=2)
        all_species=[match.group(2) for match in IC_search]
        return all_species

    def get_one_species_IC(self, species_name: str):
        # search example 2 type:
        IC_search = re.search(self.template_IC_pattern2.format(species_name), self.IC_cfg, re.IGNORECASE)
        if IC_search is None:
            raise RuntimeError(f'Couldn\'t find species named: {species_name} in config file!')
        IC_species=float(IC_search.group(2)) # 2nd group is the species value itself
        assert 1 >= IC_species >= 0
        return IC_species

    def modify_one_species_IC(self, species_name: str, new_value: float=None):
        if new_value is None: # default is just to randomly perturb
            new_value=self.get_one_species_IC(species_name)
            new_value=self.rand_perm_arg(new_value, self.spread)
        # modify example 1 type
        pattern1=self.templte_IC_pattern1.format(species_name)
        sub1=r'\1 [{}]'.format(new_value).upper() # first group is the boilerplate signature for restoration

        # if species is part of fuel
        if species_name in self.get_nonzero_IC_species(include_fuel=True, include_oxidizer=False):
            print('modify_one_species_fuel_IC(), pattern1:', pattern1, '\nsub1: ', sub1)
            self.IC_cfg=re.sub(pattern1, sub1, self.IC_cfg)
        # if species is part of oxidizer
        if species_name in self.get_nonzero_IC_species(include_fuel=False, include_oxidizer=True):
            pattern1=pattern1.replace('startValues', 'endValues')
            print('modify_one_species_fuel_IC(), pattern1:', pattern1, '\nsub1: ', sub1)
            self.IC_cfg=re.sub(pattern1, sub1, self.IC_cfg)
        
        # modify example 2 type
        pattern2=self.template_IC_pattern2.format(species_name)
        sub2=r'\1 {}'.format(new_value).upper() # first group is the boilerplate signature for restoration
        print('modify_one_species_fuel_IC(), pattern2:', pattern2, '\nsub2: ', sub2)
        self.IC_cfg=re.sub(pattern2, sub2, self.IC_cfg)

    def constrain_species_L1(self, species_names):
        IC_values = np.array([self.get_one_species_IC(name) for name in species_names])
        assert np.all(0<=IC_values), 'Yi\'s can\'t be negative!'

        IC_values /= IC_values.sum() # all Yis should sum to 1
        for name, value in zip(species_names, list(IC_values)):
            self.modify_one_species_IC(name, value)

    def __del__(self):
        # constraint modifications to be physically valid!
        fuel_species=self.get_nonzero_IC_species(True, False)
        oxidizer_species=self.get_nonzero_IC_species(False, True)
        self.constrain_species_L1(fuel_species)
        self.constrain_species_L1(oxidizer_species)

        # dump modifications to new file
        with open(self.modified_IC_cfg_path, 'w') as f:
            f.write(self.IC_cfg)

        print('dumping to {self.modified_IC_cfg_path} & perturbed.new.yaml')
        os.system(f'rm $(dirname {self.modified_IC_cfg_path})/perturbed.new.yaml 2> /dev/null')
        os.system(f'ln -s {self.modified_IC_cfg_path} $(dirname {self.modified_IC_cfg_path})/perturbed.new.yaml')
    
if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Used to Perturb Ablate ICs from the sampleDiffusionFlame.yaml template.')
    parser.add_argument('IC_cfg_path', type=str, help='The IC file to perturb, format should match that of sampleDiffusionFlame.yaml')
    parser.add_argument('--spread', type=float, default=_spread_default, help='Perturbation spread (if perturbation is enabled).')
    parser.add_argument('--IC_species_values', type=float, nargs='*', help='The IC species values to insert, should match printed order (also found in IC file).')
    parser.add_argument('--width', type=float, default=None, help='Desired width of the simulation.')
    #parser.add_argument('--print-fuel-IC', action='store_true', help='if specified it will just print the fuel ICs & exit.')
    args=parser.parse_args()

    ic_cfg=IC_Cfg(args.IC_cfg_path)
    ic_cfg.width=args.width # if we assign None (b/c width is unspecified), then it will make random perturbation of width

    species = ic_cfg.get_nonzero_IC_species()
    print('Found IC species: ', {name: ic_cfg.get_one_species_IC(name) for name in species})

    if not args.IC_species_values: # if no species values are given, then default is to perturb
        args.IC_species_values = [None]*len(species)
    
    for species_name, species_value in zip(species, args.IC_species_values):
        ic_cfg.modify_one_species_IC(species_name, species_value)
    del ic_cfg
