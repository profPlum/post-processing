import os
import pandas as pd
import numpy as np
import argparse

class Array2DF_Builder:
    def __init__(self):
        self._dfs = [] # the list of new dfs which are pending merge into master_df
        self._master_df = None # the currently build master_df (has x,y,z as indices by default)

    def add_df_array(self, df_array, columns, coords):
        """assumes same shape as other arrays & that last dimensions is the "columns" dimension """
        # order='C' means last index changes fastest for reading/writing of arrays: 
        # https://numpy.org/doc/stable/reference/generated/numpy.reshape.html#numpy.reshape
        # this is exactly what we want (since last dimension is the components dimension)
        df_array = np.asarray(df_array).reshape(-1, df_array.shape[-1], order='C').squeeze()
        df = pd.DataFrame(df_array, columns=columns)
        assert len(df)==len(coords)
        df = pd.concat([df,coords],axis=1)
        df=df.set_index(list(coords.columns)) # adds coords (df) as the multi-index for df
        self._dfs.append(df)
    
    def __str__(self):
        rep = ''
        for i, df in enumerate(self._dfs):
            rep += f'df #{i}: \n' + str(df.describe()) + '\n\n'
        return rep
    
    def compute_ablate_coordinates(self, ablate_data, dim_names = ['x','y','z']):
        """" Get coords from ablate_data (as df) & rounds to avoid machine error """
        # Coords requires manual code because it is not a "field"
        # Also ensure that the ablate_data objects are compatible!
        coords = ablate_data.compute_cell_centers(dimensions=len(dim_names))
        coords = coords.round(decimals=7)
        # round above machine epsilon so that numbers can be trusted for join!
        
        print('coords.shape: ', coords.shape)
        return pd.DataFrame(coords, columns=dim_names)
    
    def add_ablate_fields(self, ablate_data, fields: list):
        """ for adding desired fields from (possibly multiple) ablate_data objects """
        try:
            ablate_data.get_field(ablate_data.get_fields(), 1)
            raise RuntimeError("Interface has changed! Now Interwal is relevant parameter!=0")
        except IndexError: None
        
        coords = self.compute_ablate_coordinates(ablate_data)
        
        for field in fields.copy():
            # split field/aliases pairs or simply default alias=field
            ablate_field, alias = (field.split(':') + [field])[:2]
            try: 
                array, times, component_names = ablate_data.get_field(ablate_field, 0)
            except Exception as e:
                print('error! skipping: ', e)
                continue
            assert (component_names is None) == (len(array.shape)==2)
            if component_names is None: # for 'single component cases'
                component_names = [alias]
            else:
                component_names = [alias + col for col in component_names]
            self.add_df_array(array, columns=component_names, coords=coords)
            fields.remove(field)
    
    @property
    def df(self):
        if self._master_df is None: self._master_df = self._dfs.pop(0)
        if len(self._dfs): # check if we need to add new dataframes
            self._master_df = self._master_df.join(self._dfs) # does list-join across all the recorded dataframes
            self._dfs.clear()
        
        # puts x,y,z index back into columns
        return self._master_df.reset_index()

def plot_cell_groups(df):
    """ 
    Verify cell groups is working (colors are randomized for contrast).
    Cool interactive visualization!!
    """
    codes = pd.Categorical(df['group']).codes
    print(np.max(codes))
    new_codes = np.arange(np.max(codes)+1)
    np.random.shuffle(new_codes)

    df['group'] = new_codes[codes]
    import plotly.express as px

    fig = px.scatter_3d(df, x='x', y='y', z='z', color='group')
    fig.show()

def make_groups(df, n=100, dims = ['x','y','z']):
    """ 
    makes super-cell groups
    :param df: df with x,y,z coordinates (everything else optional)
    :param n: n super-cells along each axis, i.e. super-cell array is nxnxn
    :param dims: dimension names, should essentially be constant unless you want to use space-time
    """
    new_coords = df[dims] - df[dims].min() # remove negative values!
    max_ = new_coords.max().max() # complete coordinate grid must be a cube, so we use only 1 max value
    print('max_:', max_)

    #max_signed_32_int_val=2,147,483,647
    # 10^18 is as high as it can get for signed 32bit int without overflow!!
    decimal_power=9//(len(dims)-1) # <-- this is the "digit padding" used to keep the value of each coordinate from overlapping
    assert (len(dims)-1)*decimal_power<9 # -1 b/c arange goes to len(dims)-1
    assert n<=10**decimal_power, f"You're chosen n={n} would cause signed 32bit int overflow!!"
    # this must be true for hasing to be "unique"
    
    # idea is create hash by using different decimal places & sum
    powers = 10**(np.arange(len(dims), dtype='int64')*decimal_power)
    print('powers: ', powers)

    # divide by max & multiply by n
    # NOTE: the new coordinates represent the coordinates
    # of the hypothetical super-cells (i.e. in a 3d array)
    new_coords = ((new_coords/max_)*n-1e-8).astype('int64')
    groups = (powers*new_coords).sum(axis=1)
    print('new_coords: ')
    print(new_coords.describe())
    print('groups: ')
    print(groups.head())
    print(groups.describe())

    # since it is 3d and we want n across each dimension that makes n^3 super-cubes
    # (or less since technically not ever "super-group" will have members)
    assert len(np.unique(groups))<=n**len(dims)
    return groups

# Monkey Patch to Fix Matt's code locally
def AblateData_factory(files):
    from ablateData import AblateData
    if isinstance(files, str): files=[files]
    for fn in args.files:
        assert os.path.exists(fn), f"The provided Ablate HDF5 file path doesn't exist: {fn}"
    ablate_data = AblateData(files)
    if len(ablate_data.vertices)==1:
        ablate_data.vertices=ablate_data.vertices.reshape(-1,1)
    return ablate_data

# verified to work 7/23/23
def rebalance_flame(df, replace=False):
    """ 
    rebalances flame so that at least 50% of the data is flame! 
    :param replace: if true will do replacement sampling to ensure the balance is exactly 50% (if necessary will create duplicates of air)
    """
    is_flame_predicate = lambda df: np.logical_and(df['temp']>=1900, df['YiO2']<0.8)
 
    df_groups = df.groupby(['group']).median()
    df_groups['is_flame'] = is_flame_predicate(df_groups)
    df_groups= df_groups[['is_flame']]
    n_flame = df_groups['is_flame'].sum()
    #if n_flame < len(df_groups)//2:
    try: # down-sample air to match n_flame!
        df_groups = df_groups.groupby(['is_flame']).sample(n=n_flame, replace=replace).reset_index().set_index('group')
    except: None # will throw error when n_flame > n_air & replace=False

    df = df.set_index('group').join(df_groups, how='inner').reset_index()

    assert df['is_flame'].mean()>=0.5, "Rebalancing Failed!"
    if replace: assert df['is_flame'].mean()==0.5, "Rebalancing Failed!"
    return df


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='wrangle ablate data file by converting to csv.gz, collating multiple HDF5 files,'
                                    'adding UQ groups (for cell coursening), possibly rebalance based on flame.')
    parser.add_argument('--files', dest='files', type=str, required=True, nargs='+',
                        help='The path to the ablate hdf5 files containing the ablate data.'
                        'files beyond the first are only used when field cannot be found in first file')
    parser.add_argument('--fields', dest='fields', type=str,
                        help='The list of fields to map from ablate to alias in format  --field '
                             'ablate_name:alias_name e.g. --field aux_temperature:T '
                             'aux_velocity:vel', nargs='+')
    parser.add_argument('--n-cubes-per-dim', type=int, default=100, 
                        help='number of cubes touching each axis of super-cell grid array')
    parser.add_argument('--plot', action='store_true', help='whether to plot the cell groups')
    parser.add_argument('--rebalance-flame', action='store_true', help='Whether to rebalance the data flame so that 50% is flame')
    args = parser.parse_args()
    #for fn in args.files:
    #    assert os.path.exists(fn), f"The provided Ablate HDF5 file path doesn't exist: {fn}"
    print('fields: ', args.fields)
    
    ablate_data = [AblateData_factory(fn) for fn in args.files]

    #n_dims = ablate_data.vertices.shape[1]
    
    df_builder = Array2DF_Builder()
    for ablate_data_obj in ablate_data:
        df_builder.add_ablate_fields(ablate_data_obj, args.fields)
    if len(args.fields)>0:
        raise RuntimeError('Missing requested fields from Ablate HDF5 file(s): ' + ', '.join(args.fields))
    print(df_builder)
    df = df_builder.df

    dim_names=['x','y','z']
    groups = make_groups(df, n = args.n_cubes_per_dim, dims=dim_names)
    df_builder.add_df_array(groups, columns=['group'], coords=df[dim_names])
    df = df_builder.df # optional, worth it?
    #df['group'] = groups

    rebalanced_str=''
    if args.rebalance_flame:
        df = rebalance_flame(df)
        rebalanced_str = '-rebalanced'

    if args.plot: plot_cell_groups(df)
        
    df.to_csv(f"{os.path.splitext(args.files[0])[0]}{rebalanced_str}.csv.gz", index=False)

