from chrestData import *
import pandas as pd
import os
from ablateData import AblateData

class Array2DF_Builder:
    def __init__(self):
        self.dfs = []
        self._df = None

    def add_df_array(self, df_array, columns):
        """assumes same shape as other arrays & that last dimensions is the "columns" dimension """
        # order='C' means last index changes fastest for reading/writing of arrays: 
        # https://numpy.org/doc/stable/reference/generated/numpy.reshape.html#numpy.reshape
        # this is exactly what we want (since last dimension is the components dimension)
        df_array = df_array.reshape(-1, df_array.shape[-1], order='C')
        self.dfs.append(pd.DataFrame(df_array, columns=columns))
        self._df = None
    
    def __str__(self):
        rep = ''
        for i, df in enumerate(self.dfs):
            rep += f'df #{i}: \n' + str(df.describe()) + '\n\n'
        return rep

    def add_chrest_fields(self, chrest_data, fields: list): # TODO: remove chrest_data arg!!
        for field in fields:
            array, times, component_names = chrest_data.get_field(field) # this is literally the only time we use chrest_data!
            assert (component_names is None) == (len(array.shape)==4) # sanity check in case interface changes...
            if component_names is None: # for 'single component cases'
                array = np.expand_dims(array, -1)
                component_names = [field]
            else:
                component_names = [field + col for col in component_names]
            self.add_df_array(array, columns=component_names)

    @property
    def df(self):
        if self._df is None:
            self._df = pd.concat(self.dfs, axis=1)
        return self._df

def plot_cell_groups(df):
    """ 
    Verify cell groups is working (colors are randomized for contrast).
    Cool interactive visualization!!
    """
    codes = pd.Categorical(df['group']).codes
    print(np.max(codes))
    new_codes = np.arange(np.max(codes)+1)
    np.random.shuffle(new_codes)

    df['group'] = new_codes[codes] #pd.Categorical(df['group']).codes # turn group codes into regular indices
    import plotly.express as px

    fig = px.scatter_3d(df, x='x', y='y', z='z', color='group')
    fig.show()

def make_groups(df, n=100):
    """ 
    makes super-cell groups
    :param df: df with x,y,z coordinates (everything else optional)
    :param n: n super-cells along each axis, i.e. super-cell array is nxnxn
    """
    dims = ['x','y','z']
    df[dims] -= df[dims].min() # remove negative values!
    max_ = df[dims].max().max() # complete coordinate grid must be a cube, so we use only 1 max value
    print(max_)

    #max_signed_32_int_val=2,147,483,647
    # 10^18 is as high as it can get for signed 32bit int without overflow!!
    decimal_power=9//(len(dims)-1) # <-- this is the "digit padding" used to keep the value of each coordinate from overlapping
    assert (len(dims)-1)*decimal_power<9 # -1 b/c arange goes to len(dims)-1
    assert n<=10**decimal_power, f"You're chosen n={n} would cause signed 32bit int overflow!!"
    # this must be true for hasing to be "unique"
    
    # idea is create hash by using different decimal places & sum
    powers = 10**(np.arange(len(dims), dtype='int64')*decimal_power)
    print(powers)

    # divide by max & multiply by n
    # NOTE: the new coordinates represent the coordinates
    # of the hypothetical super-cells (i.e. in a 3d array)
    new_coords = ((df[dims]/max_)*n-1e-8).astype('int64')
    groups = (powers*new_coords).sum(axis=1)
    print(new_coords.describe())
    print(groups.head())

    # since it is 3d and we want n across each dimension that makes n^3 super-cubes
    # (or less since technically not ever "super-group" will have members)
    assert len(np.unique(groups))<=n**len(dims)
    return groups

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='wrangle chrest data file by adding UQ groups (for cell coursening)')
    parser.add_argument('--file', dest='file', type=str, required=True,
                        help='The path to the ablate hdf5 file containing the ablate data.')
    parser.add_argument('--fields', dest='fields', type=str,
                        help='The list of fields to map from ablate to chrest in format  --field '
                             'ablate_name:chrest_name e.g. --field aux_temperature:temperature '
                             'aux_velocity:vel', nargs='+', default=['souener', 'zmix', 'Yi', 'souspec'])
    parser.add_argument('--n-cubes-per-dim', type=int, default=100, 
                        help='number of cubes touching each axis of super-cell grid array')
    parser.add_argument('--plot', action='store_true', help='whether to plot the cell groups')
    args = parser.parse_args()
    print('unprocessed fields: ', args.fields)
    args.fields = [field.split(':')[-1] for field in args.fields] 
    print('post-processed fields: ', args.fields)
    # NOTE: our script passes in fields in this format: aux_temperature:temperature,
    # which we are supporting for conveince however we already applied the alias in 
    # ablateData.py, so now we'd just extract the alias part & use that
    

    # this is some example code for chest file post processing
    chrest_data = ChrestData(args.file) # TODO: remove ChrestData dependence!
    df_builder = Array2DF_Builder()

    df_builder.add_chrest_fields(chrest_data, args.fields) # TODO: remove ChrestData dependence!
    print(df_builder)

    # This requires manual code because it is not a "field"
    coords = chrest_data.get_coordinates()
    print(coords.shape)
    df_builder.add_df_array(coords, columns=['x','y','z'])
    df = df_builder.df

    groups = make_groups(df, n = args.n_cubes_per_dim)
    df['group'] = groups
    if args.plot: plot_cell_groups(df)
    df.to_csv(f"{os.path.dirname(args.file)}/chrest_data.csv.gz", index=False)
