from chrestData import *
import pandas as pd
import os

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

    def add_chrest_fields(self, chrest_data, fields: list):
        for field in fields:
            array, times, component_names = chrest_data.get_field(field)
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
    xyz = ['x','y','z']
    df[xyz] -= df[xyz].min() # remove negative values!
    m = df[xyz].max().max() # must be a cube, so we use only 1 max value

    decimal_power=4
    assert n<10**decimal_power # this must be true for hasing to be "unique"
    print(m)
    
    # idea is create hash by using different decimal places & sum
    powers = 10**(np.arange(3, dtype='int64')*decimal_power)
    print(powers)

    # divide by max
    new_coords = ((df[xyz]/m)*n-1e-15).astype('int64')
    print(new_coords.describe())
    groups = (powers*new_coords).sum(axis=1)

    # since it is 3d and we want n across each dimension that makes n^3 super-cubes
    #assert len(np.unique(groups))==n**3
    return groups

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='wrangle chrest data file by adding UQ groups (for cell coursening)')
    parser.add_argument('--file', dest='file', type=str, required=True,
                        help='The path to the ablate hdf5 file containing the ablate data.')
    parser.add_argument('--fields', dest='fields', type=str,
                        help='The list of fields to map from ablate to chrest in format  --field '
                             'ablate_name:chrest_name e.g. --field aux_temperature:temperature '
                             'aux_velocity:vel', nargs='+', default=['souener', 'zmix', 'Yi', 'souspec'])
    args = parser.parse_args()
    print('unprocessed fields: ', args.fields)
    args.fields = [field.split(':')[-1] for field in args.fields] 
    print('post-processed fields: ', args.fields)
    # NOTE: our script passes in fields in this format: aux_temperature:temperature,
    # which we are supporting for conveince however we already applied the alias in 
    # ablateData.py, so now we'd just extract the alias part & use that
    

    # this is some example code for chest file post processing
    chrest_data = ChrestData(args.file)
    df_builder = Array2DF_Builder()

    df_builder.add_chrest_fields(chrest_data, args.fields)
    print(df_builder)

    # This requires manual code because it is not a "field"
    coords = chrest_data.get_coordinates()
    print(coords.shape)
    df_builder.add_df_array(coords, columns=['x','y','z'])
    df = df_builder.df

    groups = make_groups(df, n = 100)
    df['group'] = groups
    df.to_csv(f"{os.path.dirname(args.file)}/chrest_data.csv.gz", index=False)
