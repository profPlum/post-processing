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
    hdf5_file = sys.argv[1] # first arg is the file
    #hdf5_file = '/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/flowField_mixtureFraction.00000.chrest/flowField_mixtureFraction.00000.chrest.00000.hdf5'

    # this is some example code for chest file post processing
    chrest_data = ChrestData(hdf5_file)
    df_builder = Array2DF_Builder()

    fields_names = ['souener', 'zmix', 'Yi', 'souspec']
    df_builder.add_chrest_fields(chrest_data, fields_names)
    print(df_builder)

    # This requires manual code because it is not a "field"
    coords = chrest_data.get_coordinates()
    print(coords.shape)
    df_builder.add_df_array(coords, columns=['x','y','z'])
    df = df_builder.df

    groups = make_groups(df, n = 100)
    df['group'] = groups
    df.to_csv(f"{os.path.dirname(hdf5_file)}/chrest_data.csv", index=False)