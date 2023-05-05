import pandas as pd
import numpy as np
import os

fn = '/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/flowField_mixtureFraction.00001.chrest/chrest_data.csv'
#fn = f'{os.environ["HOME"]}/Downloads/chrest.csv'
df = pd.read_csv(fn)

import plotly.express as px

#fig = px.scatter_3d(df, x='x', y='y', z='z', color='group')
#fig.show()

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

    fig = px.scatter_3d(df, x='x', y='y', z='z', color='group', title='cell_groups_for_coursening')
    fig.show()

plot_cell_groups(df)