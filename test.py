import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
import pylab
import struct
import matplotlib.cm as cm
import conversions as co
import pandas as pd
import sqlite3

def read_hdf5(path, p_type=1):
    """types: 1=DM, 2=disk, 3=bulge"""
    groups = ["Header", "PartType1", "PartType2", "PartType3"]
    head = {}
    
    f = h5py.File(path, "r")
    h = f["Header"]
    keys = list(h.attrs)
    for key in keys:
        head[key] = h.attrs[key]
    parts = f[groups[p_type]]
    columns = ["x","y","z", "vx", "vy", "vz"]
    df = pd.DataFrame(
        np.concatenate((parts["Coordinates"], parts["Velocities"]), axis=1), 
            columns=columns, index=parts["ParticleIDs"][:])
    df.index.name = "particleID"
    f.close()
    
    return head, df

def df_center(df):
    """Centers a data frame's x,y,z keys by subtracting the median."""
    idxs = ['x', 'y', 'z']
    for idx in idxs:
        df[idx] -= df[idx].median()
        
def df_polar(df):
    """Adds in r and phi coordinates, as well as their velocities.
    Phi is the physics spherical phi, i.e. the polar theta.
    """
    df['r'] = np.sqrt(df['x']**2 + df['y']**2)
    df['phi'] = np.arctan2(df['y'], df['x'])
    df['vr'] = (df['x']*df['vx'] + df['y']*df['vy']) / df['r']
    df['vphi'] = (df['x']*df['vy'] - df['y']*df['vx']) / (df['r']**2)


# In[5]:


def df_filter(df, key, low=None, high=None, f=None):
    """Filters a dataframe by key, greater than low, less than high,
    optionally applying function f to the data comparison.
    """
    if low is None and high is None:
        print("No filtering done")
        return df
    elif low is not None and high is not None:
        if f is None:
            return df[(df[key] > low) & (df[key] < high)]
        else:
            return df[(f(df[key]) > low) & (f(df[key]) < high)]
    elif low is not None:
        if f is None:
            return df[df[key] > low]
        else:
            return df[f(df[key]) > low]
    elif high is not None:
        if f is None:
            return df[df[key] < high]
        else:
            return df[f(df[key]) < high]
    else:
        print("Nani?")
        return df


path = "data/snap_582.hdf5"
head, df = read_hdf5(path, p_type=1)
df_center(df)
df_polar(df)


conn = sqlite3.connect("data/particles.db")
df.to_sql("darkmatter", conn, chunksize=1000, if_exists='replace')
conn.close()
