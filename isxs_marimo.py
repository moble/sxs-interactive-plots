import math
import numpy as np
import csv
import plotly.graph_objects as go
import marimo as mo
from IPython.display import display
import plotly.io as pio
pio.renderers.default = 'iframe'

G = 6.67430e-11
c = 3e8
M = 6.563e+31
r = 3.086e24

def load_data():
    """
    loads the data of the included strains
    returns array of lm modes, array of h ids, nested arrays of htilde and frequencies, and array of metadatas
    """
    file = np.load("marimo_data.npz", allow_pickle=True)
    data = file['arr_0']
    hlm = data[0]
    h_id_list = []
    strain_data = []
    metadata_list = []
    for strain in data[1:]:
        h_id_list.append(strain[0])
        strain_data.append(strain[1:3])
        metadata_list.append(strain[3])
    
    return hlm, h_id_list, strain_data, metadata_list

def load_index(h_id_list , h_id):
    """
    returns the id of corresponding strain
    """
    return find_index(h_id_list, h_id)


def load_plots(strain_data, h_id_index):
    """
    returns the htilde strain and frequency arrays of corresponding h_id
    """
    strain = strain_data[h_id_index]
    frequencies = strain[0]
    htilde = strain[1]
    return frequencies, htilde
    
def find_index(data, value):
    try:
        idx = data.index(value)
    except ValueError:
        print(f"{value} is not a supported strain")
    return idx


def iplt_lm(xStrain, yStrain, hlm, Mass, Distance):
    mscale = (Mass*1.989e+30)/M
    dscale = (Distance*3.086e+22)/r
    f_x_strain = xStrain/mscale
    f_y_strain = mscale**(3/2) * yStrain / dscale
    
    fig = go.Figure()
    for i in range(len(f_x_strain)):
        fig.add_trace(go.Scatter(x=f_x_strain[i], y=f_y_strain[i],
                         line=dict(color='royalblue', width=1),
                         name=f"{hlm[i]}"))
        fig.add_annotation(x=math.log10(f_x_strain[i][0]), y=math.log10(f_y_strain[i][0]),
            text=f"{hlm[i]}",
            showarrow=True,
            xshift=-5,
            font=dict(
                color="mediumvioletred"
            )
        )
        #print(i)

    # Edit the layout
    fig.update_layout(
        title=dict(
            text='CE vs aLIGO comparison'
        ),
        xaxis=dict(
            range=[np.log10(3), np.log10(5e3)],
            title=dict(
                 text='frequency'
            )
        ),
        xaxis_type="log",
        yaxis=dict(
            range=[np.log10(1e-25), np.log10(5e-21)],
            title=dict(
                text='amplitude'
             ),
        ),
        yaxis_type="log"
    )
    return fig