"""
Dash app displaying correlations across genomic and trait features
"""
from base64 import b64encode
import io

from dash import Dash, dcc, html, Input, Output
import pandas as pd
import plotly.express as px

FEATURES = ['coding_density',
 'gc_percentage',
 'genome_size',
 'all_ratio_acidic_pis',
 'all_mean_protein_length',
 'all_mean_pi',
 'all_mean_gravy',
 'all_mean_zc',
 'all_mean_nh2o',
 'all_mean_f_ivywrel',
 'all_weighted_mean_f_ivywrel',
 'all_weighted_mean_zc',
 'all_weighted_mean_nh2o',
 'all_weighted_mean_gravy',
 'extracellular_soluble_ratio_acidic_pis',
 'extracellular_soluble_mean_protein_length',
 'extracellular_soluble_mean_pi',
 'extracellular_soluble_mean_gravy',
 'extracellular_soluble_mean_zc',
 'extracellular_soluble_mean_nh2o',
 'extracellular_soluble_mean_f_ivywrel',
 'extracellular_soluble_weighted_mean_f_ivywrel',
 'extracellular_soluble_weighted_mean_zc',
 'extracellular_soluble_weighted_mean_nh2o',
 'extracellular_soluble_weighted_mean_gravy',
 ]

TARGETS = [
 'salinity_min', 
 'salinity_midpoint',
 'salinity_max',
 'temperature_min',
 'temperature_optimum',
 'temperature_max',
 'ph_min',
 'ph_optimum',
 'ph_max'
 ]

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets)

df = pd.read_csv('./data/feature_table.tsv.gz', sep='\t', compression='gzip')

features = [col for col in FEATURES if col in df.columns.tolist()]
targets = [col for col in TARGETS if col in df.columns.tolist()]

corr = df[features + targets].corr(method='spearman').loc[targets, features]
fig_heatmap = px.imshow(corr, text_auto=False, color_continuous_scale='RdBu_r')
fig_heatmap.update_layout(template='simple_white')

app.layout = html.Div(
    [
        dcc.Markdown("""
        # Correlations between Genomic Features and Traits

        ###### Plots:

        - **Scatter plot**: You may investigate the correlation between genomic features (x-axis) and traits (y-axis)
        using the toggle buttons to select x or y variables and the slider to restrict ranges. 
        - **Heatmap**: Spearman correlation coefficients are encoded as colors. A Spearman correlation does not need to
        be a linear relationship; it simply means that as one variable increases, the other variable increases (positive
        correlation) or decreases (negative correlation).

        ###### How to analyze this data: 
        
        _Remember that correlations may not be causal._ 
        
        Traits like halophily, acidophily, 
        alkaliphily, and (hyper)thermophily should be obvious in this dataset - but alas are only in small groups of organisms. 
        Most organisms prefer average conditions. Together, this reflect a sampling bias towards moderate conditions 
        and extreme conditions. For model training, it is preferable to have correlations with a gradation across conditions
        so we can properly interpolate between moderate and extreme conditions. Sometimes, correlations may appear weak, but 
        machine learning methods can still make use of correlations within segments of the data. Correlations may be strongly affected by outliers or
        other abnormal features in the data, while portions of the data have promise. Use the heatmap to observe correlation strength
        at a glance, then inspect scatter plots in more detail for quality of correlation.

        Technical error includes: assigning strains to different genomes within the same species, variable measurements of traits
        across laboratories, sampling bias.

        (i): To download an image, use the menu found by floating your mouse in the top right corner of either plot.
        """, style={'family' : 'helvetica'}),

        html.Div(className='ten rows', children=[
            html.Div(className='six columns', children=[
                dcc.Graph(id="scatter-plot", style={'height': '60vh'})
                ]),
            html.Div(className='six columns', children=[
                dcc.Graph(figure=fig_heatmap, id='heatmap', style={'height': '40vh'})
                ]),
        ]),

        html.Div(className='row', children=[
            dcc.Markdown("""
            **Scatter plot**:
            
            Select x variable (genome feature):""", style={'family' : 'helvetica'}),
            dcc.RadioItems(id='x-variable',
                options=features,
                value='all_mean_f_ivywrel',
                inline=True,
                ),
            dcc.Markdown("\nSelect y variable (microbial trait):", style={'family' : 'helvetica'}),
            dcc.RadioItems(id='y-variable',
                options=targets,
                value='temperature_optimum',
                inline=True,
                )
        ]),

        dcc.Markdown("""
        ###### Additional information:

        **Microbial traits** are measured directly from strains. Min and max refer to minimum and maximum reported values
        are expected to be the least accurate traits because different labs report different ranges. Optimum is the reported
        optimum growth condition and midpoint refers to average of the minimum and maximum. Temperature is in Celsius and 
        salinity is in % weight by volume sodium chloride.

        **Genome features** are measured from genomes and will be used to predict traits. Features based on proteins are divided into
        two categories: "all" means all proteins, and "extracellular soluble" means only proteins predicted to be extracellular and soluble
        (i.e. not membrane bound), which should reflect selection by extracellular conditions

        - Genome size: total basepairs
        - Coding density: percent of basepairs within coding regions
        - GC percentage: percent of basepairs that are G or C
        - pI: protein isoelectric point (ratio acidic refers to proteins < pI 7 vs. proteins > pI 7)
        - Fraction IVYWREL: fraction of amino acids that are IVYWREL, residues that lend more stability
        - GRAVY: grand average of hydropathy, which correlates with intermembrane regions
        - Zc: average carbon oxidation state of proteins
        - nH2O: stoichiometric hydration state (nH2O) of proteins

        """, style={'family' : 'helvetica'}),

    ]
)

@app.callback(
    Output("scatter-plot", "figure"),
    Input("x-variable", "value"),
    Input("y-variable", "value"),
)
def update_chart(x, y):
    mask = df[~df[x].isnull() & ~df[y].isnull()]
    fig = px.scatter(mask, x=x, y=y, title="n={} Genomes/Strains".format(len(mask)))
    fig.update_traces(marker=dict(size=8,
                                color='black',
                                line=dict(width=0, color='black')),
                    selector=dict(mode='markers'))
    fig.update_layout(template='simple_white')
    return fig

if __name__ == "__main__":
    app.run_server(debug=True)