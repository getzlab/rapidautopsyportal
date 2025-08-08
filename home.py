import dash
from dash import html, dcc, register_page, callback, Input, Output
import dash_bootstrap_components as dbc
import matplotlib
from dash import html, dcc, Input, Output
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import math
import plotly.express as px
import plotly.io as pio
import dash_bootstrap_components as dbc
import plotly.colors as pc
from utilities import treatment_process_centered, tableprocess
from dash import register_page
from dash import callback, Input, Output, dcc



register_page(__name__, path="/")

#loading in the data
treats=pd.DataFrame(treatment_process_centered('Juric_Rapid_Autopsy_MASTER-treatments1.csv', 'Juric_Rapid_Autopsy_MASTER-participants.txt', 'samples.6i2fn.tsv'))
participants=pd.read_csv('Juric_Rapid_Autopsy_MASTER-participants.txt', sep='\t')
samples=pd.read_csv('samples.6i2fn.tsv', sep='\t')

site_to_organ = {
    'BREAST (C50)': 'breast',
    'BRONCHUS AND LUNG (C34)': 'lung',
    'COLON (C18)': 'colon',
    'CORPUS UTERI (C54)': 'uterus',
    'ESOPHAGUS (C15)': 'esophagus',
    'GALLBLADDER (C23)': 'gall bladder',
    'KIDNEY (C64)': 'kidney',
    'LIVER AND INTRAHEPATIC BILE DUCTS (C22)': 'liver',
    'LYMPH NODES (C77)': 'lymph node',
    'OTHER AND UNSPECIFIED PARTS OF BILIARY TRACT (C24)':'bilary tract',
    'OVARY (C56)': 'ovary',
    'PANCREAS (C25)': 'pancreas',
    'SKIN (C44)': 'skin',
    'STOMACH (C16)': 'stomach',
    'THYROID GLAND (C73)': 'thyroid gland',
    'PROSTATE GLAND (C61)': 'prostate gland',
}
treats['organs'] = treats['site'].map(site_to_organ).fillna('Fail')

biospecdata=pd.DataFrame(tableprocess('Juric_Rapid_Autopsy_MASTER-biospecimens.txt', 'Juric_Rapid_Autopsy_MASTER-participants.txt'))
biospecdata['organs'] = biospecdata['site'].map(site_to_organ).fillna('Fail')

dropdown_options = [
    {'label': 'Race', 'value': 'race'},
    {'label': 'Gender', 'value': 'gender'},
    {'label': 'Analyte Type', 'value': 'analyte_type'},
    {'label': 'Experimental Strategy', 'value': 'experimental_strategy'}
]

#Layout for the page
layout = dbc.Container([
    dbc.Row([
        dbc.Col(
            html.H1("Welcome to the Autopsy Data Dashboard", 
                    className="text-center text-primary mb-4")
        )
    ]),
    dbc.Row([
        dbc.Col(
            html.H4("Explore cohort-wide, and drug or patient-specific insights.", 
                    className="text-secondary text-center mb-4")
        )
    ]),
    dbc.Row([
        dbc.Col(
            dbc.Button("Go to Patient Page", href="/patients", color="primary", size="lg", className="d-block mx-auto")
        )
    ]),
    dbc.Row([
        dbc.Col(
            dbc.Button("Go to Cohort Page", href="/cohort", color="primary", size="lg", className="d-block mx-auto")
        )
    ]),
    dbc.Row([
        dbc.Col(
            dbc.Button("Go to Drug Page", href="/drug", color="primary", size="lg", className="d-block mx-auto")
        )
    ]),
    dbc.Row([
        dbc.Col([
            html.P("Categories:"),
            dcc.Dropdown(
                id='names',
                options=dropdown_options,
                value='race',
                clearable=False
            ),
        ], width=4) 
    ]),
        dbc.Col(
            dcc.Graph(id="piechart"),
            width=12  
        ),
    dcc.Graph(id="barchart", config={'staticPlot': True}),
], fluid=True, style={'padding': '60px'})

#callback for both the pie chart and the barchart
@callback(
    [Output("piechart", "figure"),
     Output("barchart", "figure")],
    Input("names", "value"))

#creating the pie chart and the bar chart
def generate_chart(names):
    names=names.strip().lower()
    if names in participants.columns:
        df = participants
    elif names in samples.columns:
        df = samples
    else:
        df = pd.DataFrame(columns=[names])
 
    count_series = df[names].value_counts().reset_index()
    count_series.columns = [names, 'count']
    palette = px.colors.qualitative.Dark24
    fig = px.pie(count_series, values='count', names=names, hole=.3, color_discrete_sequence=palette)
    fig.update_layout(title=f'Distribution of {names.title()}')
    
    organ_counts = biospecdata['organs'].value_counts().reset_index()
    organ_counts.columns = ['organs', 'count']
    organ_counts = organ_counts.sort_values(by='count', ascending=False)
    bar_fig = px.bar(organ_counts, x='organs', y='count', labels={'organs': 'Cancer Type', 'count':'Biospecimen Count'})
    bar_fig.update_layout(title=f'Biospecimens by Cancer Type', showlegend=False,  xaxis_title='Cancer Type',
        yaxis_title='Biospecimen Count')
  
    bar_fig.update_traces(hoverinfo='none', selector=dict(type='bar'),
                          text=organ_counts['count'], textposition='outside' , marker_color='#ff48a5')

    return fig, bar_fig


