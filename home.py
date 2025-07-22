import dash
from dash import html, dcc, register_page
import dash_bootstrap_components as dbc
import matplotlib
from dash import html, dcc, Input, Output
import plotly.graph_objects as go
from PIL import Image
import base64
import io
import matplotlib.pyplot as plt
import pyanatomogram
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import math
import matplotlib.patches as mpatches
import matplotlib.colors as mpcolor
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import plotly.express as px
import plotly.io as pio
import dash_bootstrap_components as dbc
from dash import Dash, dash_table
import plotly.colors as pc
from utilities import treatment_process_centered, tableprocess
from dash import register_page
from dash import callback, Input, Output, dcc



register_page(__name__, path="/")

treats=pd.DataFrame(treatment_process_centered('Juric_Rapid_Autopsy_MASTER-treatments1.csv', 'Juric_Rapid_Autopsy_MASTER-participants.txt', 'samples.6i2fn.tsv'))

layout = dbc.Container([
    dbc.Row([
        dbc.Col(
            html.H1("Welcome to the Autopsy Data Dashboard", 
                    className="text-center text-primary mb-4")
        )
    ]),
    dbc.Row([
        dbc.Col(
            html.H4("Explore cohort-wide and patient-specific insights.", 
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
        dbc.Col([
            html.P("Categories:"),
            dcc.Dropdown(
                id='names',
                options=[{'label': x.capitalize(), 'value': x} for x in ['race', 'analyte', 'gender','strategy']],
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

@callback(
    [Output("piechart", "figure"),
     Output("barchart", "figure")],
    Input("names", "value"))

def generate_chart(names):
    df=treats
    count_series = df[names].value_counts().reset_index()
    count_series.columns = [names, 'count']
    palette = px.colors.qualitative.Dark24
    fig = px.pie(count_series, values='count', names=names, hole=.3, color_discrete_sequence=palette)
    fig.update_layout(title=f'Distribution of {names}')
    
    bar_fig = px.bar(treats, x='year', y='ulp_tumor', labels={'year': 'Collection Year', 'ulp_tumor': 'ULP Tumor Fraction'})
    bar_fig.update_layout(title=f'Bar Chart of ULP Tumor Fraction', showlegend=False,  xaxis_title='Year of Diagnosis',
        yaxis_title='ULP Tumor Fraction',)
  
    bar_fig.update_traces(hoverinfo='none', selector=dict(type='bar'))

    return fig, bar_fig


