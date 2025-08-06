import dash
import matplotlib
from dash import html, dcc, Input, Output, callback, register_page, Dash, dash_table
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
import plotly.colors as pc
from utilities import treatment_process_centered, tableprocess,  treatment_process_centeredbio, make_figure_with_background, tissue_to_site






tumordata=pd.DataFrame(treatment_process_centered('Juric_Rapid_Autopsy_MASTER-treatments1.csv', 'Juric_Rapid_Autopsy_MASTER-participants.txt', 'samples.6i2fn.tsv'))
treatments=pd.read_csv('Juric_Rapid_Autopsy_MASTER-treatments1.csv', sep=',')

dropdownoptions={'Abemaciclib', 'Fulvestrant', 'Abiraterone Acetate', 'Afatinib Dimaleate', 'AFP-C332',
                 'Alectinib', 'Alpelisib', 'Binimetinib', 'ALX-148', 'Amivantamab', 'Anastrozole', 'Leuprolide Acetate', 'Trastuzumab', 'Atezolizumab',
                 'Cobimetinib', 'Vermurafenib', 'Emactuzumab', 'Avelumab', 'Axitinib', 'Pembrolizumab', 'Batoprotafib', 'Bevacizumab', 'Cyclophosphamide',
                 'Doxorubicin', 'Hydrochloride Liposome', 'Everolimus', 'Fluorouracil', 'Irinotecan Hydrochloride', 'Leucovorin Calcium', 'Oxaliplatin', 'Nivolumab',
                 'Bicalutamide', 'Encorafenib','Brigatinib', 'Brilanestrant', 'Buthionine Sulfoximine','Cabozantinib', 'Capecitabine', 
                 'Epirubicin', 'Exemestane', 'Capmatinib', 'CAR-T', 'Carboplatin', 'Gemcitabine Hydrochlroide', 'Paclitaxel', 'Temsirolimus', 'Cemiplimab', 'Ceritinib',
                 'Cetuximab', 'Cisplatin', 'CLR-457', 'Cobicistat','Vemurafenib', 'Crizotinib', 'Docetaxel', 'Pegfilgrastism',
                 'Dabrafenib Mesylate', 'Trametinib','Ridaforolimus', 'Dasatinib', 'Datopotamab','Debio-1347','Denosumab', 'Dinaciclib', 'Veliparib', 'DKN-01',
                 'Dostarlimab', 'Durvalumab', 'Entinostate', 'Entrectinib', 'Enzalutamide', 'Eribulin Mesylate', 'Pertuzumab', 'Goserelin Acetate',
                 'Erlotinib', 'Etoposide', 'Ribociclib', 'Palbociclib', 'Vinorelbine Tartrate', 'Zoledronic Acid', 'Floxuridine', 'Folinic Acid',
                 'Neratinib Maleate', 'Futibatinib', 'Ganetespib', 'Gefitinib', 'Rituximab', 'Golvatinib', 'Lenvatinib', 'H3B-6545', 'Halotestin', 'Hydroxychloroquine',
                 'Ibrutinib', 'Imatinib Mesylate', 'Imiquimod', 'IMMU-132', 'Inavolisib', 'Interleukin 2', 'Ipilimumab', 'Ladiratuzumab vedotin', 'Lamivudine',
                 'Lapatinib Ditosylate', 'Letrozole', 'Pamidronate Disodium', 'Tamoxifen Citrate', 'LFA-102', 'Lisocabtagene maraleucel', 'Lorlatinib',
                 'Losatuxizumab Vedotin', 'LOXO', 'LSZ102', 'Luminespib', 'Lurbinectedin', 'LY-3214996', 'MAGE TCR', 'Margetuzimab', 'Megestrol Acetate',
                 'Methotrexate', 'Mevociclib', 'MK-2206', 'Mobocertinib', 'Navitoclax', 'Nazartinib', 'Olaparib', 'Olutasidenib', 'Osimertinib',
                 'Ramucirumab', 'Panitumumab', 'Ziv-Aflibercept', 'Pemetrexed', 'Ponatinib', 'Pralsetinib', 'Raloxifene Hydrochloride','Regorafenib', 'Repotrectinib',
                 'Rezatapopt', 'RLY', 'Rociletinib', 'Roferon-A', 'Sapanisertib', 'SAR', 'Savolitinib', 'Selpercatinib', 'Serabelisib', 'Seribantumab',
                 'Seviteronel', 'Sitravatinib', 'Sorafenib', 'Sotrastaurin Acetate', 'SRN', 'Subasumstat', 'Taletrectinib', 'Talimogene Laherparepvec',
                 'Taselisib', 'Telisotuzumab', 'Tisagenlecleucel', 'Tislelizumab', 'Toremifene', 'Tucatinib', 'Trifluridine and Tipiracil Hydrochloride', 'Ulixertinib',
                 'Vandetanib', 'WDVAX', 'XMT', 'Zotizalkib'}


register_page(__name__, path="/drug")

layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.H2("Drug-Specific Data", className="text-center text-primary mb-4"))
    ]),
    dbc.Row([ 
        dbc.Col(
            dcc.Dropdown(
                id='drug-dropdown',
                placeholder='Select Drug',
                options=[{'label': drug, 'value': drug} for drug in sorted(dropdownoptions)],
                style={'backgroundColor': 'transparent', 'color': 'mediumslateblue', 'border': '1px solid white'}
            ),
        width=6,
        style={'backgroundColor': 'transparent'}
        )
    ]),
    dbc.Row([
        dbc.Col(
            dash_table.DataTable(
                id='drugtable',
                columns=[],
                data=[],
                style_table={'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'minWidth': '100px', 'width': '150px', 'maxWidth': '200px'},
                page_size=20
            ),
            width=6,
            style={'backgroundColor': 'transparent'}
        )
    ])
])

@callback(
    Output('drugtable', 'columns'),
    Output('drugtable', 'data'),
    Input('drug-dropdown', 'value')
)
def create_drug_table(selected_drug):
    if selected_drug is None:
        return [], []
    

    filtered_df = treatments[treatments['drugs'].str.contains(selected_drug, case=False, na=False)]
    
    displaycolumns = ['participant_id', 'drugs']
    filtered_df = filtered_df[displaycolumns]

    custom_column = {
        'participant_id': 'Patient',
        'drugs': 'Drug'
    }

    columns = [{'name': custom_column[col], 'id': col} for col in filtered_df.columns]
    data = filtered_df.to_dict('records')

    return columns, data
