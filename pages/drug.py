import dash
import matplotlib
from dash import html, dcc, Input, Output, callback, register_page, Dash, dash_table
import plotly.graph_objects as go
from PIL import Image
import base64
import io
import matplotlib.pyplot as plt
import pyanatomogram
import dash_ag_grid as dag
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import math
import plotly.express as px
import plotly.io as pio
import dash_bootstrap_components as dbc
import plotly.colors as pc
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utilities import treatment_process_centered, treatmentandpatient
import lifelines
from lifelines import KaplanMeierFitter



#loading in the data
tumordata=pd.DataFrame(treatment_process_centered('Juric_Rapid_Autopsy_MASTER-treatments1.csv', 'Juric_Rapid_Autopsy_MASTER-participants.txt', 'samples.6i2fn.tsv'))
treatments=pd.read_csv('Juric_Rapid_Autopsy_MASTER-treatments1.csv', sep=',')

#HARDCODED DRUG OPTIONS
dropdownoptions={'Abemaciclib', 'Fulvestrant', 'Abiraterone Acetate', 'Afatinib Dimaleate', 'AFP-C332',
                 'Alectinib', 'Alpelisib', 'Binimetinib', 'ALX-148', 'Amivantamab', 'Anastrozole', 'Leuprolide Acetate', 'Trastuzumab', 'Atezolizumab',
                 'Cobimetinib', 'Vermurafenib', 'Emactuzumab', 'Avelumab', 'Axitinib', 'Pembrolizumab', 'Batoprotafib', 'Bevacizumab', 'Cyclophosphamide',
                 'Doxorubicin', 'Hydrochloride Liposome', 'Everolimus', 'Fluorouracil', 'Irinotecan Hydrochloride', 'Leucovorin Calcium', 'Oxaliplatin', 'Nivolumab',
                 'Bicalutamide', 'Encorafenib','Brigatinib', 'Brilanestrant', 'Buthionine Sulfoximine','Cabozantinib', 'Capecitabine', 
                 'Epirubicin', 'Exemestane', 'Capmatinib', 'CAR-T', 'Carboplatin', 'Gemcitabine Hydrochloride', 'Paclitaxel', 'Temsirolimus', 'Cemiplimab', 'Ceritinib',
                 'Cetuximab', 'Cisplatin', 'CLR-457', 'Cobicistat','Vemurafenib', 'Crizotinib', 'Docetaxel', 'Pegfilgrastism',
                 'Dabrafenib Mesylate', 'Trametinib','Ridaforolimus', 'Dasatinib', 'Datopotamab','Debio-1347','Denosumab', 'Dinaciclib', 'Veliparib', 'DKN-01',
                 'Dostarlimab', 'Durvalumab', 'Entinostat', 'Entrectinib', 'Enzalutamide', 'Eribulin Mesylate', 'Pertuzumab', 'Goserelin Acetate',
                 'Erlotinib', 'Etoposide', 'Ribociclib', 'Palbociclib', 'Vinorelbine Tartrate', 'Zoledronic Acid', 'Floxuridine', 'Folinic Acid',
                 'Futibatinib', 'Ganetespib', 'Gefitinib', 'Rituximab', 'Golvatinib', 'Lenvatinib', 'H3B-6545', 'Halotestin', 'Hydroxychloroquine',
                 'Ibrutinib', 'Imatinib Mesylate', 'Imiquimod', 'IMMU-132', 'Inavolisib', 'Interleukin 2', 'Ipilimumab', 'Ladiratuzumab vedotin', 'Lamivudine',
                 'Lapatinib Ditosylate', 'Letrozole', 'Pamidronate Disodium', 'Tamoxifen Citrate', 'LFA-102', 'Lisocabtagene Maraleucel', 'Lorlatinib',
                 'Losatuxizumab Vedotin', 'LOXO', 'LSZ102', 'Luminespib', 'Lurbinectedin', 'LY-3214996', 'MAGE TCR', 'Margetuximab', 'Megestrol Acetate',
                 'Methotrexate', 'Mevociclib', 'MK-2206', 'Mobocertinib', 'Navitoclax', 'Nazartinib', 'Olaparib', 'Olutasidenib', 'Osimertinib',
                 'Ramucirumab', 'Panitumumab', 'Ziv-Aflibercept', 'Pemetrexed', 'Ponatinib', 'Pralsetinib', 'Raloxifene Hydrochloride','Regorafenib', 'Repotrectinib',
                 'Rezatapopt', 'RLY', 'Rociletinib', 'Roferon-A', 'Sapanisertib', 'SAR', 'Savolitinib', 'Selpercatinib', 'Serabelisib', 'Seribantumab',
                 'Seviteronel', 'Sitravatinib', 'Sorafenib', 'Sotrastaurin Acetate', 'SRN', 'Subasumstat', 'Taletrectinib', 'Talimogene Laherparepvec',
                 'Taselisib', 'Telisotuzumab', 'Tisagenlecleucel', 'Tislelizumab', 'Toremifene', 'Tucatinib', 'Trifluridine and Tipiracil Hydrochloride', 'Ulixertinib',
                 'Vandetanib', 'WDVAX', 'XMT', 'Zotizalkib', 'Neratinib Maleate'}

#testing whether it should be flagged as having both pre and post treatments
def testfunction(df):
    df = df.copy()
    df['flagged'] = 'No' 

    for index, row in df.iterrows():
        current_patient = row['participant_id']
        current_start = row['start_date_dfd']
        current_stop = row['stop_date_dfd']

        # All treatments for this patient
        patient_df = df[df['participant_id'] == current_patient]

        # Check if any treatments are before and after the current one
        pre_exists = any(patient_df['start_date_dfd'] < current_start)
        post_exists = any(patient_df['stop_date_dfd'] > current_stop)

        if pre_exists and post_exists:
            df.at[index, 'flagged'] = 'Yes'

    return df

    

register_page(__name__, path="/drug")

#page layout
layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.H2("Drug-Specific Data", className="text-center text-primary mb-4"))
    ]),
    dbc.Row([
        dbc.Col(html.H2("Select a Drug to View Information", className="text-secondary mb-3"))
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
            dag.AgGrid(
                id='drugtable',
                columnDefs=[],  
                rowData=[],     
                className="ag-theme-quartz",
                style={'height': '500px', 'width': '100%'},
                defaultColDef={
                    "resizable": True,
                    "sortable": True,
                    "filter": True
                }  
            ),
            width=6,
            style={'backgroundColor': 'transparent'}
        ),
        dbc.Col(
            dbc.Spinner(
                dcc.Graph(id='drugsurvival'),
                color='info',
                type='grow',
                fullscreen=False,
                size="sm"
            ),
            width=6,
            style={'backgroundColor':'transparent'}
        )
    ]),
    dcc.Graph(id="flaggedcountbar", config={'staticPlot': True}),
    
],fluid=True,  style={'backgroundColor': 'transparent'})

#callback for drug table
@callback(
    (Output('drugtable', 'columnDefs'),
    Output('drugtable', 'rowData'),
    Output('flaggedcountbar', 'figure'),
    ),
    Input('drug-dropdown', 'value')
)

#creating drug table, as well as static bar graph
def create_drug_table(selected_drug):
    if selected_drug is None:
        return [], [], go.Figure()
    



    flagged_df = testfunction(treatments)

    table_df = flagged_df[flagged_df['drugs'].str.contains(selected_drug, case=False, na=False)]
    displaycolumns = ['participant_id', 'drugs', 'treatment_regimen','notes', 'flagged']
    table_df = table_df[displaycolumns]

    custom_column = {
        'participant_id': 'Patient',
        'drugs': 'Drug',
        'flagged': 'Flag',
        'notes': 'Notes',
        'treatment_regimen': 'Treatment Regimen'
    }

    columnDefs = [{'headerName': custom_column[col], 'field': col} for col in table_df.columns]
    rowData = table_df.to_dict('records')


    flagged_counts = (
        flagged_df[flagged_df['flagged'] == 'Yes']
        .groupby('drugs')
        .size()
        .reset_index(name='count')
        .sort_values(by='count', ascending=False)
        .head(10)
    )

    bar_fig = px.bar(
        flagged_counts,
        x='drugs',
        y='count',
        labels={'drugs': 'Drug', 'count': 'Number of Flagged Patients'},
        title='Flagged Patients by Drug'
    )

    bar_fig.update_layout(
        showlegend=False,
        xaxis_title='Drug',
        yaxis_title='Number of Flagged Patients',
        xaxis_tickangle=45,
        margin=dict(t=100,b=200)
    )

    bar_fig.update_traces(
        text=flagged_counts['count'],
        textposition='auto',
        marker_color="#ea5a41"
    )

    return columnDefs, rowData, bar_fig

#callback for survival curve
@callback(
    Output('drugsurvival', 'figure'),
    Input('drug-dropdown', 'value')
)

#survival curve plot
def update_survival(selected_drug):
    if not selected_drug:
        return go.Figure(layout={"title": "Select a drug to see survival curve"})
    
    treatments2=pd.DataFrame(treatmentandpatient('Juric_Rapid_Autopsy_MASTER-treatments1.csv', 'Juric_Rapid_Autopsy_MASTER-participants.txt'))
    df = testfunction(treatments2)
    df = df[df['drug'].str.contains(selected_drug, case=False, na=False)]
    
    df['dead'] = np.where(df['vitalstatus'] == 'dead', 1, 0)
        
    df = df.dropna(subset=['followup', 'dead'])

    df['followup'] = df['followup'].astype(int)
    df['dead'] = df['dead'].astype(int)
    
    
    fig = go.Figure()
    color_list = px.colors.qualitative.Set2
    category_list = df['flagged'].dropna().unique()
    color_map = {cat: color_list[i % len(color_list)] for i, cat in enumerate(category_list)}

    
    for i, cat in enumerate(category_list):
        group = df[df['flagged'] == cat]
        if group.empty:
            continue

        kmf = KaplanMeierFitter()
        kmf.fit(group['followup'], event_observed=group['dead'], label=cat)

        fig.add_trace(go.Scatter(
            x=kmf.survival_function_.index,
            y=kmf.survival_function_[kmf._label],
            mode='lines',
            name=cat,
            line=dict(color=color_map[cat])
        ))
      


    fig.update_layout(
        title=f'Survival Curves by Treatment Category in {selected_drug}',
        xaxis_title='Follow-up Time (days)',
        yaxis_title='Survival Probability',
        template='plotly_white',
        legend_title='Treatment Category'
    )
 
    return fig
