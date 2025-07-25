#Necessary Imports
import dash
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
from utilities import treatment_process_centered, tableprocess,  treatment_process_centeredbio, make_figure_with_background, tissue_to_site
from dash import register_page
from dash import callback, Input, Output
import lifelines
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from dash import register_page
from dash import callback, Input, Output

#Generates the Y position for the swimmers plot, puts the patients with the drug higher
def generate_patientYpos(patient_DF, selected_drug):
    all_patients = patient_DF['patient'].unique()
    patients_with_drug = patient_DF[patient_DF['drug'] == selected_drug]['patient'].unique()
    ordered_patients = list(patients_with_drug) + [p for p in all_patients if p not in patients_with_drug]
    y_pos = {patient: len(ordered_patients) - i for i, patient in enumerate(ordered_patients)}

    return y_pos


#centers the timeline based on the start date of the selected drug, leaves the other drugs as they were previously shifted
def shift_timeline_by_drug(df, drug):
    patient_offsets = (
        df[df['drug'] == drug]
        .groupby('patient')['tx_start']
        .min()
        .to_dict()
    )

    def compute_shift(row):
        offset = patient_offsets.get(row['patient'])
        if pd.isna(row['tx_start']):
            return pd.NA, pd.NA
        if offset is None:
            return row['tx_start'], row['tx_end']
        return row['tx_start'] - offset, row['tx_end'] - offset

    df[['tx_start_rel', 'tx_end_rel']] = df.apply(
        lambda row: pd.Series(compute_shift(row), index=['tx_start_rel', 'tx_end_rel']),
        axis=1
    )
    return df

#annotations of where all the different biospecimens appear on the map
annotations = {
    'thyroid gland': (214, 140),
    'esophagus': (223, 175),
    'lymph node': (240, 122),
    'lung': (200, 190),
    'left heart': (215, 104),
    'atrial heart':(220,204),
    'colon': (185, 365),
    'trans colon':(185, 360),
    'sigmoid colon':(240,380),
    'liver': (181, 275),
    'breast': (260, 215),
    'stomach': (230, 290),
    'uterus': (220, 370),
    'ovary': (247, 380),
    'pancreas': (262, 272),
    'gall bladder': (205, 295),
    'kidney cortex': (188, 305),
    'kidney med': (170,300),
    'skin': (170, 453),
    'other':(100,10),
    'muscle':(175, 460),
    'pleura': (200, 190),
    'peritoneum': (200, 190),
    'pericardium':(225, 195),
    'blood':(355,285),
    'cavity':(215,235),
    'spleen':(245,290),
    'spinal cord':(220, 250),
    'small intestine':(205,350),
    'diaphragm':(225,225),
    'chest':(240,190),
    'bone':(100,270),
    'bladder':(220,395),
    'coronary':(222,199),
    'adrenal gland':(185,300),
    'adipose':(240,280),
    'prostate gland':(219,395),
}

 #dictonary for the different stop values and all the markers associated with them
stop_marker_styles = {
        'Progression and relapse': dict(color='red', symbol='circle', size=10),
        'Completion of standard course': dict(color='black', symbol='triangle-down', size=10),
        'Toxicity': dict(color='green', symbol='triangle-down', size=10),
        'Other': dict(color='gray', symbol='circle', size=10)
    }
#making it so all the different forms of the same cohort get placed under the same category label in the drop down
morph_to_simple = {
    '8160/3 CHOLANGIOCARCINOMA': 'Chloangiocarcinoma',
    '(8521/3 INFILTRATING DUCTULAR CARCINOMA)': 'Infiltrating Ductular Carcinoma',
    '(8522/3 INFILTRATING DUCT AND LOBULAR CARCINOMA': 'Infiltrating Duct and Lobular Carcinoma',
    '(8520/3 INFILTRATING LOBULAR CARCINOMA':'Infiltrating Lobular Carcinoma',
    'ENDOMETRIOID ADENOCARCINOMA':'Endometrioid Adenocarinoma',
    'NON-SMALL CELL CARCINOMA': 'Non-Small Cell Carcinoma',
    'Metastatic Melanoma': 'Metastatic Melanoma',
    'Non-Small Cell Lung Cancer':'Non-Small Cell Lung Cancer',
    'Non small cell lung cancer': 'Non-Small Cell Lung Cancer',
    'NSCLC/Adenocarcinoma':'Non-Small Cell Lung Cancer and Adenocarcinoma',
    'NSCLC/poorly-differentiated carcinoma':'Non-Small Cell Lung Cancer and Carcinoma',
    'NSCLC/Squamous cell carcinoma': 'Non-Small Cell Lung Cancer and Squamous Cell Carcinoma',
    'Squamous Cell': 'Squamous Cell Carcinoma',
    'Metaplastic carcinoma':'Metaplastic carcinoma',
    'Adenocarcinoma':'Adenocarcinoma',
    'SIGNET RING CELL CARCINOMA':'Signet Ring Cell Carcinoma',
    'Hepatocellular carcinoma ':'Hepatocellular',
    'CRIBRIFORM CARCINOMA':'Cribiform Carcinoma',
    'MALIGNANT MELANOMA':'Melanoma',
    'Melanoma':'Melanoma',
    'Mullerian carcinoma':'Mullerian Carcinoma',
    'MEDULLARY THYROID CARCINOMA':'Medullary Thyroid Carcinoma',
    'Follicular Lymphoma ':'Follicular Lymphoma ',
    'Renal Cell Carcinoma':'Renal Cell Carcinoma',
    '8201/2':'Cribiform Carcinoma',
    'cholangiocarcinoma': 'Chloangiocarcinoma',
    'Malignant lymphoma':'Malignant Lymphoma',
    'moderately differentiated invasive adenocarcinoma':'Adenocarcinoma',
    'METASTATIC POORLY DIFFERENTIATED CARCINOMA':'Carcinoma',
    'Breast Cancer and NSCLC':'Non-Small Cell Lung Cancer and Breast Cancer ',
    'squamous cell':'Squamous Cell Carcinoma'

    
    
    
}

#mapping function to help creat this new variable
def map_site(text):
    if pd.isnull(text):
        return 'Unknown'
    for keyword, label in morph_to_simple.items():
        if keyword.lower() in text.lower(): 
            return label
    return 'Unknown'

#importing the data and doing some basic data management and creating some new variables
treats=pd.DataFrame(treatment_process_centeredbio('Juric_Rapid_Autopsy_MASTER-treatments1.csv', 'Juric_Rapid_Autopsy_MASTER-participants.txt', 'Juric_Rapid_Autopsy_MASTER-biospecimens.txt'))
treats['tx_end'] = pd.to_numeric(treats['tx_end'], errors='coerce')
xmax = math.ceil(treats['tx_end'].max())
treats['biospeclocation'] = treats['tissue_site'].map(tissue_to_site).fillna('Fail')
biospecdata=pd.DataFrame(tableprocess('Juric_Rapid_Autopsy_MASTER-biospecimens.txt', 'Juric_Rapid_Autopsy_MASTER-participants.txt'))
biospecdata['biospeclocation'] = biospecdata['tissue_site'].map(tissue_to_site).fillna('Fail')
treats['simplemorph'] = treats['morph'].apply(map_site).fillna('Unknown')
biospecdata['simplemorph'] = biospecdata['morph'].apply(map_site).fillna('Unknown')
treats['simplemorph'] = treats['simplemorph'].astype(str)

#creating the background image, making sure it is in the right form, this is where the io, and base 64 is necessary

buffer = io.BytesIO()
encoded_image = base64.b64encode(buffer.getvalue()).decode()

#registering this page as the cohort page in the dash app
register_page(__name__, path="/cohort")

#creating the layout of the page
layout = dbc.Container([
    dcc.Location(id='url', refresh=False),
    dbc.Row([
        dbc.Col(html.H2("Cohort Data", className="text-center text-primary mb-4"))
    ]),
    dbc.Row([ 
        dbc.Col(
            [
            dcc.Dropdown(
                id='cohort-dropdown',
                options=[{'label': morph, 'value': morph} for morph in sorted(treats['simplemorph'].dropna().unique())],
                placeholder="Select a Cohort",
                style={'backgroundColor': 'transparent', 'color': 'mediumslateblue', 'border': '1px solid white'}
            ),
            dcc.Dropdown(
                id='drug-dropdown',
                placeholder='Select Drug',
                style={'backgroundColor': 'transparent', 'color': 'mediumslateblue', 'border': '1px solid white'}
            )
        ],
        width=6,
        style={'backgroundColor': 'transparent'}
        )
    ]),
    dbc.Row([
        dbc.Col(
            dbc.Spinner(
                dcc.Graph(id='timeline-plot-with'),
                color='info',
                type='grow',
                fullscreen=False,
                size="sm"
            ),
            width=6,
            style={'backgroundColor': 'transparent'} 
        ),
        dbc.Col(
            dbc.Spinner(
                dcc.Graph(id='timeline-plot-without'),
                color='info',
                type='grow',
                fullscreen=False,
                size="sm"
            ),
            width=6,
            style={'backgroundColor': 'transparent'} 
        ),
    ]),
    dbc.Row([
        dbc.Col(
            dbc.Spinner(
                dcc.Graph(id='cohort-anatomy-plot'),
                color='info',
                type='grow',
                fullscreen=False,
                size="sm"
            ),
            width=6,
            style={'backgroundColor':'transparent'}
            
        ),
        dbc.Col(
            dash_table.DataTable(
                id='tableinfo2',
                columns=[],
                data=[],
                style_table={'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'minWidth': '100px', 'width': '150px', 'maxWidth': '200px'},
                page_size=20
            ),
            width=6,
            style={'backgroundColor': 'transparent'}
        )
    ], className="mb-4"),
    dbc.Row([   
        dbc.Col(
            dbc.Spinner(
                dcc.Graph(id='survival'),
                color='info',
                type='grow',
                fullscreen=False,
                size="sm"
            ),
            width=6,
            style={'backgroundColor':'transparent'}
        )
    ])
],fluid=True,  style={'backgroundColor': 'transparent'})


#callback for the drug selection drodown
@callback(
    Output('drug-dropdown', 'options'),
    Output('drug-dropdown', 'value'),
    Input('cohort-dropdown', 'value')
)
#this helps us decide what the drug options are going to be based on what drugs are in the selected cohort
def update_drug_options(selected_morph):
    if selected_morph is None:
        return [], None
    filtered = treats[treats['simplemorph'] == selected_morph]
    drugs = sorted(filtered['drug'].dropna().unique())
    options = [{'label': d, 'value': d} for d in drugs]

    value = drugs[0] if drugs else None
    return options, value

#this is the callback for both the with drug and without drug swimmer plot, takes inputs from the cohort and drug dropdowns
@callback(
    Output('timeline-plot-with', 'figure'),
    Output('timeline-plot-without', 'figure'),
    Input('cohort-dropdown', 'value'),
    Input('drug-dropdown', 'value')
)
#updates the swimmer plots based on the selected cohort and drug, refrences shifting the drug by timeline
def update_swimmer_figure(selected_morph, selected_drug):
    if selected_morph is None or selected_drug is None:
        return go.Figure(), go.Figure()
    
    filtered = treats[treats['simplemorph'] == selected_morph].copy()
    filtered = shift_timeline_by_drug(filtered, selected_drug)

    patients_with_drug = filtered[filtered['drug'] == selected_drug]['patient'].unique()

    
    with_drug = filtered[filtered['patient'].isin(patients_with_drug)]
    without_drug = filtered[~filtered['patient'].isin(patients_with_drug)]

    def create_figure(subset, title_suffix):
        patient_pos = generate_patientYpos(subset, selected_drug)
        unique_drugs = sorted(subset['drug'].dropna().unique())
        color_list = px.colors.qualitative.Alphabet
        while len(color_list) < len(unique_drugs):
            color_list += color_list
        drug_color_map = {drug: color_list[i % len(color_list)] for i, drug in enumerate(unique_drugs)}
        fig = go.Figure()
        plottable = subset.dropna(subset=['tx_start_rel', 'tx_end_rel'])
        added_drug_legends = set()

        for _, row in plottable.iterrows():
            drug = row['drug']
            color = drug_color_map.get(drug, 'lightgray')
            show_legend = drug not in added_drug_legends
            if show_legend:
                added_drug_legends.add(drug)
            patient_y = patient_pos[row['patient']]
            fig.add_trace(go.Scatter(
                x=[row['tx_start_rel'], row['tx_end_rel']],
                y=[patient_y, patient_y],
                mode='lines',
                line=dict(color=color, width=10),
                hoverinfo='text',
                name=drug,
                text=f"Patient: {row['patient']}<br>Drug: {drug}<br>Start: {row['tx_start_rel']}<br>End: {row['tx_end_rel']}",
                showlegend=show_legend
            ))

        for _, row in plottable.iterrows():
            stop = row['stop']
            if stop in stop_marker_styles:
                style = stop_marker_styles[stop]
                patient_y = patient_pos[row['patient']]
                fig.add_trace(go.Scatter(
                    x=[row['tx_end_rel']],
                    y=[patient_y],
                    mode='markers',
                    marker=dict(color=style['color'], symbol=style['symbol'], size=style['size']),
                    name=stop,
                    hoverinfo='text',
                    text=f"Stop Reason: {row['stop']}<br>Morph: {row['morph']}",
                    showlegend=False
                ))

        fig.update_layout(
            title=f'Treatment Timeline {title_suffix}',
            xaxis_title='Days from Drug Start',
            yaxis_title='Patient',
            yaxis=dict(
                type='category',
                categoryorder='array',
                categoryarray=list(reversed(list(patient_pos.keys())))
            ),
            height=max(600, len(patient_pos) * 30),
            margin=dict(l=120)
        )
        return fig

    fig_with = create_figure(with_drug, f'with {selected_drug}')
    fig_without = create_figure(without_drug, f'without {selected_drug}')

    return fig_with, fig_without

#cohort anatomy plot callback
@callback(
    Output('cohort-anatomy-plot', 'figure'),
    Input('cohort-dropdown', 'value')
)
def update_anatomy_figure(selected_morph):
    if not selected_morph:
        return make_figure_with_background(encoded_image)  

    patient_info = treats[treats['simplemorph'] == selected_morph]
    if patient_info.empty:
        return make_figure_with_background(encoded_image)

    gender = patient_info.iloc[0]['gender'].lower()

    fig, ax = plt.subplots(figsize=(6, 10))

    cohort_organs = patient_info['biospeclocation'].dropna().unique().tolist()

    cohort_biospecimens = biospecdata[biospecdata['simplemorph'] == selected_morph]
    
    anatomogram = pyanatomogram.Anatomogram('homo_sapiens.female')
    anatomogram.to_matplotlib(ax=ax)

    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)

    img = Image.open(buf)
    w, h = img.size
    img_str = base64.b64encode(buf.getvalue()).decode()

    fig_out = go.Figure()
    fig_out.add_layout_image(
        dict(
            source=f"data:image/png;base64,{img_str}",
            xref="x",
            yref="y",
            x=0,
            y=h,
            sizex=w,
            sizey=h,
            sizing="stretch",
            layer="below"
        )
    )

   
    organ_counts = {}
    blood_counts={}
  
    for _, row in cohort_biospecimens.iterrows():
        organ = row['biospeclocation']
        if row['ogmaterial'] == 'blood':
            blood_counts[organ] = blood_counts.get(organ, 0) + 1
        else:
            organ_counts[organ] = organ_counts.get(organ, 0) + 1
    unique_organs = list(organ_counts.keys())
    organ_colors = {organ: color for organ, color in zip(unique_organs, pc.qualitative.Plotly)}

  
    for organ, (x, y_pixel) in annotations.items():
        if organ in organ_counts:
            count = organ_counts[organ]
            fig_out.add_trace(go.Scatter(
                x=[x],
                y=[h - y_pixel],  # invert y axis for image coords
                mode='markers',
                marker=dict(
                    size=5 + count* 0.2,
                    color=organ_colors.get(organ, 'lightpink'),
                    opacity=0.6,
                    line=dict(width=1, color='darkgrey')
                ),
                name=organ,
                text=f"{organ}: {count} biospecimens",
                hoverinfo='text',
                customdata=[organ]
            ))
            
    total_blood = sum(blood_counts.values())

    if total_blood > 0:
        fig_out.add_trace(go.Scatter(
            x=[w * 0.2],  
            y=[h * 0.8],
            mode='markers',
            marker=dict(size=8 + total_blood * 2, color='red', symbol='diamond', opacity=0.6),
            name="Blood",
            text=f"Blood: {total_blood} biospecimens",
            hoverinfo='text',
            customdata=["blood"]
        ))

    fig_out.update_layout(
        title=f'Biospecimen Map',
        width=w,
        height=h,
        xaxis=dict(visible=False, range=[0, w]),
        yaxis=dict(visible=False, range=[0, h]),
        margin=dict(l=20, r=30, t=50, b=0),
        dragmode=False,
        clickmode='event+select',
        showlegend=False,
    )

    return fig_out

@callback(
    Output('tableinfo2', 'columns'),
    Output('tableinfo2', 'data'),
    Input('cohort-dropdown', 'value'),
    Input('cohort-anatomy-plot', 'clickData')
)
def createtable2(selected_morph, clickData):
    if selected_morph is None:
        return [], []

    cohort_data = biospecdata[biospecdata['simplemorph'] == selected_morph].copy()

    if cohort_data.empty:
        return [], []

    
    if clickData and 'points' in clickData:
        clicked_organ = clickData['points'][0]['customdata']
        cohort_data['biospeclocation'] = cohort_data['biospeclocation'].astype(str).str.lower().str.strip()
        clicked_organ = clicked_organ.lower().strip()
        cohort_data = cohort_data[cohort_data['biospeclocation'] == clicked_organ]

    displaycolumns = ['patient', 'ogmaterial', 'submaterial', 'collection', 'tissue_site', 'tissuedetails', 'main_normal']
    cohort_data = cohort_data[displaycolumns]

    custom_column = {
        'patient': 'Patient',
        'ogmaterial': 'Original Material Type',
        'submaterial': 'Submitted Material Type',
        'collection': 'Collection Date',
        'tissue_site': 'Sample Location',
        'tissuedetails': 'Tissue Site Details',
        'main_normal': 'Main Normal'
    }

    columns = [{'name': custom_column[col], 'id': col} for col in cohort_data.columns]
    data = cohort_data.to_dict('records')

    return columns, data


@callback(
    Output('survival', 'figure'),
    Input('cohort-dropdown', 'value')
)
def update_survival(selected_cohort):
    if not selected_cohort:
        return go.Figure(layout={"title": "Select a patient to see survival curve"})
    df = treats.copy()
    df['dead'] = np.where(df['vitalstatus'] == 'dead', 1, 0)

    if df['followup'].dtype == 'object':
        df['followup'] = df['followup'].str.replace(',', '').str.strip()
        df['followup'] = pd.to_numeric(df['followup'], errors='coerce')
        df['dead'] = pd.to_numeric(df['dead'], errors='coerce')
        
    df = df.dropna(subset=['followup', 'dead'])

    df['followup'] = df['followup'].astype(int)
    df['dead'] = df['dead'].astype(int)

    
    cohort_data = df[df['simplemorph'] == selected_cohort]
    
    
    fig = go.Figure()
    color_list = px.colors.qualitative.Set3
    category_list = cohort_data['categories'].dropna().unique()
    color_map = {cat: color_list[i % len(color_list)] for i, cat in enumerate(category_list)}

    
    for i, cat in enumerate(category_list):
        group = cohort_data[cohort_data['categories'] == cat]
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
        title=f'Survival Curves by Treatment Category in {selected_cohort}',
        xaxis_title='Follow-up Time (days)',
        yaxis_title='Survival Probability',
        template='plotly_white',
        legend_title='Treatment Category'
    )
 
    return fig



