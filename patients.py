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
import lifelines
from lifelines import KaplanMeierFitter
kmf=KaplanMeierFitter()
from lifelines import CoxPHFitter
from utilities import treatment_process_centered, tableprocess, treatment_process_centeredbio, make_figure_with_background
from dash import register_page
from dash import callback, Input, Output



matplotlib.use('Agg')


#generate the y position needed for the swimmer plot
def generate_patientYpos(patient_DF):
    """
    Return a dictionary with unique patients as keys, and 
    y positions as values, for plotting. Assumes the column with 
    the patient names is called "patient"
    """
    
    patient_ids=list(patient_DF['patient'])
    y_pos={}
    new_pos=0000
    
    for patient_name in patient_ids:
        if patient_name not in y_pos:
            new_pos=new_pos+1
            y_pos[patient_name]=new_pos
    
    return y_pos

#this function makes the background imaging



def add_legend_annotation(fig, labels, colors, row, col, xshift=0.1, yshift=0.05):
    base_x = 0.0 + (col - 1.2) * 0.5 + xshift
    base_y = 1.0 - (row - 0.75) * 0.25 - yshift
    for i, (label, color) in enumerate(zip(labels, colors)):
        fig.add_annotation(
            x=base_x,
            y=base_y - i * 0.04,
            text=f"<span style='color:{color}'>■</span> {label}",
            showarrow=False,
            xanchor="left",
            yanchor="top",
            font=dict(size=12)
        )
#Loading in the Data Set
tumordata=pd.DataFrame(treatment_process_centered('Juric_Rapid_Autopsy_MASTER-treatments1.csv', 'Juric_Rapid_Autopsy_MASTER-participants.txt', 'samples.6i2fn.tsv'))
biospecdata=pd.DataFrame(tableprocess('Juric_Rapid_Autopsy_MASTER-biospecimens.txt', 'Juric_Rapid_Autopsy_MASTER-participants.txt'))
tumordata['tx_end'] = pd.to_numeric(tumordata['tx_end'], errors='coerce')
xmax = math.ceil(tumordata['tx_end'].max())
patient_pos=generate_patientYpos(tumordata)
drugtype = tumordata['drug'].unique()
patients = tumordata['patient'].unique()
#Data managmenet to get the tumor primary site to be in a workable form in order to automatically plot using pyanatomogram package
site_to_organ = {
    'BREAST (C50)': 'breast',
    'COLON (C18)': 'colon',
    'BRONCHUS AND LUNG (C34)': 'lung',
    'CORPUS UTERI (C54)': 'uterus',
    'ESOPHAGUS (C15)': 'esophagus',
    'GALLBLADDER (C23)': 'gall bladder',
    'KIDNEY (C64)': 'kidney',
    'LIVER AND INTRAHEPATIC BILE DUCTS (C22)': 'liver',
    'LYMPH NODES (C77)': 'lymph node',
    'OVARY (C56)': 'ovary',
    'PANCREAS (C25)': 'pancreas',
    'SKIN (C44)': 'skin',
    'STOMACH (C16)': 'stomach',
    'THYROID GLAND (C73)': 'thyroid gland',
    'PROSTATE GLAND (C61)': 'prostate gland',
    'OTHER AND UNSPECIFIED PARTS OF BILIARY TRACT (C24)':'bilary tract'
}

tissue_to_site={
    'Lung (UBERON:0008952)':'lung',
    'Muscle - Skeletal (UBERON:0011907)': 'muscle',
    'Whole Blood (UBERON:0013756)': 'blood',
    'Uterus (UBERON:0000995)': 'uterus',
    'Thyroid (UBERON:0002046)' :'thyroid gland',
    'Thoracic Cavity (UBERON:0002224)': 'cavity',
    'Stomach (UBERON:0000945)':'stomach',
    'Spleen (UBERON:0002106)':'spleen',
    'Spinal Cord (UBERON:0002240)':'spinal cord',
    'Small Intestine - Terminal Ileum (UBERON:0001211)':'small intestine',
    'Skin of Body (UBERON:0002097)':'skin',
    'Skin - Sun Exposed (Lower leg) (UBERON:0004264)':'skin',
    'Skin - Not Sun Exposed (Suprapubic) (UBERON:0036149)':'skin',
    'Prostate (UBERON:0002367)':'prostate',
    'Pleural Effusion (UBERON:0000175)': 'pleura',
    'Pleura (UBERON:0000977)':'pleura',
    'Peritoneum (UBERON:0002358)':'peritoneum',
    'Pericardium (UBERON:0002407)':'pericardium',
    'Pancreas (UBERON:0001150)':'pancreas',
    'Other - Anatomical Entity (UBERON:0001062)':'other',
    'Lymph Node (UBERON:0000029)':'lymph node',
    'Liver (UBERON:0001114)':'liver',
    'Kidney - Medulla (UBERON:0001293)':'kidney med',
    'Kidney - Cortex (UBERON:0001225)':'kidney cortex',
    'Heart - Left Ventricle (UBERON:0006566)': 'left heart',
    'Heart - Atrial Appendage (UBERON:0006631)':'atrail heart',
    'Esophagus - Mucosa (UBERON:0006920)':'esophagus',
    'Diaphragm (UBERON:0001103)':'diaphragm',
    'Colon - Transverse (UBERON:0001157)':'trans colon',
    'Colon - Sigmoid (UBERON:0001159)':'sigmoid colon',
    'Colon (UBERON:0001155)':'colon',
    'Chest Wall (UBERON:0016435)':'chest',
    'Breast - Mammary Tissue (UBERON:0008367)':'breast',
    'Brain - Frontal Cortex (BA9) (UBERON:0009834':'brain',
    'Brain - Cortex (UBERON:0001870)':'brain',
    'Brain - Cerebellum (UBERON:0002037)':'brain',
    'Brain - Anterior cingulate cortex (BA24) (UBERON:0009835)':'brain',
    'Brain (UBERON:0000955)':'brain',
    'Bone - Bone Element (UBERON:0001474)':'bone',
    'Bladder (UBERON:0001255)':'bladder',
    'Artery - Coronary (UBERON:0001621)':'coronary',
    'Adrenal Gland (UBERON:0002369)':'adrenal gland',
    'Adipose - Subcutaneous (UBERON:0002190)':'adipose',
    'Ovary (UBERON:0000992)':'ovary',
}

tumordata['organs'] = tumordata['site'].map(site_to_organ).fillna('Fail')
biospecdata['biospeclocation'] = biospecdata['tissue_site'].map(tissue_to_site).fillna('Fail')
biospecdata['organs']=biospecdata['site'].map(site_to_organ).fillna('Fail')
#DICTONARIES
#these two dictonarys hold the organs and a possible color number for the different organs in the study
#Note-there are more organs in the package
allfemaleorgans = {
    'heart': 3, 'lung': 2, 'brain': 18, 'colon': 1, 'liver': 4, 'breast': 5, 
    'stomach': 6, 'esophagus': 7, 'ovary': 8, 'uterus': 9, 'pancreas': 10, 
    'lymph node': 11, 'gall bladder': 12, 'thyroid gland': 13, 'kidney': 14,
    'uterine cervix': 15, 'bronchus': 16, 'skin': 17
}
allmaleorgans={
    'heart': 3, 'lung': 2, 'brain': 18, 'colon': 1, 'liver': 4, 'breast': 5, 
    'stomach': 6, 'esophagus': 7, 'pancreas': 8, 
    'lymph node': 11, 'gall bladder': 12, 'thyroid gland': 13, 'kidney': 14,
    'uterine cervix': 15, 'bronchus': 16, 'skin': 10, 'prostate gland':9
}
#These two dictonarys hold a textbubble and arrow location for each gland
annotationsfemale = {
    'lung': ((0, 0), (48, 44)),
    'heart': ((0, 10), (52, 48)),
    'brain': ((25, 0), (50, 5)),
    'colon': ((0,95), (43, 85)),
    'liver': ((0, 53), (43, 65)),
    'breast': ((100, 30), (62, 50)),
    'pancreas': ((100, 62), (51, 70)),
    'esophagus': ((80, 20), (51, 40)),
    'uterus': ((95, 78), (50.5, 88.5)),
    'ovary': ((90, 88), (58, 90.75)),
    'stomach': ((93, 50), (62, 64)),
    'lymph node': ((75, 0), (56, 27)),
    'gall bladder':((0,62),(49,70)),
    'thyroid gland':((22,10),(52,31)),
    'kidney':((0,85),(45,72)),
    'skin':((0,100),(40,110))
}

annotationsmale = {
    'lung': ((0, 0), (48, 44)),
    'heart': ((0, 10), (52, 48)),
    'brain': ((25, 0), (50, 5)),
    'colon': ((0,95), (43, 85)),
    'liver': ((0, 53), (43, 65)),
    'breast': ((100, 30), (62, 50)),
    'pancreas': ((100, 62), (51, 70)),
    'esophagus': ((80, 20), (51, 40)),
    'stomach': ((93, 50), (62, 64)),
    'lymph node': ((75, 0), (56, 27)),
    'gall bladder':((0,62),(49,70)),
    'thyroid gland':((22,10),(52,31)),
    'kidney':((0,85),(45,72)),
    'skin':((0,100),(40,110)),
    'prostate gland':((0,50),(0,10))
}
#These two dictonaries have just a dot, no arrows, best for plotly interactive steps since it is not a tuple
annotations = {
    'thyroid gland': (299, 240),
    'esophagus': (308, 275),
    'lymph node': (325, 222),
    'lung': (285, 290),
    'left heart': (300, 304),
    'atrial heart':(305,304),
    'colon': (270, 465),
    'trans colon':(270, 460),
    'sigmoid colon':(275,480),
    'liver': (266, 375),
    'breast': (345, 315),
    'stomach': (315, 390),
    'uterus': (305, 470),
    'ovary': (332, 480),
    'pancreas': (347, 372),
    'gall bladder': (290, 395),
    'kidney cortex': (273, 405),
    'kidney med': (265,400),
    'skin': (255, 553),
    'other':(0,0),
    'muscle':(260, 560),
    'pleura': (285, 290),
    'peritoneum': (285, 290),
    'pericardium':(310, 295),
    'blood':(440,385),
    'cavity':(300,335),
    'spleen':(100,100),
    'spinal cord':(305, 350),
    'small intestine':(290,450),
    'diaphragm':(310,325),
    'chest':(325,290),
    'bone':(100,100),
    'bladder':(305,495),
    'coronary':(307,299),
    'adrenal gland':(270,400),
    'adipose':(100,100),
    'prostate gland':(304,495),
}


#this creates the original background image before a person has chosen a patient, defaults to female but could simply change to male if that is a need
fig,ax=plt.subplots(figsize=(6,10))

scatter_fig = go.Figure()


anatomogram2 = pyanatomogram.Anatomogram('homo_sapiens.female')
anatomogram2.highlight_tissues(allfemaleorgans, cmap='gist_heat_r')
anatomogram2.to_matplotlib(ax=ax)



for organ, ((x_text, y_text), (x_arrow, y_arrow)) in annotationsfemale.items():
    ax.annotate(
        organ.capitalize(),
        xy=(x_arrow, y_arrow),
        xytext=(x_text, y_text),
        arrowprops=dict(facecolor='black', arrowstyle='->'),
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='black', lw=1),
        fontsize=8
    )
    


buf = io.BytesIO()
plt.tight_layout()  
plt.savefig(buf, format="png")
plt.close(fig)
buf.seek(0)

img = Image.open(buf)
w, h = img.size

buffer = io.BytesIO()
img.save(buffer, format="PNG")
encoded_image = base64.b64encode(buffer.getvalue()).decode()

scatter_fig = go.Figure()
scatter_fig.add_layout_image(
    dict(
        source=f"data:image/png;base64,{encoded_image}",
        xref='x',
        yref='y',
        x=0,
        y=h,
        sizex=w,
        sizey=h,
        sizing="stretch",
        layer="below"
    )
)
for organ, (x, y) in annotations.items():
    scatter_fig.add_trace(go.Scatter(
        x=[x],
        y=[h - y],
        mode="markers",
        marker=dict(size=5, color="black"),
        name=organ,
        text=organ,
        hoverinfo="text",
        customdata=[organ]
    ))

scatter_fig.update_layout(
    showlegend=False,
    xaxis=dict(visible=False, range=[0, w]),
    yaxis=dict(visible=False, range=[0, h]),
    margin=dict(l=0, r=0, t=0, b=0),
    clickmode='event+select'
)

#Dash App Creationg occurs here



register_page(__name__, path="/patients")

layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.H2("Patient-Specific Data", className="text-center text-primary mb-4"))
    ]),
    dbc.Row([
        dbc.Col(html.H2("Select a Patient to View Information", className="text-secondary mb-3"))
    ]),
    dbc.Row([
        dbc.Col(
            dcc.Dropdown(
                id='patient-dropdown',
                options=[{'label': pid, 'value': pid} for pid in sorted(tumordata['patient'].unique())],
                placeholder="Select a patient",
                style={'backgroundColor': 'transparent','color': 'mediumslateblue', 'border': '1px solid white',}
        ),
        width=6,
        style={'backgroundColor': 'transparent'}
        )
    ]),
    dbc.Row([
        dbc.Col(dcc.Graph(id='swimmer-plot'), width=12)
    ], className="mb-4"),

    dbc.Row([
        dbc.Col(
            dash_table.DataTable(
                id='tableinfo',
                columns=[],
                data=[],
                style_table={'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'minWidth': '100px', 'width': '150px', 'maxWidth': '200px'},
                page_size=10
            ),
            width=12,
            style={'backgroundColor': 'transparent'}
        )
    ], className="mb-4"),

    dbc.Row([
        dbc.Col(
            dcc.Graph(id='anatomy-plot', figure=scatter_fig),
            width=6,
            style={'backgroundColor': 'transparent'}
        ),
        dbc.Col(
            dcc.Graph(id='organ-data-plot'),
            width=6,
            style={'backgroundColor': 'transparent'}
        )
    ])
], fluid=True,  style={'backgroundColor': 'transparent'})




@callback(
    Output('anatomy-plot', 'figure'),
    Input('patient-dropdown', 'value')
)
def update_anatomy_figure(selected_patient):
    if not selected_patient:
        return make_figure_with_background(encoded_image)

    buf = io.BytesIO()

    patient_info=tumordata[tumordata['patient']==selected_patient]
    if patient_info.empty:
        return make_figure_with_background(encoded_image)
    gender = patient_info.iloc[0]['gender'].lower()
    
    fig, ax = plt.subplots(figsize=(6, 10))
    
    if gender == 'female':
        anatomogram = pyanatomogram.Anatomogram('homo_sapiens.female')
        patient_organs = tumordata[tumordata['patient'] == selected_patient]['organs'].unique().tolist()
        highlight_organs = {organ: allfemaleorgans[organ] for organ in patient_organs if organ in allfemaleorgans}
        anatomogram.highlight_tissues(highlight_organs, cmap='twilight')
        anatomogram.to_matplotlib(ax=ax)
    else:
        anatomogram = pyanatomogram.Anatomogram('homo_sapiens.male')
        patient_organs = tumordata[tumordata['patient'] == selected_patient]['organs'].unique().tolist()
        highlight_organs = {organ: allmaleorgans[organ] for organ in patient_organs if organ in allmaleorgans}
        anatomogram.highlight_tissues(highlight_organs, cmap='twilight')
        anatomogram.to_matplotlib(ax=ax)
    
 
    plt.savefig(buf, format='png')
    plt.close(fig)
    buf.seek(0)
    img = Image.open(buf)
    w, h = img.size
    img_str = base64.b64encode(buf.getvalue()).decode()

    
    fig_out = go.Figure()

    patient_tissue = biospecdata[biospecdata['patient'] == selected_patient]['biospeclocation'].tolist()
    fig_out.add_layout_image(
        dict(
            source=f"data:image/png;base64,{img_str}",
            xref="x", yref="y",
            x=0, y=h,
            sizex=w, sizey=h,
            sizing="stretch",
            layer="below"
        )
    )
    organ_counts = {}
    for organ in patient_tissue:
        organ_counts[organ] = organ_counts.get(organ, 0) + 1
    
    unique_organs = list(organ_counts.keys())
    organ_colors = {organ: color for organ, color in zip(unique_organs, pc.qualitative.Plotly)}
  
    
    for organ, (x, y_pixel) in annotations.items():
        if organ in organ_counts:
            count = organ_counts[organ]
            fig_out.add_trace(go.Scatter(
                x=[x],
                y=[h - y_pixel],
                mode='markers',
                marker=dict(size=8 + count * 2, color=organ_colors.get(organ, 'lightpink'), opacity=0.6),
                name=organ,
                text=organ,
                hoverinfo='text',
                customdata=[organ]
            ))
  
        
        


    fig_out.update_layout(
        width=w,
        height=h,
        xaxis=dict(visible=False, range=[0, w]),
        yaxis=dict(visible=False, range=[0, h]),
        margin=dict(l=0, r=0, t=0, b=0),
        dragmode=False,
        clickmode='event+select',
        showlegend=False,
    )

    return fig_out


@callback(
    Output('organ-data-plot', 'figure'),
    Input('patient-dropdown', 'value')
)
def display_organ_data(selected_patient):
   
    patient_data = tumordata[tumordata['patient'] == selected_patient]

    if patient_data.empty:
        return go.Figure(layout=go.Layout(title=f"No data for patient {selected_patient}"))


    patient_organs = patient_data['organs'].dropna().unique()


    cohort_data = tumordata[tumordata['organs'].isin(patient_organs)].copy()
    patient_data = patient_data[patient_data['organs'].isin(patient_organs)]
    patientbiodata=biospecdata[biospecdata['patient']==selected_patient]
    cohortbiodata=biospecdata[biospecdata['organs'].isin(patient_organs)].copy()
    
    for df in [patient_data, cohort_data]:
        df['biocount'] = pd.to_numeric(df['biocount'], errors='coerce')
        df['treatcount'] = pd.to_numeric(df['treatcount'], errors='coerce')

    patient_data = patient_data.dropna(subset=['biocount', 'treatcount'])
    cohort_data = cohort_data.dropna(subset=['biocount', 'treatcount'])


    if patient_data.empty:
        return go.Figure(layout=go.Layout(title=f"No valid tumor data for patient {selected_patient}"))

    tumorspec = patient_data['tum_nor'].value_counts()
    tumor_counts = cohort_data['tum_nor'].value_counts()
    analytespec = patient_data['analyte'].value_counts()
    analyte_counts = cohort_data['analyte'].value_counts()
    experspec = patient_data['strategy'].value_counts()
    exper_counts = cohort_data['strategy'].value_counts()
    statespec=patientbiodata['tumorstate'].value_counts()
    state_counts=cohortbiodata['tumorstate'].value_counts()

    if tumorspec.empty and tumor_counts.empty:
        return go.Figure(layout=go.Layout(title=f"No tumor categories found for patient {selected_patient}"))

   
    fig = make_subplots(
        rows=4, cols=2,
        subplot_titles=[
            "Patient Tumor vs Normal", "Cohort Tumor vs Normal",
            "Patient Analyte Types", "Cohort Analyte Types",
            "Patient Strategies", "Cohort Strategies",
            "Patient Tumor States", "Cohort Tumor States"
        ],
        specs=[[{"type": "domain"}, {"type": "domain"}],
               [{"type": "domain"}, {"type": "domain"}],
               [{"type": "domain"}, {"type": "domain"}],
               [{"type": "domain"}, {"type": "domain"}]]
    )


    def add_pie_trace(data, title, row, col, colors=None):
        if not data.empty:
            fig.add_trace(go.Pie(
                labels=data.index,
                values=data.values,
                name=title,
                hovertemplate="%{label}: %{value} (%{percent})<extra></extra>",
                marker=dict(colors=colors) if colors else None
            ), row=row, col=col)

    add_pie_trace(tumorspec, 'Patient Tumor vs Normal', row=1, col=1, colors=['lavender', 'mediumslateblue'])
    add_pie_trace(tumor_counts, 'Cohort Tumor vs Normal', row=1, col=2, colors=['lavender', 'mediumslateblue'])
    add_pie_trace(analytespec, 'Patient Analyte Types', row=2, col=1, colors=['lightblue', 'blue'])
    add_pie_trace(analyte_counts, 'Cohort Analyte Types', row=2, col=2, colors=['lightblue', 'blue'])
    add_pie_trace(experspec, 'Patient Strategies', row=3, col=1, colors=['lightpink', 'palevioletred','orchid','hotpink'])
    add_pie_trace(exper_counts, 'Cohort Strategies', row=3, col=2, colors=['lightpink', 'palevioletred','orchid','hotpink'])
    add_pie_trace(statespec, "Patient Tumor States", row=4, col=1, colors=['coral', 'orangered','darkorange','tomato'])
    add_pie_trace(state_counts, "Cohort Tumor States", row=4, col=2, colors=['coral', 'orangered','darkorange','tomato'])

    add_legend_annotation(fig, tumorspec.index, ['lavender', 'mediumslateblue'], row=1, col=1)
    add_legend_annotation(fig, tumor_counts.index, ['lavender', 'mediumslateblue'], row=1, col=2)
    add_legend_annotation(fig, analytespec.index, ['lightblue', 'blue'], row=2, col=1)
    add_legend_annotation(fig, analyte_counts.index, ['lightblue', 'blue'], row=2, col=2)
    add_legend_annotation(fig, experspec.index, ['lightpink', 'palevioletred', 'orchid', 'hotpink'], row=3, col=1)
    add_legend_annotation(fig, exper_counts.index, ['lightpink', 'palevioletred', 'orchid', 'hotpink'], row=3, col=2)
    add_legend_annotation(fig, statespec.index, ['coral', 'orangered','darkorange','tomato'], row=4, col=1)
    add_legend_annotation(fig, state_counts.index, ['coral', 'orangered','darkorange','tomato'], row=4, col=2)
            
    
    organ_list = ', '.join([str(org).capitalize() for org in patient_organs])
    fig.update_layout(
        title_text=f"Tumor Data for Patient {selected_patient} and {organ_list} Cancer",
        height=750,
        width=1000,
        showlegend=False,
        paper_bgcolor='rgba(0,0,0,0)', 
        plot_bgcolor='rgba(0,0,0,0)',
        font=dict(color='black')
    )

    return fig



    # Fallback (default empty plot)
    return go.Figure(layout=go.Layout(title="Click on an organ to view tumor information"))

@callback(
    Output('swimmer-plot', 'figure'),
    Input('patient-dropdown', 'value')
)
def update_swimmer_figure(selected_patient):
    if not selected_patient:
        return go.Figure(layout=go.Layout(title="Select a patient to view treatment timeline"))

    plt.rcParams.update({'font.size': 14})
    plt.rcParams['svg.fonttype'] = 'none'
    pio.templates.default = "seaborn"

    fig = go.Figure()

    patient_data = tumordata[tumordata['patient'] == selected_patient]
    if patient_data.empty:
        return go.Figure(layout=go.Layout(title="No data found for selected patient"))

    unique_drugs = sorted(patient_data['drug'].dropna().astype(str).unique()) 
    color_list = px.colors.qualitative.Alphabet
    while len(color_list) < len(unique_drugs):
        color_list += color_list  

    drug_color_map = {drug: color_list[i % len(color_list)] for i, drug in enumerate(unique_drugs)}

    added_drug_legends = set()

    for idx, row in patient_data.iterrows():
        drug = row['drug']
        color = drug_color_map.get(drug, 'lightgray')
        show_legend = drug not in added_drug_legends
        if show_legend:
            added_drug_legends.add(drug)

        fig.add_trace(go.Scatter(
            x=[row['tx_start'], row['tx_end']],
            y=[row['patient'], row['patient']],
            mode='lines',
            line=dict(color=color, width=10),
            hoverinfo='text',
            name=drug,
            text=f"Patient: {row['patient']}<br>Drug: {row['drug']}<br>Start: {row['tx_start']}<br>End: {row['tx_end']}<br>Stop: {row['stop']}, <br>Morph: {row['morph']}",
            showlegend=show_legend
        ))

    stop_marker_styles = {
        'Progression and relapse': dict(color='red', symbol='circle', size=10),
        'Completion of standard course': dict(color='black', symbol='triangle-down', size=10),
        'Toxicity': dict(color='green', symbol='triangle-down', size=10),
        'Other': dict(color='gray', symbol='circle', size=10)
    }

    for idx, row in patient_data.iterrows():
        stop = row['stop']
        if stop in stop_marker_styles:
            style = stop_marker_styles[stop]
            fig.add_trace(go.Scatter(
                x=[row['tx_end']],
                y=[row['patient']],
                mode='markers',
                marker=dict(color=style['color'], symbol=style['symbol'], size=style['size']),
                name=stop,
                hoverinfo='text',
                text=f"Stop Reason: {row['stop']}<br>Morph: {row['morph']}",
                showlegend=False
            ))

    morph_options = patient_data['morph'].dropna().unique()
    buttons = []
    for morph in morph_options:
        visibility = [morph.lower() in str(trace.text).lower() for trace in fig.data]
        buttons.append(dict(
            label=morph,
            method='update',
            args=[{'visible': visibility},
                  {'title': f'Drug Progression in {selected_patient} for {morph}'}]
        ))

    fig.update_layout(
        title=dict(
            text=f'Drug Progression in {selected_patient}',
            x=0.5,
            xanchor='center',
            font=dict(size=20)
        ),
        xaxis_title='Time (Days)',
        yaxis_title='Patient',
        yaxis=dict(type='category'),
        legend_title='Drug / Stop Reason',
        height=400,
        hovermode='closest',
        updatemenus=[dict(
            active=0,
            buttons=buttons,
            direction="down",
            showactive=True,
            x=0.0,
            y=1,
            xanchor="left",
            yanchor="top"
        )] if buttons else []
    )

    return fig


@callback(
    Output('tableinfo', 'columns'),
    Output('tableinfo', 'data'),
    Input('patient-dropdown', 'value')
)
def createtable(selected_patient):
    if selected_patient is None:
        return [],[]

  
    patient_data = biospecdata[biospecdata['patient'] == selected_patient]
    if patient_data.empty:
        return [],[]


    displaycolumns=['patient', 'ogmaterial', 'submaterial', 'collection', 'tissue_site', 'tissuedetails']


    patient_data=patient_data[displaycolumns]
    custom_column={
        'patient':'Patient ID',
        'ogmaterial':'Orignal Material Type',
        'submaterial':'Submitted Material Type',
        'collection':'Collection Date',
        'tissue_site':'Sample Location',
        'tissuedetails':'Tissue Site Details'
    }
    columns=[{'name':custom_column[col],'id':col} for col in patient_data.columns]
    data=patient_data.to_dict('records')

    return columns,data


    
