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
from dash import register_page
from dash import callback, Input, Output, dcc

def treatment_process_centered(filename, filename2, filename3):

    df = pd.read_csv(filename, sep=',')
    df2=pd.read_csv(filename2, sep='\t')
    df3=pd.read_csv(filename3, sep='\t')
    df=pd.merge(df,df2, on='participant_id', how='inner',  suffixes=('', '_df2'))
    df=pd.merge(df,df3, on='participant_id', how='inner',  suffixes=('', '_df3'))
    patient_val = df['participant_id'].tolist()
    treatment_val = df['categories'].tolist()
    tx_start = df['start_date_dfd'].tolist()
    tx_end = df['stop_date_dfd'].tolist()
    stop_reason=df['stop_reason'].tolist()
    tumor_morph=df['tumor_morphology'].tolist()
    drug_type=df['drugs'].tolist()
    bio_count = pd.to_numeric(df['biospecimen_total_count'], errors='coerce').tolist()
    treat_count = pd.to_numeric(df['treatment_count'], errors='coerce').tolist()
    prim_site=df['tumor_primary_site'].tolist()
    gender_value=df['gender'].tolist()
    tumorornormal=df['tumor_normal'].tolist()
    analytetype=df['analyte_type'].tolist()
    experimentalstrategy=df['experimental_strategy'].tolist()
    ulpfrac=df['ulp_tumor_fraction'].tolist()
    followupdate=df['follow_up_date'].tolist()
    vital_status=df['vital_status'].tolist()
    categories=df['categories'].tolist()
    main_normal=df['main_normal'].tolist()
    race=df['race'].tolist()
    gender=df['gender'].tolist()
    stage=df['cancer_stage'].tolist()
    year=df['year_of_diagnosis']
    df = df.drop_duplicates()
    
    
    
    new_dict = {'patient': patient_val,
               'treatment': treatment_val,
               'tx_start': tx_start,
               'tx_end': tx_end,
               'stop':stop_reason,
               'morph':tumor_morph,
               'drug':drug_type,
               'biocount':bio_count,
               'treatcount':treat_count,
               'site':prim_site,
               'gender':gender_value,
               'tum_nor':tumorornormal,
               'analyte':analytetype,
               'strategy':experimentalstrategy,
               'ulp_tumor':ulpfrac,
               'followup':followupdate,
               'vitalstatus':vital_status,
               'categories':categories,
               'main_normal':main_normal,
               'race':race,
               'gender':gender,
               'stage':stage,
               'year':year
    }
    return new_dict

def treatment_process_centeredbio(filename, filename2, filename3):

    df = pd.read_csv(filename, sep=',')
    df2=pd.read_csv(filename2, sep='\t')
    df3=pd.read_csv(filename3, sep='\t')
    df=pd.merge(df,df2, on='participant_id', how='inner',  suffixes=('', '_df2'))
    df=pd.merge(df,df3, on='participant_id', how='inner',  suffixes=('', '_df3'))
    patient_val = df['participant_id'].tolist()
    treatment_val = df['categories'].tolist()
    tx_start = df['start_date_dfd'].tolist()
    tx_end = df['stop_date_dfd'].tolist()
    stop_reason=df['stop_reason'].tolist()
    tumor_morph=df['tumor_morphology'].tolist()
    drug_type=df['drugs'].tolist()
    bio_count = pd.to_numeric(df['biospecimen_total_count'], errors='coerce').tolist()
    treat_count = pd.to_numeric(df['treatment_count'], errors='coerce').tolist()
    prim_site=df['tumor_primary_site'].tolist()
    gender_value=df['gender'].tolist()
    tumorornormal=df['tumor_normal'].tolist()
    followupdate=df['follow_up_date'].tolist()
    vital_status=df['vital_status'].tolist()
    categories=df['categories'].tolist()
    tiss_site=df['tissue_site'].tolist()
    details=df['tissue_site_detail'].tolist()
    tumor_state=df['tumor_state'].tolist()
    df = df.drop_duplicates()
    
    new_dict = {'patient': patient_val,
               'treatment': treatment_val,
               'tx_start': tx_start,
               'tx_end': tx_end,
               'stop':stop_reason,
               'morph':tumor_morph,
               'drug':drug_type,
               'biocount':bio_count,
               'treatcount':treat_count,
               'site':prim_site,
               'gender':gender_value,
               'tum_nor':tumorornormal,
               'followup':followupdate,
               'vitalstatus':vital_status,
                'categories':categories,
                'tissue_site':tiss_site,
               'tissuedetails':details,
               'tumorstate':tumor_state,
    }
    return new_dict

def tableprocess(filename, filename1):
    df=pd.read_csv(filename1, sep='\t')
    df2 = pd.read_csv(filename, sep='\t')
    df=pd.merge(df,df2, on='participant_id', how='inner',  suffixes=('', '_df2'))
    patient_val = df['participant_id'].tolist()
    collectiondate=df['collection_date_dfd'].tolist()
    sub_material_type=df['submitted_material_type'].tolist()
    og_material_type=df['original_material_type'].tolist()
    prim_site=df['primary_site'].tolist()
    tiss_site=df['tissue_site'].tolist()
    details=df['tissue_site_detail'].tolist()
    tumor_state=df['tumor_state'].tolist()
    tumor_morph=df['tumor_morphology'].tolist()
    main_normal=df['main_normal'].tolist()
    
    new_dict = {'patient': patient_val,
               'site':prim_site,
               'collection':collectiondate,
               'submaterial':sub_material_type,
               'ogmaterial':og_material_type,
               'tissue_site':tiss_site,
               'tissuedetails':details,
               'tumorstate':tumor_state,
                'morph':tumor_morph,
                'main_normal':main_normal
               }
     
    return new_dict
    
def make_figure_with_background(base64_img):
    fig = go.Figure()
   
    fig.add_layout_image(
        dict(
            source=f"data:image/png;base64,{base64_img}",
            xref="x",
            yref="y",
            x=0,
            y=600, 
            sizex=400,
            sizey=600,
            sizing="stretch",
            layer="below"
        )
    )
    
    fig.update_layout(
        width=400,
        height=600,
        xaxis=dict(visible=False, range=[0, 400]),
        yaxis=dict(visible=False, range=[0, 600]),
        margin=dict(l=0, r=0, t=0, b=0),
        dragmode=False,
        clickmode='event+select',
        showlegend=False
    )
    
    return fig

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
