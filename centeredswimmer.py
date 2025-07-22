import dash
from dash import dcc, html, Output, Input
import plotly.graph_objects as go
import pandas as pd
import math
import plotly.express as px
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
    followupdate=df['follow_up_date'].tolist()
    vital_status=df['vital_status'].tolist()
    
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
               'vitalstatus':vital_status
    }
    return new_dict


treats = pd.DataFrame(treatment_process_centered('Juric_Rapid_Autopsy_MASTER-treatments1.csv', 'Juric_Rapid_Autopsy_MASTER-participants.txt', 'Juric_Rapid_Autopsy_MASTER-biospecimens.txt'))
treats['tx_end'] = pd.to_numeric(treats['tx_end'], errors='coerce')
treats['tx_start'] = pd.to_numeric(treats['tx_start'], errors='coerce')

app = dash.Dash(__name__)

all_morphs = sorted(treats['morph'].dropna().unique())

app.layout = html.Div([
    dcc.Dropdown(
        id='cohort-dropdown',
        options=[{'label': morph, 'value': morph} for morph in all_morphs],
        placeholder='Select Cohort (Morphology)'
    ),
    dcc.Dropdown(
        id='drug-dropdown',
        placeholder='Select Drug',
    ),
    dcc.Graph(id='timeline-plot-with'),
    dcc.Graph(id='timeline-plot-without')
])

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

def generate_patientYpos(patient_DF, selected_drug):
    all_patients = patient_DF['patient'].unique()
    patients_with_drug = patient_DF[patient_DF['drug'] == selected_drug]['patient'].unique()
    ordered_patients = list(patients_with_drug) + [p for p in all_patients if p not in patients_with_drug]
    y_pos = {patient: len(ordered_patients) - i for i, patient in enumerate(ordered_patients)}

    return y_pos

stop_marker_styles = {
        'Progression and relapse': dict(color='red', symbol='circle', size=10),
        'Completion of standard course': dict(color='black', symbol='triangle-down', size=10),
        'Toxicity': dict(color='green', symbol='triangle-down', size=10),
        'Other': dict(color='gray', symbol='circle', size=10)
    }
 
@app.callback(
    Output('drug-dropdown', 'options'),
    Output('drug-dropdown', 'value'),
    Input('cohort-dropdown', 'value')
)
def update_drug_options(selected_morph):
    if selected_morph is None:
        return [], None
    filtered = treats[treats['morph'] == selected_morph]
    drugs = sorted(filtered['drug'].dropna().unique())
    options = [{'label': d, 'value': d} for d in drugs]
    # default select first drug
    value = drugs[0] if drugs else None
    return options, value


@app.callback(
    Output('timeline-plot-with', 'figure'),
    Output('timeline-plot-without', 'figure'),
    Input('cohort-dropdown', 'value'),
    Input('drug-dropdown', 'value')
)
def update_figure(selected_morph, selected_drug):
    if selected_morph is None or selected_drug is None:
        return go.Figure(), go.Figure()
    
    filtered = treats[treats['morph'] == selected_morph].copy()
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

    fig_with = create_figure(with_drug, f'WITH {selected_drug}')
    fig_without = create_figure(without_drug, f'WITHOUT {selected_drug}')

    return fig_with, fig_without



if __name__ == '__main__':
    app.run(debug=True)


