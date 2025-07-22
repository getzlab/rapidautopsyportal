import dash
import dash_bootstrap_components as dbc
from dash import Dash, html, dcc

app = Dash(__name__, use_pages=True, external_stylesheets=[dbc.themes.LUX], suppress_callback_exceptions=True)

navbar=dbc.NavbarSimple(
    brand="Rapid Autopsy Portal",
    brand_href='/',
    color='primary',
    children=[
        dbc.NavItem(dbc.NavLink("Home", href='/')),
        dbc.NavItem(dbc.NavLink("Patients", href="/patients")),
        dbc.NavItem(dbc.NavLink("Cohorts", href="/cohort"))
                    
    ]
)
app.layout = dbc.Container([
    navbar,
    html.Hr(),
    dash.page_container 
], fluid=True)

if __name__ == "__main__":
    app.run(debug=True)