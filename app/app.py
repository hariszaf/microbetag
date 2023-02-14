# import dash
# from dash import dcc
# from dash import html

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# app.layout = html.Div(children=[
#   html.H1(children='Hello Dash'),

#   html.Div(children='''
#     Dash: A web application framework for Python.
#   '''),

#   dcc.Graph(
#     id='example-graph',
#     figure={
#       'data': [
#         {'x': [1, 2, 3], 'y': [4, 1, 2], 'type': 'bar', 'name': 'SF'},
#         {'x': [1, 2, 3], 'y': [2, 4, 5], 'type': 'bar', 'name': u'Montréal'},
#       ],
#       'layout': {
#         'title': 'Dash Data Visualization'
#       }
#     }
#   )
# ])

# if __name__ == '__main__':
#   app.run_server(host='0.0.0.0', port=8050)



# EXAMPLE WHERE YOU HAVE A DICTIONARY ABOVE YOUR GRAPH WITH THE NODES METADATA
import json
import dash
import dash_cytoscape as cyto
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output

app = dash.Dash(__name__)

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}


nodes = [
    {
        'data': {'id': short, 'label': label},
        'helloFriend': {'x': 20*lat, 'y': -20*long}
    }
    for short, label, long, lat in (
        ('la', 'Los Angeles', 34.03, -118.25),
        ('nyc', 'New York', 40.71, -74),
        ('to', 'Toronto', 43.65, -79.38),
        ('mtl', 'Montreal', 45.50, -73.57),
        ('van', 'Vancouver', 49.28, -123.12),
        ('chi', 'Chicago', 41.88, -87.63),
        ('bos', 'Boston', 42.36, -71.06),
        ('hou', 'Houston', 29.76, -95.37)
    )
]

edges = [
    {'data': {'source': source, 'target': target}}
    for source, target in (
        ('van', 'la'),
        ('la', 'chi'),
        ('hou', 'chi'),
        ('to', 'mtl'),
        ('mtl', 'bos'),
        ('nyc', 'bos'),
        ('to', 'hou'),
        ('to', 'nyc'),
        ('la', 'nyc'),
        ('nyc', 'bos')
    )
]


default_stylesheet = [
    {
        'selector': 'node',
        'style': {
            'background-color': '#BFD7B5',
            'label': 'data(label)'
        }
    }
]


app.layout = html.Div([
    cyto.Cytoscape(
        id='cytoscape-event-callbacks-1',
        layout={'name': 'preset'},
        elements=edges+nodes,
        stylesheet=default_stylesheet,
        style={'width': '100%', 'height': '450px'}
    ),
    html.Pre(id='cytoscape-tapNodeData-json', style=styles['pre'])
])


@app.callback(Output('cytoscape-tapNodeData-json', 'children'),
              Input('cytoscape-event-callbacks-1', 'tapNodeData'))
def displayTapNodeData(data):
    return json.dumps(data, indent=2)


if __name__ == '__main__':
  app.run_server(host='0.0.0.0', port=8050)

