#!/usr/bin/python env

# This script get a FlashWeave-like edgelist as input and returns 
# an object that can be loaded in DASH Cytoscape app 

import sys
import dash 
import dash_cytoscape as cyto  # pip install dash-cytoscape==0.2.0 or higher
from dash import html
from dash import dcc
from dash.dependencies import Output, Input



external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


# ------------------------------------------------------------------------------------------------


edge_list_file = open(sys.argv[1])

def edge_list_to_graph(edgelist):
   list_of_edges = []
   for line in edgelist: 

      if "# header" in line:

         node_names = (line.split("\t")[1]).split(",")
         node_names[-1] = node_names[-1][:-1] # remove newline symbol for the last node name

      if "#" not in line:

         line = line.split("\t")
         list_of_edges.append(line)

   # List of dictionaries
   potential_nodes = []
   counter = 0
   for node in node_names: 
      counter += 1

      new_node                  = {}
      new_node['data']          = {}
      new_node['data']['id']    = "Taxon_" + str(counter)
      new_node['data']['label'] = node
      potential_nodes.append(new_node)

   edges = []
   actual_nodes = set()
   for edge in list_of_edges:

      partner_a = edge[0]
      partner_b = edge[1]
      weight    = edge[2]

      if partner_a == partner_b: 
         print("Same source with target: ", partner_a_id)
         continue

      else:

         for node_candidate in potential_nodes:

            if node_candidate['data']['label'] == partner_a: 

               partner_a_id = node_candidate['data']['id']
            
            if node_candidate['data']['label'] == partner_b:

               partner_b_id = node_candidate['data']['id']

         new_edge         = {}
         new_edge['data'] = {}
         new_edge['data']['source'] = partner_a_id
         new_edge['data']['target'] = partner_b_id
         new_edge['data']['weight'] = weight
         new_edge['classes']: 'square'
         edges.append(new_edge)


         another_edge         = {}
         another_edge['data'] = {}
         another_edge['data']['source'] = partner_a_id
         another_edge['data']['target'] = partner_b_id
         another_edge['data']['weight'] = weight
         another_edge['classes']: 'diamond'
         edges.append(another_edge)




         actual_nodes.add(partner_a_id)
         actual_nodes.add(partner_b_id)

   nodes = []
   for node in potential_nodes: 
      if node['data']['id'] in actual_nodes:
         nodes.append(node)


   elements = nodes + edges

   return(elements)

my_elements = edge_list_to_graph(edge_list_file)
# ------------------------------------------------------------------------------------------------


# DASH CYTOSCAPE APP 


app.layout = html.Div(
   [
   html.Div(
      [
       # Here we build the dropdown menu with the options of how to show the graph
      dcc.Dropdown(
         id='dpdn',
         value='breadthfirst',
         clearable=False,
         options=[
                  {'label': name.capitalize(), 'value': name}
                  for name in ['random', 'circle', 'cose', 'concentric'] # 'breadthfirst' ,'grid',
         ]
        ),

      # Here is the actual graph
      cyto.Cytoscape(

                     id       = 'org-chart',   # 
                     minZoom  = 0.2,
                     maxZoom  = 2,
                     layout   = {
                                 'name': 'cose'
                                 },
                     elements = my_elements,
                     style     = {
                                    "height": "80em", 
                                    "width": "100%"                           
                                 },
                     stylesheet= [
                                    # Group selectors for NODES
                                    {
                                       'selector': 'node',
                                       'style': {
                                             'label': 'data(label)'
                                       }
                                    },

                                    #  Group selectors for EDGES - once this is on, the weight is displayed all the time in the middle of the edge
                                     {
                                         'selector': 'edge',
                                         'style': {
                                             'label': 'data(weight)'
                                         }
                                     },

                                     # Class selectors
                                     {
                                         'selector': '.purple',
                                         'style': {
                                             'background-color': 'purple',
                                             'line-color': 'purple'
                                         }
                                     },
                                    {
                                       'selector': '.square',
                                       'style': {
                                             'shape': 'square',
                                       }
                                    },
                                     {
                                         'selector': '.myimage',
                                         'style': {
                                             'width': 100,
                                             'height': 100,
                                             'background-image': ['./assets/sunny-and-cloud.jpg']
                                         }
                                     },
                                     {
                                         'selector': '.green',
                                         'style': {
                                             'background-color': 'green',
                                             'line-color': 'green'
                                         }
                                     },
                                    {
                                       'selector': '.diamond',
                                       'style': {
                                             'shape': 'diamond',
                                       }
                                    },
                                    # Conditional Styling
                                    # this weight class only exists within the EDGES
                                    {
                                       'selector': '[weight > 3]',
                                       'style': {
                                             'width': 20
                                       }
                                    },
                                    # *= means string contains...
                                    {
                                       'selector': '[label *= "rah"]',
                                       'style': {
                                             'background-color': '#000000',
                                       }
                                    }
                                 ]  # HERE IS WHERE STYLESHEET ENDS
               )
            ], className='eight columns'
            ),
            ], className='row'
            )



# 
@app.callback(
               Output('org-chart', 'layout'),
               Input('dpdn', 'value')
            )





# This function is about to work once the app is on 
def update_layout(layout_value):

    if layout_value == 'breadthfirst':
        return {
        'name': layout_value,
        'roots': '[id = "ed"]',
        'animate': True
        }

    else:
        return {
            'name': layout_value,
            'animate': True
        }





if __name__ == '__main__':
    app.run_server(debug=False, host="0.0.0.0", port=8080)

   



