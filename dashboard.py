import os
import numpy as np
import pandas
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from dash import Dash, dcc, html, Input, Output
import dash_table


"""

Interactive dashboard implemented in Dash for accessible evaluation of GWAS (Genome Wide Association Study) outputs.
GWAS is a widely used genomic technique which applies statistical models to evaluate the strength of associations
between the variation in phenotype (a biological trait) and the genotype (genetic variants).
In other words, GWAS is used to map genetic positions which are associated with complex traits.

The dashboard bridges GEMMA (GWAS pipeline, https://github.com/genetics-statistics/GEMMA) and Manhattan Harvester 
(GWAS peak calling tool, https://genomics.ut.ee/en/tools) to visualise the association peaks in interactive way, 
with additional feature to include parameters to filter the noise and yield strict and objective peak calling. 

"""



# path to the directory with har_processed files
pathway = "harv_processed/"

# names of chromosomes  for the studied species
chromosomes_list = ("Chromosome 1", "Chromosome 2", "Chromosome 3", "Chromosome 4", "Chromosome 5")


def read_harv_processed(file):
    """
        Inputs:
    file: harv_processed file from a given pathway.

        Outputs:
    averages: the list generated from the first line of harv_processed
    file (average LODs of top 5-10 peaks).

    df (pandas data frame) containing the GEMMA (GWAS pipeline) and harvester
    combined information about SNPs and parameters in their allocated peaks.
    """
    with open(file) as file:
        line = file.readline().strip()
        averages = line.split("\t")
        df = pandas.read_csv(file, sep="\t")
        return df, averages



# make list of all harv_processed files in the given directory
files = [f for f in os.listdir(pathway) if os.path.isfile(os.path.join(pathway, f))  and "harv_processed" in f]


# app layout written with Dash, renders in HTML.
app = Dash(__name__)
app.layout = html.Div([

    html.H1("Semi-Natural Intelligence Peak calling (SNIP calling)", style = {'text-align':'center', "font-family":"Arial"}),
    html.H2("Interactive dashboard for visualising GWAS output", style = {'text-align':'center', "font-family":"Arial"}),


        html.Br(),
        html.Div(className="row", children = [

        html.Div([html.Label("Select file")], style = {"width":'48%', "font-family":"Arial", 'margin-left': "14%"}) ,
           html.Div([html.Label("Peaks")], style = {"width":'25%', "font-family":"Arial", 'margin-left': "-22.5%"}) ,
           html.Div([html.Label("Factor")],style = {"width":'6%', "font-family":"Arial", 'margin-left': "-16.15%"})  ,
           html.Div([html.Label("Min SNPs")], style = {"width":'6%', "font-family":"Arial", 'margin-left': "3.6%"}) ,
           html.Div([html.Label("Max spacing")], style = {"width":'6%', "font-family":"Arial", 'margin-left': "3.5%"}),
           html.Div([html.Label("Min vbal")], style = {"width":'6%', "font-family":"Arial", 'margin-left': "3.4%"})
        ], style = {"display" : "flex", "margin-top" : "1%" }),

        html.Div( className="row",children = [

        # html.Button('Next', id='next_button'),                          ###### to add

        html.Div(dcc.Dropdown(id='filename',
                     options=[
                         {'label': i, 'value': i} for i in files
                     ],
                     multi=False,
                     style = {"width":"80%", "margin-left":"0%"},
                     value = files[0]

                     ), style = {"width":"35%", "margin-left":"14%"}),

        html.Div(dcc.Dropdown(id="Top avg",
                     options=[
                              {"label":"5","value":0},
                              {"label":"6","value":1},
                              {"label":"7","value":2},
                              {"label":"8","value":3},
                              {"label":"9","value":4},
                              {"label":"10","value":5}
                              ],

                     multi = False,
                     value= 0,
                     ), style = {"width":"6%", "margin-left":"-9.7%", "margin-right":"-0%"}),

        dcc.Input(id="factor",
                  type = "number",
                  style = {"width":"6%",'margin-left': "3%"},
                  value= 1.33
                     ),

        dcc.Input(id="peak_count",
                  type = "number",
                  style={"width":"6%",'margin-left': "3%"},
                  # size = "2",
                  value = 50,
                  ),

        dcc.Input(id="max_spacing",
                  type = "number",
                  style={"width":"6%",'margin-left': "3%"},
                  value = 20000,
                  ),
        dcc.Input(id="min_vbal",
                  type = "number",
                  style={"width":"6%",'margin-left': "3%"},
                  value = 0.1,
                  )
        ], style = dict(display = "flex",margin_bottom = "-10%")),


    dcc.Graph(id='man_plot', figure = {}, style={'width': '100%', 'height': '75vh', "margin_top":"-70%"}),

    html.Br(),
    html.H3("Parameters of called SNPs",
            style = {'text-align':'center', "font-family":"Arial"}),

    dash_table.DataTable(id="called_table",
                         style_cell={
                             'overflow': 'hidden',
                             'textOverflow': 'ellipsis',
                             'maxWidth': 0,
                              "font-family":"Arial",
                              'text-align':'center'
                         },
                         style_table={"width":"95%",
                                      "margin-left":"2.5%"}
                         ),
    html.Div(id= "output_container", children =[]),
])

def handle_inputs(input, new_type = int, default_value = 0):
    """
    function which handles the app dcc.input values and returns them in desired type/format.
    Handles dcc.input of type number. If input is not given, default value is returned.
    """
    if input is not None and len(str(input)) > 0:
         return new_type(input)
    else:
        return default_value

# ------------------------------------------------------------------------------
# App callback that connects the graph/table with input Dash components
@app.callback(
    [Output(component_id='output_container', component_property='children'),
     Output(component_id='man_plot', component_property='figure'),
     Output(component_id="called_table",component_property="data"),
     Output(component_id="called_table", component_property = "columns")
     ],
    [Input(component_id="filename",component_property="value"),
    Input(component_id='Top avg', component_property='value'),
    Input(component_id='factor',component_property='value'),
    Input(component_id='peak_count',component_property='value'),
    Input(component_id='max_spacing',component_property='value'),
    Input(component_id='min_vbal',component_property='value')
    ])




def update_graph(filen,tavg,fact,peak_c,max_spac, min_vb):

    """
    Callback function which takes the app inputs and returns the app outputs.

        Inputs
    filen: selected harv_processed file showing GWAS results for the given splice-site.
    tavg: number of peaks used to compute average value in the noise threshold formula.
    fact: factor to multiply the tavg value and calculate the noise threshold   (noise = tag*fact)
    peak_c: minimal total number of SNPs with LOD over 3 in a called peak.
    max_spac: Maximum spacing (average distance between SNPs in a peak) in a called peak.
    min_vb: minimal vbal (vertical balance) in a called peak.

        Outputs
    container: Optional text to display the selected inputs.
    fig: interactie plotly figure for the Manhattan plot of the selected GWAS output.
    data: table with example output data having the called SNPs and their attributes.
    columns: columns for the table.
    """

    if type(filen) != list:
        df_trim,noise_borders = read_harv_processed(pathway +filen)
    else:
        df_trim,noise_borders = read_harv_processed(pathway + filen[-1])

    df_trim["LOD"] = np.log10(df_trim["p_wald"]) * (-1)

    bonferroni = np.log10( 0.05 / 10709466) * (-1)

    dff = df_trim.copy()  # store a copy to enable resetting manipulation of parameters

    fact = handle_inputs(fact,float)

    noise_border = np.log10(float(noise_borders[tavg])) * (-1) *fact

    peak_c = handle_inputs(peak_c,int)
    max_spac = handle_inputs(max_spac,int)
    min_vb = handle_inputs(min_vb,float)


    fig = make_subplots(rows = 1,
                        cols = 5,
                        shared_yaxes = True,
                        horizontal_spacing = 0.01,
                        subplot_titles=chromosomes_list)

    table_list = []

    for chromosome in range(1,len(chromosomes_list) +1):
        df_chr = df_trim[df_trim["chr"] == chromosome]

        fig.add_trace(go.Scatter(x = df_chr["ps"],
                                 y = df_chr["LOD"],
                                 hovertemplate=
                                 "<br>pos: %{x}<br>" +
                                 "<br>LOD: %{y}<br>"
                                 ,
                                 name='not called',
                                 mode='markers',
                                 marker=dict(
                                     color='gray',
                                     size=5,
                                     )
                                 ),
                            row = 1,
                            col = chromosome
                        )

        dff = df_chr[df_chr["GQS"] != "None"]
        dff["GQS"] = dff["GQS"].astype(float)
        dff = dff[dff["GQS"] > 3]

        fig.add_trace(go.Scatter(x=dff["ps"],
                                 y=dff["LOD"],
                                 customdata = dff[["GQS", "spacing","count","monot","vbal1"]].values.tolist(),
                                 hovertemplate =
                                 "<br>pos: %{x}<br>" +
                                 "<br>LOD: %{y}<br>" +
                                 "<br>GQS: %{customdata[0]}<br>" +
                                 "<br>spac: %{customdata[1]}<br>" +
                                 "<br>count: %{customdata[2]}<br>" +
                                 "<br>monot: %{customdata[3]}<br>" +
                                 "<br>vbal1: %{customdata[4]}<br>",

                                 # color_discrete_map = COLOUR_MAP,

                                 mode='markers',
                                 name='harvester called',
                                 marker=dict(
                                     color=dff["GQS"],
                                     colorscale='Viridis_r',
                                     cmin=2,
                                     cmax=5,
                                     colorbar=dict(x=- 0.11,
                                                   title="GQS",
                                                   tickvals=[2, 3, 4, 5],
                                                   ticktext=["2", "3", "4", "5"]

                                                   ),
                                     size=5,
                                 )
                                 ),
                      row=1,
                      col=chromosome
                      )

        dff["count"] = dff["count"].astype(float)
        dff["spacing"] = dff["spacing"].astype(float)
        dff["vbal1"] = dff["vbal1"].astype(float)

        dff = dff[(dff["LOD"] > noise_border) &
                  (dff["LOD"] > bonferroni) &
                  (dff["count"] > peak_c) &
                  (dff["spacing"] < max_spac) &
                  (dff["vbal1"] > min_vb)]


        table_list.append(dff)
        fig.add_trace(go.Scatter(x= dff["ps"],
                                 y=dff["LOD"],

                                 customdata=dff[["GQS", "spacing","count","monot","vbal1"]].values.tolist(),
                                 hovertemplate=
                                 "<br>pos: %{x}<br>" +
                                 "<br>LOD: %{y}<br>" +
                                 "<br>GQS: %{customdata[0]}<br>" +
                                 "<br>spac: %{customdata[1]}<br>" +
                                 "<br>count: %{customdata[2]}<br>" +
                                 "<br>monot: %{customdata[3]}<br>" +
                                 "<br>vbal1: %{customdata[4]}<br>",

                                 mode='markers',
                                 name='called',
                                 marker=dict(
                                     color='red',
                                     size=5,)
                                 ),
                      row=1,
                      col=chromosome
                      )

        fig.add_trace(go.Scatter(x=df_chr["ps"], y = [noise_border for _ in range(len(df_chr))],
                                 name='noise',
                                 mode= "lines",
                                 line=dict(color='royalblue', width=2, dash='dash')),
                        row=1,
                        col=chromosome)

        fig.add_trace(go.Scatter(x=df_chr["ps"], y=[bonferroni for _ in range(len(df_chr))],
                                 mode = "lines",
                                 name='bonferroni',
                                 line=dict(color='red', width=2, dash='dash')),
                      row=1,
                      col=chromosome)


    fig.update_layout(showlegend=False, yaxis_title="-log10(p-value)",)
    container = ""
    data = pandas.concat(table_list)
    columns = [{'name': col, 'id': col} for col in data.columns]
    data = data.to_dict(orient='records')

    return  container, fig,data,columns

if __name__ == '__main__':
    app.run_server(debug=True)
