import escher,escher.urls,json,os
from IPython.display import HTML


def show_map(sol,map_loc,color=0):
    ''' Returns an escher Builder object for solution 'sol', map 'map_loc' and the supplied color scheme.
        sol:        the solution object containing the simulation results.
        map_loc:    filename of the map json
        color:      color scheme to use
    '''

    # read css
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    with open(script_dir+'/my_escher_stylesheet.css', 'r') as myfile:
        my_css=myfile.read().replace('\n', ' ')

    if color == 0: # grey to yellow to red
        colors = [{'type': 'min', 'color': '#cccccc', 'size': 5},
                  {'type': 'mean', 'color': '#FFFF00', 'size': 10},
                  {'type': 'max', 'color': '#FF0000', 'size': 20}]
    else:
        print('Color scheme not defined!')
        return

    if type(sol) != dict:
        try:
            d = sol.fluxes # it is a cobrapy solution object
        except:
            print('An empty solution was passed.')
            d = {}
    else:
        d = sol # shorthand

    d2 = {} # for some reason my types (from float 64 to float) were not updating and with a new dictionary they do
    for key in list(d.keys()): # remove output like this: 1.653e-15
        d2[key] = round(float(d[key]),6)


    network = escher.Builder(map_json=map_loc, reaction_data=d2, 
                       reaction_styles=['color', 'size', 'abs', 'text'],
                       # change the default colors, blue to purple to red
                       reaction_scale=colors,
                       hide_secondary_metabolites=False,secondary_metabolite_radius=10,
                       highlight_missing=True,embedded_css=my_css) 
    return network