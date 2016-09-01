
# coding: utf-8

# # Simulations reproducing the (Schuchmann, 2014) analysis on the C. ljungdhahlii GEMM

# In[15]:

import cobra
from cobra.flux_analysis import parsimonious
from cobra import Model, Reaction, Metabolite

import pandas as pd
pd.set_option('display.max_colwidth', -1)

import re # reg. expr.
import traceback

# For Escher visualization
import escher, escher.urls
import json
import os
from IPython.display import HTML

# define function that overlays simulation results on network map
def show_map(sol):
    import escher,escher.urls,json,os
    from IPython.display import HTML
    '''Takes 1 argument, the solution object containing the simulation results.'''
    for key in sol.x_dict.keys(): # remove output like this: 1.653e-15
        sol.x_dict[key] = round(sol.x_dict[key],3)
    network = escher.Builder(map_json='./escher/escher_map_c_ljungdahlii.json',reaction_data=sol.x_dict, 
                       reaction_styles=['color', 'size', 'abs', 'text'],
                       # change the default colors, blue to purple to red
                       reaction_scale=[{'type': 'min', 'color': '#cccccc', 'size': 2},
                                       {'type': 'mean', 'color': '#0000dd', 'size': 10},
                                       {'type': 'max', 'color': '#ff0000', 'size': 30}],
                       hide_secondary_metabolites=False,secondary_metabolite_radius=10,highlight_missing=True) 
    return network

M = cobra.io.read_sbml_model("./models/c_ljungdahlii_nagarajan_2013_update.xml")


# ## List of blocked reactions in the model

# In[16]:

model = M.copy()
d = {}
for rxn in model.reactions:
    if rxn.lower_bound == 0 and rxn.upper_bound == 0:
        d[rxn.id] = [rxn.name,rxn.reaction]
        
df = pd.DataFrame.from_dict(d)
df = df.transpose()
df.columns = ['Name','Reaction']
df


# ### The RNF reaction

# In[17]:

M.reactions.RNF.build_reaction_string


# ## Reproduce some growth numbers from Table 1 in (Nagarajan, 2013)
# This should return 0.034 and 0.06

# In[3]:

model = M.copy()
model.objective = model.reactions.BIOMASS_Cl_DSM_WT_46p666M1
model.reactions.ATPM.lower_bound = 0.45
model.reactions.EX_ac_e.lower_bound = 0 # do not force acetate flux
model.reactions.EX_co2_e.lower_bound = -10
model.reactions.EX_h2_e.lower_bound = -20
print model.optimize().f

model.reactions.EX_co2_e.lower_bound = 0
model.reactions.EX_h2_e.lower_bound = 0
model.reactions.EX_co_e.lower_bound = -20
print model.optimize().f


# ## Gibbs free energy from nothing checks
# Check that in no combination of enzyme alternatives there can be energy production from nothing

# In[4]:

results = [[],[],[]] # this will be filled with real results
for MTHFR in range(2):
    for MTHFD in range(2):
        for FDH in range(3):    
            for HYD in range(2):
                pfbaSol = []
                model = M.copy()
                
                # remove carbon and H2 from medium
                model.reactions.EX_co2_e.lower_bound = 0
                model.reactions.EX_h2_e.lower_bound = 0

                model.reactions.EX_ac_e.lower_bound = 0 # do not force acetate flux

                if FDH == 0:
                    pass # keep the default FDH7
                elif FDH == 1:
                    model.reactions.FDH7.upper_bound = 0; model.reactions.FDH7.lower_bound = 0; 
                    model.reactions.FDHH2.lower_bound = 0; model.reactions.FDHH2.upper_bound = 0
                    model.reactions.FDHFDNADPH.lower_bound = -1000; model.reactions.FDHFDNADPH.upper_bound = 1000
                else:
                    model.reactions.FDH7.upper_bound = 0; model.reactions.FDH7.lower_bound = 0; 
                    model.reactions.FDHH2.lower_bound = -1000; model.reactions.FDHH2.upper_bound = 1000
                    model.reactions.FDHFDNADPH.lower_bound = 0; model.reactions.FDHFDNADPH.upper_bound = 0 

                if HYD == 0:
                    model.reactions.HYDFDN.lower_bound = -1000; model.reactions.HYDFDN.upper_bound = 0 # reversed flux! The Fd + NADH hydrogenase
                    model.reactions.HYDFDN2r.lower_bound = 0; model.reactions.HYDFDN2r.upper_bound = 0 # the Fd + NADPH hydrogenase
                else:
                    model.reactions.HYDFDN.lower_bound = 0; model.reactions.HYDFDN.upper_bound = 0 # reversed flux! The Fd + NADH hydrogenase
                    model.reactions.HYDFDN2r.lower_bound = -1000; model.reactions.HYDFDN2r.upper_bound = 1000 # the Fd + NADPH hydrogenase
                    
                if MTHFD == 0:
                    model.reactions.MTHFD.upper_bound = 0; model.reactions.MTHFD_alt.upper_bound = 1000
                else:
                    model.reactions.MTHFD.upper_bound = 1000; model.reactions.MTHFD_alt.upper_bound = 0
                    
                if MTHFR == 0:    
                    model.reactions.MTHFR5.upper_bound = 0; model.reactions.MTHFR5_alt.upper_bound = 1000
                else: 
                    model.reactions.MTHFR5.upper_bound = 1000; model.reactions.MTHFR5_alt.upper_bound = 0
                try:
                    pfbaSol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
                    results[FDH].append(str(round(abs(pfbaSol.f),3)) + ' (' + str(round(pfbaSol.x_dict['FRNDPR2r_1'],2)) + ')' )
                    b = show_map(pfbaSol)
                    b.save_html('./predictions/fig_3_entry_row_'+str(FDH)+'_col_'+str(HYD)+'_MTHFD_'+str(MTHFD)+'_MTHFR_'+str(MTHFR)+'.html',overwrite=True)
                except:
                    results[FDH].append('NP')          

iterables = [['NADH', '2 NADH + Fd'], ['NADH', 'NADPH'],['Fd + NADH','Fd + NADPH']]
index = pd.MultiIndex.from_product(iterables, names=['MTHFR', 'MTHFD','FDH \ HYD'])
pd.DataFrame(results, columns=index,index=['Fd','FD+NADPH','H2'])


# ## Gibbs free energy from the full network: every alternative allowed
# this only works if the alternatives are made irreversible

# In[5]:

model = M.copy()
model.reactions.MTHFD.upper_bound = 1000
model.reactions.MTHFR5.upper_bound = 1000
model.reactions.HYDFDN.lower_bound = -1000
model.reactions.FDHFDNADPH.upper_bound = 1000
model.reactions.FDHH2.upper_bound = 1000
model.reactions.FDH7.lower_bound = 0

model.reactions.EX_co2_e.lower_bound = 0
model.reactions.EX_h2_e.lower_bound = 0
model.reactions.EX_ac_e.lower_bound = 0


# ## Attempt at gear-shifting

# In[6]:

model = M.copy()
model.reactions.MTHFD.lower_bound = -1000; model.reactions.MTHFD.upper_bound = 0 # defined in reverse
model.reactions.MTHFD_alt.lower_bound = 0; model.reactions.MTHFD_alt.upper_bound = 1000
model.reactions.MTHFR5.lower_bound = 0; model.reactions.MTHFR5.upper_bound = 1000
model.reactions.MTHFR5_alt.lower_bound = 0; model.reactions.MTHFR5_alt.upper_bound = 1000
model.reactions.HYDFDN2r.lower_bound = 0; model.reactions.HYDFDN2r.upper_bound = 1000
model.reactions.HYDFDN.lower_bound = -1000; model.reactions.HYDFDN.upper_bound = 0 # defined in reverse
model.reactions.FDHFDNADPH.upper_bound = 1000
model.reactions.FDHH2.upper_bound = 1000
model.reactions.FDH7.lower_bound = 0; model.reactions.FDH7.upper_bound = 1000

model.reactions.EX_co2_e.lower_bound = -2
model.reactions.EX_h2_e.lower_bound = -4
model.reactions.EX_ac_e.lower_bound = 1

model.reactions.ATPM.lower_bound = 0;

model.objective = model.reactions.ATPM

mother_model = model.copy()


# Net zero ATP

# In[7]:

model = mother_model.copy()
model.reactions.ATPM.upper_bound = 0.91
sol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
b = show_map(sol)
b.display_in_notebook()


# ## Visualize the flux distribution when growing on CO2 vs glucose

# In[8]:

model = M.copy()

# using the Schuchmann WLP
model.reactions.MTHFD.lower_bound,model.reactions.MTHFD.upper_bound = (0,0)
model.reactions.MTHFR5.lower_bound,model.reactions.MTHFR5.upper_bound = (0,0)
model.reactions.MTHFD_alt.lower_bound,model.reactions.MTHFD_alt.upper_bound = (-1000,1000)
model.reactions.MTHFR5_alt.lower_bound,model.reactions.MTHFR5_alt.upper_bound = (-1000,1000)

# the bad hydrogenase
model.reactions.HYDFDN2r.lower_bound = 0; model.reactions.HYDFDN2r.upper_bound = 0
model.reactions.HYDFDN.lower_bound = -1000; model.reactions.HYDFDN.upper_bound = 1000

model.objective = model.reactions.ATPM
model.reactions.ATPM.lower_bound = 0
model.reactions.EX_ac_e.lower_bound = 1
model.reactions.EX_co2_e.lower_bound = -2
model.reactions.EX_h2_e.lower_bound = -4

pfbaSol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
b = show_map(pfbaSol)
b.save_html('./predictions/acetate_CO2+H2_schuchmann_FDH7.html',overwrite=True)
b.display_in_notebook()


# ## Reproduce all 6 situations from Fig. 3 of (Schuchmann, 2014)
# Including extended analysis by varying the MTHFR5 and MTHFD

# In[9]:

results = [[],[],[]] # this will be filled with real results
for MTHFR in range(2):
    for MTHFD in range(2):
        for FDH in range(3):    
            for HYD in range(2):
                pfbaSol = []
                model = M.copy()

                if FDH == 0:
                    pass # keep the default FDH7
                elif FDH == 1:
                    model.reactions.FDH7.upper_bound = 0; model.reactions.FDH7.lower_bound = 0; 
                    model.reactions.FDHH2.lower_bound = 0; model.reactions.FDHH2.upper_bound = 0
                    model.reactions.FDHFDNADPH.lower_bound = -1000; model.reactions.FDHFDNADPH.upper_bound = 1000
                else:
                    model.reactions.FDH7.upper_bound = 0; model.reactions.FDH7.lower_bound = 0; 
                    model.reactions.FDHH2.lower_bound = -1000; model.reactions.FDHH2.upper_bound = 1000
                    model.reactions.FDHFDNADPH.lower_bound = 0; model.reactions.FDHFDNADPH.upper_bound = 0 

                if HYD == 0:
                    model.reactions.HYDFDN.lower_bound = -1000; model.reactions.HYDFDN.upper_bound = 0 # reversed flux! The Fd + NADH hydrogenase
                    model.reactions.HYDFDN2r.lower_bound = 0; model.reactions.HYDFDN2r.upper_bound = 0 # the Fd + NADPH hydrogenase
                else:
                    model.reactions.HYDFDN.lower_bound = 0; model.reactions.HYDFDN.upper_bound = 0 # reversed flux! The Fd + NADH hydrogenase
                    model.reactions.HYDFDN2r.lower_bound = -1000; model.reactions.HYDFDN2r.upper_bound = 1000 # the Fd + NADPH hydrogenase
                    
                if MTHFD == 0:
                    model.reactions.MTHFD.upper_bound = 0; model.reactions.MTHFD_alt.upper_bound = 1000
                else:
                    model.reactions.MTHFD.upper_bound = 1000; model.reactions.MTHFD_alt.upper_bound = 0
                    
                if MTHFR == 0:    
                    model.reactions.MTHFR5.upper_bound = 0; model.reactions.MTHFR5_alt.upper_bound = 1000
                else: 
                    model.reactions.MTHFR5.upper_bound = 1000; model.reactions.MTHFR5_alt.upper_bound = 0
                try:
                    pfbaSol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
                    results[FDH].append(str(round(abs(pfbaSol.f),3)) + ' (' + str(round(pfbaSol.x_dict['FRNDPR2r_1'],2)) + ')' )
                    b = show_map(pfbaSol)
                    b.save_html('./predictions/fig_3_entry_row_'+str(FDH)+'_col_'+str(HYD)+'_MTHFD_'+str(MTHFD)+'_MTHFR_'+str(MTHFR)+'.html',overwrite=True)
                except:
                    #traceback.print_exc()
                    results[FDH].append('NP')          

iterables = [['NADH', '2 NADH + Fd'], ['NADH', 'NADPH'],['Fd + NADH','Fd + NADPH']]
index = pd.MultiIndex.from_product(iterables, names=['MTHFR', 'MTHFD','FDH \ HYD'])
pd.DataFrame(results, columns=index,index=['Fd','FD+NADPH','H2'])


# ## Show that the entire Wood-Ljungdahlii pathway and the hydrogenases are needed to produce acetate
# The exception is the formate dehydrogenase.

# In[10]:

list_of_knockouts = ['FDH7','FTHFLi','MTHFC',                    'MTHFD','MTHFR5','METR','CODH_ACS']
results = {}
for rxn in list_of_knockouts:
    model = M.copy()
    model.reactions.EX_ac_e.lower_bound = 0
    model.objective = model.reactions.EX_ac_e
    #model.reactions.POR_2.lower_bound = 0
    model.reactions.get_by_id(rxn).lower_bound = 0; model.reactions.get_by_id(rxn).upper_bound = 0
    results[rxn] = round(abs(model.optimize().f),4)
    
print results
pd.DataFrame(results.values(),index=results.keys(),columns=['Acetate flux'])


# Now the hydrogenases

# In[11]:

model = M.copy()
model.reactions.EX_ac_e.lower_bound = 0
model.objective = model.reactions.EX_ac_e
results = {}

# no hydrogenases
model.reactions.HYDFDN2r.lower_bound = 0; model.reactions.HYDFDN2r.upper_bound = 0;
model.reactions.HYD2.lower_bound = 0; model.reactions.HYD2.upper_bound = 0
results['No hydrogenase'] = round(abs(model.optimize().f),4)

model.reactions.HYDFDN2r.lower_bound = 0; model.reactions.HYDFDN2r.upper_bound = 0;
results['Only HYD2'] = round(abs(model.optimize().f),4)

model.reactions.HYD2.lower_bound = 0; model.reactions.HYD2.upper_bound = 0;
model.reactions.HYDFDN.lower_bound = -1000; model.reactions.HYDFDN2r.upper_bound = 0;
results['only HYDFDN'] = round(abs(model.optimize().f),4)

model.reactions.HYD2.lower_bound = 0; model.reactions.HYD2.upper_bound = 0;
model.reactions.HYDFDN.lower_bound = 0; model.reactions.HYDFDN2r.upper_bound = 0;
model.reactions.HYDFDN2r.lower_bound = -1000; model.reactions.HYDFDN2r.upper_bound = 1000;
results['only HYDFDN2r'] = round(abs(model.optimize().f),4)

model.reactions.HYD2.lower_bound = 0; model.reactions.HYD2.upper_bound = 0;
model.reactions.HYDFDN.lower_bound = 0; model.reactions.HYDFDN2r.upper_bound = 0;
model.reactions.HYDFDN2r.lower_bound = 0; model.reactions.HYDFDN2r.upper_bound = 0;
model.reactions.HYDFDi.lower_bound = -1000; model.reactions.HYDFDi.upper_bound = 0;
results['only HYDFDi'] = round(abs(model.optimize().f),4)

pd.DataFrame(results.values(),index=results.keys(),columns=['Acetate flux'])


# In[12]:

model = M.copy()
print model.reactions.HYD2.lower_bound, model.reactions.HYD2.upper_bound,model.reactions.HYD2.reaction
print model.reactions.HYDFDN.lower_bound, model.reactions.HYDFDN.upper_bound,model.reactions.HYDFDN.reaction
print model.reactions.HYDFDN2r.lower_bound, model.reactions.HYDFDN2r.upper_bound,model.reactions.HYDFDN2r.reaction
print model.reactions.HYDFDi.lower_bound, model.reactions.HYDFDi.upper_bound,model.reactions.HYDFDi.reaction


# ## Gear-shifting: step-wise increase of the maintenance flux

# In[13]:

mother_model = M.copy()

# make sure all options are on the table

model.reactions.ATPM.lower_bound,model.reactions.ATPM.upper_bound, = (0,0); 

pfbaSol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
b = show_map(pfbaSol)
b.save_html('./predictions/growth_CO2+H2.html',overwrite=True)
b.display_in_notebook()


# ## Effect of the redox on BHB yieldd

# In[14]:

results = [[],[],[]] # this will be filled with real results
for MTHFR in range(2):
    for MTHFD in range(2):
        for FDH in range(3):    
            for HYD in range(2):
                pfbaSol = []
                model = M.copy()
                
                model.reactions.ATPM.lower_bound = 0
                model.reactions.EX_ac_e.lower_bound = 0
                
                model.reactions.EX_co2_e.lower_bound = 0
                model.reactions.EX_co_e.lower_bound = -2
                model.reactions.EX_h2_e.lower_bound = -4
                
                model.objective = model.reactions.DM_3hbcoa_c

                if FDH == 0:
                    pass # keep the default FDH7
                elif FDH == 1:
                    model.reactions.FDH7.upper_bound = 0; model.reactions.FDH7.lower_bound = 0; 
                    model.reactions.FDHH2.lower_bound = 0; model.reactions.FDHH2.upper_bound = 0
                    model.reactions.FDHFDNADPH.lower_bound = -1000; model.reactions.FDHFDNADPH.upper_bound = 1000
                else:
                    model.reactions.FDH7.upper_bound = 0; model.reactions.FDH7.lower_bound = 0; 
                    model.reactions.FDHH2.lower_bound = -1000; model.reactions.FDHH2.upper_bound = 1000
                    model.reactions.FDHFDNADPH.lower_bound = 0; model.reactions.FDHFDNADPH.upper_bound = 0 

                if HYD == 0:
                    model.reactions.HYDFDN.lower_bound = -1000; model.reactions.HYDFDN.upper_bound = 0 # reversed flux! The Fd + NADH hydrogenase
                    model.reactions.HYDFDN2r.lower_bound = 0; model.reactions.HYDFDN2r.upper_bound = 0 # the Fd + NADPH hydrogenase
                else:
                    model.reactions.HYDFDN.lower_bound = 0; model.reactions.HYDFDN.upper_bound = 0 # reversed flux! The Fd + NADH hydrogenase
                    model.reactions.HYDFDN2r.lower_bound = -1000; model.reactions.HYDFDN2r.upper_bound = 1000 # the Fd + NADPH hydrogenase
                    
                if MTHFD == 0:
                    model.reactions.MTHFD.upper_bound = 0; model.reactions.MTHFD_alt.upper_bound = 1000
                else:
                    model.reactions.MTHFD.upper_bound = 1000; model.reactions.MTHFD_alt.upper_bound = 0
                    
                if MTHFR == 0:    
                    model.reactions.MTHFR5.upper_bound = 0; model.reactions.MTHFR5_alt.upper_bound = 1000
                else: 
                    model.reactions.MTHFR5.upper_bound = 1000; model.reactions.MTHFR5_alt.upper_bound = 0
                try:
                    pfbaSol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
                    results[FDH].append(str(round(abs(pfbaSol.f),3)) + ' (' + str(round(pfbaSol.x_dict['FRNDPR2r_1'],2)) + ')' )
                    b = show_map(pfbaSol)
                    b.save_html('./predictions/3hbcoa_synth_row_'+str(FDH)+'_col_'+str(HYD)+'_MTHFD_'+str(MTHFD)+'_MTHFR_'+str(MTHFR)+'.html',overwrite=True)
                except:
                    results[FDH].append('NP')          

iterables = [['NADH', '2 NADH + Fd'], ['NADH', 'NADPH'],['Fd + NADH','Fd + NADPH']]
index = pd.MultiIndex.from_product(iterables, names=['MTHFR', 'MTHFD','FDH \ HYD'])
pd.DataFrame(results, columns=index,index=['Fd','FD+NADPH','H2'])


# In[ ]:




# In[ ]:



