
# coding: utf-8

# # Maps for when the living gets tough: Maneuvering through a hostile energy landscape
# ##### Thierry D.G.A Mondeel, Samrina Rehman, Yanfei Zhang, Malkhey Verma, Peter Dürre, Matteo Barberis and Hans V. Westerhoff
# 
# Notebook by: Thierry Mondeel
# 
# ---
# 
# This jupyter notebook (http://jupyter.org/) contains all code (and a bit extra) to reproduce the analysis in the conference paper for http://www.fosbe2016.ovgu.de/. Credit goes to the http://mybinder.org/ project which provides the server and software to make this notebook available to you. 
# 
# This notebook is a demo attempt at achieving computational reproducibility: https://doi.org/10.1371/journal.pcbi.1003285.
# 
# The story depends heavily on
# - the publication of a genome-wide metabolic map of Clostridium ljungdahlii (Nagarajan,2013) https://doi.org/10.1186/1475-2859-12-118 
# - the review paper by Schuchmann and Müller https://doi.org/10.1038/nrmicro3365
# 
# #### Abstract
# With genome sequencing of thousands of organisms, a scaffold has become available for data integration:  molecular information can now be organized by attaching it to the genes and their gene-expression products.  It is however, the genome that is selfish not the gene, making it necessary to organize the information into maps that enable functional interpretation of the fitness of the genome. Using flux balance analysis one can calculate the theoretical capabilities of the living organism. Here we examine whether according to this genome organized information, organisms such as the ones present when life on Earth began, are able to assimilate the Gibbs energy and carbon that life needs for its reproduction and maintenance, from a relatively poor Gibbs-energy environment. We shall address how Clostridium ljungdahlii may use at least two special features and one special pathway to this end: gear-shifting, electron bifurcation and the Wood-Ljungdahl pathway. Additionally, we examined whether the C. ljungdahlii map can also help solve the problem of waste management. We find that there is a definite effect of the choices of redox equivalents in the Wood-Ljungdahl pathway and the hydrogenase on the yield of interesting products like hydroxybutyrate. We provide a drawing of a subset of the metabolic network that may be utilized to project flux distributions onto by the community in future works. Furthermore, we make all the code leading to the results discussed here publicly available for the benefit of future work. 

# ### How to interact with this document
# This Jupyter notebook (http://jupyter.org/) contains text cells (like this one), input and output cells. Input cells will contain code to do something in Python and (usually) perform a simulation. The output cells contain the result, i.e. a table or an image. 
# 
# The idea of this document is that you can:
# - Check the code by reading it
# - Run the code by going to "cell > run all" or by selecting one input cell and clicking the ">|" button in the toolbar at the top of the screen.
# - See the results (in the output cells
# - Edit the code (if you know how to write code) and see the results of your changes

# ### Introduction
# In this work, we aim to compare the analysis by Schuchmann and Müller with the predictions of maximal ATP synthesis emanating from the genome-wide metabolic reconstruction of the model acetogen Clostridium ljungdahlii through flux-balance analysis (FBA) (Orth, 2010). Some analysis on the effect of redox equivalents on growth and product synthesis was already present in (Nagarajan, 2013). Specifically, it was shown that the genome-wide map predicts the possibility of growth on CO2/H2 and CO and the effect of various options in redox equivalents were analyzed under the knockout of acetate kinase. However, we hope to extend that analysis here by including various alternative reactions that were considered in the treatment by Schuchmann and Müller. Specifically, we will investigate alternatives in the electron donors/acceptors for various enzymes and their effect on ATP yield coupled to acetogenesis. Additionally, we will focus on the importance of the Wood-Ljungdahl pathway as opposed to single enzymes, the need for electron bifurcation and the Nfn complex, the concept of gear-shifting, the requirement of low gear and advantages of high gear operation, and how much product yield might be attained when engineering C. Ljungdahlii with two additional genes for producing poly-hydroxybutyrate (PHB) under various redox alterations.
# 
# We started from the model from (Nagarajan,2013) downloaded from http://bigg.ucsd.edu/models/iHN637 
# Then we simply added the reactions considered in Schuchmann and Müller that were not present yet in the model.
# 
# In all simulations unless stated otherwise the objective is the ATP maintenance reactions, we allow CO2/H2 uptake in a (2/4) ratio and the output flux of acetate is forced to be 1. 

# ### On flux balance analysis
# For all simulation we use flux balance analysis (Orth,2010): http://dx.doi.org/10.1038/nbt.1614 with help from COBRApy (Ebrahim, 2013): https://doi.org/10.1186/1752-0509-7-74 
# 
# Briefly, this technique concerns the following linear programming problem:
# maximize or 
# $$\text{minimize } Z=c^T v, \text{such that for all } k:$$
# $$Sv=0$$ 
# $$\alpha_k \leq v_k \leq \beta_k$$
# 
# where $S$ is the stoichiometric matrix for the metabolites, $v$ is the vector of fluxes through all reactions including exchange reactions with the environment of the system considered, $c$ is a vector of weights generating the linear combination of fluxes that make up the objective function $Z$ and $\alpha$ and $\beta$ are the vectors of lower and upper bounds on these fluxes.
# 
# A flux distribution returned by FBA is therefore such that all metabolites are produced and consumed in equal amounts, the flux boundaries are accommodated and the flux distribution maximizes (or minimizes) a linear combination of fluxes in the model. 

# ### Set up the python environment
# Nothing interesting here, just execute this. This cell of code loads the required python modules and the updated (nagarajan,2013) model.
# 
# 

# In[1]:

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
# The following reactions were blocked in the original (Nagarajan,2013) model and still blocked in our model, FYI.

# In[2]:

model = M.copy()
d = {}
for rxn in model.reactions:
    if rxn.lower_bound == 0 and rxn.upper_bound == 0:
        d[rxn.id] = [rxn.name,rxn.reaction]
        
df = pd.DataFrame.from_dict(d)
df = df.transpose()
df.columns = ['Name','Reaction']
df


# ## Make sure there is no Gibbs free energy being generated from nothing
# Check that in no combination of enzyme alternatives there can be energy production (ATP) from nothing. 
# By "from nothing" we mean that we allow nothing to enter the cell from the outside (i.e. no glucose uptake). 
# 
# The table generated should contain only zeros, indicating zero ATP synthesis.

# In[3]:

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
                    model.reactions.FDH7.lower_bound = 0; model.reactions.FDH7.upper_bound = 0; 
                    model.reactions.FDHH2.lower_bound = 0; model.reactions.FDHH2.upper_bound = 0
                    model.reactions.FDHFDNADPH.lower_bound = -1000; model.reactions.FDHFDNADPH.upper_bound = 1000
                else:
                    model.reactions.FDH7.lower_bound = 0; model.reactions.FDH7.upper_bound = 0; 
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


# ## Visualize the flux distribution when performing acetogenesis
# Using https://escher.github.io/ developed by Zak King we developed a focused network drawing of the C. ljungdahlii genome-wide metabolic map containing a subset of the reactions of interest here. This allows us to visualize the predicted flux distributions.
# 
# For simplicity we only draw: glycolysis, the Wood-Ljungdahl pathway, the branched TCA cycle, the synthesis pathways of acetate and butanol and the main exchange reactions (input/output) of interest.
# 
# Here we visualize the case where we allow uptake of CO2 and H2 (at 2 and 4 units respectively) and force 1 unit of acetate to be produced. We define the metabolic network to have the exact reactions Schuchmann and Müller considered in the first cell of the table in Figure 3: the Fd dependent formate dehydrogenase together with the NADH dependent hydrogenase. 
# 
# The objective function is set to the ATP maintenance reaction to predict the maximal amount of ATP that may be coupled to this acetogenesis process. According to Schuchmann and Müller this should return 0 for this network configuration.

# In[4]:

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
print 'ATP flux coupled to acetogenesis:',pfbaSol.x_dict['ATPM']
b = show_map(pfbaSol)
b.save_html('./predictions/acetate_CO2+H2_schuchmann_FDH7.html',overwrite=True)
b.display_in_notebook()


# ## Show that the entire Wood-Ljungdahlii pathway and the hydrogenases are needed to produce acetate
# To show that the entire Wood-Ljungdahl pathway needs to be present we sequentially knock out one step in the pathway and ask the model to produce acetate from CO2/H2. 
# If the result is zero, acetate cannot be produced. If it is non-zero acetate may still be synthesized and the considered reaction is not predicted to be essential.
# 
# From the resulting table we see that almost all reactions are essential: the exception is the formate dehydrogenase since it may be generated from pyruvate (see the network drawing)

# In[5]:

list_of_knockouts = ['FDH7','FTHFLi','MTHFC',                    'MTHFD','MTHFR5','METR','CODH_ACS']
results = {}
for rxn in list_of_knockouts:
    model = M.copy()
    model.reactions.EX_ac_e.lower_bound = 0
    model.objective = model.reactions.EX_ac_e
    model.reactions.get_by_id(rxn).lower_bound = 0; model.reactions.get_by_id(rxn).upper_bound = 0
    results[rxn] = round(abs(model.optimize().f),4)
    
pd.DataFrame(results.values(),index=results.keys(),columns=['Acetate flux'])


# ### Show also that the hydrogenases are essential
# Using the same approach as for the WLP, we check whether the hydrogenases are essential for the acetogenesis process.
# 
# The ljungdahlii GEMM actually contains 4 hydrogenase reactions. The result table clearly shows that at least one hydrogenase must be present and that acetogenesis is possible with either the HYDFDN2r, HYDFDN or the HYDFDi.

# In[6]:

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


# ## Reproduce all 6 situations from Fig. 3 of (Schuchmann, 2014)
# Here we consider a similar analysis (but performed with FBA on the metabolic map instead of on paper on a small reaction network) as in Figure 3 in the publication by Schuchmann and Müller.
# 
# We extended the considered scenarios of Figure 3 by also considering alternatives in the MTHFD and MTHFR reactions since Schuchmann and Müller consider different versions of these enzymes than the model by Nagarajan et al. 
# 
# In brackets we show the flux through the Nfn complex to show its contribution to each flux pattern. The set of Nfn fluxes shows that first of all it is not essential and second it may function at different flux levels and in different directions.

# In[7]:

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


# ### Effect of the redox on BHB yield
# Since the enzymes affect ATP yield coupled to acetogenesis, do they also affect product yield? Here we perform the same simulations as above but using beta-hydroxybutyrate synthesis as the objective function

# In[8]:

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



