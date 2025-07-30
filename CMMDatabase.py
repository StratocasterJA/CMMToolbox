import sys

# STRAINS
eco_met_glc ={    
    'N': ['Glucose','Acetate'], # Phenotype names
    'S': ['GLC','ACE'], # Substrate 
    'R': [1,1], # Phenotype Regulation dimensionless
    'type': 'inhibition', # Reguraltion type
    'v': [0.7,0.3], # Maximum specific growth rate 1/h
    'Ks': [0.25,0.01], # Half-saturation constant g/L
    'Yx': [0.5,0.3],  # Yield coefficient g biomass/g substrate
    'P': ['ACE',None], # Product
    'Yp': [0.3,0.0], # Yield coefficient g product/g substrate
    'I': ['ACE','GLC'],  # Inhibition constant g/L
    'Ki': [sys.maxsize,0.25]  # Product inhibition constant g/L
    }

# Reactor initial parameters
eco_bionet_std = {
    'V': 1.0,  # Volume of the bioreactor L
    'pH': 7.0,  # pH level
    'T': 37.0,  # Temperature in Celsius
    'DOT': 0.005,  # Dissolved oxygen g/L
    'GLC': 5.0,   # Carbon source concentration g/L
    'ACE': 0.0,   # Product source concentration g/L
    'kla': 0.1, # Volumetric mass transfer coefficient 1/h
    }

# Operation parameters
bionet_batchstd= {
    'F_in': 0.0,  # Inflow rate of carbon source L/h
    'M_in': ['GLC'], # Metabolite names in the inflow
    'C_in': [0.0],  # Concentration of M in the inflow g/L
    'F_out': 0.0,  # Outflow rate L/h
    'M_pl': ['ACE'],  # Pulsed metabolites names
    'C_pl': [0.0],  # Concentration of M in the pulse g/L
    'V_pl': 0.00,  # Volume of the pulse L
    'PVA': 0,  #  Pulse Volume Addition, 1 for no correction volume added, 0 for correction volume not added due to outflow adjusting to pulse (height of the outflow point and max outflow > to inflow + pulse volume)
    }
bionet_cc10std= {
    'F_in': 0.1,  # Inflow rate of carbon source L/h
    'M_in': ['GLC'], # Metabolite names in the inflow
    'C_in': [5.0],  # Concentration of M in the inflow g/L
    'F_out': 0.1,  # Outflow rate L/h
    'M_pl': ['ACE'],  # Pulsed metabolites names
    'C_pl': [0.5],  # Concentration of M in the pulse g/L
    'V_pl': 0.00,  # Volume of the pulse L
    'PVA': 0,  #  Pulse Volume Addition, 1 for no correction volume added, 0 for correction volume not added due to outflow adjusting to pulse (height of the outflow point and max outflow > to inflow + pulse volume)
    }



# Controllers 