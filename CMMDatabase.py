import sys

#### STRAINS
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

whic_met_cmm= {
    'N': ['Glucose','Acetate','Xylose','Arabinose'], # Phenotype names
    'S': ['GLC','ACE','XYL','ARA'], # Substrate
    'R': [1,1,1,1], # Phenotype Regulation dimensionless
    'type': 'inhibition', # Regulation type
    'v': [0.7,0.3,0.3,0.5], # Maximum specific growth rate 1/h
    'Ks': [0.25,0.01,0.01,0.01], # Half-saturation constant g/L
    'Yx': [0.5,0.3,0.3,0.4],  # Yield coefficient g biomass/g substrate
    'P': ['ACE',None,None,None], # Product
    'Yp': [0.3,0.0,0.0,0.0], # Yield coefficient g product/g substrate
    'I': ['ACE','GLC','GLC','GLC'],  # Inhibition constant g/L
    'Ki': [sys.maxsize,0.25,0.25,0.25]  # Product inhibition constant g/L
    }

scer_met_cmm= {
    'N': ['Glucose','Acetate','Ethanol'], # Phenotype names
    'S': ['GLC','ACE','ETH'], # Substrate
    'R': [1,1,1], # Phenotype Regulation dimensionless
    'type': 'inhibition', # Regulation type
    'v': [0.5,0.1,0.2], # Maximum specific growth rate 1/h
    'Ks': [0.25,0.01,0.01], # Half-saturation constant g/L
    'Yx': [0.5,0.3,0.2],  # Yield coefficient g biomass/g substrate
    'P': ['ETH',None,None], # Product
    'Yp': [0.1,0.0,0.0], # Yield coefficient g product/g substrate
    'I': ['ACE','GLC','GLC'],  # Inhibition constant g/L
    'Ki': [sys.maxsize,1,0.5]  # Product inhibition constant g/L
}

#### Reactor initial parameters
bionet_std = {
    'V': 1.0,  # Volume of the bioreactor L
    'pH': 7.0,  # pH level
    'T': 37.0,  # Temperature in Celsius
    'DOT': 0.005,  # Dissolved oxygen g/L
    'GLC': 5.0,   # Carbon source concentration g/L
    'ACE': 0.0,   # Product source concentration g/L
    'kla': 0.1, # Volumetric mass transfer coefficient 1/h
    }

bionet_cmm1= {
    'V': 1.0,  # Volume of the bioreactor L
    'pH': 7.0,  # pH level
    'T': 37.0,  # Temperature in Celsius
    'DOT': 0.005,  # Dissolved oxygen g/L
    'GLC': 0.0,   # Carbon source concentration g/L
    'XYL': 5.0,   # Xylose concentration g/L
    'ARA': 0.0,   # Arabinose concentration g/L
    'ACE': 0.0,   # Product source concentration g/L
    'ETH': 0.0,   # Ethanol concentration g/L
    'kla': 0.1, # Volumetric mass transfer coefficient 1/h
}

bionet_cmm2= {
    'V': 1.0,  # Volume of the bioreactor L
    'pH': 7.0,  # pH level
    'T': 37.0,  # Temperature in Celsius
    'DOT': 0.005,  # Dissolved oxygen g/L
    'GLC': 5.0,   # Carbon source concentration g/L
    'XYL': 0.0,   # Xylose concentration g/L
    'ARA': 0.0,   # Arabinose concentration g/L
    'ACE': 0.0,   # Product source concentration g/L
    'ETH': 0.0,   # Ethanol concentration g/L
    'kla': 0.1, # Volumetric mass transfer coefficient 1/h
}

bionet_cmm3= {
    'V': 1.0,  # Volume of the bioreactor L
    'pH': 7.0,  # pH level
    'T': 37.0,  # Temperature in Celsius
    'DOT': 0.005,  # Dissolved oxygen g/L
    'GLC': 0.25,   # Carbon source concentration g/L
    'XYL': 0.25,   # Xylose concentration g/L
    'ARA': 0.0,   # Arabinose concentration g/L
    'ACE': 0.0,   # Product source concentration g/L
    'ETH': 0.0,   # Ethanol concentration g/L
    'kla': 0.1, # Volumetric mass transfer coefficient 1/h
}

#### Operation parameters
batch_std= {
    'F_in': 0.0,  # Inflow rate of carbon source L/h
    'M_in': ['GLC'], # Metabolite names in the inflow
    'C_in': [0.0],  # Concentration of M in the inflow g/L
    'F_out': 0.0,  # Outflow rate L/h
    'M_pl': ['ACE'],  # Pulsed metabolites names
    'C_pl': [0.0],  # Concentration of M in the pulse g/L
    'V_pl': 0.00,  # Volume of the pulse L
    'PVA': 0,  #  Pulse Volume Addition, 1 for no correction volume added, 0 for correction volume not added due to outflow adjusting to pulse (height of the outflow point and max outflow > to inflow + pulse volume)
    }

chemostat_dr10std= {
    'F_in': 0.1,  # Inflow rate of carbon source L/h
    'M_in': ['GLC'], # Metabolite names in the inflow
    'C_in': [5.0],  # Concentration of M in the inflow g/L
    'F_out': 0.1,  # Outflow rate L/h
    'M_pl': ['ACE'],  # Pulsed metabolites names
    'C_pl': [0.5],  # Concentration of M in the pulse g/L
    'V_pl': 0.00,  # Volume of the pulse L
    'PVA': 0,  #  Pulse Volume Addition, 1 for no correction volume added, 0 for correction volume not added due to outflow adjusting to pulse (height of the outflow point and max outflow > to inflow + pulse volume)
    }

chemostat_dr10glc= {
    'F_in': 0.1,  # Inflow rate of carbon source L
    'M_in': ['GLC'], # Metabolite names in the inflow
    'C_in': [5.0],  # Concentration of M in the inflow g/L
    'F_out': 0.1,  # Outflow rate L/h
    'M_pl': ['GLC'],  # Pulsed metabolites names
    'C_pl': [50],  # Concentration of M in the pulse
    'V_pl': 0.00,  # Volume of the pulse L
    'PVA': 0,  #  Pulse Volume Addition, 1 for no
    # correction volume added, 0 for correction volume not added due to outflow adjusting to pulse (height of the outflow point and max outflow > to inflow + pulse volume)
    }

chemostat_dr10xyl= {
    'F_in': 0.1,  # Inflow rate of carbon source L
    'M_in': ['XYL'], # Metabolite names in the inflow
    'C_in': [5.0],  # Concentration of M in the inflow g/L
    'F_out': 0.1,  # Outflow rate L/h
    'M_pl': ['XYL'],  # Pulsed metabolites names
    'C_pl': [50],  # Concentration of M in the pulse
    'V_pl': 0.00,  # Volume of the pulse L
    'PVA': 0,  #  Pulse Volume Addition, 1 for no
    # correction volume added, 0 for correction volume not added due to outflow adjusting to
    # pulse (height of the outflow point and max outflow > to inflow + pulse volume)
    }

chemostat_dr10cascade= {
    'F_in': 0.1,  # Inflow rate of carbon source L
    'M_in': ['GLC', 'XYL','ARA','ACE','ETH'], # Metabolite names in the inflow
    'C_in': [0,0,0,0,0],  # Concentration of M in the inflow g/L
    'F_out': 0.1,  # Outflow rate L/h
    'M_pl': ['GLC','XYL'],  # Pulsed metabolites names
    'C_pl': [2.5,2.5],  # Concentration of M in the pulse g/L
    'V_pl': 0.001,  # Volume of the pulse L (100 mL/h for D=0.1 at dt=0.01h)
    'PVA': 0,  #  Pulse Volume Addition, 1 for no correction volume added, 0 for correction volume not added due to outflow adjusting to pulse (height of the outflow point and max outflow > to inflow + pulse volume)
}

# Controllers 