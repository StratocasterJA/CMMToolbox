import sys
import numpy as np

class operation:
    def __init__(self, ID, parameters,dt=0.01,regulation=None):
        self.name = ID
        self.parameters = parameters 
        self.dt = dt
        self.validate_parameters()
        self.regulation = regulation
    
    def validate_parameters(self):
        if not self.parameters:
            raise ValueError("Operation parameters are missing.")
        if 'F_in' not in self.parameters or 'F_out' not in self.parameters:
            raise ValueError("Inflow and outflow rates must be defined.")
        if self.parameters['F_in'] < 0 or self.parameters['F_out'] < 0:
            raise ValueError("Inflow and outflow rates must be positive.")

        if 'M_in' not in self.parameters or 'C_in' not in self.parameters:
            raise ValueError("Metabolite names and concentrations in the inflow must be defined.")
        if len(self.parameters['M_in']) != len(self.parameters['C_in']):
            raise ValueError("Mismatch in number of metabolite names and concentrations in the inflow.")
        if 'M_pl' not in self.parameters or 'C_pl' not in self.parameters or 'V_pl' not in self.parameters:
            raise Warning("Pulsed metabolites parameters are missing, assuming no pulsed metabolites.")
            self.parameters['M_pl'] = [None]
            self.parameters['C_pl'] = [0.0]
        if len(self.parameters['M_pl']) != len(self.parameters['C_pl']) :
            raise ValueError("Mismatch in number of pulsed metabolites names and concentrations.")
        if 'PVA' not in self.parameters:
            raise Warning("Pulse Volume Addition is missing, assuming no volume added by pulse.")
            self.parameters['PVA'] = 0
        if self.dt <= 0:
            raise ValueError("Time step must be positive.")
        if self.dt == 0:
            raise ValueError("Time step cannot be zero.")
        if self.dt >= 0.5:
            raise Warning("Time step is too large, consider reducing it for more accurate results.")        
        self.parameters['F_in'] = float(self.parameters['F_in'])
        self.parameters['F_out'] = float(self.parameters['F_out'])
        self.parameters['C_in'] = np.array(self.parameters['C_in'], dtype=float)
        self.parameters['C_pl'] = np.array(self.parameters['C_pl'], dtype=float)
        self.parameters['V_pl'] = float(self.parameters['V_pl'])
        self.parameters['PVA']  = float(self.parameters['PVA'])

class bioreactor:
    def __init__(self, biologicals, parameters, id='Reactor_1', t=0):
        self.id = id
        self.strains = 0
        self.time = t
        self.biologicals = biologicals
        self.parameters = parameters
        self.check_biologicals()
        self.check_parameters() 
        self.metabolites = list(self.parameters.keys())
        for item in ['V', 'pH', 'T', 'DOT','kla']:
            if item in self.metabolites:
                self.metabolites.remove(item)
        for metabolite in self.metabolites:
            self.parameters[metabolite] = float(self.parameters[metabolite])
    
    def check_biologicals(self):
        def check_metabolism(met):
            if not met:
                raise ValueError("biologicals components are missing.")
            if not met['S']:
                raise ValueError("Substrate components are missing in biologicals.")
            if not met['N']:
                for i in range(len(met['S'])):
                    met['N'].append(str(i+1))
            if not met['R']:
                for i in range(len(met['N'])):
                    met['R'].append(1)
            elif len(met['R']) != len(met['N']):
                raise Warning("Mismatch in number of phenotypes and regulation variables, no regulation will be used, setting values to 1.")
            met['R'] = np.array(met['R'], dtype=float)
            if not met['type']:
                met['type'] = 'default'  # Default metabolism type
            if not met['v']:
                raise ValueError("Maximum specific growth rates are missing in biologicals.")
            elif len(met['v']) != len(met['N']):
                raise ValueError("Mismatch in number of maximum specific growth rates and phenotypes.")
            met['v'] = np.array(met['v'], dtype=float)
            if not met['Ks']:
                raise ValueError("Half-saturation constants are missing in biologicals.")
            elif len(met['Ks']) != len(met['N']):
                raise ValueError("Mismatch in number of half-saturation constants and phenotypes.")
            met['Ks'] = np.array(met['Ks'], dtype=float)
            if not met['Yx']:
                raise ValueError("Yield coefficients are missing in biologicals.")
            elif len(met['Yx']) != len(met['N']):
                raise ValueError("Mismatch in number of yield coefficients and phenotypes.")
            met['Yx'] = np.array(met['Yx'], dtype=float)
            if not met['P']:
                print("No product defined in biologicals, assuming None for all phenotypes.")
                for i in range(len(met['N'])):
                    met['P'].append(None)
            elif len(met['P']) != len(met['N']):
                raise ValueError("Mismatch in number of products and phenotypes.")
            if not met['Yp']:
                print("No yield coefficients for product defined in biologicals, assuming 0 for all phenotypes.")
                met['Yp'] = [0.0] * len(met['N'])
            elif len(met['Yp']) != len(met['N']):
                raise ValueError("Mismatch in number of yield coefficients for product and phenotypes.")
            met['Yp'] = np.array(met['Yp'], dtype=float)
            if not met['I']:
                print("No inhibition constants defined in biologicals, assuming None for all phenotypes.")
                met['I'] = [None] * len(met['N'])
            elif len(met['I']) != len(met['N']):
                raise ValueError("Mismatch in number of inhibition constants and phenotypes.")
            if not met['Ki']:
                print("No product inhibition constants defined in biologicals, assuming sys.maxsize for all phenotypes.")
                met['Ki'] = [sys.maxsize] * len(met['N'])
            elif len(met['Ki']) != len(met['N']):
                raise ValueError("Mismatch in number of product inhibition constants and phenotypes.")
            met['Ki'] = np.array(met['Ki'], dtype=float)
            return True
        if not self.biologicals:
            raise ValueError("biologicals components are missing.")
        self.strains= len(self.biologicals['strains'])
        if self.strains == 0:
            raise ValueError("No strains defined in biologicals components.")
        if not self.biologicals['met']:
            raise ValueError("Metabolism type is missing in biologicals components.")
        if len(self.biologicals['met']) != self.strains:
            raise ValueError("Mismatch in number of strains and metabolism types in biologicals components.")
        r= [False] * self.strains
        for i in range(self.strains):
            r[i]=check_metabolism(self.biologicals['met'][i])
            if not r[i]:
                raise ValueError(f"Invalid metabolism type for strain {self.biologicals['strains'][i]}.")
        if 'X' not in self.biologicals:
            raise ValueError("Biomass concentration is missing in biologicals.")
        if len(self.biologicals['X']) != self.strains:
            raise ValueError("Mismatch in number of biomass concentrations and strains.")
        self.biologicals['X'] = np.array(self.biologicals['X'], dtype=float)
    
    def check_parameters(self):
        if not self.parameters:
            raise ValueError("parameters components are missing.")
        if 'V' not in self.parameters:
            raise Warning("Volume of the bioreactor is missing in parameters. setting to default 1.0 L.")
            self.parameters['V'] = 1.0
        if self.parameters['V'] <= 0:
            raise Warning("Volume of the bioreactor must be positive. setting to default 1.0 L.")
            self.parameters['V'] = 1.0
        self.parameters['V'] = float(self.parameters['V'])
        if 'pH' not in self.parameters:
            raise Warning("pH level is missing in parameters. setting to default 7.0.")
            self.parameters['pH'] = 7.0
        if self.parameters['pH'] < 0 or self.parameters['pH'] > 14:
            raise ValueError("pH level must be between 0 and 14.")    
        self.parameters['pH'] = float(self.parameters['pH'])
        if 'T' not in self.parameters:
            raise Warning("Temperature is missing in parameters. setting to default 37.0 Celsius.")
            self.parameters['T'] = 37.0
        if self.parameters['T'] < 0:
            raise Warning("Temperature is negative.")
        if self.parameters['T'] > 100:
            raise Warning("Temperature is above 100 Celsius.")
        self.parameters['T'] = float(self.parameters['T'])
        if 'DOT' not in self.parameters:
            raise Warning("Dissolved oxygen is missing in parameters. setting to default 8.0 mg/L. 100%% saturation at 37 Celsius.")
            self.parameters['DOT'] = 0.008
        if self.parameters['DOT'] < 0:
            raise ValueError("Dissolved oxygen concentration cannot be negative.")
        if self.parameters['DOT'] > 0.008:
            raise ValueError("Dissolved oxygen concentration cannot exceed 0.008 g/L at 37 Celsius. make sure to adjust units accordingly.")
        self.parameters['DOT'] = float(self.parameters['DOT'])
        if 'kla' not in self.parameters:
            raise Warning("Volumetric mass transfer coefficient is missing in parameters. Setting to default 0.1 1/h.")
            self.parameters['kla'] = 0.1
        self.parameters['kla'] = float(self.parameters['kla'])
        for i in range(self.strains):
            for substrate in self.biologicals['met'][i]['S']:
                if substrate not in self.parameters:
                    raise Warning(f"Concentration for substrate {substrate} is missing in parameters.seting to 0.")
                    self.parameters[substrate] = 0.0
                if self.parameters[substrate] < 0:
                    raise Warning(f"Concentration for substrate {substrate} cannot be negative. setting to 0.")
                    self.parameters[substrate] = 0.0
            for product in self.biologicals['met'][i]['P']:
                if product not in self.parameters:
                    if product is None:
                        continue
                    else:
                        raise Warning(f"Concentration for product {product} is missing in parameters.seting to 0.")
                        self.parameters[product] = 0.0
                if self.parameters[product] < 0:
                    raise Warning(f"Concentration for product {product} cannot be negative. setting to 0.")
                    self.parameters[product] = 0.0
        return True

    def print_info(self):
        print(f"Bioreactor ID: {self.id}")
        print(f"Bioreactor state at time {self.time} hours.")
        print(f"biologicals Components: {self.biologicals}")
        print(f"parameters Components: {self.parameters}")
    
    def operate_dt(self,operation,verbose=False):
    # Update bioreactor volume based on inflow and outflow rates
    # VOLUME
        dV=+ operation.parameters['F_in'] * operation.dt - operation.parameters['F_out'] * operation.dt + operation.parameters['V_pl']* operation.parameters['PVA']
        Vi=self.parameters['V']
        self.parameters['V'] = self.parameters['V'] + dV
    # Outflow on MASS effects
        for metabolite in self.metabolites:
            self.parameters[metabolite] += -self.parameters[metabolite] * operation.parameters['F_out'] * operation.dt

        for biomass in range(len(self.biologicals['X'])):
            print(f"Biomass {biomass} concentration before outflow: {self.biologicals['X'][biomass]}")
            self.biologicals['X'][biomass] += - self.biologicals['X'][biomass] * operation.parameters['F_out'] * operation.dt
            print(- self.biologicals['X'][biomass] * operation.parameters['F_out'] * operation.dt)
    # Inflow on MASS effects
        for metabolite in operation.parameters['M_in']:
            if metabolite not in self.parameters:
                self.parameters[metabolite] = 0.0 + operation.parameters['C_in'][operation.parameters['M_in'].index(metabolite)] * operation.parameters['F_in'] * operation.dt
            else:
                self.parameters[metabolite] += operation.parameters['C_in'][operation.parameters['M_in'].index(metabolite)] * operation.parameters['F_in'] * operation.dt
    # Pulses on MASS effects
        for metabolite in operation.parameters['M_pl']:
            if metabolite not in self.parameters:
                self.parameters[metabolite] = 0.0 + operation.parameters['C_pl'][operation.parameters['M_pl'].index(metabolite)] * operation.parameters['V_pl'] / self.parameters['V']
            else:
                self.parameters[metabolite] += operation.parameters['C_pl'][operation.parameters['M_pl'].index(metabolite)] * operation.parameters['V_pl'] / self.parameters['V']
            if operation.parameters['PVA'] == 1:
                self.parameters[metabolite] = self.parameters[metabolite] * Vi / self.parameters['V']
                for biomass in range(len(self.biologicals['X'])):
                    self.biologicals['X'][biomass] = self.biologicals['X'][biomass] * Vi / self.parameters['V']
    
    # BIOMASS METABOLISM EFFECTS
    #def metabolism_dt(self,dt,type='default'):
        Vi=self.parameters['V']
        totalphenotypes = 0
        for biomass in range(len(self.biologicals['X'])):
            totalphenotypes += len(self.biologicals['met'][biomass]['N'])
        mu = np.zeros(totalphenotypes)
        qs = np.zeros(totalphenotypes)
        qp = np.zeros(totalphenotypes) 
        ph_index = 0  
        for biomass in range(len(self.biologicals['X'])):
            Xi=self.biologicals['X'][biomass]
            for phenotype in range(len(self.biologicals['met'][biomass]['N'])):
                # Simplify the access to the parameters
                r=self.biologicals['met'][biomass]['R'][phenotype]
                v=self.biologicals['met'][biomass]['v'][phenotype]
                ks=self.biologicals['met'][biomass]['Ks'][phenotype]
                ki=self.biologicals['met'][biomass]['Ki'][phenotype]
                Yx=self.biologicals['met'][biomass]['Yx'][phenotype]
                Yp=self.biologicals['met'][biomass]['Yp'][phenotype]
                S=self.biologicals['met'][biomass]['S'][phenotype]
                s=self.parameters[S]
               
                if self.biologicals['met'][biomass]['I'][phenotype] is None:
                    i=0.0
                    I= None
                    iflag=False
                else:
                    I=self.biologicals['met'][biomass]['I'][phenotype]
                    i=self.parameters[I]
                    iflag=True
                if self.biologicals['met'][biomass]['P'][phenotype] is not None:
                    P=self.biologicals['met'][biomass]['P'][phenotype]
                    p=self.parameters[P]
                    pflag=True
                else:
                    P=None
                    p=0.0
                    pflag=False
                # Calculate the specific growth rate,
                if self.biologicals['met'][biomass]['type'] == 'None':
                    # Default metabolism dynamics
                    mu[ph_index] = v * s / (ks + s)
                    ph_index += 1
                elif self.biologicals['met'][biomass]['type'] == 'inhibition':
                    # Default metabolism dynamics
                    mu[ph_index] = v * s / (ks + s)*(1/(1 + (i / ki)))
                    ph_index += 1
                elif self.biologicals['met'][biomass]['type'] == 'expression_regulation':
                    # External regulation dynamics
                    mu[ph_index] = r * v * s / (ks + s)*(1/(1 + (i / ki)))
                    ph_index += 1
                else:
                    raise ValueError(f"Unknown metabolic regulation type: {self.biologicals['met'][biomass]['type']}. Currently supported types are 'None', 'inhibition', and 'expression_regulation'.")
                # Outputting Biologicals feedback to the bioreactor
                So= self.parameters[S]
                self.parameters[S] -= mu[ph_index-1] / Yx * Xi * operation.dt
                if self.parameters[S] < 0:
                    self.parameters[S] = 0.0
                    dS= So 
                else:
                    dS= mu[ph_index-1] / Yx * Xi * operation.dt
                if pflag:
                    self.parameters[P] += Yp*dS
                    if self.parameters[P] < 0:
                        self.parameters[P] = 0.0
                print(self.biologicals['X'][biomass])
                self.biologicals['X'][biomass] += Yx*dS
                if self.biologicals['X'][biomass] < 0:
                    self.biologicals['X'][biomass] = 0.0
                print(self.biologicals['X'][biomass])

        #return self, mu, qs, qp

        self.time += operation.dt
        if verbose:
            print(f"Bioreactor state updated to {self.time} hours.")
            return dV, Yx*dS/operation.dt/Xi, dS/operation.dt/Xi, Yp*dS/operation.dt/Xi

        

                
   
