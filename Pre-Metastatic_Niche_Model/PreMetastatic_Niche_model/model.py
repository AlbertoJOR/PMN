
from pickle import TRUE
from re import T
from mesa import Model
from mesa.space import SingleGrid
from mesa.datacollection import DataCollector
from mesa.time import RandomActivation
from .agents import *
import csv
from math import floor
# Moléculas
M_IL_12       = 100
M_fibronectin = 101
M_TGFB        = 102
M_S100        = 103
M_IL6         = 104
M_VEGF        = 105


class PreMetastaticNicheModel(Model):
    """
        A model for the formation and dynamics of the Pre-Metastatic Niche
    """
    
    height = 20
    width = 20    

    def __init__(self, width, height, R_BMDC, R_CD8, R_M1, R_N1, R_MCDEV,R_MDE, R_PDL1,R_ODEV, R_MPCDE, R_CTC,
    R_IL_12, R_fibronectin, R_TGFB, R_S100, R_IL6, R_VEGF, Stromal_init, Fibroblast_init):        
        
        self.schedule = RandomActivation(self)
        self.grid = SingleGrid(width, height, torus=True)
        self.current_id =0
        self.height = height
        self.width = width
        self.total_agents = 0
        self.max_agents = (width*height)/2


        # Num of molecules in the whole space
        self.num_molecules ={M_IL_12       : 0,
                            M_fibronectin : 0,
                            M_TGFB        : 0,
                            M_S100        : 0,
                            M_IL6         : 0,
                            M_VEGF        : 0} 
        
        

        # Consumption production rates
        self.mol_rates =   {M_IL_12       : R_IL_12,
                            M_fibronectin : R_fibronectin,
                            M_TGFB        : R_TGFB,
                            M_S100        : R_S100,
                            M_IL6         : R_IL6,
                            M_VEGF        : R_VEGF} 
        # Limit production 
        self.mol_max_lim = {M_IL_12       : 1000,
                            M_fibronectin : 1000,
                            M_TGFB        : 1000,
                            M_S100        : 1000,
                            M_IL6         : 1000,
                            M_VEGF        : 1000} 
        # Limit for biological effect
        self.mol_bio_limit = {M_IL_12       : 200,
                            M_fibronectin : 200,
                            M_TGFB        : 200,
                            M_S100        : 200,
                            M_IL6         : 200,
                            M_VEGF        : 200} 
        self.N_vasc_per=0
        self.N_rec_CD11b=0
        self.N_immunosupres=0

        self.pro_vasc_per=1
        self.pro_rec_CD11b=1
        self.pro_immunosupres=1

        self.th_vasc_per=200
        self.th_rec_CD11b=200
        self.th_immunosupres=200

        #Ising
        self.ising = 0

        #recruitment rate
        self.recuit_rates={T_BMDC : R_BMDC
        ,T_CD8 : R_CD8
        ,T_M1 : R_M1
        ,T_N1 : R_N1
        # Misma tasa para M2
        ,T_M2 : R_M1
        ,T_N2 : R_N1
        ,E_MCDEV : R_MCDEV
        ,E_MDE : R_MDE
        ,E_PD_L1 :R_PDL1
        ,E_ODEV : R_ODEV
        ,E_MPCDE :R_MPCDE, 
        T_CTC: R_CTC}
        #Acumulators ( when the rates are fractions)
        self.acumulators= {T_BMDC : 0
                            ,T_CD8 : 0 
                            ,T_M1 :0 
                            ,T_N1 : 0 
                            ,E_MCDEV : 0
                            ,E_MDE : 0 
                            ,E_PD_L1 :0
                            ,E_ODEV : 0
                            ,E_MPCDE :0
                            ,T_CTC : 0
                            ,T_M2 : 0
                            ,T_N2 : 0}
        self.acumulators_molecules= {T_BMDC : 0
                            ,T_CD8 : 0 
                            ,T_M1 :0 
                            ,T_N1 : 0 
                            ,E_MCDEV : 0
                            ,E_MDE : 0 
                            ,E_PD_L1 :0
                            ,E_ODEV : 0
                            ,E_MPCDE :0
                            , T_CTC : 0
                            ,T_M2 : 0
                            ,T_N2 : 0}
        

        # Agents

        self.num_agents = {T_BMDC : 0
                            ,T_NK: 0
                            ,T_BCELL: 0
                            ,T_CD4 : 0
                            ,T_CD8: 0
                            ,T_TREG: 0
                            ,T_M1 : 0
                            ,T_N1 : 0
                            ,T_M2 : 0
                            ,T_N2 : 0
                            ,T_STROMAL: 0
                            ,T_FIBROBLAST:0
                            ,T_CTC: 0
                            ,T_CD11:0
                            ,E_MCDEV:0
                            ,E_MDE:0
                            ,E_MPCDE:0
                            ,E_ODEV:0
                            ,E_PD_L1: 0} 
        self.cont_BMDC = 0
        # Banderas
        self.PMN_var= False
        self.ctc_col = False
        self.vasc_per = False
        self.BMDC_REC = False
        self.inmunosup = False
    
        
        
        

        self.datacollector = DataCollector(
            {
                "BMDC": lambda m: m.num_agents[T_BMDC],
                "FIBROBLAST": lambda m: m.num_agents[T_FIBROBLAST],
                "B": lambda m: m.num_agents[T_BCELL],                
                "T_CD4": lambda m: m.num_agents[T_CD4],                
                "T_CD8": lambda m: m.num_agents[T_CD8],
                "M1": lambda m: m.num_agents[T_M1],
                "N1": lambda m: m.num_agents[T_N1],
                "M2": lambda m: m.num_agents[T_M2],
                "N2": lambda m: m.num_agents[T_N2],
                "CD11": lambda m: m.num_agents[T_CD11],
                "CTC": lambda m: m.num_agents[T_CTC],
                "PDL1": lambda m: m.num_agents[E_PD_L1],
                "MCDEV": lambda m: m.num_agents[E_MCDEV],
                "MDE": lambda m: m.num_agents[E_MDE],
                "ODEV": lambda m: m.num_agents[E_ODEV],
                "MPCDE": lambda m: m.num_agents[E_MPCDE],
                "ISING" :lambda m: m.ising,
                "IL_12" : lambda m: m.num_molecules[M_IL_12],
                "fibronectin": lambda m: m.num_molecules[M_fibronectin],
                "TGFB": lambda m: m.num_molecules[M_TGFB],
                "S100": lambda m: m.num_molecules[M_S100],
                "IL6" : lambda m: m.num_molecules[M_IL6],
                "VEGF" : lambda m: m.num_molecules[M_VEGF],
                
            }
        )

        self.datacollector.collect(self)

        #CSV
        self.csv_name = "result_cell_abundance.csv"
        csv_row = ["step","ISING","BMDC", "NK", "B", "T_CD4", "T_CD8", "Treg", "M1", "N1", "M2", "N2","Stromal","Fibroblast", 
         "CTC","CD11","MCDEV", "MDE","MPCDE","ODEV", "PD_L1", "IL_12", "fibronectin", "TGFB", "S100", "IL6", "VEGF"]
        with open(self.csv_name,'w',newline='') as file:
            writer = csv.writer(file)
            writer.writerow(csv_row)

        #Initial agents
        self.create_agents(T_STROMAL, Stromal_init)
        self.create_agents(T_FIBROBLAST, Fibroblast_init)
        
        self.step_counter = 0
        self.running = True

    def step(self):
        """Acciona las interacciones que ocurren dentro del sistema
        y alamacena los datos
        """
        # Collect Data 
        self.total_agents = 0
        # Total de agentes en el sistema 
        for k, v in self.num_agents.items():
            self.total_agents += v
        # Cuenta los agentes individualmente para el csv
        csv_row = []
        csv_row.append(self.step_counter)
        csv_row.append(self.ising)
        self.step_counter += 1
        for key, value in self.num_agents.items():
            csv_row.append(value)
        for k, val in self.num_molecules.items():
            csv_row.append(val)
        self.cont_BMDC += self.num_agents[T_BMDC]


        # activa las interacciones de todos los agentes 
        self.schedule.step()
    


        #Recruit new agents
        
        if(self.max_agents>self.total_agents):
            # Entrada al sistema
            self.recruit_cells(T_BMDC)
            self.recruit_cells(T_CD8)
            self.recruit_cells(T_CTC)
            self.recruit_cells(T_M1)
            self.recruit_cells(T_N1)
            self.recruit_cells(E_PD_L1)
            # self.recruit_cells(T_NK)
            # Reclutamiento por moléculas.
            self.molecule_recruitment(T_M2, M_TGFB)
            self.molecule_recruitment(T_BMDC, M_VEGF)
            self.molecule_recruitment(T_CTC, M_S100)
            self.molecule_recruitment(T_M1, M_fibronectin)
            # Entrada de moléculas
            self.recruit_cells(E_MCDEV)
            self.recruit_cells(E_MDE)
            self.recruit_cells(E_MPCDE)
            self.recruit_cells(E_ODEV)
        # Simula la entrada de moléculas por parte del tumor primario
        self.increment_molecules(M_VEGF)
        self.increment_molecules(M_IL6)
        self.increment_molecules(M_TGFB)
        self.increment_molecules(M_S100)

        # Simula la producción de moleculas por los agentes dentro del sistema
        self.agent_molecule_production(T_STROMAL, M_fibronectin)
        self.agent_molecule_production(T_STROMAL, M_TGFB)
        self.agent_molecule_production(T_BMDC, M_fibronectin)

        # Simula la expresión de moléculas dentro del sistema
        self.molecule_molecule_production(M_TGFB, M_S100)
        self.molecule_molecule_production(M_VEGF, M_S100)

        # Evalúa las banderas
        self.PMN_flags()

            
            
                 
        
        # Recolecta datos para ser presentados en las gráficas del servidor
        # Escribe los datos en el csv
        self.datacollector.collect(self)
        with open(self.csv_name,'a',newline='') as file:
            writer = csv.writer(file)
            writer.writerow(csv_row)
        
        
    def create_agents(self, type , n):
        """ Crea una cantidad n de agentes de un tipo determinado

        Args:
            type (int): Indica que tipo de agentes se quieren crear.
            n (int): La cantidad de agentes a crear
        """
        for i in range (n):
            pos = self.grid.find_empty()
            energy = self.random.randrange(100)
            agent = recruit_new_cell(type, self.next_id(),pos,self,True, energy)
            self.grid.position_agent(agent, pos)
            self.schedule.add(agent)
    def protumoral(self):
        """Aumenta la energía del hamiltoniano de ising por cada interacción 
        protumoral que ocurra
        """
        self.ising +=1
    def antitumoral(self):
        """Disminuye la energía del hamiltoniano de ising por cada interacción
        antitumoral
        """
        self.ising -=1
    def recruit_cells (self, type_cell):
        """Recluta nuevos agentes al sistema, dependiendo del tipo de célula que
        a su vez determinará la tasa a la que se reclutan"""
        n_reclutados = 0
        if(self.acumulators[type_cell]>=1):
            n_reclutados = floor(self.acumulators[type_cell])
            self.acumulators[type_cell] -= n_reclutados
            self.create_agents(type_cell,n_reclutados)
        self.acumulators[type_cell] += self.recuit_rates[type_cell]

    def increment_molecules (self, type_molecule):
        """Incrementa la cantidad de moléculas dentro del sistema dependiendo 
        de la tasa del la molécula."""
        if (self.mol_max_lim[type_molecule]>self.num_molecules[type_molecule]):
            self.num_molecules[type_molecule] += self.mol_rates[type_molecule]

    def agent_molecule_production (self, type_cell, type_molecule):
        """Producción de moléculas por parte de los agentes, el total de agentes se multiplica
        por la tasa de producción."""
        self.num_molecules[type_molecule] += self.mol_rates[type_molecule] * self.num_agents[type_cell]

    def molecule_recruitment (self, type_cell, type_molecule):
        """Reclutamiento de agentes por moléculas similar al método de recruit_cells
        solo que se disminuye el número de moleculas a cierta tasa por agente reclutado"""
        n_reclutados = 0
        if(self.num_molecules[type_molecule]>self.mol_bio_limit[type_molecule]):
            if(self.acumulators_molecules[type_cell]>=1 ):
                n_reclutados = floor(self.acumulators_molecules[type_cell])
                self.acumulators_molecules[type_cell] -= n_reclutados
                self.num_molecules[type_molecule] -= n_reclutados
                self.create_agents(type_cell,n_reclutados)
            self.acumulators_molecules[type_cell] += self.mol_rates[type_molecule]

    def molecule_molecule_production(self, type1, type2, stoichiometry1 = 1,  stoichiometry2 = 1):
        """Expresion de moléculas por parte de otra molécula"""
        n_new_mol = stoichiometry1*self.mol_rates[type1]
        if (self.mol_max_lim[type2]>self.num_molecules[type2]):
            self.num_molecules[type1] -= n_new_mol
            self.num_molecules[type2] += n_new_mol* stoichiometry2
        
    def PMN_flags(self):
        """Revisa si hay nicho premetastásico"""
        if(self.cont_BMDC > 200):
            self.BMDC_REC = True
        if(self.N_immunosupres>=self.th_immunosupres):
            self.inmunosup = True
        if (self.N_vasc_per>= self.th_vasc_per):
            self.vasc_per = True
        if(self.BMDC_REC and self.inmunosup and self.vasc_per):
            self.ctc_col = True
            self.PMN_var = TRUE



