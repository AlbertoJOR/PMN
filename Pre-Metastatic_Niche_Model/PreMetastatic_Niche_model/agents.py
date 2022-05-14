

from mesa import Agent
import random


# Cell Types id
T_BMDC = 1
T_NK = 6
T_BCELL = 7
T_CD8 = 8
T_CD4 = 9
T_TREG = 10
T_M1 = 11
T_N1 = 12
T_M2 = 13
T_N2 = 14
T_STROMAL = 15
T_FIBROBLAST =16
T_CTC = 17 
T_CD11 = 24

# Exosomes types id
E_PD_L1= 19
E_MDE = 20 
E_MPCDE = 21
E_ODEV = 22
E_MCDEV = 23



def recruit_new_cell(typecell,unique_id, pos, model, moore, energy=100):
    """Crea un nuevo agentes de un tipo determinato typecell

    Args:
        typecell (int):
        unique_id (int):
        pos ((int, int)):
        model ((model)):
        moore (boolean): 
        energy (int):

    Returns:
        _type_: _description_
    """
    if(typecell == 1):
        return BMDC(unique_id, pos, model, moore, energy)
    elif(typecell == 6):
        return NK(unique_id, pos, model, moore, energy)
    elif(typecell == 7):
        return BCell(unique_id, pos, model, moore, energy)
    elif(typecell == 8):
        return CD8(unique_id, pos, model, moore, energy)
    elif(typecell == 9):
        return CD4(unique_id, pos, model, moore, energy)
    elif(typecell == 10):
        return Treg(unique_id, pos, model, moore, energy)
    elif(typecell == 11):
        return M1(unique_id, pos, model, moore, energy)
    elif(typecell == 12):
        return N1(unique_id, pos, model, moore, energy)
    elif(typecell == 13):
        return M2(unique_id, pos, model, moore, energy)
    elif(typecell == 14):
        return N2(unique_id, pos, model, moore, energy)
    elif(typecell == 15):
        return StromalC(unique_id, pos, model, moore, energy)
    elif(typecell == 16):
        return Fibroblast(unique_id, pos, model, moore, energy)
    elif(typecell == 17):
        return CTC(unique_id, pos, model, moore, energy)
    elif(typecell == 19):
        return PD_L1(unique_id, pos, model, moore, energy)
    elif(typecell == E_MDE):
        return MDE(unique_id, pos, model, moore, energy)
    elif(typecell == E_MPCDE):
        return MPCDE(unique_id, pos, model, moore, energy)
    elif(typecell == E_ODEV):
        return ODEV(unique_id, pos, model, moore, energy)
    elif(typecell == E_MCDEV):
        return MCDEV(unique_id, pos, model, moore, energy)
    elif(typecell == T_CD11):
        return CD11(unique_id, pos, model, moore, energy)
    else:
        return None


class Exosome(Agent):
    """ Clase principal de un exosoma.

    Args:
        Agent (class): Hereda de la clase agente de mesa
    """
    def __init__(self, unique_id, pos, model, moore, energy=None) -> None:
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 18
        self.model.num_agents[self.type] +=1
        
    def empty_spaces(self):
        """Encuentra los espacios vacios al rededor del agente.

        Returns:
            list: espacios vacios al rededor
        """
        list_neigh = []
        x , y = self.pos
        for i in range (-1, 2):
            for j in range(-1,2):
                a = (i+x) %self.model.width
                b = (j+y) %self.model.height
                list_neigh.append((a,b))
        list_neigh.remove((x,y))
        for neighbor in self.model.grid.neighbor_iter(self.pos):
            list_neigh.remove(neighbor.pos)
        return list_neigh
    def random_move(self):
        """ realiza un movimiento aleatorio en la malla del modelo
        """
        empty = self.empty_spaces()
        if len(empty)>0:
            new_position = random.choice(empty)
            #print(new_position)
            self.model.grid.move_agent(self, new_position)
    def destroy_self(self):
        """Elimina a si mismo, representa la muerte del agente"""
        self.model.num_agents[self.type] -=1
        self.model.schedule.remove(self)
        self.model.grid.remove_agent(self)
        
    def add_new_cell(self, agent, position):
        """crea una nueva célula"""
        self.model.grid.place_agent(agent, position)
        self.model.schedule.add(agent)
            
    def step(self):
        """movimiento e interacción del agente"""
        self.random_move()
    
    def change_energy(self):
        """cambia su propia energia, temporizador del tiempo del agente 
        dentro del sistema.
        """
        # cellmates = self.model.grid.get_cell_list_contents([self.pos])
        self.energy -= 1
      
        

class PD_L1(Exosome):
    """_Representa a la clase de PD_L1 que interactua con macrófagos y neutrófilos
    """
    def __init__(self, unique_id, pos, model, moore, energy=None) -> None:
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 19
        self.model.num_agents[self.type] +=1
        self.model.protumoral()
    def step(self):
        """Movimiento aleratorio del agente hasta que se le acabe la energía
        """
        
        if self.energy > 0:
            self.random_move()
            self.change_energy()
        else :
            # death
            self.destroy_self()

class MDE(Exosome):
    def __init__(self, unique_id, pos, model, moore, energy=None) -> None:
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 20
        self.model.num_agents[self.type] +=1
        self.model.protumoral()
    def step(self):
        
        if self.energy > 0:
            self.random_move()
            self.change_energy()
            # Vascular permeability increase
            
        else :
            # death
            if(self.model.N_vasc_per< self.model.th_vasc_per):
                self.model.N_vasc_per += self.model.pro_vasc_per
            self.destroy_self()

        

class MPCDE(Exosome):
    def __init__(self, unique_id, pos, model, moore, energy=None) -> None:
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 21
        self.model.num_agents[self.type] +=1
        self.model.protumoral()
    def step(self):
        
        if self.energy > 0:
            self.random_move()
            self.change_energy()
            
        else :
            # death
            if(self.model.N_rec_CD11b< self.model.th_rec_CD11b):
                self.model.N_rec_CD11b += self.model.pro_rec_CD11b
            self.destroy_self()


class MCDEV(Exosome):
    def __init__(self, unique_id, pos, model, moore, energy=None) -> None:
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 22
        self.model.num_agents[self.type] +=1
        self.model.protumoral()
    def step(self):
        
        if self.energy > 0:
            self.random_move()
            self.change_energy()
            
        else :
            # death
            if(self.model.N_vasc_per< self.model.th_vasc_per):
                self.model.N_vasc_per += self.model.pro_vasc_per
            self.destroy_self()

class ODEV(Exosome):
    def __init__(self, unique_id, pos, model, moore, energy=None) -> None:
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 23
        self.model.num_agents[self.type] +=1
        self.model.protumoral()
    def step(self):
        
        if self.energy > 0:
            self.random_move()
            self.change_energy()
            
        else :
            # death
            if(self.model.N_immunosupres< self.model.th_immunosupres):
                self.model.N_immunosupres += self.model.pro_immunosupres
            self.destroy_self()




    
        
class StromalC(Agent):
    """Agente que representa a la células estromales"""
    def __init__(self, unique_id, pos, model, moore, energy=None) -> None:
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 15
        self.model.num_agents[self.type] +=1
        #print("hi: ",self.unique_id, "\n")
    
    def step(self):
        """No realiza nada, su comportamiento está definido en model.py"""
        i = 0
        i += 1

        
class Fibroblast(Agent):
    """Clase que representa a los fibrobastos"""
    def __init__(self, unique_id, pos, model, moore, energy=None) -> None:
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 16
        self.model.num_agents[self.type] +=1
        self.protumoral = False
    def step(self):
        """No realiza nada, su comportamiento está definido en model.py"""
        i = 0
        i += 1
class CTC(Agent):
    """Célula cancerígena circulante"""
    def __init__(self, unique_id, pos, model, moore, energy=None) -> None:
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 17
        self.model.num_agents[self.type] +=1
        self.move = True
        self.colonize = False
        self.model.protumoral()
    def destroy_self(self):
        self.model.num_agents[self.type] -=1
        self.model.schedule.remove(self)
        self.model.grid.remove_agent(self)
    
    def change_energy(self):
        # cellmates = self.model.grid.get_cell_list_contents([self.pos])
        self.energy -= 1
    def empty_spaces(self):
        list_neigh = []
        x , y = self.pos
        for i in range (-1, 2):
            for j in range(-1,2):
                a = (i+x) %self.model.width
                b = (j+y) %self.model.height
                list_neigh.append((a,b))
        list_neigh.remove((x,y))
        for neighbor in self.model.grid.neighbor_iter(self.pos):
            list_neigh.remove(neighbor.pos)
        return list_neigh
    def random_move(self):
        empty = self.empty_spaces()
        if len(empty)>0:
            new_position = random.choice(empty)
            #print(new_position)
            self.model.grid.move_agent(self, new_position)
            
    def step(self):
        """consiste en cuatro fases:
         1. Circula en el sistema
         2. Si encuentra donde adherirce se detienen y permanece latente
         3. Coloniza el sistema
         4. Muere si se le acaba la energía
        """
        
        if self.energy > 0:
            if self.move:
                self.random_move()
                for neighbor in self.model.grid.neighbor_iter(self.pos):
                    if neighbor.type == T_STROMAL or neighbor.type == T_FIBROBLAST: # Interacción with the Stromal Cells.
                        
                        self.move = False
                    elif neighbor.type == T_BMDC:
                        if neighbor.move == False:
                            self.move == False
                    elif neighbor.type == T_CTC:
                        if neighbor.move == False:
                            self.move == False
                self.change_energy()
            # Fibronectin production when the cell is attached
            else:
                if(self.colonize):
                    empty_list = self.empty_spaces()
                    if len(empty_list)>0:
            
                        new_position = random.choice(empty_list)
                        energy =self.model.random.randrange(100)
                        newCancerCell = CTC(self.model.next_id(),new_position,self.model,True,energy)
                        newCancerCell.move =False
                        newCancerCell.colonize = True

                        self.model.grid.place_agent(newCancerCell, new_position)
                        self.model.schedule.add(newCancerCell)
                    self.change_energy()
                else:
                    if(self.model.ctc_col):
                        self.colonize = True
                    self.change_energy()
                # Remember the production of the other molecules. Proinflamatory
                # Angiogenic 
            
        else :
            # death
            self.destroy_self()
        
      


class BMC(Agent):
    """Bone Marrow Cell agent."""

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.move = True
        self.type = 0
        self.model.num_agents[self.type] +=1


    def change_energy(self):
        """Reduce su energía en cada step"""
        # cellmates = self.model.grid.get_cell_list_contents([self.pos])
        self.energy -= 1
    def kill_cell(self, agent):
        """Elimina un agente si se encuentra con el"""
        self.model.num_agents[agent.type] -=1
        self.model.schedule.remove(agent)
        self.model.grid.remove_agent(agent)
        
    def destroy_self(self):
        """muerte celular propia"""
        self.model.num_agents[self.type] -=1
        self.model.schedule.remove(self)
        self.model.grid.remove_agent(self)
        
    def add_new_cell(self, agent, position):
        """Crea un nuevo agente (reclutamiento)"""
        self.model.grid.place_agent(agent, position)
        self.model.schedule.add(agent)
    def empty_spaces(self):
        list_neigh = []
        x , y = self.pos
        for i in range (-1, 2):
            for j in range(-1,2):
                a = (i+x) %self.model.width
                b = (j+y) %self.model.height
                list_neigh.append((a,b))
        list_neigh.remove(self.pos)
        for neighbor in self.model.grid.neighbor_iter(self.pos):
            list_neigh.remove(neighbor.pos)
        return list_neigh
    def random_move(self):
        empty = self.empty_spaces()
        if len(empty)>0:
            new_position = random.choice(empty)
            #print(new_position)
            self.model.grid.move_agent(self, new_position)
            
    def step(self):
        if self.energy >0:
            self.random_move()
        else:
            self.destroy_self()
        self.change_energy()
      
            
            
###################################################################            
            
            
            
class BMDC(BMC):
    """Bone Marrow Derived Cell agent."""

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.type = 1       
        self.energy = energy
        self.move = True
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.model.num_agents[self.type] +=1 
        self.model.protumoral()
            
    def step(self):
        """Step 3 fases:
        1. Circula a travez del medio
        2. Se adhiere
        3. Muere"""
        if self.energy > 0:
            if self.move:
                self.random_move()
                for neighbor in self.model.grid.neighbor_iter(self.pos):
                    if neighbor.type == T_STROMAL or neighbor.type == T_FIBROBLAST:
                        # Interacción with the Stromal Cells
                        # thredshold of fibronectin for the BMDC to attach 
                        self.move = False
                    elif neighbor.type == T_BMDC:
                        if neighbor.move == False:
                            self.move == False
            # Fibronectin production when the cell is attached
            # else:
            #     # Fibronectin production when the cell is attached
            #     if(self.model.N_fibronectin< self.model.th_fibronectin):
            #         self.model.N_fibronectin += self.model.pro_fibronectin
            # # Remember the production of the other molecules. Proinflamatory
            # Angiogenic
        else :
            # death
            self.destroy_self()
        self.change_energy()

        
        
        
class HSC(BMC):
    """Hematopoietic Stem Cell agent."""

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 2
        self.model.num_agents[self.type] +=1             
            
        
 
### ====================================================================================
        
        
class MyeloidCell(HSC):
    """ Myeloid progenitor Cell agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 3
        self.model.num_agents[self.type] +=1            
            
    
    def become_protumoral(self, new_type, pos):
        """Cambia de M1 a M2 o de N1 a N2"""
        self.destroy_self()
        new_cell_type = recruit_new_cell(new_type, self.model.next_id(),self.pos, self.model, self.moore, self.energy)
        self.add_new_cell(new_cell_type,pos)    
         
        

### -------------------------------------------------------------------------------


class LympohidCell(HSC):
    """ Lympohid progenitor Cell agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 4 
        self.model.num_agents[self.type] +=1
            
    def sense_protumoral(self): 
        """La T_CD8 puede ser ihbibida si se encuentra con los siguientes agentes"""  
        count = 0 
        n_t = -1
        for neighbor in self.model.grid.neighbor_iter(self.pos):
            n_t = neighbor.type
            if n_t== T_N2 or n_t == T_M2 or n_t == E_PD_L1:
                count += 1
        return count 
   
### ====================================================================================   
        
        
        
class MDSC(MyeloidCell):
    """ Myeloid-Derived Suppressor agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.TLR = True
        self.type = 5
        self.model.num_agents[self.type] +=1            

    
            
      
      
  
#########################################################################################################################################        
      
 
            
class NK(LympohidCell):
    """ Natural Killer Cell agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.nk = True
        self.type = 6
        self.model.num_agents[self.type] +=1   
        self.model.antitumoral()           
            
    def step(self):  
        """Elimina CTC o cancer"""
        if self.energy >0:
            self.random_move()  
            for neigh in self.model.grid.neighbor_iter(self.pos):
                if neigh.type == T_CTC:
                    # Kill the circuling cancer Cell
                    self.kill_cell(neigh)
                    self.model.antitumoral()   
        else:
            self.destroy_self()
        self.change_energy()      
            

class BCell(LympohidCell):
    """ B Cell agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.b = True
        self.type = 7
        self.model.num_agents[self.type] +=1   
        self.model.antitumoral()                  



class CD8(LympohidCell):
    """ T Cell agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.t = True
        self.type = 8
        self.model.num_agents[self.type] +=1 
        self.model.antitumoral()            
     
            
    def step(self):   
        """Elimina CTC o cancer pero puede ser inhibida"""
        if self.energy >0:
            self.random_move()  
            inhibition_agents = self.sense_protumoral()
            for neigh in self.model.grid.neighbor_iter(self.pos):
                # Kill the cancer cell if the TCD8 is not inhibited
                if neigh.type == T_CTC and inhibition_agents == 0:
                    self.kill_cell(neigh)
                    self.model.antitumoral() 
                if neigh.type == T_CTC and inhibition_agents >0:
                    if(self.model.N_immunosupres< self.model.th_immunosupres):
                        self.model.N_immunosupres += self.model.pro_immunosupres

        else:
            self.destroy_self()
        self.change_energy()
        
        
class CD4(LympohidCell):
    """ T Cell agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.t = True
        self.type = 9
        self.model.num_agents[self.type] +=1
        self.model.antitumoral()           
    
class Treg(LympohidCell):
    """ T Cell agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.t = True
        self.type = 10
        self.model.num_agents[self.type] +=1
        self.model.antitumoral() 
            


### ====================================================================================   

class CD11(MyeloidCell):
    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.type = 24
        self.model.num_agents[self.type] +=1
               
            
    def step(self):        
        if self.energy >0:
            self.random_move()
        else:
            self.destroy_self()
        self.change_energy()     


class M1(MyeloidCell):
    """ Type 1 Macrophage agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.m1 = True
        self.type = 11
        self.model.num_agents[self.type] +=1
        self.model.antitumoral() 
               
            
    def step(self):  
        """ELimna CTC y puede convertise en protumoral"""  
        if self.energy >0:
            self.random_move()
            for neigh in self.model.grid.neighbor_iter(self.pos):
                if neigh.type == E_PD_L1:
                    self.kill_cell(neigh)
                    self.become_protumoral(T_M2, self.pos)
                    break
                if neigh.type == T_CTC:
                    self.kill_cell(neigh)
                    self.model.antitumoral() 
        else:
            self.destroy_self()
        self.change_energy()         
        



class N1(MyeloidCell):
    """ Type 1 Neutrophil agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.n1 = True
        self.type = 12
        self.model.num_agents[self.type] +=1
        self.model.antitumoral() 
              
            
    def step(self): 
        """ELimna CTC y puede convertise en protumoral"""  
        if self.energy >0:
            self.random_move()
            for neigh in self.model.grid.neighbor_iter(self.pos):
                if neigh.type == E_PD_L1:
                    self.kill_cell(neigh)
                    self.become_protumoral(T_N2, self.pos)
                    break
                if neigh.type == T_CTC:
                    self.kill_cell(neigh)
                    self.model.antitumoral() 
        else:
            self.destroy_self()
        self.change_energy()           
        



### -------------------------------------------------------------------------------



class M2(M1):
    """ Type 2 Macrophage agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.pd_l1 = True
        self.type = 13
        self.model.num_agents[self.type] +=1
        self.model.protumoral() 
               
            
    def step(self): 
        "Inhibe a las CD8"       
        if self.energy >0:
            self.random_move()
        else:
            self.destroy_self()
        self.change_energy()    



class N2(N1):
    """ Type 2 Neutrophil agent."""   

    def __init__(self, unique_id, pos, model, moore, energy=None):
        self.unique_id = unique_id
        self.pos =  pos
        self.model = model
        self.moore = moore
        self.energy = energy
        self.pd_l1 = True
        self.type = 14
        self.model.num_agents[self.type] +=1
        self.model.protumoral() 
             
            
    def step(self):     
        "Inhibe a las CD8"   
        if self.energy >0:
            self.random_move()

        else:
            self.destroy_self()
        self.change_energy()    


