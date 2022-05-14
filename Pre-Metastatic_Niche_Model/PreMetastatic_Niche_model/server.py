# Alberto Josué Ortiz Rosales
# 18 feb 2022


from PreMetastatic_Niche_model.agents import E_MCDEV, E_MDE, E_MPCDE, E_ODEV
from mesa.visualization.ModularVisualization import ModularServer
from mesa.visualization.modules import CanvasGrid, ChartModule, TextElement
from mesa.visualization.UserParam import UserSettableParameter

from  .model import PreMetastaticNicheModel

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

# Exosomes types id
E_PD_L1= 19

def PMN_draw(agent):
    """
    Portrayal Method for canvas
    """
    if agent is None:
        return
    portrayal = {"Shape": "rect", "w":1, "h": 1, "Filled": "true", "Layer": 0}

    
    
    if agent.type == T_BMDC:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/BMDC.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "BMDC"
           
    if agent.type == T_BCELL:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/BCELL.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "B"
    if agent.type == T_CD4:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/CD4.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "T_CD4"
    if agent.type == T_CD8:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/CD8.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "T_CD8"
    if agent.type == T_M1:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/M1.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "M1"
    if agent.type == T_N1:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/N1.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "N1"
    if agent.type == T_M2:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/M2.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "M2"
    if agent.type == T_N2:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/N2.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "N2"
    if agent.type == T_CTC:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/CTC_.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "CTC"  
    if agent.type == T_STROMAL:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/STRO_.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "Stromal"  
    if agent.type == T_FIBROBLAST:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/FIB_.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "Fibrobast"
    if agent.type == E_PD_L1:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/PDL1_.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "PD_L1"      
    
    if agent.type == E_ODEV:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/EXO.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "ODEV" 
    if agent.type == E_MDE:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/EV.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "MDE" 
    if agent.type == E_MCDEV:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/EXO.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "MCDEV" 
    if agent.type == E_MPCDE:
        portrayal["Shape"] = "PreMetastatic_Niche_model/img/EV.png"
        portrayal["scale"] = 1.5
        portrayal["type"] = "MPCDE"                                  
    else:
        portrayal["Color"] = ["#000000", "#000000"]
        portrayal["stroke_color"] = "#000000"
    return portrayal

#size of the grid

height_canvas = 100
width_canvas = 100
canvas_element = CanvasGrid(PMN_draw, height_canvas, width_canvas, 500, 500)


e_chart = ChartModule(
    [{"Label": "BMDC", "Color": "#641E16"}, 
    {"Label": "B", "Color": "#58D68D"},
    {"Label": "T_CD8", "Color": "#AF601A"},
    {"Label": "T_CD4", "Color": "#28B463"},
    {"Label": "M1", "Color": "#2980B9"},
    {"Label": "N1", "Color": "#138D75"},
    {"Label": "M2", "Color": "#A93226"},
    {"Label": "N2", "Color": "#8E44AD"},
    {"Label": "CTC", "Color":  "#221B01"},
    {"Label": "PDL1", "Color":  "#D35400"},]
)
i_chart = ChartModule(
    [
    {"Label": "ISING", "Color":  "#D35400"},]
)

m_chart = ChartModule(
    [
    {"Label": "IL_12", "Color":  "#D35400"},
    {"Label": "fibronectin", "Color":  "#221B01"},
    {"Label": "TGFB", "Color":  "#D35400"},
    {"Label": "S100", "Color": "#2980B9"},
    {"Label": "IL6", "Color": "#138D75"},
    {"Label": "VEGF", "Color": "#8E44AD"},]
)
step_slider = 10
model_params = {
    
    "height": height_canvas,
    "width": width_canvas,
    
    
    "R_CTC":UserSettableParameter("slider", "CTC Recruitment rate", 0.5, 0, 1, 0.1),
    "R_BMDC": UserSettableParameter("slider", "BMDC Recruitment rate",  0.5, 0, 1, 0.1),
    "R_CD8": UserSettableParameter("slider", "TCD8 Rec1ruitment rate",  0.5, 0, 1, 0.1),
    "R_M1": UserSettableParameter("slider", "M1 Recruitment rate",  0.5, 0, 1, 0.1),
    "R_N1": UserSettableParameter("slider", "N1 Recruitment rate",  0.5, 0, 1, 0.1),
    "R_MCDEV": UserSettableParameter("slider", "MCDEV Recruitment rate",  0.5, 0, 1, 0.1),
    "R_MDE": UserSettableParameter("slider", "MDE Recruitment rate",  0.5, 0, 1, 0.1),
    "R_PDL1": UserSettableParameter("slider", "PDL1 Recruitment rate",  0.5, 0, 1, 0.1),
    "R_ODEV" :UserSettableParameter("slider", "ODEV Recruitment rate",  0.5, 0, 1, 0.1),
    "R_MPCDE":UserSettableParameter("slider", "MPCDEA Recruitment rate",  0.5, 0, 1, 0.1),

    "R_IL_12":UserSettableParameter("slider", "IL_12 production rate",  0.5, 0, 1, 0.1),
    "R_fibronectin":UserSettableParameter("slider", "fibronectin production  rate",  0.5, 0, 1, 0.1),
    "R_TGFB":UserSettableParameter("slider", "TGFB production  rate",  0.5, 0, 1, 0.1),
    "R_S100":UserSettableParameter("slider", "S100 production  rate",  0.5, 0, 1, 0.1), 
    "R_IL6" :UserSettableParameter("slider", "IL_6 production  rate",  0.5, 0, 1, 0.1),
    "R_VEGF" :UserSettableParameter("slider", "VEGF production  rate",  0.5, 0, 1, 0.1),
    "Stromal_init" :UserSettableParameter("slider", "Stromal init value",  10, 0, 100, step_slider),
    "Fibroblast_init" :UserSettableParameter("slider", "Fibroblast init value",  10, 0, 100, step_slider),

}

server = ModularServer(
    PreMetastaticNicheModel, [canvas_element,  e_chart, i_chart, m_chart], "PMN", model_params
)
