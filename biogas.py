''' 
BIOGAS toolkit
Author: Matteo Vegezzi - Politecnico di Milano 

It calculates the parameters for the biogas homesystem  of one household given the Energy demand and the selected substrate
'''

import math

# INPUT DATA

#E_demand = 18 * 6 # MJ daily household Energy demand
#LHV_bio = 33 # MJ Lower Heating Value biogas
#MM_bio = 22 # kg/kmol Molecular Mass biogas
#P = 101.325 # kPa Pressure
#T = 25 + 273.15 # K Temperature
#GY = 0.35 # Gas Yield
#TS = 0.45 # Total Solids 
#rho_sub = 760 # kg/m3 mass density of the substrate
#COD = 0.5 # kg_org/kg Chemical Oxygen Demand Substrate
#HRT = 15 # days Hydraulic Retention Time
#OLR = 1.2 # kg_org/m3 day Organic Loading Rate

# SOLUTION 
def run(LHV_fw, m_fw, LHV_bio, MM_bio, GY, TS, rho_sub, COD, HRT, OLR):
    
    
    P = 101.325 # kPa Pressure
    T = 25 + 273.15 # K Temperature
    
    E_demand = float(m_fw) * float(LHV_fw)
    m_bio = E_demand / float(LHV_bio)
    R = 8.314 / float(MM_bio)
    rho_bio = P / (R * T)
    V_bio = m_bio / rho_bio
    m_sub = V_bio / float(GY) / (float(TS)/100)
    V_sub = m_sub / float(rho_sub)
    m_org = m_sub * float(COD)
    
    V_HRT = V_sub * float(HRT)
    V_OLR = m_org / float(OLR)
    
    V_dig = max(V_HRT,V_OLR)
    V_dig = min(V_dig, 12)
    
    V_sub = math.ceil(V_sub)
    V_dig = math.ceil(V_dig)

    return V_dig, V_sub