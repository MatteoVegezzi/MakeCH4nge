''' 
OPTIMIZER toolkit
Author: Matteo Vegezzi - Politecnico di Milano 
'''

# LIBRARIES


from gurobipy import *
import numpy as np
import math
import csv

def run(filename,HH, m_fw, LHV_fw, eta_fw, stima, ramp, m_lpg, LHV_lpg, eta_lpg, P_lpg, C_lpg, c_lpg, e_lpg, t_change, eta_ics, P_ics, C_ics, c_ics, e_ics, t_ignition, t_collect_ics, MM_bio, LHV_bio, eta_bio, GY, TS, rho_sub, COD, HRT, OLR, V_dig, C_bio, c_bio, e_bio, t_collect_bio, min_lpg, max_lpg, min_ics, max_ics, min_bio, max_bio, MAX_inv, MAX_e, MAX_h, investment_is_a_constrain, obb1, obb2, obb3, W1, W2, W3):
    
    #--------------------TD Constants
    P = 101.325 # kPa Atmosphere Pressure
    R = 8.314 #kJ/kmol K 
    T = 25 # °C Ambient Temperature
    T = T + 273.15 # K 
    
    #System
    
    HH = int(HH) # number of households
    
    # ---------------------OUTPUT DEMAND 
    
    m_fw = float(m_fw) # kg of firewood daily burned by an household
    LHV_fw = float(LHV_fw) # MJ/kg 
    eta_fw = float(eta_fw) # yield three- stones technology
    e_fw = 1.519 # kgCO2/kg 
    t_fw = 8/6 # h/kg 
    
    # ---------------------INPUT 
    
    life = 15 # years useful life
    i_f = 0.061
    i_r = 0.067
    CRF = 0
   
    for n in range(life):
        CRF = CRF + (( 1 + i_f ) ** n) / (( 1 + i_r ) ** n)
     

        

    E_demand = 296489.5166 * 15

    
        
        
    # elif ramp == True:
        
    #     with open(filename, newline="", encoding="ISO-8859-1") as filecsv:
            
    #         lettore = csv.reader(filecsv,delimiter=",")
    #         matrice = list(lettore)
    #         j = 0
    #         E_demand = np.zeros(1441)
    #         for i in range(1,1438):
    #             E_demand[j] = matrice[i][1]
    #             j = j + 1
            
    #         E_demand = (sum(E_demand) * 365 * life)/pow(10,6)
    #         m_fw = E_demand / (LHV_fw * eta_fw * HH * 365 * life)
            
            
    
            
    
    #LPG 
    
    #Technical 
    LHV_lpg = float(LHV_lpg) # MJ/ kg https://www.claverton-energy.com/wordpress/wp-content/uploads/2012/08/the_energy_and_fuel_data_sheet1.pdf
    m_lpg_cylinder = float(m_lpg) #kg 
    eta_lpg = float(eta_lpg) # yield lpg technology
    
    #Economical 
    C_lpg = float(C_lpg) # KSH/cylinder https://totalenergies.ke/products/totalenergies-gas/totalenergies-gas-prices-cylinders-accessories
    c_lpg = float(c_lpg) # KSH/kg 
    #Environmental 
    e_lpg = float(e_lpg) #kgCO2/kg emission factor
    #Social
    P_lpg = float(P_lpg) # kW power by lpg stove
    t_change = float(t_change) #min 
    t_refill = 60
    #ICS
    
    #Technical
    eta_ics = float(eta_ics)# yield ICS technology
    max_fuel = m_fw * eta_ics / eta_fw
    #Economical 
    C_ics = float(C_ics) #KSH/stove
    c_ics = float(c_ics) #KSH/kg
    #Environmental 
    e_ics= float(e_ics) #kgCO2/kg emission factor
    #Social
    t_ignition = float(t_ignition) # min
    P_ics = float(P_ics) # kW power by ics stove
    t_collect = float(t_collect_ics) # min/kg
    
    #BIOGAS 
    
    #Technical 
    MM_bio = float(MM_bio) # kg/kmol
    LHV_bio = float(LHV_bio) # MJ/kg
    GY = float(GY) # m3/kg Gas Yield biogas 
    max_gas = m_fw * LHV_fw * eta_fw/ (LHV_bio * eta_lpg ) # kg/day max biogas outlet for one household
    rho_bio = P/(R*T/MM_bio)
   # V_dig = float(V_dig)
    #Economical 
    C_bio = float(C_bio) # KSH/m3
    c_bio = float(c_bio)#KSH/kg
    #Environmental 
    e_bio= float(e_bio) #kgCO2/kg emission factor
    #Social
    P_bio = 5 #kW 
    TS = float(TS)/100 # % Total Solids
    t_collect_wheel = float(t_collect_bio) #min 
    rho_sub = float(rho_sub) # kg/m3 Substrate density
    
    #--------------------CONSTRAINS
    
    weighted_model = False
    max_cost = float(MAX_inv)# KSH Investment Cost
    max_emission = float(MAX_e) # ton of CO2/ year
    max_hours = float(MAX_h) # h/day HH
    
    max_cylinders = float(max_lpg)
    max_ics = float(max_ics)
    max_biogas = float(max_bio) # m3 of biodigester
    
    min_cylinders = float(min_lpg)
    min_ics = float(min_ics)
    min_biogas = float(min_bio)
    
    # weight definition
    
    W1 = float(W1)
    W2 = float(W2)
    W3 = float(W3)
    
    #priority definition
    
    if obb1 == "Costo Totale":
        p1 = 1
        if obb2 == "Impatto Ambientale":
            p2 = 2
            p3 = 3
        else:
            p2 = 3
            Waux = W2
            W3 = W2
            W3 = Waux
            p3 = 2
            
    elif obb1 == "Impatto Ambientale":
        p1 = 2
        Waux = W2
        W2 = W1
        if obb2 == "Costo Totale":
            p2 = 1
            p3 = 3
            W1 = Waux
        else:
            p2 = 3
            p3 = 1
            W1 = W3
            W3 = Waux
    else:
        p1 = 3 
        Waux = W3
        W3 = W1
        if obb2 == "Impatto Ambientale":
            p2 = 2 
            p3 = 1
            W1 = Waux
        else:
            p2 = 1 
            p3 = 2
            W1 = W2
            W2 = Waux
            
    
    #Energy Demand Constrain coefficients

    
    d = LHV_lpg * eta_lpg
    e = LHV_fw * eta_ics
    f = LHV_bio * eta_lpg
    
    #Cost function coefficients
    
    if not(investment_is_a_constrain) == True :   
        a1 = c_lpg
        b1 = c_ics
        c1 = c_bio
            
    else:
        a1 = C_lpg / (m_lpg_cylinder * 365 * life) + c_lpg
        b1 = C_ics / (max_fuel * 365 * life) + c_ics
        c1 = V_dig * C_bio / (max_gas * 365 * life) + c_bio 
        
    #Emission function coefficients
    a2 = e_lpg
    b2 = e_ics
    c2 = e_bio
        
    #Time function coefficients 
        
    a3 = (LHV_lpg * eta_lpg * 1000)/(P_lpg*60) + t_change / (m_lpg_cylinder * 365 * life) + t_refill * (1/ m_lpg_cylinder - 1 / (m_lpg_cylinder * 365 * life))
    b3 = eta_ics / (m_fw * eta_fw) * 3 * t_ignition + (LHV_fw * eta_ics * 1000)/(P_ics * 60) + t_collect
    c3 = (LHV_bio * eta_lpg * 1000)/(P_bio * 60) +  t_collect_wheel / (GY * rho_bio * rho_sub * 0.1)
    
    """
    WEIGHTED MODEL:
        All objective functions are normalised and summed with user-selected relative weights
    """
        
        # Multi-Objective weighted MODEL 
    if weighted_model == True:
        
        # Normalizing factors
    
        Na = a1 + b1 + c1
        Nb = a2 + b2 + c2
        Nc = a3 + b3 + c3
        
        m = Model()
        x1 = m.addVar(vtype = GRB.INTEGER, name = "x1")
        x2 = m.addVar(vtype = GRB.INTEGER, name = "x2")
        x3 = m.addVar(vtype = GRB.INTEGER, name = "x3")
        
        A = W1*a1/Na + W2*a2/Nb + W3*a3/Nc
        B = W1*b1/Na + W2*b2/Nb + W3*b3/Nc
        C = W1*c1/Na + W2*c2/Nb + W3*c3/Nc
        
            
        m.setObjective(A*x1 + B*x2 + C*x3, GRB.MINIMIZE)
        
        # positivity constrains
        constr1 = m.addConstr(x1 >= 0)
        constr2 = m.addConstr(x2 >= 0)
        constr3 = m.addConstr(x3 >= 0)
        # mix constrains
        constr4 = m.addConstr(x1 / (m_lpg_cylinder * 365 * life) >= min_cylinders)
        constr5 = m.addConstr(x2 * eta_ics / (m_fw * eta_fw) >= min_ics)
        constr6 = m.addConstr(x3 / (max_gas * 365 * life) >= min_biogas)
        
        constr7 = m.addConstr(x1 / (m_lpg_cylinder * 365 * life) <= max_cylinders)
        constr8 = m.addConstr(x2 * eta_ics / (m_fw * eta_fw) <= max_ics)
        constr9 = m.addConstr(x3 / (max_gas * 365 * life) <= max_biogas)
        
        #objective functions constrains
        constr10 = m.addConstr( x1 * a1 + x2 * b1 + x3 * c1 <= max_cost)
        constr11 = m.addConstr((x1 * a2 + x2 * b2 + x3 * c2)/1000 <= max_emission)
        constr12 = m.addConstr((x1 * a3 + x2 * b3 + x3 * c3)/(365*HH*60) <= max_hours)
        
        # Satisfy Energy demand constrain
        constr13 = m.addConstr(d*x1 + e*x2 + f*x3 >= E_demand)
            
        m.optimize()
            
        m.printAttr("X")
    else:
        
        """
        PRIORITISATION MODEL:
        Objective functions are minimised according to a user-defined priority classification
        """
        m = Model()
        x1 = m.addVar(vtype = GRB.INTEGER, name = "x1")
        x2 = m.addVar(vtype = GRB.INTEGER, name = "x2")
        x3 = m.addVar(vtype = GRB.INTEGER, name = "x3")
                
        m.setObjectiveN(a1*x1 + b1*x2 + c1*x3, 0, p1, W1)
        m.setObjectiveN(a2*x1 + b2*x2 + c2*x3, 1, p2, W2)
        m.setObjectiveN(a3*x1 + b3*x2 + c3*x3, 2, p3, W3)
                
        # positivity constrains
        constr1 = m.addConstr(x1 >= 0)
        constr2 = m.addConstr(x2 >= 0)
        constr3 = m.addConstr(x3 >= 0)
        # mix constrains
        constr4 = m.addConstr(x1 / (m_lpg_cylinder * 365 * life) >= min_cylinders)
        constr5 = m.addConstr(x2 / (max_fuel * 365 * life) >= min_ics)
        constr6 = m.addConstr((x3  * V_dig)/ (max_gas * 365 * life) >= min_biogas)
        
        constr7 = m.addConstr(x1 / (m_lpg_cylinder * 365 *life)  <= max_cylinders)
        constr8 = m.addConstr(x2 / (max_fuel * 365 * life) <= max_ics)
        constr9 = m.addConstr((x3  * V_dig)/ (max_gas * 365 * life) <= max_biogas)
        
        #objective functions constrains
        constr10 = m.addConstr( x1 * a1 + x2 * b1 + x3 * c1 <= max_cost)
        constr11 = m.addConstr((x1 * a2 + x2 * b2 + x3 * c2)/1000 <= max_emission)
        constr12 = m.addConstr((x1 * a3 + x2 * b3 + x3 * c3)/(365*HH*60) <= max_hours)
        
        # Satisfy Energy demand constrain
        constr13 = m.addConstr(d*x1 + e*x2 + f*x3 >= E_demand)
                    
        m.optimize()
                    
        m.printAttr("X")
        
    X1 = x1.X
    X2 = x2.X
    X3 = x3.X
    
    ''' Data elaboration '''
    
    # Quantities for Energy source
    n_lpg = math.ceil(X1 / (m_lpg_cylinder * life)) 
    n_ics = math.ceil(X2 / (max_fuel * 365 * life))
    n_bio = math.ceil((X3  * V_dig)/ (max_gas * 365 * life))
    
    
    E_lpg = d * X1
    E_ics = e * X2
    E_bio = f * X3
    
    # Energy share per Energy source
    share_lpg = abs(round(E_lpg/E_demand * 100, 2))
    share_ics = abs(round(E_ics/E_demand * 100, 2))
    share_bio = abs(round(E_bio/E_demand * 100, 2))
    
    # Energy Demand [MJ]
    E_demand = E_demand/(life * 365)
    E_demand = math.ceil(E_demand)
    
    #Investment cost and annual variable cost
    
    Io =  C_lpg * n_lpg + C_ics * n_ics + C_bio * n_bio
    Io = math.ceil(Io)
    c_es = c_lpg * X1 / life + c_ics * X2 / life + c_bio * X3 / life
    c_es = math.ceil(c_es)
    
    #Emission ante/post intervention
    
    e_ante = m_fw * HH * 365 * e_fw / 1000 # ton/year
    e_post  = (X1 * a2 + X2 * b2 + X3 * c2) / (1000 * life) # ton/year
    
    e_ante = math.ceil(e_ante)
    e_post = math.ceil(e_post)
    
    #Hours/day ante/post intervention
    
    h_ante = m_fw * t_fw
    h_post = (X1 * a3 + X2 * b3 + X3 * c3)/(365*HH*60*life)
    
    h_ante = math.ceil(h_ante)
    h_post = math.ceil(h_post)

    
    return n_lpg, n_ics, n_bio, share_lpg, share_ics, share_bio, E_demand, Io, c_es, e_ante, e_post, h_ante, h_post
    
    
    
    


    
    