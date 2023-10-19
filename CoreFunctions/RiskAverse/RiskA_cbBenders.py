#Custom functions
from CoreFunctions.RiskAverse.RA_subproblem import InitSubproblem
from CoreFunctions.masterproblem import InitMasterProblem
from CoreFunctions.extended_new import InitExtended

#General libraries
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import time


#For constraints (1b)
def T_x_sum_1b(data,T_l,T_r,x,i,k_prime,k):
    T_arr = []
    x_arr = []
    for j in data['upstream'][i]:
        if (j,i,k_prime) in data['valid_arcs']:
            T_arr.append(T_l[j,i,k_prime])
            x_arr.append(x[j,i,k_prime])
    for j in data['downstream'][i]:
        if (i,j,k) in data['valid_arcs']:
            T_arr.append(T_r[i,j,k])
            x_arr.append(x[i,j,k])
    return gp.LinExpr(T_arr,x_arr)

#For constraints (1d)
def T_x_sum_1d(data,T,x,i):
    T_arr = []
    x_arr = []
    for j in data['downstream'][i]:
        for k in data['K']:
            if (i,j,k) in data['valid_arcs']:
                T_arr.append(T[i,j,k])
                x_arr.append(x[i,j,k])
    return gp.LinExpr(T_arr,x_arr)

#For constraints (1e)
def T_x_sum(data,T,x,i,k):
    T_arr = []
    x_arr = []
    for j in data['upstream'][i]:
        if (j,i,k) in data['valid_arcs']:
            T_arr.append(T[j,i,k])
            x_arr.append(x[j,i,k])
    return gp.LinExpr(T_arr,x_arr)

#General multiplication for any constraints (RHS)
def pi_h(pi,h):
    return np.dot(np.array([*pi.values()]),np.array([*h.values()]))


# For constraints related to the time components (2a)-(2e) replacing (1f)
def pi_T(pi,T,xi,domain):
    return np.multiply([pi[i,j,k,xi] for (i,j,k) in domain],[*T.values()])
def pi_T_y(pi,T,xi,y,domain):
    return  gp.LinExpr(pi_T(pi,T,xi,domain),[y[i,j,k] for (i,j,k) in domain])

#For risk-averse constraints
def pi_Tik(pi,T,xi,domain):
    return np.multiply([pi[i,k,xi] for (i,k) in domain],[*T.values()])
def pi_T_eta_u(pi,T,xi,eta_u,domain):
    return  gp.LinExpr(pi_Tik(pi,T,xi,domain),[eta_u[i,k] for (i,k) in domain])

def cb_benders(model,where):
    if where == GRB.Callback.POLLING:
            # Ignore polling callback
        pass
    if where == GRB.Callback.MIPSOL:
        #Taking solution from MP and resolving subproblem
        x = model.cbGetSolution(model._x)
        y = model.cbGetSolution(model._y)
        eta_u = model.cbGetSolution(model._eta_u)
        eta_v = model.cbGetSolution(model._eta_v)
        theta_sol = model.cbGetSolution(model._theta)
        mpobj = model.cbGet(GRB.Callback.MIPSOL_OBJBST)

        model._sub.modify_RHS(first_stage=[x,y,eta_u,eta_v])
        model._sub.solve_subproblem(silent=True)
        #print(f'Tomo {time.time()-start} segundos para resolver subproblem!')
        
        subobj = model._sub.m.ObjVal
        #theta = model._theta_sol
        theta = (1/len(model._Omega))*sum(theta_sol.values())
        U_B = mpobj-theta+subobj
        L_B = mpobj
        gap = (U_B-L_B)/U_B
        # print(gap,model._iterations)
        if model._compare_optimal == True:
            model._convergence_stats.append({'Upper bound':U_B,'Lower bound':L_B,'Optimal':model._OptObj})
        else:
            model._convergence_stats.append({'Upper bound':U_B,'Lower bound':L_B})
        
        start = time.time()

        if (model._sub.m.status == GRB.Status.OPTIMAL):
            if (gap > model._set_gap) & (model._iterations < 1000):
                #Getting information for (1b)
                model._pi_x_1b = {}
                for (i,k_prime,k,xi) in model._i_kprime_k_xi:
                    sub_cstr_x_1b = model._sub.m.getConstrByName(f'Flow_balance_{i}_{k_prime}_{k}_{xi}')
                    model._pi_x_1b[i,k_prime,k,xi] = sub_cstr_x_1b.Pi
                
                #Getting information for cuts (1c)
                model._pi_1c = {}
                #Getting information for cuts (1d)
                model._pi_x_1d = {}
                for i in model._V:
                    if i not in model._v_subsets['Retail']:
                        for xi in model._Omega:
                            sub_cstr_x_1d = model._sub.m.getConstrByName(f'Total_flow_capacity_{i}_{xi}')
                            model._pi_x_1d[i,xi] = sub_cstr_x_1d.Pi

                #Getting information for cuts for (1e)
                model._pi_x= {}
                
                model._i_k_xi = []
                for (i,k,xi),d in model._d.items():
                    sub_cstr_x = model._sub.m.getConstrByName(f'Demand_penalty_{i}_{k}_{xi}')
                    model._pi_x[i,k,xi] = sub_cstr_x.Pi
                    model._i_k_xi.append((i,k,xi))
                            
            
                #Lead-time cuts (2a),(2c),(2d)
                model._pi_y_1 ={}
                model._pi_y_2 = {}
                model._pi_y_3 = {}

                #Risk-averse cuts
                model._pi_eta_v = {}


                for (i,j,k) in model._valid_arcs:
                    for xi in model._Omega:
                        # individual Capacity cut (1c)
                        sub_cstr_x_1c = model._sub.m.getConstrByName(f'Capacity_{i}_{j}_{k}_{xi}')
                        model._pi_1c[i,j,k,xi] = sub_cstr_x_1c.Pi

                        #Lead time cuts (2a)
                        sub_cstr_y_1 = model._sub.m.getConstrByName(f'o_y_linear1_{i}_{j}_{k}_{xi}')
                        model._pi_y_1[i,j,k,xi] = sub_cstr_y_1.Pi
                        # Lead time (2c) 
                        sub_cstr_y_2 = model._sub.m.getConstrByName(f'o_y_linear3_{i}_{j}_{k}_{xi}')
                        model._pi_y_2[i,j,k,xi] = sub_cstr_y_2.Pi
                        # Lead time (2d)
                        sub_cstr_y_3 = model._sub.m.getConstrByName(f'o_y_linear4_{i}_{j}_{k}_{xi}')
                        model._pi_y_3[i,j,k,xi] = sub_cstr_y_3.Pi

                        #Risk averse cuts
                        sub_cstr_eta_v = model._sub.m.getConstrByName(f'CVaR(v)_{i}_{j}_{k}_{xi}')
                        model._pi_eta_v[i,j,k,xi] = sub_cstr_eta_v.Pi

                #Risk-averse cuts
                model._pi_eta_u = {}
                
                for (i,k) in model._active_i_k:
                    for xi in model._Omega:
                        #Risk-averse cuts
                        sub_cstr_eta_u = model._sub.m.getConstrByName(f'CVaR(u)_{i}_{k}_{xi}')
                        model._pi_eta_u[i,k,xi] = sub_cstr_eta_u.Pi

                pi_h_x = pi_h(model._pi_x,model._h_x)
                pi_h_y_1 = pi_h(model._pi_y_1,model._h_y_1)
                pi_h_y_2 = pi_h(model._pi_y_2,model._h_y_2)
                pi_h_y_3 = pi_h(model._pi_y_3,model._h_y_3)
                pi_h_1c = pi_h(model._pi_1c,model._h_1c)
                pi_h_x_1b = pi_h(model._pi_x_1b,model._h_x_1b)
                pi_h_x_1d = pi_h(model._pi_x_1d,model._h_x_1d)
                pi_h_eta_u = pi_h(model._pi_eta_u,model._h_eta_u)
                pi_h_eta_v = pi_h(model._pi_eta_v,model._h_eta_v)

                #Adding single cut
                model.cbLazy(
                            gp.quicksum((1/len(model._Omega))*model._theta[xi] for xi in model._Omega) 
                            #+gp.quicksum((1/len(model._Omega))*gp.quicksum(model._pi_x[i,k,xxi]*T_x_sum(data,model._T_x,model._x,i,k) for (i,k,xxi) in model._i_k_xi if xxi == xi) for xi in model._Omega)
                            +gp.quicksum(model._pi_x[i,k,xi]*T_x_sum(model._data,model._T_x,model._x,i,k) for (i,k,xi) in model._i_k_xi)
                            +gp.quicksum((1/len(model._Omega))*pi_T_y(model._pi_y_1,model._T_y_1,xi,model._y,model._valid_arcs) for xi in model._Omega)
                            +gp.quicksum((1/len(model._Omega))*pi_T_y(model._pi_y_2,model._T_y_2,xi,model._y,model._valid_arcs) for xi in model._Omega)
                            +gp.quicksum((1/len(model._Omega))*pi_T_y(model._pi_y_3,model._T_y_3,xi,model._y,model._valid_arcs) for xi in model._Omega)
                            +gp.quicksum((1/len(model._Omega))*pi_T_y(model._pi_1c,model._T_y_1c,xi,model._y,model._valid_arcs) for xi in model._Omega)
                            +gp.quicksum((1/len(model._Omega))*pi_T_y(model._pi_1c,model._T_x_1c,xi,model._x,model._valid_arcs) for xi in model._Omega)
                            +gp.quicksum((1/len(model._Omega))*pi_T_eta_u(model._pi_eta_u,model._T_eta_u,xi,model._eta_u,model._active_i_k) for xi in model._Omega)
                            + gp.quicksum(model._pi_x_1b[i,k_prime,k,xi]*T_x_sum_1b(model._data,model._T_x_1b_l,model._T_x_1b_r,model._x,i,k_prime,k) for (i,k_prime,k,xi) in model._i_kprime_k_xi)
                            +gp.quicksum(model._pi_x_1d[i,xi]*T_x_sum_1d(model._data,model._T_x_1d,model._x,i) for (i,xi) in model._i_xi)
                            +gp.quicksum(pi_T_y(model._pi_eta_v,model._T_eta_v,xi,model._y,model._valid_arcs) for xi in model._Omega)
                            >= 
                            (1/len(model._Omega))*pi_h_y_1
                            +(1/len(model._Omega))*pi_h_y_2
                            +(1/len(model._Omega))*pi_h_y_3
                            +(1/len(model._Omega))*pi_h_1c
                            +(1/len(model._Omega))*pi_h_eta_u
                            +pi_h_x_1b
                            +pi_h_x_1d
                            +pi_h_x
                            +pi_h_eta_v
                            #+gp.quicksum((1/len(mp.m._Omega))*gp.quicksum(mp.m._pi_x[i,k,xxi]*mp.m._h_x[i,k,xxi] for (i,k,xxi) in mp.m._i_k_xi if xxi == xi) for xi in mp.m._Omega)
                            )
            #print(f'Adding cut took {time.time()-start} seconds!')
            model.update()
            model._iterations += 1

def solve_using_benders_decomposition(data,compare_optimal=False,set_gap=0.01,c_e=1.5,alpha=[0.7,0.7,0.7,0.7]):
        
        #Callback ends!!
        
        mp = InitMasterProblem(cb=True)
        #mp.register_callback(mp.m._callback)
        mp.load_data(data)
        mp.set_masterproblem(risk_averse=True)
        mp.m._set_gap = set_gap
        mp.m._iterations = 0
        mp.m._same_solution = 0
        #mp.solve_masterproblem(silent=True)

        mp.m._x_sols,mp.m._y_sols = [],[]
        mp.m._eta_u_sols,mp.m._eta_v_sols = [],[]
        
        U_B = np.infty
        L_B = -np.infty
        gap = np.inf

        mp.m._compare_optimal = compare_optimal
        #BEnders decomposition
        mp.m._leadtime,mp.m._linearized1,mp.m._fixedPenalty,mp.m._linearized2 = True,True,False,True
        extended = InitExtended()
        extended.load_data(data)
        extended.set_extended(risk_averse=True,alpha=alpha,
                              leadtime=mp.m._leadtime,linearized1=mp.m._linearized1,
                              fixedPenalty=mp.m._fixedPenalty,linearized2=mp.m._linearized2,c_e=c_e)
        extended_time = 0
        if compare_optimal == True:
            start = time.time()
            mp.m._leadtime,mp.m._linearized1,mp.m._fixedPenalty,mp.m._linearized2 = True,True,False,True
            ext = InitExtended()
            ext.load_data(data)
            ext.set_extended(risk_averse=True,alpha=alpha,
                                leadtime=mp.m._leadtime,linearized1=mp.m._linearized1,
                                fixedPenalty=mp.m._fixedPenalty,linearized2=mp.m._linearized2,c_e=c_e)
            ext.solve_extended(silent=True)
            mp.m._x_opt = ext.m._x_sol
            mp.m._y_opt = ext.m._y_sol
            mp.m._OptObj = ext.m.ObjVal
            extended_time = time.time()-start
            print(f'Extended formulation took {extended_time} seconds')
        extended.m.update()
        print('Benders decomposition begins')
        #Forming data
        master_start = time.time()
        start = time.time()

        #flow balance (1b)
        mp.m._i_kprime_k_xi = []
        mp.m._h_x_1b,mp.m._T_x_1b_l,mp.m._T_x_1b_r = {},{},{}
        for i in mp.m._V:
            if i not in mp.m._v_subsets['Retail']:
                for ii,k in mp.m._directed_output:
                    if ii == i:
                        if i in mp.m._v_subsets['Dist']:
                            k_primes = [k]
                            rs = [1]
                        else:
                            conversion_df_i_k = mp.m._prod_conv[mp.m._prod_conv['downstream']==k]
                            rs = list(conversion_df_i_k['Amount'])
                            k_primes = list(conversion_df_i_k['upstream'])

                        for k_prime,r in zip(k_primes,rs):
                            for xi in mp.m._Omega:
                                mp.m._i_kprime_k_xi.append((i,k_prime,k,xi))
                                mp.m._ext_cstr_x_1b = extended.m.getConstrByName(f'Flow_balance_{i}_{k_prime}_{k}_{xi}')
                                mp.m._h_x_1b[i,k_prime,k,xi] = mp.m._ext_cstr_x_1b.RHS

                                for j in extended.m._upstream[i]:
                                    if (j,i,k_prime) in data['valid_arcs']:
                                        mp.m._T_x_1b_l[j,i,k_prime] = extended.m.getCoeff(mp.m._ext_cstr_x_1b,extended.m._x[j,i,k_prime])
                                for j in data['downstream'][i]:
                                    if (i,j,k) in data['valid_arcs']:
                                        mp.m._T_x_1b_r[i,j,k] = extended.m.getCoeff(mp.m._ext_cstr_x_1b,extended.m._x[i,j,k])

        

        #mixed flow capacity (1d)    
        mp.m._i_xi = []
        mp.m._h_x_1d,mp.m._T_x_1d = {},{}  
        for i in mp.m._V:
            if i not in mp.m._v_subsets['Retail']:
                for xi in mp.m._Omega:
                    mp.m._i_xi.append((i,xi))
                    mp.m._ext_cstr_x_1d = extended.m.getConstrByName(f'Total_flow_capacity_{i}_{xi}')
                    mp.m._h_x_1d[i,xi] = mp.m._ext_cstr_x_1d.RHS
                    for j in extended.m._downstream[i]:
                        for k in extended.m._K:
                            if (i,j,k) in extended.m._valid_arcs:
                                mp.m._T_x_1d[i,j,k] = extended.m.getCoeff(mp.m._ext_cstr_x_1d,extended.m._x[i,j,k])


        #Unmet demand computation (1e)
        mp.m._h_x,mp.m._T_x = {},{}
        for (i,k,xi),d in extended.m._d.items():
            mp.m._ext_cstr_x = extended.m.getConstrByName(f'Demand_penalty_{i}_{k}_{xi}')
            mp.m._h_x[i,k,xi] = mp.m._ext_cstr_x.RHS
            for j in extended.m._upstream[i]:
                if (j,i,k) in extended.m._valid_arcs:
                    mp.m._T_x[j,i,k] = extended.m.getCoeff(mp.m._ext_cstr_x, extended.m._x[j,i,k])
        
        # individual flow capacity (1c)
        mp.m._h_1c,mp.m._T_x_1c,mp.m._T_y_1c = {},{},{}
        #Lead-time cuts (2a),(2c),(2d)
        mp.m._h_y_1,mp.m._h_y_2,mp.m._h_y_3 = {},{},{}
        mp.m._T_y_1, mp.m._T_y_2,mp.m._T_y_3 = {},{},{}
        #Risk averse cuts
        mp.m._h_eta_v = {}
        mp.m._T_eta_v = {}
        for (i,j,k) in extended.m._valid_arcs:
            for xi in extended.m._Omega:
                #Individual flow capacity (1c)
                mp.m._ext_cst_1c = extended.m.getConstrByName(f'Capacity_{i}_{j}_{k}_{xi}')
                mp.m._h_1c[i,j,k,xi] = mp.m._ext_cst_1c.RHS
                mp.m._T_x_1c[i,j,k] = extended.m.getCoeff(mp.m._ext_cst_1c,extended.m._x[i,j,k])
                mp.m._T_y_1c[i,j,k] = extended.m.getCoeff(mp.m._ext_cst_1c,extended.m._y[i,j,k])

                #Lead-time cuts info
                mp.m._ext_cstr_y_1 = extended.m.getConstrByName(f'o_y_linear1_{i}_{j}_{k}_{xi}')
                mp.m._h_y_1[i,j,k,xi] = mp.m._ext_cstr_y_1.RHS
                mp.m._T_y_1[i,j,k] = extended.m.getCoeff(mp.m._ext_cstr_y_1, extended.m._y[i,j,k])
                
                mp.m._ext_cstr_y_2 = extended.m.getConstrByName(f'o_y_linear3_{i}_{j}_{k}_{xi}')
                mp.m._h_y_2[i,j,k,xi] = mp.m._ext_cstr_y_2.RHS
                mp.m._T_y_2[i,j,k] = extended.m.getCoeff(mp.m._ext_cstr_y_2, extended.m._y[i,j,k])
                
                mp.m._ext_cstr_y_3 = extended.m.getConstrByName(f'o_y_linear4_{i}_{j}_{k}_{xi}')
                mp.m._h_y_3[i,j,k,xi] = mp.m._ext_cstr_y_3.RHS
                mp.m._T_y_3[i,j,k] = extended.m.getCoeff(mp.m._ext_cstr_y_3, extended.m._y[i,j,k])

                #Risk-averse cuts
                mp.m._ext_cstr_eta_v = extended.m.getConstrByName(f'CVaR(v)_{i}_{j}_{k}_{xi}')
                mp.m._h_eta_v[i,j,k,xi] = mp.m._ext_cstr_eta_v.RHS
                mp.m._T_eta_v[i,j,k] = extended.m.getCoeff(mp.m._ext_cstr_eta_v,extended.m._eta_v[i,j,k])

        mp.m._h_eta_u = {}
        mp.m._T_eta_u = {}
        for (i,k) in mp.m._active_i_k:
            for xi in mp.m._Omega:
                mp.m._ext_cstr_eta_u = extended.m.getConstrByName(f'CVaR(u)_{i}_{k}_{xi}')
                mp.m._h_eta_u[i,k,xi] = mp.m._ext_cstr_eta_u.RHS
                mp.m._T_eta_u[i,k] = extended.m.getCoeff(mp.m._ext_cstr_eta_u,extended.m._eta_u[i,k]) 
        
        #print(f'Tomo {time.time()-start} segundos para obtener coeficientes!')
        #print('--------------')
        mp.m._sub = InitSubproblem()
        mp.m._sub.load_data(data)
        initial_x = {}
        initial_y = {}
        initial_eta_u = {}
        initial_eta_v = {}
        for (i,j,k) in mp.m._valid_arcs:
            initial_x[i,j,k] = 1000
            initial_y[i,j,k] = 1
            initial_eta_v[i,j,k] = 0
        for (i,k) in mp.m._active_i_k:
            initial_eta_u[i,k] = 0

        mp.m._sub.set_subproblem(
                            first_stage=[initial_x,initial_y,initial_eta_u,initial_eta_v],
                            leadtime=mp.m._leadtime,linearized1=mp.m._linearized1,
                            fixedPenalty=mp.m._fixedPenalty,linearized2=mp.m._linearized2,
                            subgradient=False,c_e=c_e,alpha=alpha)
        
        mp.register_callback(cb_benders)
        mp.solve_masterproblem(silent=True)
        mp.m._gap = np.infty

        #FOr graphing purposes adding last iteration of algorithm
        mpobj = mp.m.ObjVal

        mp.m._sub.modify_RHS(first_stage=[mp.m._x_sol,mp.m._y_sol,mp.m._eta_u_sol,mp.m._eta_v_sol])
        mp.m._sub.solve_subproblem(silent=True)
        #print(f'Tomo {time.time()-start} segundos para resolver subproblem!')
        
        subobj = mp.m._sub.m.ObjVal
        #theta = model._theta_sol
        theta = (1/len(mp.m._Omega))*sum(mp.m._theta_sol.values())
        U_B = mpobj-theta+subobj
        L_B = mpobj
        gap = (U_B-L_B)/U_B
        if mp.m._compare_optimal == True:
            mp.m._convergence_stats.append({'Upper bound':U_B,'Lower bound':L_B,'Optimal':mp.m._OptObj})
        else:
            mp.m._convergence_stats.append({'Upper bound':U_B,'Lower bound':L_B})
        
        benders_time = time.time()-master_start
        #print(f'Benders took {benders_time} seconds')
        return mp,extended_time,benders_time


        