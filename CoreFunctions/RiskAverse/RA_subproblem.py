import gurobipy as gp
from gurobipy import GRB
import time


class InitSubproblem:
    def __init__(self):
        self.m = gp.Model()
        self.m.setParam("Threads", 8)
        #self.m.setParam("MIPFocus",2)
        

    def load_data(self,data):
        #preloading relevant parameters to model structure
        self.m._Omega,self.m._K,self.m._V,self.m._valid_arcs,self.m._active_i_k,self.m._h = data['Omega'],data['K'],data['V'],data['valid_arcs'],data['active_i_k'],data['h']
        self.m._c_f = data['c'] #Flow cost parameter
        self.m._rho_d = data['rho_d']
        self.m._r = data['r']
        self.m._q_ind,self.m._d,self.m._q_mix = data['q_ind'],data['d'],data['q_mix']
        self.m._v_subsets = data['v_subsets']
        self.m._l,self.m._t,self.m._rho_l_u,self.m._rho_l_f = data['l'],data['t'],data['unitLatePenalty'],data['fixLatePenalty']
        self.m._active_prods = data['active_prods']
        self.m._bigM = data['bigM']


        #Parameters for OneShot Formulation\
        self.m._upstream,self.m._downstream = data['upstream'],data['downstream'] #Upstream and downstream sets for each entity
        self.m._directed_output,self.m._prod_conv = data['directed_output'],data['prod_conv']

        #PArameters for risk-averse formulation
        self.m._lambda = data['lambda']

        #Precomputed constraint indices
        self.m._cst_3c,self.m._cst_3d_1,self.m._cst_3d_2 = data['cst_3c'],data['cst_3d_1'],data['cst_3d_2']
 

    def set_subproblem(self,type=1,mu=0.5,c_e=1.5,first_stage=list,leadtime=True,linearized1=True,fixedPenalty=True,linearized2=True,out_of_sample=False,subgradient=False,alpha=[0.7,0.7,0.7,0.7]):
        self.leadtime = leadtime
        self.fixedPenalty = fixedPenalty 
        self.m._subgradient = subgradient
        #First stage solution parsing
        self.m._x = first_stage[0]
        self.m._y = first_stage[1]
        self.m._eta_u = first_stage[2]
        self.m._eta_v = first_stage[3]

        #SEtting risk attitudes
        alpha_map = {'Part':alpha[0],'Manuf':alpha[1],'Dist':alpha[2],'Retail':alpha[3]}
        alpha_tier = {}
        for i in self.m._V:
            if i in self.m._v_subsets['Part']:
                alpha_tier[i] = alpha_map['Part']
            elif i in self.m._v_subsets['Manuf']:
                alpha_tier[i] = alpha_map['Manuf']
            elif i in self.m._v_subsets['Dist']:
                alpha_tier[i] = alpha_map['Dist']
            elif i in self.m._v_subsets['Retail']:
                alpha_tier[i] = alpha_map['Retail']
            

        ###Setting decision variables dictionaries
        start = time.time()
        #Risk neutral variables
        #Emergency flows
        self.m._z = {}
        #Unmet demand variables
        self.m._u = {}
        #Arrival times
        self.m._a,self.m._o,self.m._v,self.m._m1 = {},{},{},{}

        self.m._alpha_u = alpha_tier
        self.m._alpha_v = alpha_tier
        self.m._s_u = self.m.addVars(self.m._active_i_k,self.m._Omega,name='s(u)',lb=0,ub=GRB.INFINITY)
        self.m._s_v = self.m.addVars(self.m._valid_arcs,self.m._Omega,name='s(v)',lb=0,ub=GRB.INFINITY)

        ###Registering the gurobi variables
        # self.m._s_u = self.m.addVars(self.m._active_i_k,self.m._Omega,name='s(u)',lb=0,ub=GRB.INFINITY)
        # self.m._s_w = self.m.addVars(self.m._valid_arcs,self.m._Omega,name='s(w)',lb=0,ub=GRB.INFINITY)
        # self.m._u = self.m.addVars(self.m._active_i_k,self.m._Omega, name='unmetDemand_penalty') #demand penalty term
        # self.m._a = self.m.addVars(self.m._valid_arcs,self.m._Omega, name='arrival_time') #arrival time variable
        # self.m._o = self.m.addVars(self.m._active_i_k,self.m._Omega, name='start_time') #activie start time variable
        # self.m._w = self.m.addVars(self.m._valid_arcs,self.m._Omega, name='Total_flow_penalty') #term to compute the whole lateness penalty
        # self.m._z = self.m.addVars(self.m._valid_arcs,self.m._Omega, name='deliveryLateness') #delivery lateness
        # self.m._m1 = self.m.addVars(self.m._valid_arcs,self.m._Omega, name='o_y_product')
        for xi in self.m._Omega:
            for (i,k) in self.m._active_i_k:
                self.m._o[i,k,xi] = self.m.addVar(name=f'start_time_{i}_{k}_{xi}') #earliest start time variable
                self.m._u[i,k,xi] = self.m.addVar(name='unmetDemand_penalty') #unmet demand variable

            for (i,j,k) in self.m._valid_arcs:
                self.m._z[i,j,k,xi] = self.m.addVar(name=f'emergency_flow_{i}_{j}_{k}_{xi}')
                self.m._a[i,j,k,xi] = self.m.addVar(name=f'arrival_time_{i}_{j}_{k}_{xi}') #arrival time variable
                
                self.m._m1[i,j,k,xi] = self.m.addVar(name=f'o_y_product_{i}_{j}_{k}_{xi}')
                self.m._v[i,j,k,xi] = self.m.addVar(name=f'delivery_lateness_{i}_{j}_{k}_{xi}')

        #print(f'variable setting took {time.time()-start} seconds')
        ###Objective Function
        if out_of_sample == True:
            self.m._Omega = [0]
        
        self.m.setObjective((1/len(self.m._Omega))*gp.quicksum(c_e*self.m._c_f[i,j,k]*self.m._z[i,j,k,xi] for (i,j,k) in self.m._valid_arcs for xi in self.m._Omega)
                                    +   (
                                        gp.quicksum((1-self.m._lambda[i,k])*((1/((1-self.m._alpha_u[i])*len(self.m._Omega)))*gp.quicksum(self.m._s_u[i,k,xi] for xi in self.m._Omega)) for i,k in self.m._active_i_k)
                                        + gp.quicksum(self.m._lambda[i,k]*((1/((1-self.m._alpha_v[i])*len(self.m._Omega)))*gp.quicksum(self.m._s_v[i,j,k,xi] for xi in self.m._Omega)) for i,j,k in self.m._valid_arcs)
                                        )
                            ,GRB.MINIMIZE)
        #print(f'Objective setting took {time.time()-start} seconds')
            

        #Constraints (1b), flow balance considering multi-product structure
        for i in self.m._V:
            if i not in self.m._v_subsets['Retail']:
                for ii,k in self.m._directed_output:
                        if ii == i:
                            if i in self.m._v_subsets['Dist']:
                                k_primes = [k]
                                rs = [1]
                            else:
                                conversion_df_i_k = self.m._prod_conv[self.m._prod_conv['downstream']==k]
                                rs = list(conversion_df_i_k['Amount'])
                                k_primes = list(conversion_df_i_k['upstream'])

                            for k_prime,r in zip(k_primes,rs):
                                for xi in self.m._Omega:
                                    self.m.addConstr((1/r)*gp.quicksum(self.m._z[j,i,k_prime,xi] for j in self.m._upstream[i] if (j,i,k_prime) in self.m._valid_arcs) 
                                                    - gp.quicksum(self.m._z[i,j,k,xi] for j in self.m._downstream[i] if (i,j,k) in self.m._valid_arcs)
                                                    >= -(1/r)*gp.quicksum(self.m._x[j,i,k_prime] for j in self.m._upstream[i] if (j,i,k_prime) in self.m._valid_arcs)
                                                       +gp.quicksum(self.m._x[i,j,k] for j in self.m._downstream[i] if (i,j,k) in self.m._valid_arcs)
                                                    ,name=f'Flow_balance_{i}_{k_prime}_{k}_{xi}')
            
        #Constraints (1c) for flow capacity
        for (i,j,k) in self.m._valid_arcs:
            for xi in self.m._Omega:
                self.m.addConstr(self.m._z[i,j,k,xi] <= self.m._q_ind[i,j,k]*self.m._y[i,j,k] - self.m._x[i,j,k],name=f'Capacity_{i}_{j}_{k}_{xi}')
        #start = time.time()
        
        #Constraints (1d) for total mixed-flow capacity
        for i in self.m._V:
            if i not in self.m._v_subsets['Retail']:
                for xi in self.m._Omega:
                    self.m.addConstr(gp.quicksum(self.m._z[i,j,k,xi] for j in self.m._downstream[i] for k in self.m._K if (i,j,k) in self.m._valid_arcs)
                                                <= self.m._q_mix[i] - gp.quicksum(self.m._x[i,j,k] for j in self.m._downstream[i] for k in self.m._K if (i,j,k) in self.m._valid_arcs),name=f'Total_flow_capacity_{i}_{xi}')
        
        ### Demand penalty constraints (1e)
        for (i,k,xi),d in self.m._d.items():
            self.m.addConstr((self.m._u[i,k,xi] + gp.quicksum(self.m._z[j,i,k,xi] for j in self.m._upstream[i] if (i,j,k) in self.m._valid_arcs) >= d - gp.quicksum(self.m._x[j,i,k] for j in self.m._upstream[i] if (i,j,k) in self.m._valid_arcs) ),name=f'Demand_penalty_{i}_{k}_{xi}')
        #print(f'2b Constraints setting took {time.time()-start} seconds')

        #Unmet demand risk-averse constraints
        for (i,k) in self.m._active_i_k:
            for xi in self.m._Omega:
                self.m.addConstr(self.m._s_u[i,k,xi] - (self.m._rho_d[i,k]*self.m._u[i,k,xi])  >= - self.m._eta_u[i,k],name=f'CVaR(u)_{i}_{k}_{xi}')

        
        for xi in self.m._Omega:
            for (i,j,k) in self.m._valid_arcs: 
                #Constraints (3e)
                if (j,k) in self.m._active_i_k:
                    
                    coefs = [1,-1]
                    varbs = [self.m._a[i,j,k,xi],self.m._v[i,j,k,xi]]
                    self.m.addConstr((gp.LinExpr(coefs,varbs) <= self.m._t[j,k]),name=f'FlowLatenessComputation_{i}_{j}_{k}_{xi}')
            
                #Block of linearized contraints replacing (3b)
                
                #self.m.addConstr((self.m._m1[i,j,k,xi] <= self.m._bigM*self.m._y[i,j,k]),name=f'o_y_linear1_{i}_{j}_{k}_{xi}')
                #self.m.addConstr((self.m._m1[i,j,k,xi] - self.m._bigM*self.m._y[i,j,k] <= 0 ),name=f'o_y_linear1_{i}_{j}_{k}_{xi}')
                coefs = [1]
                varbs = [self.m._m1[i,j,k,xi]]
                self.m.addConstr((gp.LinExpr(coefs,varbs) <= self.m._bigM[xi]*self.m._y[i,j,k]),name=f'o_y_linear1_{i}_{j}_{k}_{xi}')

                #self.m.addConstr((self.m._m1[i,j,k,xi] <= self.m._o[i,k,xi]),name=f'o_y_linear2_{i}_{j}_{k}_{xi}')
                #self.m.addConstr((self.m._m1[i,j,k,xi] - self.m._o[i,k,xi] <= 0 ),name=f'o_y_linear2_{i}_{j}_{k}_{xi}')
                coefs = [1,-1]
                varbs = [self.m._m1[i,j,k,xi],self.m._o[i,k,xi]]
                self.m.addConstr((gp.LinExpr(coefs,varbs) <= 0 ),name=f'o_y_linear2_{i}_{j}_{k}_{xi}')

                #self.m.addConstr((self.m._m1[i,j,k,xi] >= self.m._o[i,k,xi] - ((1-self.m._y[i,j,k])*self.m._bigM)),name=f'o_y_linear3_{i}_{j}_{k}_{xi}')
                #self.m.addConstr((self.m._m1[i,j,k,xi] - self.m._o[i,k,xi] >=  -self.m._bigM+self.m._bigM*self.m._y[i,j,k]),name=f'o_y_linear3_{i}_{j}_{k}_{xi}')
                coefs = [1,-1]
                varbs = [self.m._m1[i,j,k,xi],self.m._o[i,k,xi]]
                self.m.addConstr((gp.LinExpr(coefs,varbs) >=  -self.m._bigM[xi]+self.m._bigM[xi]*self.m._y[i,j,k]),name=f'o_y_linear3_{i}_{j}_{k}_{xi}')

                #self.m.addConstr((self.m._a[i,j,k,xi] >= self.m._l[i,j,k,xi]*self.m._y[i,j,k] + self.m._m1[i,j,k,xi]),name=f'o_y_linear4_{i}_{j}_{k}_{xi}')
                #self.m.addConstr((self.m._a[i,j,k,xi] - self.m._m1[i,j,k,xi] >= self.m._l[i,j,k,xi]*self.m._y[i,j,k] ),name=f'o_y_linear4_{i}_{j}_{k}_{xi}')
                coefs = [1,-1]
                varbs = [self.m._a[i,j,k,xi],self.m._m1[i,j,k,xi]]
                self.m.addConstr((gp.LinExpr(coefs,varbs) >= self.m._l[i,j,k,xi]*self.m._y[i,j,k]),name=f'o_y_linear4_{i}_{j}_{k}_{xi}')
                
                #lateness risk averse constraints
                coefs = [1,-self.m._rho_l_u[i,j,k]]
                varbs = [self.m._s_v[i,j,k,xi],self.m._v[i,j,k,xi]]
                self.m.addConstr((gp.LinExpr(coefs,varbs) >= -self.m._eta_v[i,j,k] ),name=f'CVaR(v)_{i}_{j}_{k}_{xi}')

                #Constraints (3f)
                #self.m.addConstr((self.m._w[i,j,k,xi] >= self.m._rho_l_u[i,j,k] * self.m._z[i,j,k,xi]),name=f'FlowLatenessPenalty_{i}_{j}_{k}_{xi}')
                #self.m.addConstr((self.m._w[i,j,k,xi] - self.m._rho_l_u[i,j,k] * self.m._z[i,j,k,xi] >= 0 ),name=f'FlowLatenessPenalty_{i}_{j}_{k}_{xi}')
                #coefs = [1,-self.m._rho_l_u[i,j,k]]
                #varbs = [self.m._w[i,j,k,xi],self.m._z[i,j,k,xi]]
                #self.m.addConstr((gp.LinExpr(coefs,varbs) >= 0 ),name=f'FlowLatenessPenalty_{i}_{j}_{k}_{xi}')

        #print(f'One-loop Constraints setting took {time.time()-start} seconds')

        
        #Constraints (3c)
        for i,k,xi in self.m._cst_3c:
            self.m.addConstr(self.m._o[i,k,xi] == 0 ,name=f'Last_tier_start_time_{i}_{k}_{xi}')
        start = time.time()
        #Constraints (3d)1
        for i,j,k,xi in self.m._cst_3d_1:
           #self.m.addConstr((self.m._o[j,k,xi] >= self.m._a[i,j,k,xi]),name=f'Initial_processing_time_non_transformation_{i}_{j}_{k}_{xi}')
           self.m.addConstr((self.m._o[j,k,xi] - self.m._a[i,j,k,xi] >= 0 ),name=f'Initial_processing_time_non_transformation_{i}_{j}_{k}_{xi}')
        for i,j,k,k_sub,xi in self.m._cst_3d_2:
           #self.m.addConstr((self.m._o[j,k,xi] >= self.m._a[i,j,k_sub,xi]),name=f'Initial_processing_time_transformation_{i}_{j}_{k_sub}_{xi}')
           self.m.addConstr((self.m._o[j,k,xi] - self.m._a[i,j,k_sub,xi] >= 0 ),name=f'Initial_processing_time_transformation_{i}_{j}_{k_sub}_{xi}')
        #print(f'Complex? Constraints setting took {time.time()-start} seconds')
        self.m.update()

    def modify_RHS(self,first_stage=None):
        #updating first-stage varibales
        self.m._x = first_stage[0]
        self.m._y = first_stage[1]
        self.m._eta_u = first_stage[2]
        self.m._eta_v = first_stage[3]

        
        n_RHS = []
        cstrs = []
        #Updating constraints 1b
        for i in self.m._V:
            if i not in self.m._v_subsets['Retail']:
                for ii,k in self.m._directed_output:
                        if ii == i:
                            if i in self.m._v_subsets['Dist']:
                                k_primes = [k]
                                rs = [1]
                            else:
                                conversion_df_i_k = self.m._prod_conv[self.m._prod_conv['downstream']==k]
                                rs = list(conversion_df_i_k['Amount'])
                                k_primes = list(conversion_df_i_k['upstream'])

                            for k_prime,r in zip(k_primes,rs):
                                for xi in self.m._Omega:
                                    sum_x_prime = 0
                                    for j in self.m._upstream[i]:
                                        if (j,i,k_prime) in self.m._valid_arcs:
                                            sum_x_prime += self.m._x[j,i,k_prime]
                                    sum_x = 0
                                    for j in self.m._downstream[i]: 
                                        if (i,j,k) in self.m._valid_arcs:
                                            sum_x += self.m._x[i,j,k]
                                    n_RHS.append((-(1/r)*sum_x_prime)+sum_x)
                                    cstrs.append(self.m.getConstrByName(f'Flow_balance_{i}_{k_prime}_{k}_{xi}'))
        
        #Constraints (1d) for total mixed-flow capacity
        for i in self.m._V:
            if i not in self.m._v_subsets['Retail']:
                for xi in self.m._Omega:
                    sum_x = 0
                    for j in self.m._downstream[i]: 
                        for k in self.m._K: 
                            if (i,j,k) in self.m._valid_arcs:
                                sum_x += self.m._x[i,j,k]
                    n_RHS.append(self.m._q_mix[i]-sum_x)
                    cstrs.append(self.m.getConstrByName(f'Total_flow_capacity_{i}_{xi}'))
                    
        #Updating unmet demand constraints ()
        for (i,k,xi) in self.m._d.keys():
            sum_x = 0
            for j in self.m._upstream[i]:
                if (i,j,k) in self.m._valid_arcs:
                    sum_x += self.m._x[j,i,k]
            d_realization = self.m._d[i,k,xi]
            n_RHS.append(d_realization-sum_x)
            cstrs.append(self.m.getConstrByName(f'Demand_penalty_{i}_{k}_{xi}'))

        #Main loop for constraints
        for xi in self.m._Omega:
            for (i,k) in self.m._active_i_k:
                #Unmet demand risk-averse constraints
                n_RHS.append(-self.m._eta_u[i,k])
                cstrs.append(self.m.getConstrByName(f'CVaR(u)_{i}_{k}_{xi}'))

            for (i,j,k) in self.m._valid_arcs:
                
                #Lateness risk-averse constraints
                n_RHS.append(-self.m._eta_v[i,j,k])
                cstrs.append(self.m.getConstrByName(f'CVaR(v)_{i}_{j}_{k}_{xi}'))

                #Flow capacity constraints (1c)
                n_RHS.append((self.m._q_ind[i,j,k]*self.m._y[i,j,k])-self.m._x[i,j,k])
                cstrs.append(self.m.getConstrByName(f'Capacity_{i}_{j}_{k}_{xi}'))

                #Linearized constraints replacing (3b)_1
                n_RHS.append(self.m._bigM[xi]*self.m._y[i,j,k])
                cstrs.append(self.m.getConstrByName(f'o_y_linear1_{i}_{j}_{k}_{xi}'))

                #Linearized constraints replacing (3b)_3
                n_RHS.append(-self.m._bigM[xi]+self.m._bigM[xi]*self.m._y[i,j,k])
                cstrs.append(self.m.getConstrByName(f'o_y_linear3_{i}_{j}_{k}_{xi}'))

                #Linearized constraints replacing (3b)_4
                n_RHS.append(self.m._l[i,j,k,xi]*self.m._y[i,j,k])
                cstrs.append(self.m.getConstrByName(f'o_y_linear4_{i}_{j}_{k}_{xi}'))
        
        #updating all RHS
        self.m.setAttr('RHS',cstrs,n_RHS)
        self.m.update()


    def solve_subproblem(self,silent,first_stage=None):
        if self.m._subgradient == True:
            #Fixing the value based on first-stage solutions
            for idx in self.m._x:
                #self.m._x[idx].setAttr('UB',first_stage[0][idx])
                #self.m._x[idx].setAttr('LB',first_stage[0][idx])
                self.m.addConstr(self.m._x[idx] <= first_stage[0][idx])
                self.m.addConstr(self.m._x[idx] >= first_stage[0][idx])
                
                #self.m._y[idx].setAttr('UB',first_stage[1][idx]) 
                #self.m._y[idx].setAttr('LB',first_stage[1][idx])
                self.m.addConstr(self.m._y[idx] <= first_stage[1][idx])
                self.m.addConstr(self.m._y[idx] >= first_stage[1][idx])
            
        sc = gp.StatusConstClass
        stat_dict = {sc.__dict__[k]: k for k in sc.__dict__.keys() if k[0] >= 'A' and k[0] <= 'Z'}
        if silent == False:
            self.m.Params.LogToConsole = 1
        else:
            self.m.Params.LogToConsole = 0
        self.m.optimize()
        stat = self.m.Status
        if stat_dict[stat] == 'OPTIMAL': 
            self.m._z_sol = self.m.getAttr('x',self.m._z)
            self.m._u_sol = self.m.getAttr('x',self.m._u)

            if self.leadtime == True:
                self.m._a_sol = self.m.getAttr('x',self.m._a)
                self.m._o_sol = self.m.getAttr('x',self.m._o)
                self.m._v_sol = self.m.getAttr('x',self.m._v)
                

        
        
        
