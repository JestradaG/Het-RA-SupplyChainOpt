import gurobipy as gp
from gurobipy import GRB
import pandas as pd


class InitExtended:
    def __init__(self):
        self.m = gp.Model()
        self.m.Params.TimeLimit = 60*10
        self.m.reset()
        

    def load_data(self,data):
        #preloading relevant parameters to model structure
        self.m._d_bars = data['demand_bar']
        self.m._Omega,self.m._K,self.m._V,self.m._valid_arcs,self.m._active_i_k,self.m._h = data['Omega'],data['K'],data['V'],data['valid_arcs'],data['active_i_k'],data['h']
        self.m._rho_d = data['rho_d']
        self.m._r = data['r']
        self.m._q_ind,self.m._d,self.m._q_mix = data['q_ind'],data['d'],data['q_mix']
        self.m._v_subsets = data['v_subsets']
        self.m._l,self.m._t,self.m._rho_l_u,self.m._rho_l_f = data['l'],data['t'],data['unitLatePenalty'],data['fixLatePenalty']
        self.m._active_prods = data['active_prods']
        self.m._bigM = data['bigM']
        self.m._lambda = data['lambda']
        self.m._depth = data['depth']

        #Parameters for OneShot Formulation
        self.m._c_a = data['cc'] #contract of arc cost parameter
        self.m._c_f = data['c'] #Flow cost parameter
        self.m._upstream,self.m._downstream = data['upstream'],data['downstream'] #Upstream and downstream sets for each entity
        self.m._directed_output,self.m._prod_conv = data['directed_output'],data['prod_conv']

            
 

    def set_extended(self,risk_averse=False,c_e=1.5,leadtime=True,linearized1=True,fixedPenalty=True,linearized2=True,alpha=[0.7,0.7,0.7,0.7]):
        self.leadtime = leadtime
        self.fixedPenalty = fixedPenalty
        ##Flow Variables
        self.m._x = self.m.addVars(self.m._valid_arcs, name=f'flow',lb=0,ub=GRB.INFINITY) #flow i,j,k
        self.m._y = self.m.addVars(self.m._valid_arcs,name='contract',vtype=GRB.BINARY) #contract between entities
        #Emergency flow
        self.m._z = self.m.addVars(self.m._valid_arcs,self.m._Omega,name='emergency_flow',lb=0,ub=GRB.INFINITY) #emergency flows
        ##Unmet demand
        self.m._u = self.m.addVars(self.m._active_i_k,self.m._Omega, name='penalty_demand',lb=0,ub=GRB.INFINITY) #demand penalty term

        if leadtime == True:
            self.m._a = self.m.addVars(self.m._valid_arcs,self.m._Omega, name='arrival_time',lb=0,ub=GRB.INFINITY) #arrival time variable
            self.m._o = self.m.addVars(self.m._active_i_k,self.m._Omega, name='start_time',lb=0,ub=GRB.INFINITY) #activie start time variable
            #self.m._w = self.m.addVars(self.m._valid_arcs,self.m._Omega, name='Total_flow_penalty',lb=0,ub=GRB.INFINITY) #term to compute the whole lateness penalty
            
            if fixedPenalty == False:
                self.m._v = self.m.addVars(self.m._valid_arcs,self.m._Omega, name='deliveryLateness',lb=0,ub=GRB.INFINITY) #delivery lateness
            else:
                self.m._v = self.m.addVars(self.m._valid_arcs,self.m._Omega, name='deliveryLatenessSwitch',vtype=GRB.BINARY) #switch to check if flow is late or not

            if linearized1 == True:
                    self.m._m1 = self.m.addVars(self.m._valid_arcs,self.m._Omega, name='o_y_product',lb=0,ub=GRB.INFINITY)
            
            if linearized2 == True:
                self.m._m2 = self.m.addVars(self.m._valid_arcs,self.m._Omega, name='v_a_product',lb=0,ub=GRB.INFINITY)

        if risk_averse == True:
            #SEtting risk attitudes
            #alpha_map = {'Part':alpha[0],'Manuf':alpha[1],'Dist':alpha[2],'Retail':alpha[3]}
            alpha_map = {'Retail':alpha[0],3:alpha[1],2:alpha[2],1:alpha[3]}
            alpha_tier = {}

            v_level = {}
            for (i,k),level in pd.DataFrame(self.m._depth).iterrows():
                #print(i,k,level.values[0])
                v_level[i] = level.values[0]

            for i in self.m._V:
                if i in self.m._v_subsets['Retail']:
                    alpha_tier[i] = alpha_map['Retail']
                else: 
                    alpha_tier[i] = alpha_map[v_level[i]]
            self.m._alpha_u = alpha_tier
            self.m._alpha_v = alpha_tier
            self.m._mu = 0.5

            self.m._eta_u = self.m.addVars(self.m._active_i_k, name='VaR(u)',lb=0,ub=GRB.INFINITY)
            self.m._eta_v = self.m.addVars(self.m._valid_arcs,name='VaR(v)',lb=0,ub=GRB.INFINITY)

            self.m._s_u = self.m.addVars(self.m._active_i_k,self.m._Omega,name='s(u)',lb=0,ub=GRB.INFINITY)
            self.m._s_v = self.m.addVars(self.m._valid_arcs,self.m._Omega,name='s(v)',lb=0,ub=GRB.INFINITY)
        

        ###Objective Function
        if leadtime == True:
            if risk_averse == True:
                self.m.setObjective((gp.quicksum(self.m._c_f[i,j,k]*self.m._x[i,j,k] + self.m._c_a[i,j,k]*self.m._y[i,j,k] for i,j,k in self.m._valid_arcs))         
                                    + (1/len(self.m._Omega))*gp.quicksum(c_e*self.m._c_f[i,j,k]*self.m._z[i,j,k,xi] for (i,j,k) in self.m._valid_arcs for xi in self.m._Omega)
                                    +   (
                                        gp.quicksum((1-self.m._lambda[i,k])*(self.m._eta_u[i,k]+(1/((1-self.m._alpha_u[i])*len(self.m._Omega)))*gp.quicksum(self.m._s_u[i,k,xi] for xi in self.m._Omega)) for i,k in self.m._active_i_k)
                                        + gp.quicksum(self.m._lambda[i,k]*(self.m._eta_v[i,j,k]+(1/((1-self.m._alpha_v[i])*len(self.m._Omega)))*gp.quicksum(self.m._s_v[i,j,k,xi] for xi in self.m._Omega)) for i,j,k in self.m._valid_arcs)
                                        )
                            ,GRB.MINIMIZE)
            else:
                self.m.setObjective(gp.quicksum(self.m._c_f[i,j,k]*self.m._x[i,j,k] + self.m._c_a[i,j,k]*self.m._y[i,j,k] for i,j,k in self.m._valid_arcs) +         
                                    ((1/len(self.m._Omega))*(
                                    gp.quicksum(self.m._rho_d[i, k] * self.m._u[i, k, xi] for i,k in self.m._active_i_k for xi in self.m._Omega) + 
                                    gp.quicksum( c_e*self.m._c_f[i,j,k]*self.m._z[i,j,k,xi] + self.m._rho_l_u[i,j,k]*self.m._v[i,j,k,xi] for (i,j,k) in self.m._valid_arcs for xi in self.m._Omega))) 
                            ,GRB.MINIMIZE)
        else:
            self.m.setObjective(gp.quicksum(self.m._c_f[i,j,k]*self.m._x[i,j,k] + self.m._c_a[i,j,k]*self.m._y[i,j,k] for i,j,k in self.m._valid_arcs)                
                    +((1/len(self.m._Omega))*(gp.quicksum(self.m._rho_d[i, k] * self.m._u[i, k, xi] for i,k in self.m._active_i_k for xi in self.m._Omega)+
                                              gp.quicksum( c_e*self.m._c_f[i,j,k]*self.m._z[i,j,k,xi] for (i,j,k) in self.m._valid_arcs for xi in self.m._Omega))) 
                    ,GRB.MINIMIZE)


        #First Stage Constraints
        #Constraints (1b), customer demand satisfaction
        for (i,k),d_bar in self.m._d_bars.items():
            self.m.addConstr((gp.quicksum(self.m._x[j,i,k] for j in self.m._upstream[i] if (j,i,k) in self.m._valid_arcs) >= d_bar),f'Customer_{i}_Product_{k}_demand_satisfaction')
        
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
                                    self.m.addConstr((1/r)*gp.quicksum(self.m._x[j,i,k_prime] + self.m._z[j,i,k_prime,xi] for j in self.m._upstream[i] if (j,i,k_prime) in self.m._valid_arcs) 
                                                    + gp.quicksum(-self.m._x[i,j,k] - self.m._z[i,j,k,xi] for j in self.m._downstream[i] if (i,j,k) in self.m._valid_arcs)
                                                    >= 0
                                                    ,name=f'Flow_balance_{i}_{k_prime}_{k}_{xi}')

        #Constraints (1c) for flow capacity
        for (i,j,k) in self.m._valid_arcs:
            for xi in self.m._Omega:
                self.m.addConstr(self.m._x[i,j,k] + self.m._z[i,j,k,xi] - self.m._q_ind[i,j,k]*self.m._y[i,j,k] <= 0 ,name=f'Capacity_{i}_{j}_{k}_{xi}')

        #Constraints (1d) for total mixed-flow capacity
        for i in self.m._V:
            if i not in self.m._v_subsets['Retail']:
                for xi in self.m._Omega:
                    self.m.addConstr(gp.quicksum(self.m._x[i,j,k] + self.m._z[i,j,k,xi] for j in self.m._downstream[i] for k in self.m._K if (i,j,k) in self.m._valid_arcs)
                                                <= self.m._q_mix[i],name=f'Total_flow_capacity_{i}_{xi}')
            
        ### Demand penalty constraints (1e)
        for (i,k,xi),d in self.m._d.items():
            self.m.addConstr((self.m._u[i,k,xi] + gp.quicksum(self.m._x[j,i,k]+self.m._z[j,i,k,xi] for j in self.m._upstream[i] if (j,i,k) in self.m._valid_arcs) >= d ),name=f'Demand_penalty_{i}_{k}_{xi}')

        if risk_averse == True:
            #Unmet demand risk-averse constraints
            for (i,k) in self.m._active_i_k:
                for xi in self.m._Omega:
                    self.m.addConstr(self.m._s_u[i,k,xi] - (self.m._rho_d[i,k]*self.m._u[i,k,xi]) + self.m._eta_u[i,k] >= 0,name=f'CVaR(u)_{i}_{k}_{xi}')

        if leadtime == True:
            if risk_averse == True:
                #lateness risk-averse constraints
                for (i,j,k) in self.m._valid_arcs: 
                    for xi in self.m._Omega:
                        self.m.addConstr(self.m._s_v[i,j,k,xi] - (self.m._rho_l_u[i,j,k]*self.m._v[i,j,k,xi]) + self.m._eta_v[i,j,k]  >= 0,name=f'CVaR(v)_{i}_{j}_{k}_{xi}')
            #Constraints (1f)
            if linearized1 == True:
                #Block of linearized contraints
                for (i,j,k) in self.m._valid_arcs:
                    for xi in self.m._Omega:
                        self.m.addConstr((self.m._m1[i,j,k,xi] - self.m._bigM[xi]*self.m._y[i,j,k] <= 0 ),name=f'o_y_linear1_{i}_{j}_{k}_{xi}')
                        self.m.addConstr((self.m._m1[i,j,k,xi] - self.m._o[i,k,xi] <= 0 ),name=f'o_y_linear2_{i}_{j}_{k}_{xi}')
                        self.m.addConstr((self.m._m1[i,j,k,xi] - self.m._o[i,k,xi] - self.m._bigM[xi]*self.m._y[i,j,k]) >=  - self.m._bigM[xi],name=f'o_y_linear3_{i}_{j}_{k}_{xi}')
                        self.m.addConstr((self.m._a[i,j,k,xi] - self.m._l[i,j,k,xi]*self.m._y[i,j,k] - self.m._m1[i,j,k,xi] >= 0 ),name=f'o_y_linear4_{i}_{j}_{k}_{xi}')
            else:
                self.m.addConstrs((self.m._a[i,j,k,xi] >= (self.m._l[i,j,k,xi]+self.m._o[i,k,xi])*self.m._y[i,j,k] for (i,j,k) in self.m._valid_arcs for xi in self.m._Omega),name='Arrival time')

            #Constraints (1h)-(1g)

            for e_type,i_s in self.m._v_subsets.items():
                if e_type == 'Part':
                    for i in i_s:
                        for k in self.m._K:
                            if (i,k) in self.m._active_i_k:
                                for xi in self.m._Omega:
                                    self.m.addConstr((self.m._o[i,k,xi] == 0),name=f'Last_tier_start_time_{i}_{k}_{xi}')

            for i,j,k in self.m._valid_arcs:
                if ((j in self.m._v_subsets['Dist']) |(j in self.m._v_subsets['Retail'])):
                    if (j,k) in self.m._active_i_k:
                        for xi in self.m._Omega:
                            self.m.addConstr((self.m._o[j,k,xi] - self.m._a[i,j,k,xi] >= 0 ),name=f'Initial_processing_time_non_transformation_{i}_{j}_{k}_{xi}')
                else:
                    for k in self.m._active_prods[j]:
                        for k_sub in self.m._prod_conv[self.m._prod_conv['downstream']==k]['upstream']:
                            if (i,j,k_sub) in self.m._valid_arcs:
                                if (j,k) in self.m._active_i_k:
                                    for xi in self.m._Omega:
                                        self.m.addConstr((self.m._o[j,k,xi]-self.m._a[i,j,k_sub,xi] >= 0),name=f'Initial_processing_time_transformation_{i}_{j}_{k_sub}_{xi}')

            if fixedPenalty == False: 
                #Constraints (3e)
                self.m.addConstrs((self.m._a[i,j,k,xi]  -  self.m._v[i,j,k,xi] <= self.m._t[j,k] for i,j,k in self.m._valid_arcs if (j,k) in self.m._active_i_k for xi in self.m._Omega),name='FlowLatenessComputation')
                #Constraints (3f)
                #self.m.addConstrs((self.m._w[i,j,k,xi] >= self.m._rho_l_u[i,j,k] * self.m._z[i,j,k,xi] for i,j,k in self.m._valid_arcs for xi in self.m._Omega),name='FlowLatenessPenalty')
            else:
                #Constraints (5a)
                self.m.addConstrs((self.m._a[i,j,k,xi] <= self.m._t[j,k]+self.m._bigM[xi]*self.m._v[i,j,k,xi] for i,j,k in self.m._valid_arcs if (j,k) in self.m._active_i_k for xi in self.m._Omega),name='FlowLatenessCheck')
                #Constraints (5b)
                self.m.addConstrs((self.m._w[i,j,k,xi] >= self.m._v[i,j,k,xi]*(self.m._rho_l_f[i,j,k]+self.m._rho_l_u[i,j,k]*(self.m._a[i,j,k,xi]-self.m._t[j,k])) for i,j,k in self.m._valid_arcs if (j,k) in self.m._active_i_k for xi in self.m._Omega),name='TotalLAtenessPenaltyComputation')


    def solve_extended(self,silent):
        sc = gp.StatusConstClass
        stat_dict = {sc.__dict__[k]: k for k in sc.__dict__.keys() if k[0] >= 'A' and k[0] <= 'Z'}
        if silent == False:
            self.m.Params.LogToConsole = 1
        else:
            self.m.Params.LogToConsole = 0
        self.m.optimize()
        self.m._x_sol = self.m.getAttr('x',self.m._x)
        self.m._y_sol = self.m.getAttr('x',self.m._y)
        self.m._z_sol = self.m.getAttr('x',self.m._z)
        self.m._u_sol = self.m.getAttr('x',self.m._u)

        if self.leadtime == True:
            self.m._a_sol = self.m.getAttr('x',self.m._a)
            self.m._o_sol = self.m.getAttr('x',self.m._o)
            self.m._v_sol = self.m.getAttr('x',self.m._v)
            if self.fixedPenalty == True:
                self.m._v_sol = self.m.getAttr('x',self.m._v)
            else:
                self.m._z_sol = self.m.getAttr('x',self.m._z)

        
        
