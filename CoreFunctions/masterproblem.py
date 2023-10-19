import gurobipy as gp
from gurobipy import GRB
import numpy as np
from typing import Callable


#Master Problem
class InitMasterProblem:
    def __init__(self,id='',cb=False):
        self.m = gp.Model(f'MasterProblem{id}')
        self.m.Params.TimeLimit = 60*30
        self.m._callback = cb
        if self.m._callback == True:
            #gp.setParam('Heuristics', 0)
            self.m._cb: Callable = None
            self.m.Params.lazyConstraints = 1
            self.m._iterations = 0
        self.m._convergence_stats = []
        self.m._MPinst = []
        self.m._SubPinst = []
        self.m._real_runtime = 0

    def register_callback(self, cb: Callable):
        self.m._cb = cb

    def load_data(self,data):
        self.m._data = data
        self.m._d_bars = data['demand_bar']
        self.m._Omega,self.m._K,self.m._V,self.m._valid_arcs,self.m._active_i_k,self.m._h = data['Omega'],data['K'],data['V'],data['valid_arcs'],data['active_i_k'],data['h']
        self.m._rho_d = data['rho_d']
        self.m._r = data['r']
        self.m._q_ind,self.m._d,self.m._q_mix = data['q_ind'],data['d'],data['q_mix']
        self.m._v_subsets = data['v_subsets']
        self.m._l,self.m._t,self.m._unitLatePenalty,self.m._fixLatePenalty = data['l'],data['t'],data['unitLatePenalty'],data['fixLatePenalty']

        #Parameters for OneShot Formulation
        self.m._c_a = data['cc'] #contract of arc cost parameter
        self.m._c_f = data['c'] #Flow cost parameter
        self.m._upstream,self.m._downstream = data['upstream'],data['downstream'] #Upstream and downstream sets for each entity
        self.m._directed_output,self.m._prod_conv = data['directed_output'],data['prod_conv']
        
        #Risk-averse parameters
        self.m._lambda = data['lambda']
        
            
    def set_masterproblem(self,risk_averse=False,mu=0.5):
        self.m._risk_averse = risk_averse
        if self.m._risk_averse == True:
            self.m._mu = mu

        self.m._x = self.m.addVars(self.m._valid_arcs, name=f'flow',lb=0,ub=GRB.INFINITY) #flow i,j,k
        self.m._y = self.m.addVars(self.m._valid_arcs,name='contract',vtype=GRB.BINARY) #contract between entities
        
        if risk_averse == True:
            self.m._eta_u = self.m.addVars(self.m._active_i_k, name='VaR(u)',lb=0,ub=GRB.INFINITY)
            self.m._eta_v = self.m.addVars(self.m._valid_arcs,name='VaR(v)',lb=0,ub=GRB.INFINITY)
        #self.m._theta = self.m.addVar(name='SecondStageApproximation',lb=0,ub=GRB.INFINITY)
        self.m._theta = self.m.addVars(self.m._Omega,name='SecondStageApproximation',lb=0,ub=GRB.INFINITY) #Approximation to second stage optimal objective
        #Objective function (1a)
        if self.m._risk_averse == True:
            self.m.setObjective((gp.quicksum(self.m._c_f[i,j,k]*self.m._x[i,j,k] + self.m._c_a[i,j,k]*self.m._y[i,j,k] for i,j,k in self.m._valid_arcs))
                                + 
                                        gp.quicksum((1-self.m._lambda[i,k])*self.m._eta_u[i,k] for i,k in self.m._active_i_k)
                                        + gp.quicksum(self.m._lambda[i,k]*self.m._eta_v[i,j,k] for i,j,k in self.m._valid_arcs)
                                        
                                + (1/len(self.m._Omega))*gp.quicksum(self.m._theta[xi] for xi in self.m._Omega)  
                                #+ self.m._theta
                                ,GRB.MINIMIZE)
        else:
            self.m.setObjective(gp.quicksum(self.m._c_f[i,j,k]*self.m._x[i,j,k] + self.m._c_a[i,j,k]*self.m._y[i,j,k] for i,j,k in self.m._valid_arcs)
                            + (1/len(self.m._Omega))*gp.quicksum(self.m._theta[xi] for xi in self.m._Omega)  
                            #+ self.m._theta
                            ,GRB.MINIMIZE)
    
        ##First stage constraints

        #Constraints (1b), customer demand satisfaction
        for (i,k),d_bar in self.m._d_bars.items():
            self.m.addConstr((gp.quicksum(self.m._x[j,i,k] for j in self.m._upstream[i] if (j,i,k) in self.m._valid_arcs) >= d_bar),f'Customer_{i}_Product_{k}_demand_satisfaction')

        #Constraints (1c), flow balance considering multi-product structure
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
                                self.m.addConstr((1/r)*gp.quicksum(self.m._x[j,i,k_prime] for j in self.m._upstream[i] if (j,i,k_prime) in self.m._valid_arcs) 
                                                 - gp.quicksum(self.m._x[i,j,k] for j in self.m._downstream[i] if (i,j,k) in self.m._valid_arcs)
                                                >= 0
                                                ,name=f'Flow_balance_{i}_{k_prime}_{k}')

        #Constraints (1d) for individual flow capacity
        for (i,j,k) in self.m._valid_arcs:
            self.m.addConstr(self.m._x[i,j,k] - self.m._q_ind[i,j,k]*self.m._y[i,j,k] <= 0,name=f'Arc_Capacity_{i}_{j}_{k}')
        
        #Constraints (1e) for total mixed-flow capacity
        for i in self.m._V:
            if i not in self.m._v_subsets['Retail']:
                self.m.addConstr(gp.quicksum(self.m._x[i,j,k] for j in self.m._downstream[i] for k in self.m._K if (i,j,k) in self.m._valid_arcs)
                                            <= self.m._q_mix[i],name=f'Total_flow_capacity_{i}')
    
    def check_solution(self,x,y):
        for (i,j,k) in self.m._valid_arcs:
            self.m.addConstr(self.m._x[i,j,k] == x[i,j,k])
            self.m.addConstr(self.m._y[i,j,k] == y[i,j,k])
        self.m.update()
        
    def solve_masterproblem(self,silent):
        sc = gp.StatusConstClass
        stat_dict = {sc.__dict__[k]: k for k in sc.__dict__.keys() if k[0] >= 'A' and k[0] <= 'Z'}
        if silent == False:
            self.m.Params.LogToConsole = 1
        else:
            self.m.Params.LogToConsole = 0
        
        if self.m._callback == True:
            self.m.optimize(self.m._cb)
        else:
            self.m.optimize()
        self.m._x_sol = self.m.getAttr('x',self.m._x)
        self.m._y_sol = self.m.getAttr('x',self.m._y)
        #self.m._theta_sol = self.m._theta.X
        self.m._theta_sol = self.m.getAttr('x',self.m._theta)
        if self.m._risk_averse == True:
            self.m._eta_u_sol = self.m.getAttr('x',self.m._eta_u)
            self.m._eta_v_sol = self.m.getAttr('x',self.m._eta_v)
            
