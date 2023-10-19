import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
from collections.abc import MutableMapping
from scipy import stats
import copy

def flatten(d, parent_key='', sep='_'):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def sample_log_norm(mean,sd,n_samples):
    samples = stats.lognorm.rvs(s=sd,loc=mean,size=n_samples)
    return samples

def sample_norm(mean,sd,n_samples):
    samples = stats.norm.rvs(scale=sd,loc=mean,size=n_samples)
    return samples

def parse_data(file_name,num_scenarios,decomp=False,detDemand=False,detTime=False,out_of_sample=True,disruptions=None): #I would need to run this everytime to change the master iteration counter
    #reading excel file with input data
    demand = pd.read_excel(file_name,sheet_name='Demand',index_col=[0,1])
    demand['AgentType'] = 'Retail'
    demand['ProductionCost'] = 0
    demand['ProductionLineCost'] = 0
    demand['Level']=4
    demand['ProdLevel']=4
    demand['ProductionCapacity']=0
    market = pd.read_excel(file_name,sheet_name='Market',index_col=[0,1])
    market['Demand'] = 0
    product_struct = pd.read_excel(file_name,sheet_name='ProductStructure',index_col=[0,1])
    #first_stage = pd.read_excel(file_name,sheet_name='FirstStageTest',index_col=[0,1,2]).to_dict()['val']

    cols = ['Level','AgentType','Demand','HoldingCost','DemandPenalty','lambda','InventoryPenalty','DueDate','FixedLateness','UnitLateness','ProductionCost','ProductionLineCost','ProdLevel','ProductionCapacity']
    # concatenating the market and the demand
    all_entity_df = pd.concat([market[cols], demand], axis=0)

    #Initilizing sets
    V = set(all_entity_df.index.get_level_values('AgentName')) #set of entities (nodes)
    K = set(all_entity_df.index.get_level_values('ProductType')) #set of product types

    #Getting 
    v_level = {}
    for (i,k),level in pd.DataFrame(all_entity_df['Level']).iterrows():
        #print(i,k,level.values[0])
        v_level[i] = level.values[0]

    #Direct output i_k of involvement of entities and products
    directed_output = list(all_entity_df.index)

    #Active products
    active_prods = {}
    for i in V:
        active_prods[i] = []
        for k in K:
            if (i,k) in directed_output:
                active_prods[i].append(k)

    #Subset of type of agent V^p,V^d,V^t
    V_type = all_entity_df['AgentType'].to_dict()
    v_subsets = {}
    v_subsets['Part'] = []
    v_subsets['Manuf'] = []
    v_subsets['Dist'] = []
    v_subsets['Retail'] = []
    for k,i in V_type.items():
        v_subsets[i].append(k[0])
    for type,lst in v_subsets.items():
        v_subsets[type] = list(set(v_subsets[type]))

    V_type = {}
    for type,ents in v_subsets.items():
        for ent in ents:
            V_type[ent] = type
    #Generating product structure
    prd_str = product_struct.to_dict()['Amount']
    prod_structure = {}
    prod_fwd = {}
    all_prods = []
    all_prods_fwd = []
    for idx,val in prd_str.items():
        if val > 0:
            if idx[0] not in all_prods_fwd:
                all_prods_fwd.append(idx[0])
                prod_fwd[idx[0]] = {}
            if idx[1] not in all_prods:
                all_prods.append(idx[1])
                prod_structure[idx[1]] = {}
            prod_fwd[idx[0]][idx[1]] = val
            prod_structure[idx[1]][idx[0]] = val

    #Initializing demands with those from clients (final products)
    demand_bar = demand.to_dict()['D_AVG']
    #all_demands = demand.to_dict()['Demand'] #Demand vector

    #Encoding active i_k
    active_i_k = []
    for i,prods in active_prods.items():
        for k in prods:
            if i in v_subsets['Manuf']:
                for kk in prod_structure[k].keys():
                    active_i_k.append((i,kk))
            active_i_k.append((i,k))
    active_i_k = list(set(active_i_k))

    #Encoding arcs information (only valid ones)
    #sent_direction = {'Part':'Manuf','Manuf':'Dist','Dist':'Retail'}
    #sent_direction = {'Part':'Manuf','Manuf':'Retail'}
    arc_subsets = {}
    for type,entities in v_subsets.items():
        arc_subsets[type] = {}
        prod_list = []
        if type == 'Retail':
            continue
        else:
            for i in entities:
                for j in V:
                    for k in active_prods[i]:
                        if (i,k) in active_i_k:
                            if (j,k) in active_i_k:
                                if k not in prod_list:
                                    arc_subsets[type][k] = []
                                    prod_list.append(k)
                                arc_subsets[type][k].append((i,j,k))
    valid_arcs = []
    for idx,lst in flatten(arc_subsets).items():
        for arc in lst:
            valid_arcs.append(arc)
    valid_arcs = list(set(valid_arcs))

    upstream = {}
    downstream = {}
    for i in V:
        upstream[i] = []
        downstream[i] = []
    for i_j_k in valid_arcs:
        i = i_j_k[0]
        j = i_j_k[1]
        k = i_j_k[2]
        if i not in upstream[j]:
            upstream[j].append(i)
        if j not in downstream[i]:
            downstream[i].append(j)


    ### First stage parameter
    ccs = market['ContractCost'].to_dict()
    #set for optimality cut
    #So that I can access the entities related to some products
    reverse_active_prods = {}
    for entity,prods in active_prods.items():
        for prod in prods:
            if prod not in list(reverse_active_prods.keys()):
                reverse_active_prods[prod] = []
            reverse_active_prods[prod].append(entity)
    #dictionary with relation of products and tiers
    tier_prod = {}
    for prod,entities in reverse_active_prods.items():
        for tier,entities_tier in v_subsets.items():
            if tier != 'Retail':
                intersect = set(entities).intersection(set(entities_tier))
                if len(intersect) > 0:
                    tier_prod[prod,tier] = list(intersect)
    
    entering_prods = {}
    for i,j,k in valid_arcs:
        if j not in entering_prods.keys():
            entering_prods[j] = []
        entering_prods[j].append(k)
    for key in entering_prods.keys():
        entering_prods[key] = list(set(entering_prods[key]))

    ##parameters
    qs = market['ArcCapacity'].to_dict()
    qbars = market['TransportCapacity'].to_dict()
    fs = market['ArcCost'].to_dict()
    cs = market['TransportCost'].to_dict()
    es = all_entity_df['ProductionCost'].to_dict()
    phis = all_entity_df['ProductionLineCost'].to_dict()
    p_bars = all_entity_df['ProductionCapacity'].to_dict()
    p_caps = all_entity_df['ProdLevel'].to_dict()
    hs = all_entity_df['HoldingCost'].to_dict()
    fixlatepens = all_entity_df['FixedLateness'].to_dict()
    unitlatepens = all_entity_df['UnitLateness'].to_dict()
    rho_ds = all_entity_df['DemandPenalty'].to_dict()
    rho_Is = all_entity_df['InventoryPenalty'].to_dict()
    ts = all_entity_df['DueDate'].to_dict()

    #parsing lambdas
    lambdas = all_entity_df['lambda'].to_dict()

    #Network flow parameters
    q_mix = {} #Link mixed-product capacity
    q_ind = {} #Link single product capacity

    cc = {} #Contract cost for first stage decision
    f = {} #Cost for using arc in the period
    c = {} #Transportation cost
    e = {} #Production cost
    r = {} #product conversion rates
    for k in K:
        for kk in K:
            r[k,kk] = 0
    for k_kk,val in flatten(prod_structure).items():
        prod_rel = k_kk.split('_')
        k = prod_rel[0]
        kk = prod_rel[1] 
        r[k,kk] = val

    prd_str = pd.DataFrame(r,index=[0]).T.reset_index()
    prd_str = prd_str[prd_str[0]>0]
    #Production parameters
    phi = {} #production line opening fixed cost
    p_cap = {} #production rate
    p_bar = {} #Production run length capacity

    #Demand/Inventory
    h =  {} #Unit holding cost
    I_0 = {} #intial inventory
    I_s = {} #safety inventory level

    #Inventory/demand satisfaction penalties
    rho_I = {} #Unmet demand penalty
    rho_d = {} #Unmet inventory penalty

    t = {} #Due date of the entity
    fixLatePenalty = {} #Fixed lateness penalty
    unitLatePenalty = {} #unit lateness penalty

    ##Parsing parameters
    for (i,j,k) in valid_arcs:
        
        q_mix[i] = qbars[i,k]
        q_ind[(i,j,k)] = qs[i,k]
        cc[(i,j,k)] = ccs[i,k]
        f[(i,j,k)] = fs[i,k]
        c[(i,j,k)] = cs[i,k]
        fixLatePenalty[(i,j,k)] = fixlatepens[i,k]
        unitLatePenalty[(i,j,k)] = unitlatepens[i,k]

    for i_k,val in phis.items():
        phi[i_k[0]] = val

    disp_info = list(hs.keys()) 
    for i,k in active_i_k:
        if (i,k) in disp_info:
            lambdas[i,k] = lambdas[i,k]
            rho_I[i,k] = rho_Is[i,k]
            rho_d[i,k] = rho_ds[i,k]
            e[i,k] = es[i,k]
            h[i,k] = hs[i,k]
            p_cap[i,k] = p_caps[i,k]
            p_bar[i] = p_bars[i,k]
            t[i,k] = ts[i,k]
        else:
            p_cap[i,k] = 0
            p_bar[i] = 0
            for ii,kk in disp_info:
                lambdas[i,k] = lambdas[ii,kk]
                t[i,k] = ts[ii,kk]
                rho_I[i,k] = rho_Is[ii,kk]
                rho_d[i,k] = rho_ds[ii,kk]
                e[i,k] = es[ii,kk]
                h[i,k] = hs[ii,kk]
        I_0[i,k] = 0
        I_s[i,k] = 0 

    #Generating scenarios for uncertain parameters
    
    #Hard-coding the number of scenarios
    scenario_num = num_scenarios # Hardcoding the number of scenarios to consider in the SAA
    omega = range(scenario_num)
    #Parsing distribution data for the lead times
    leadtime_df = market[['LeadTime','LT_SD']].to_dict() #As of now we just consider log-norm
    demand_df = demand[['D_AVG','D_SD']].to_dict()
    d ={}
    l = {}
    for (i,j,k) in valid_arcs:
        if detTime == False:
            if 'All' in disruptions.keys():
                mean = leadtime_df['LeadTime'][i,k]*disruptions['All']
                sd = float(mean*leadtime_df['LT_SD'][i,k])
                l_realizations = sample_log_norm(mean=mean,sd=sd,n_samples=len(omega))
                for xi in omega:
                    l[i,j,k,xi] = max(1,np.floor(l_realizations[xi]))
            elif i in disruptions.keys():
                mean = leadtime_df['LeadTime'][i,k]*disruptions[i]
                sd = float(mean*leadtime_df['LT_SD'][i,k])
                l_realizations = sample_log_norm(mean=mean,sd=sd,n_samples=len(omega))
                for xi in omega:
                    l[i,j,k,xi] = max(1,np.floor(l_realizations[xi]))
            else:
                l_realizations = np.ones(len(omega))*leadtime_df['LeadTime'][i,k]
                for xi in omega:
                    l[i,j,k,xi] = max(1,np.floor(l_realizations[xi]))
        else:
            l_realizations = np.ones(len(omega))*leadtime_df['LeadTime'][i,k]
            for xi in omega:
                    l[i,j,k,xi] = max(1,np.floor(l_realizations[xi]))

        
        #l[i,j,k,xi] = l_realization
    for (i,k) in demand.index:
        if detDemand == False:
            if 'All' in disruptions.keys():
                mean = demand_df['D_AVG'][i,k]*disruptions['All']
                sd = float(mean*demand_df['D_SD'][i,k])
                d_realizations = sample_norm(mean=mean,sd=sd,n_samples=len(omega))
                for xi in omega:
                    d[i,k,xi] = max(0,np.floor(d_realizations[xi]))
            elif i in disruptions.keys():
                mean = demand_df['D_AVG'][i,k]*disruptions[i]
                sd = float(mean*demand_df['D_SD'][i,k])
                d_realizations = sample_norm(mean=mean,sd=sd,n_samples=len(omega))
                for xi in omega:
                    d[i,k,xi] = max(0,np.floor(d_realizations[xi]))
            else:
                mean = demand_df['D_AVG'][i,k]
                sd = float(mean*demand_df['D_SD'][i,k])
                d_realizations = sample_norm(mean=mean,sd=sd,n_samples=len(omega))
                for xi in omega:
                    d[i,k,xi] = max(0,np.floor(d_realizations[xi]))
        #d[i,k,xi] = d_realization

    #Computing index of constraints that will be reused for subproblem building

    #For constraints 3c
    cst_3c = []
    for e_type,i_s in v_subsets.items():
            if e_type == 'Part':
                for i in i_s:
                    for k in K:
                        if (i,k) in active_i_k:
                            for xi in omega:
                                cst_3c.append((i,k,xi))

    prod_conv = product_struct.reset_index()

    cst_3d_1,cst_3d_2 = [],[]
    for i,j,k in valid_arcs:
            if ((j in v_subsets['Dist']) |(j in v_subsets['Retail'])):
                if (j,k) in active_i_k:
                    for xi in omega:
                        cst_3d_1.append((i,j,k,xi))
            else:
                for k in active_prods[j]:
                    for k_sub in prod_conv[prod_conv['downstream']==k]['upstream']:
                        if (i,j,k_sub) in valid_arcs:
                            if (j,k) in active_i_k:
                                for xi in omega:
                                    cst_3d_2.append((i,j,k,k_sub,xi))

    #Big-m computation:
    tier_lists = {}
    tiers = market['Level'].unique()
    for tier in tiers:
        tier_lists[tier] = list(market['Level'][market['Level']==tier].reset_index()['AgentName'].unique())

    bigM = {}
    for xi in omega:
        bigM[xi] = 0
    tiers = market['Level'].unique()
    for xi in omega:
        for tier in tiers:    
            temp = []
            for (i,j,k) in valid_arcs:
                if i in tier_lists[tier]:
                    temp.append(l[i,j,k,xi])
            bigM[xi] += max(temp)


    data_dict = {'cc':cc,'demand_bar':demand_bar,'tier_prod':tier_prod,'reverse_active_prods':reverse_active_prods,'entering_prods':entering_prods,'prd_str':prd_str, #First stage parameters
                'depth':all_entity_df['Level'],'K':K,'V':V,'Omega':omega,'valid_arcs':valid_arcs,
                'active_i_k':active_i_k,'c':c,'h':h,'e':e,'f':f,'phi':phi,'rho_I':rho_I,
                'rho_d':rho_d,'r':r,'p_cap':p_cap,'I_0':I_0,'q_ind':q_ind,'q_mix':q_mix,'p_bar':p_bar,
                'd':d,'I_s':I_s,'V_type':V_type,'v_subsets':v_subsets,'active_prods':active_prods,'bigM':bigM,
                'prod_conv':prod_conv,'l':l,'t':t,'unitLatePenalty':unitLatePenalty,'fixLatePenalty':fixLatePenalty,
                'upstream':upstream,'downstream':downstream,'directed_output':directed_output,
                'lambda':lambdas,
                'cst_3c':cst_3c,'cst_3d_1':cst_3d_1,'cst_3d_2':cst_3d_2,
                'v_level':v_level}
    if decomp == False:
        return data_dict
    
    elif decomp == True:
        data_instances = []
        
        for xi in omega:
            l_temp = {}
            for (i,j,k) in valid_arcs:
                l_temp[i,j,k,0] = l[i,j,k,xi]
            
            d_temp = {}
            for (i,k) in demand.index:
                d_temp[i,k,0] = d[i,k,xi]
                


            data_dict_temp = {'num_scenarios':num_scenarios,
                    'id':xi,'cc':cc,'demand_bar':demand_bar,'tier_prod':tier_prod,'reverse_active_prods':reverse_active_prods,'entering_prods':entering_prods,'prd_str':prd_str, #First stage parameters
                    'depth':market['Level'],'K':K,'V':V,'Omega':range(1),'valid_arcs':valid_arcs,
                    'active_i_k':active_i_k,'c':c,'h':h,'e':e,'f':f,'phi':phi,'rho_I':rho_I,
                    'rho_d':rho_d,'r':r,'p_cap':p_cap,'I_0':I_0,'q_ind':q_ind,'q_mix':q_mix,'p_bar':p_bar,
                    'd':d_temp,'I_s':I_s,'V_type':V_type,'v_subsets':v_subsets,'active_prods':active_prods,'bigM':bigM,
                    'prod_conv':prod_conv,'l':l_temp,'t':t,'unitLatePenalty':unitLatePenalty,'fixLatePenalty':fixLatePenalty,
                    'upstream':upstream,'downstream':downstream,'directed_output':directed_output,
                'lambda':lambdas}
            data_instances.append(data_dict_temp)
        
        return data_instances,data_dict
    

def inject_pipe_id(data_instances,master_iteration,pipe):
    
    data_insts = copy.deepcopy(data_instances)
    for data_instance in data_insts:
        data_instance['id'] = (master_iteration,data_instance["id"])
        data_instance['lambda_k']=pipe[0], 
        data_instance['anticipatory_penalties']=pipe[1],
        data_instance['S_x']=pipe[2],
        data_instance['cuts_x_1']=pipe[3],
        data_instance['cuts_x_0']=pipe[4],
        data_instance['cuts_g_1']=pipe[5],
        data_instance['cuts_g_0']=pipe[6]

    return data_insts