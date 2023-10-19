import pandas as pd
import numpy as np

def solve_Heuristic(data,criterion):
    #setting up variables and parameters for use
    #Checking the demand to satisfy for the customers
    x = {}
    y = {}
    for i,j,k in data['valid_arcs']:
        x[i,j,k] = 0
        y[i,j,k] = 0

    #building capacity dataframe

    caps = {}
    for (i),cap in data['q_mix'].items():
        caps[i] = cap
    capacities = pd.DataFrame(caps,index=[0]).T
    capacities.columns = ['cap']


    leadTimes = {}
    for (i,j,k,_),l in data['l'].items():
        if (i,j,k) not in leadTimes.keys():
            leadTimes[i,j,k] = 0
        leadTimes[i,j,k] += l
    for (i,j,k) in data['valid_arcs']:
        leadTimes[i,j,k] = leadTimes[i,j,k]/len(data['Omega'])


    #Building flow cost dataframe
    flow_cost = pd.DataFrame(data['c'],index=[0]).T.reset_index()
    flow_cost.columns = ['i','j','k','cost']

    info_df = pd.DataFrame([data['c'],data['q_ind'],leadTimes],index=['Cost','Capacity','LeadTime']).T
    info_df.index = pd.MultiIndex.from_tuples(info_df.index,names=['i','j','k'])
    info_df.reset_index(inplace=True)

    demand2satisfy = {} #general dictionary for target demand
    demands = {} #general dictionary for generated demands

    #initialization of demands with retail
    ret_dem = {}
    for (i,k,_),d in data['d'].items():
        if (i,k) not in ret_dem.keys():
            ret_dem[i,k] = 0
        ret_dem[i,k] += np.ceil(d/len(data['Omega']))
    demand2satisfy['Retail'] = pd.DataFrame(ret_dem,index=[0]).T
    demand2satisfy['Retail'][0] = demand2satisfy['Retail'][0]

    #Getting information on the tiers so that we can greedily fill approaximate demand
    levels = list(data['depth'].unique())
    levels.sort(reverse=True)
    levels = levels[1:]
    levels.insert(0,'Retail')
    #print(levels)

    for tier_idx,tier in enumerate(levels):
        if tier != 'Retail':
            demand2satisfy[tier] = pd.DataFrame(demands[tier],index=[0]).T
        #Filling the demand for the retailers (with the last tier)
        while demand2satisfy[tier][0].sum() > 0:
            (j,k) = demand2satisfy[tier][0].idxmax() #Starting to fill demand of entity with the most
            while demand2satisfy[tier][0][j,k] > 0: #Working to fill current entities demand
                possible_flows = info_df[(info_df['j'] == j)&(info_df['k']==k)] #Filtering those entities which can serve the entity
                
                #Selecting which entity to choose in our greedy way (depending on the criterion)
                if criterion == 'Capacity':
                    idx = possible_flows[criterion].idxmax()
                else:
                    idx = possible_flows[criterion].idxmin()
                best_flow = possible_flows.loc[idx]
                i,j,k = best_flow['i'],best_flow['j'],best_flow['k'] #Getting scpecific information of the best entity
                total_flow = min(capacities.cap[i],demand2satisfy[tier][0][j,k])
                #assigning value to decision variables
                x[i,j,k] += total_flow 
                y[i,j,k] = 1
                demand2satisfy[tier][0][j,k] -= total_flow 
                #updating capacities
                capacities.cap[i] -= total_flow
                info_df.loc[info_df[(info_df['i']==i) & (info_df['k']==k)].index,'Capacity'] -= total_flow

                ##Spawning demands from flow that just tooks place (demand from the top demand) taking product structure in consideration
                if tier_idx < len(levels)-1:
                    if levels[tier_idx+1] not in demands.keys():
                        demands[levels[tier_idx+1]] = {}
                    prod_conv = data['prod_conv'][data['prod_conv']['downstream']==k]
                    for idx,row_data in prod_conv.iterrows():
                        if (i,row_data['upstream']) not in demands[levels[tier_idx+1]].keys():
                            demands[levels[tier_idx+1]][i,row_data['upstream']] = 0
                        demands[levels[tier_idx+1]][i,row_data['upstream']] += row_data['Amount']*total_flow

                    #Removing emptied capacities
                    #display(info_df)
                    for (a) in capacities[capacities['cap']==0].index:
                        capacities.drop(index=(a),inplace=True)
                        rows_idx = info_df[(info_df['i']==a)].index
                        info_df.drop(index=rows_idx,inplace=True)
    return x,y
    
    