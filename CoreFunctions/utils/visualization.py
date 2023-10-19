import pandas as pd
##Plotting chatacteritics
import matplotlib.pyplot as plt
import matplotlib
import networkx as nx
import numpy as np
from matplotlib.ticker import PercentFormatter
import seaborn as sns

def generateViz(data,x,y,z,emergency=False,save=False,file_name='',figsize=(9,8),legend_loc='B',flow=True):
    BIGGER_SIZE = 13
    plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE-3)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE-3)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE-3)    # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE-3)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    NodeColorMap = {'Manuf': '#6283A3', 'Part': '#7EA6E0', 'Retail': '#8F7775'}
    NodeTypeLabels = {'Manuf': 'Assembly', 'Part': 'Supplier', 'Retail': 'Customer'}
    linestyle_map = {'Vehicle1':':','Vehicle2':':','PowerTrain1':':','Dashboard':':','PowerTrain2':':','Transmission':':'}
    ArcColorMap = {'INFT$_{1}$':'tab:red','INFT$_{2}$':'tab:red','INFT$_{3A}$':'tab:red','INFT$_{3B}$':'tab:red',
                   'BTN':'brown','SWT':'tab:blue','RAD':'green','NAV':'black','CNTR':'gold','WR':'blue',
                   'WRN$_{1}$':'gray','WRN$_{2}$':'gray','WRN$_{3A}$':'gray','WRN$_{3B}$':'gray',
                   'SCR$_{1}$':'black','SCR$_{2A}$':'black','SCR$_{2B}$':'black',
                   'CHP$_{1}$':'orange','CHP$_{2}$':'orange'}
    #product_type = {}
    ArcTypeLabels = {'Vehicle1':'Vehicle$_{1}$','Dashboard':'Dashboard','Vehicle2':'Vehicle$_{2}$','PowerTrain1':'Powertrain$_{1}$','Transmission':'Transmission','PowerTrain2':'Powertrain$_{2}$'}

    flow_arcs = []
    active_arcs = []
    active_arcs_width = []
    active_arcs_color = []
    active_emergency = []
    active_emergency_color = []
    node_color = []

    #Correcting valid arcs (from the numerical instability of translating my model without too much modification)
    plot_arcs = []
    for (i,j,k) in data['valid_arcs']:
        if (i == j)|(data['v_level'][i] == data['v_level'][j]):
            pass
        else:
            plot_arcs.append((i,j,k))
    #suplementary arcs to try to improve the visualization
    for i in data['V']:
        if 'WRN' in i:
            for j in data['V']:
                if data['v_level'][j] == 1:
                    if j not in data['upstream'][i]:
                        plot_arcs.append((i,j,k))

    for (i,j,k) in plot_arcs:
        if (i,j) not in flow_arcs:
            flow_arcs.append((i,j))
        if (i,j,k) in data['valid_arcs']:
            if flow == True:
                if y[i,j,k] > 0.5:
                    active_arcs.append((i,j))
                    active_arcs_width.append(x[i,j,k]/300)
                    active_arcs_color.append(ArcColorMap[k])
                    if emergency == True:
                        for xi in data['Omega']:
                            if z[i,j,k,xi] > 0:
                                active_emergency.append((i,j))
                                active_emergency_color.append(ArcColorMap[k])
    for i in data['V']:
        node_color.append(NodeColorMap[data['V_type'][i]])

    f = plt.figure(figsize=(9,9))
    ax = f.add_subplot(1,1,1)

    G = nx.DiGraph(rankdir='LR',seed=2)
    G.add_nodes_from(data['V'])
    G.add_edges_from(flow_arcs)
    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
    nodes = nx.draw_networkx_nodes(G, pos,node_size=800,node_color='white',ax=ax,linewidths=2)
    # Set edge color to red
    nodes.set_edgecolor(node_color)
    # Uncomment this if you want your labels
    nx.draw_networkx_labels(G, pos,font_size=BIGGER_SIZE-5,ax=ax)
    # #Adding layout structure
    nx.draw_networkx_edges(
            G,
            pos,
            edgelist=flow_arcs,
            connectionstyle='arc3',
            edge_color='white',ax=ax
            )
    if emergency == True:
        #adding emergency arcs
        nx.draw_networkx_edges(
                G,
                pos,
                arrows=False,
                style=':',
                edgelist=active_emergency,
                connectionstyle='arc3',
                alpha=0.7,
                edge_color=active_emergency_color,ax=ax,
                width=5
                )
    if flow == True:
        nx.draw_networkx_edges(
                G,
                pos,
                edgelist=active_arcs,
                connectionstyle='arc3',
                edge_color=active_arcs_color,
                arrows=True,
                arrowsize=20,
                arrowstyle='-|>',
                width=active_arcs_width,
                style='-',
                min_target_margin=20,ax=ax
                )
    # #plt.savefig(f'TASE_SCLayout.pdf',bbox_inches='tight')
    for label in NodeColorMap:
        ax.scatter(None,None,color=NodeColorMap[label],label=NodeTypeLabels[label],facecolors='none',s=200,lw=3)




        
    #legend1 = plt.legend(bbox_to_anchor=(0.9, 1.35),ncol=2,borderaxespad=0,title='Hola',title_fontsize=10)
    LabArcColorMap = {'Infotainment':'tab:red','Button':'brown','Switch':'tab:blue',
                'Radio':'green','Navigation': 'black','Connector': 'gold',
                'Wire':'blue','Wiring assembly':'gray','Touchscreen':'black',
                'Chip':'orange'}
    if flow == True:
        ax2 = ax.twinx()
        for label in LabArcColorMap:
            ax2.plot([0],[0],color=LabArcColorMap[label],label=label,lw=2)
        ax2.get_yaxis().set_visible(False)
        ax2.legend(bbox_to_anchor=(0.47, 0.99),ncol=2,borderaxespad=0,title='Product type',title_fontsize=10) #top

    if emergency == True:
        ax3 = ax.twinx()
        labelmap = {'-':'Planned',':':'Emergent'}
        for type in ['-',':']:
            ax3.plot([0],[0],color='k',linestyle=type,label=labelmap[type],lw=4)
        ax3.get_yaxis().set_visible(False)
        ax3.legend(bbox_to_anchor=(0.97, 1.085),ncol=2,borderaxespad=0,title='Arc (flow) type',title_fontsize=10,borderpad=0.65)

    handles, labels = ax.get_legend_handles_labels()
    # specify order
    order = [1, 0, 2]

    # pass handle & labels lists along with order as below
    #plt.legend()
    ax.legend([handles[i] for i in order], [labels[i] for i in order],bbox_to_anchor=(0.58, 1.085),ncol=4,borderaxespad=0,title='Entity type',title_fontsize=10,borderpad=0.65)
    # if all_nodes == False:
    #     if legend_loc == 'B':
    #         ax2.legend(bbox_to_anchor=(0.99, 0.15),ncol=3,borderaxespad=0,title='Product type',title_fontsize=10) #bottom
    #     else:
    #         

    if save == True:
        plt.savefig(f'{file_name}.pdf',bbox_inches='tight')
    return f

def ConvergencePlot(convergence_stats,save=False,file_name='',figsize=(8,6)):
    BIGGER_SIZE = 17
    plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    pd.DataFrame(convergence_stats).plot(figsize=figsize,linewidth=4,style=['-','--',':'])
    plt.xlabel('Iterations')
    plt.ylabel('Objective')
    plt.legend(loc='lower right')
    if save == True:
        plt.savefig(f'{file_name}.pdf',bbox_inches='tight')
    return 

def compareDemand_hist(approaches,approach_labels,data,save=True,file_name='',figsize=(8,6),bw=1):
    BIGGER_SIZE = 17
    plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    #approaches = [optimal_insample,lt_second_stage_sols,cst_second_stage_sols,cap_second_stage_sols]
    #approach_labels = ['Optimization','Lead-time','Cost','Capacity']
    #approaches = approaches[0:2]
    tot_dem = sum(data['d'].values())
    scenario_dem = {}
    for xi in data['Omega']:
        scenario_dem[xi] = 0
        for (i,k,scen),d in data['d'].items():
            if scen == xi:
                scenario_dem[xi] += d
    unmet_demand = []
    for idx_app in range(len(approaches)):
        solutions = approaches[idx_app]
        
        solution_approach = approach_labels[idx_app]
        for (i,k,xi),d in solutions['u'].items():
            if i in data['v_subsets']['Retail']:
                tot_unmet = sum(solutions["u"].values())
                #print(d,tot_unmet)
                unmet_demand.append({'Client':i,'Product':k,'Scenario':xi,'Approach':solution_approach,'UnmetDemand':d/scenario_dem[xi]})#/(tot_dem/len(solutions['Omega']))})#np.round(d/tot_dem,8)}) #np.round(d/tot_dem*100,2)
    ax = plt.figure(figsize=figsize).add_subplot(111)
    df = pd.DataFrame(unmet_demand)[['Approach','UnmetDemand']]
    sns.set_context('paper')
    colors = ['#F58E62','#6EBAF5','#A8943B','#F55546']
    sns.violinplot(x='Approach',y='UnmetDemand',data=df,ax=ax, linewidth=2, width=0.95,bw=bw,cut=0,palette=colors)
    plt.xlabel(f'(% of total demand not satisfied) \n Approach')
    plt.ylabel(r'% of unmet demand of infotainment')
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    xticks = list(range(len(approaches)))
    x_ticks_label = []
    for i,approach in enumerate(approaches):
        x_ticks_label.append(f'({np.round(sum(approach["u"].values())/tot_dem*100,2)}%) \n {approach_labels[i]}')
    plt.xticks(xticks,x_ticks_label,rotation=0)
    if save == True:
        plt.savefig(f'{file_name}.pdf',bbox_inches='tight')
    
    return

def compareLateness_hist(approaches=[],approach_labels=[],data={},save=True,file_name='',figsize=(8,6),final_day=10,bbox=(0.9, 1.18)):
    BIGGER_SIZE = 17
    plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    #approaches = [optimal_insample,lt_second_stage_sols,cst_second_stage_sols,cap_second_stage_sols]
    #approach_labels = ['Optimization','Lead-time','Cost','Capacity']
    client_lateness = []
    unmet_dems = []
    for idx_app in range(len(approaches)):
        solutions = approaches[idx_app]
        solution_approach = approach_labels[idx_app]
        for (i,k,xi),d in data['d'].items():
            all_flows = 0
            for j in data['upstream'][i]:
                if (j,i,k) in data['valid_arcs']:
                    all_flows += solutions['x'][j,i,k]+solutions['z'][j,i,k,xi]
            realized_flows = min(d,all_flows)

            for j in data['upstream'][i]:
                if (j,i,k) in data['valid_arcs']:
                    client_lateness.append({'Client':i,'Product':k,'Approach':solution_approach,'Lateness':int(solutions['v'][j,i,k,xi]),'Flow':int(realized_flows)})
        unmet_dems.append(sum(solutions['u'].values()))
    k_data = pd.DataFrame(client_lateness)[['Lateness','Flow','Approach']].groupby(['Lateness','Approach']).sum().unstack().fillna(0)
    k_data_idx_exist = list(k_data.index)
    for i in range(0,final_day):
        if (i in k_data_idx_exist) == False:
            #print(i)
            dummyrow = pd.DataFrame(k_data.iloc[-1])
            dummyrow[k_data.index[-1]] = 0
            dummyrow = dummyrow.T
            dummyrow.index = [i]
            k_data = pd.concat([k_data, dummyrow], ignore_index=False)
            #k_data = k_data.append(dummyrow)
    k_data = k_data.loc[0:final_day]
    k_data.sort_index(inplace=True)
    k_data.columns = [' '.join(col).split()[1] for col in k_data.columns.values]

    #display(k_data)
    #display()
    #k_data = k_data.append(pd.DataFrame([unmet_dems],index=['Unmet'],columns=k_data.columns))
    k_data = pd.concat([k_data,pd.DataFrame([unmet_dems],index=['Unmet'],columns=k_data.columns)],ignore_index=False)
    total_data = sum(data['d'].values())
    #Normalizing data
    for c in k_data.columns: 
        col = k_data[c]
        all_flow = col.sum()
        for r in col.index:
            k_data.loc[r,c] = k_data.loc[r,c]/all_flow

    colors = ['#F58E62','#6EBAF5','#A8943B','#F55546']
    ax = plt.figure(figsize=figsize).add_subplot(111)
                    
    k_data.plot.bar(ax=ax,stacked=False,color=colors,width=0.8)
    
    n_bars = len(list(k_data.columns))
    bars = ax.patches
    hatches = ['x','+','O','-','.','|','*','o']
    patterns = hatches[0:n_bars]# set hatch patterns in the correct order
    hatches = []  # list for hatches in the order of the bars
    for h in patterns:  # loop over patterns to create bar-ordered hatches
        for i in range(int(len(bars) / len(patterns))):
            hatches.append(h)
    for bar, hatch in zip(bars, hatches):  # loop over bars and hatches to set hatches in correct order
        bar.set_hatch(hatch)
    plt.xlabel(f'Days of lateness to customers')
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.ylabel(r'% of delivered units of infotainments')
    plt.legend(
            bbox_to_anchor=bbox,ncol=n_bars,
            borderaxespad=0,title='Approach',title_fontsize=BIGGER_SIZE)
    plt.xticks(rotation=0)
    if save == True:
        plt.savefig(f'{file_name}.pdf',bbox_inches='tight')
    
    return ax

def compare_costs(data,approaches,approach_labels,c_e,save=False,file_name=''):
    BIGGER_SIZE = 17
    plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    setup_arr = []
    plan_arr = []
    emergency_arr = []
    app_labs = []
    costs_arr = {}
    demand_arr = []
    lateness_arr = []


    for idx_app in range(len(approaches)):
        approach = approaches[idx_app]
        solution_approach = approach_labels[idx_app]

        #Computing unmet demand penalty
        demand_penalty = 0
        for (i,k,xi),unmet_demand in approach['u'].items():
            demand_penalty += unmet_demand*data['rho_d'][i,k]
        demand_penalty = demand_penalty/(len(data['Omega'])*10000)
        demand_arr.append(demand_penalty)

        #Computing unmet demand penalty
        lateness_penalty = 0
        for (i,j,k,xi),lateness in approach['v'].items():
            if (approach['x'][i,j,k]+approach['z'][i,j,k,xi]) > 0:
                lateness_penalty += lateness*data['unitLatePenalty'][i,j,k]
        lateness_penalty = lateness_penalty/(len(data['Omega'])*10000)

        lateness_arr.append(lateness_penalty)

        #Computing costs

        #Setup costs (y var)
        total_setup = 0
        for (i,j,k),setup in approach['y'].items():
            total_setup += setup*data['cc'][i,j,k]
        costs_arr[('Setup',solution_approach)] = total_setup

        #Planned flow costs (x var)
        total_plan = 0
        for (i,j,k),plan in approach['x'].items():
            total_plan += plan*data['c'][i,j,k]
        costs_arr[('Plan',solution_approach)] = total_plan

        #Average emergency flow costs (z var)
        total_emergency = 0
        for (i,j,k,xi),emergency in approach['z'].items():
            total_emergency += emergency*c_e*data['c'][i,j,k]
        total_emergency = total_emergency/len(data['Omega'])
        costs_arr[('Emergency',solution_approach)] = total_emergency

        tot_cost = np.round((total_setup+total_plan+total_emergency)/1000,1)

        setup_arr.append(total_setup)
        plan_arr.append(total_plan)
        emergency_arr.append(total_emergency)


        app_labs.append(f"[{int(demand_penalty)} | {int(lateness_penalty)}]\n ({tot_cost}) \n {solution_approach}")

    cost_labels = []
    cost_types = ['Setup','Plan','Emergency']
    for cost_type in cost_types:
        for solution_approach in approach_labels:
            cost_labels.append(costs_arr[cost_type,solution_approach])

    app_costs = {'Setup':setup_arr,'Planned flow':plan_arr,'Average emergent flow':emergency_arr}

    norm_costs = app_costs.copy()
    for lab,costs in norm_costs.items():
        max_cost = max(costs)
        direct_norm_costs = []
        for cost in costs:
            direct_norm_costs.append(cost/max_cost)
        norm_costs[lab] = direct_norm_costs


    width = 0.7
    ax = plt.figure(figsize=(9,6)).add_subplot(111)
    bottom = np.zeros(3)
    df_norm_costs = pd.DataFrame(norm_costs)
    plot = df_norm_costs.plot.bar(ax=ax,stacked=False,width=0.9)

    n_bars = len(list(df_norm_costs.columns))
    bars = ax.patches
    hatches = ['*','x','+','O','-','.','|','*','o']
    patterns = hatches[0:n_bars]# set hatch patterns in the correct order
    hatches = []  # list for hatches in the order of the bars
    for h in patterns:  # loop over patterns to create bar-ordered hatches
        for i in range(int(len(bars) / len(patterns))):
            hatches.append(h)
    count = 0
    colors = ['#F58E62','#6EBAF5','#A8943B','#F55546','#F58E62','#6EBAF5','#A8943B','#F55546','#F58E62','#6EBAF5','#A8943B','#F55546']
    #colors = ['#31C7DE','#8ECDE8','#42A1FF','#31C7DE','#8ECDE8','#42A1FF','#31C7DE','#8ECDE8','#42A1FF']
    for bar, hatch in zip(bars, hatches):  # loop over bars and hatches to set hatches in correct order
        plot.annotate(f'[{format(cost_labels[count]/10000,".1f")}]',
                    (bar.get_x()+bar.get_width()/2,
                    bar.get_height()),ha='center',va='center',
                    size=BIGGER_SIZE-2,xytext=(0,8),textcoords='offset points'
                    )
        
        bar.set_facecolor(colors[count])
        bar.set_hatch(hatch)
        count += 1

    
    


    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.ylabel(r'Normalized cost [Absolute cost $\$ 10^{4}$]')
    plt.xlabel("[Unmet demand penalty $\$ 10^{4}$ | Lateness penalty $\$ 10^{4}$]\n (Total cost $\$ 10^{4}$) \n Approach")
    ax.legend().set_visible(False)
    #             bbox_to_anchor=(1.05,1.18),ncol=3,
    #             borderaxespad=0,title='Cost type',title_fontsize=BIGGER_SIZE)
    xticks = list(range(len(approaches)))
    plt.xticks(xticks,app_labs,rotation=0)
    
    ax3 = ax.twinx()
    labelmap = {'Setup':'*','Planned flow':'x','Average emergent flow':'+'}
    for label,hatch in labelmap.items():
        ax3.bar([0],[0],color='None',label=label)
    ax3.get_yaxis().set_visible(False)
    bars = ax3.patches
    hatches = ['*','x','+']
    for bar,hatch in zip(bars,hatches):
        bar.set_hatch(hatch)
    ax3.legend(bbox_to_anchor=(1.05, 1.18),ncol=3,borderaxespad=0,title='Cost type',title_fontsize=BIGGER_SIZE)

    if save == True:
        plt.savefig(f'{file_name}.pdf',bbox_inches='tight')
    plt.show()
    