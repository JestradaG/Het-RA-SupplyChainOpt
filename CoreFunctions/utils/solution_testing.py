from CoreFunctions.RiskNeutral.RN_subproblem import InitSubproblem
import pickle

def test_sample(data,x,y,file_name):
    #solving one by one for stastitical check
    # full_obj_val = 0
    # for scenario in range(num_scenarios):
    #     sub = InitSubproblem()
    #     sub.load_data(data=data_instances[scenario])
    #     sub.set_subproblem(first_stage=[x,y],leadtime=True,fixedPenalty=False,out_of_sample=True)
    #     sub.solve_subproblem(silent=True)
    #     full_obj_val += sub.m.ObjVal
    # print(full_obj_val/num_scenarios)
    
    #solving all in one to recover all solutions
    sub = InitSubproblem()
    sub.load_data(data=data)
    sub.set_subproblem(first_stage=[x,y],leadtime=True,fixedPenalty=False,out_of_sample=False)
    sub.solve_subproblem(silent=True)
    
    #Saving solutions
    sample_sols = {}
    sample_sols['x'] = x
    sample_sols['y'] = y
    sample_sols['a'] = sub.m._a_sol
    sample_sols['o'] = sub.m._o_sol
    sample_sols['z'] = sub.m._z_sol
    sample_sols['u'] = sub.m._u_sol
    sample_sols['v'] = sub.m._v_sol
    sample_sols['Obj'] = sub.m.ObjVal
    sample_sols['Omega'] = sub.m._Omega
    #with open(f'Experiments/{file_name}.pickle', 'wb') as handle:
    #    pickle.dump(sample_sols, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    return sample_sols
