import bioscrape as bs
from bioscrape.types import Model
from bioscrape.simulator import py_simulate_model
import scipy.stats as stats
import numpy as np
import pylab as plt
#Create a Model
species = ["X"]
reactions = [(["X"], [], "massaction", {"k":"d1"}), ([], ["X"], "massaction", {"k":"k1"})]
k1 = 10.0
d1 = .2
params = [("k1", k1), ("d1", d1)]
initial_condition = {"X":0}
M = Model(species = species, reactions = reactions, parameters = params, initial_condition_dict = initial_condition)
M.py_initialize()

N = 10 #Number of flow distributions taken
time_samples = [2.0,5.0,10.0,20.0] #time samples at which data is taken
nS = len(time_samples)
noise_std = 0.05 #Standar deviation of the guassian noise added onto the measurements

#Generate Distributions
R = [] #Results as Pandas Dataframes
Data = [] #Results will become a numpy array
MultipleInitialConditions = True #Different initial conditions for each trajectory?
N_simulations = 100
meas = 1 # Number of measurements 
Data = np.zeros( (nS,N,N_simulations,meas))
X0_list = [] #multiple initial conditions will be saved for inference
for nS_i in range(nS):
    t = time_samples[nS_i]
    for n in range(N):
        if MultipleInitialConditions:
            initial_condition = {"X": np.random.randint(0, 100)}
            X0_list.append(initial_condition)
        M.set_species(initial_condition)
        for i in range(N_simulations):
            r = py_simulate_model(np.array([0,t]), Model = M, stochastic = True)
            R.append(r["X"][1])
            noisy_data = r["X"][1] + np.random.normal(loc = 0, scale = noise_std, size = 1)
            Data[nS_i,n,i] = noisy_data
        
#Create LogLikelihoodFunction
from bioscrape.inference import BulkData
from bioscrape.inference import StochasticStatesLikelihood as SSLL
from bioscrape.inference import FlowData


# Create a Data Objects
DataFlow = FlowData(np.array(time_samples), Data, ["X"], N)

#Ceate Likelihodo objects:
N_simulations = 100 #Number of simulations per sample to compare to
# norm_order = 1 #(integer) Which norm to use: 1-Norm, 2-norm, etc.
moment_order = 2 # (integer) Upto what norm do you want to include for cost computation. 

#If there are multiple initial conditions in a data-set, should correspond to multiple initial conditions for inference.
#Note len(X0_list) must be equal to the number of trajectories N
if MultipleInitialConditions:    
    # Either use norms or moments for cost computation
#     LL_stoch = SSLL(model = M, init_state = X0_list, data = DataFlow, N_simulations = N_simulations, norm_order = norm_order)
    LL_flow = SSLL(model = M, init_state = X0_list, data = DataFlow, N_simulations = N_simulations, moment_order = moment_order)

    

#Multiple samples with a single initial only require a single initial condition.
else:
    # Either use norms or moments for cost computation
#     LL_stoch = SSLL(model = M, init_state = X0_list, data = DataFlow, N_simulations = N_simulations, norm_order = norm_order)
    LL_flow = SSLL(model = M, init_state = initial_condition, data = DataFlow, N_simulations = N_simulations, moment_order = moment_order)
print("LL_flow.N_simulations", LL_flow.py_get_N_simulations())


npoints = 15
d_list = np.logspace(-12, 12, base = 2, num = npoints)
k_list = np.logspace(-12, 12, base = 2, num = npoints)
HM_flow = np.zeros((len(d_list), len(k_list)))

for di in range(len(d_list)):
    print("di=", di, end = "...ki =")
    for ki in range(len(k_list)):
        print(ki, end = " ")
        LL_flow.set_init_params({"d1":d_list[di], "k1":k_list[ki]})
        
        vs = LL_flow.py_log_likelihood()
        
        HM_flow[di, ki] = -vs

plt.figure(figsize = (18, 8))
plt.title("Stochastic States Log Likelihood of a Non-Identifiable Manifold\n k ="+str(k1)+" d="+str(d1))
plt.xlabel("log k")
plt.xticks(range(npoints), [round(np.log(k), 1) for k in k_list])
plt.ylabel("log d")
plt.yticks(range(npoints), [round(np.log(d), 1) for d in d_list])
cb = plt.pcolor(np.log(HM_flow))
_ = plt.colorbar(cb)