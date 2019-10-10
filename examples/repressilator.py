#matplotlib inline
import matplotlib.pyplot as plt
import matplotlib as mpl

#%config InlineBackend.figure_f.ormats=['svg']

mpl.rc('axes', prop_cycle=(mpl.cycler('color', ['r', 'k', 'b','g','y','m','c']) ))

mpl.rc('xtick', labelsize=12) 
mpl.rc('ytick', labelsize=12)

import numpy as np
from bioscrape.types import Model
from bioscrape.simulator import DeterministicSimulator, SSASimulator
from bioscrape.simulator import ModelCSimInterface

import bioscrape
m = bioscrape.types.read_model_from_sbml('models/repressilator_sbml.xml')
s = bioscrape.simulator.ModelCSimInterface(m)
s.py_prep_deterministic_simulation()
s.py_set_initial_time(0)

sim = bioscrape.simulator.DeterministicSimulator()
timepoints = np.linspace(0,1000,10000)
result = sim.py_simulate(s,timepoints)

plt.plot(timepoints,result.py_get_result())
plt.legend(m.get_species_list())
plt.title('Repressilator Model')
plt.xlabel('Time')
plt.ylabel('Amount')
plt.show()