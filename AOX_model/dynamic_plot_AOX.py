import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time
from matplotlib.lines import Line2D


###
# Define constants
###


### all parameters taken from cerevisiae


#Transcription:


# Kos M, Tollervey D. Yeast pre-rRNA processing and modification occur cotranscriptionally. Mol Cell. 2010 Mar 26 37(6):809-20. doi: 10.1016/j.molcel.2010.02.024. abstract and p.811 left column 4th paragraph and p.816 table 1 & p.817 right column 5th paragraph --> 40 nt/s
# we have 481 nucleatides for hemoglobin, so we calculate how much hemoglobin is transcribed per second


ktx = 8.3e-2     #M/s maximum transcription rate

# Bonven B, Gulløv K. Peptide chain elongation rate and ribosomal activity in Saccharomyces cerevisiae as a function of the growth rate. Mol Gen Genet. 1979 Feb 26 170(2):225-30
# numer of amino acids is 481 --> (481/3)/7.5 (7.5 is "medium" value from paper for AA elongation per second) --> translation per second 


ktl =  4.6e-2               #M/s maximum translation constant


#Brandon Ho et al., Comparative analysis of protein abundance studies to quantify the Saccharomyces cerevisiae proteome, bioRxiv preprint first posted online Feb. 2, 2017
# they specified typical half-life at 32 min, so we took 30 minutes=1800 seconds to calculate decay per second

deg_mRNA = 1.7e-4        #/s degredation constant of mRNA

# Christiano R, Nagaraj N, Fröhlich F, Walther TC. Global proteome turnover analyses of the Yeasts S. cerevisiae and S. pombe. Cell Rep. 2014 Dec 11 9(5):1959-65. doi: 10.1016/j.celrep.2014.10.065. Supplemental Information p.S12 table S4
# 
# average and median half life is 43 min according to : Belle A, Tanay A, Bitincka L, Shamir R, O'Shea EK. Quantification of protein half-lives in the budding yeast proteome. Proc Natl Acad Sci U S A. 2006 Aug 29 103(35):13004-9 p.13004 right column 4th paragraph

# as hemoglobin is big and fulfils important function, it could have longer half life, like that of the ones with high half life defined in paper by christiano et al.
# so 5hrs = 18000s are taken ### we could look up the half life of our particular hemoglobin further

deg_Protein = 1.67e-5      #/s degredation constant of Protein

###
#Define ODE
###

def ODEs(initial_conditions , t):
    #variables = list of concentrations, so here, [mRNA , Protein]. t = time
    mRNA = initial_conditions[0] 
    Protein = initial_conditions[1]  #
    hill_coefficient = 1.539 
    K = 200 #nM
    hill = 200**hill_coefficient/(K**hill_coefficient+200**hill_coefficient) #nM we can vary TF and so indirectly methanol here
    
# arbitrary numbers, try so that concentration is just a little above Kd (steep curve will make it big fast)

    coef_repr = 100
    K_repressor = 50
    conc_repr = 0
    repressor = K_repressor**coef_repr/(K_repressor**coef_repr+conc_repr**coef_repr)

    leakiness = 0.0000001

   
    dTF_dt = ktx - deg_mRNA*mRNA

    # RNA

    dmRNA_dt =      leakiness + (1-leakiness)*repressor*ktx*hill - deg_mRNA*mRNA
    
    # Protein

    dProtein_dt =   ktl*mRNA - deg_Protein*Protein  

    return [dmRNA_dt, dProtein_dt] 


#####
#Solving the ODEs
#####
t0 = 0              #Initial time
t1 = 360000          #Final time
total =  10000   #Number of time steps (larger the better)

initial_conditions = [0.0, 0.0]        #set the initial values for [mRNA] and [Protein]
t = sp.linspace(t0,t1,total)                       #set the array of time values to integrate over

solution = odeint(ODEs , initial_conditions , t) #Produces an 2d array of solutions
                                                 #for each variable wrt time

mRNA = solution[:,0]    #Index all values in first column
Protein = solution[:,1] #Index all values in second column

simulation_len = len(mRNA)

#####
#Plot the data
#####

#Set the parameters for the figure   (arbitrary values, varry as you like)
params = {
    'axes.labelsize':10,
    'font.size':15,
    'legend.fontsize':10,
    'xtick.labelsize':8,
    'ytick.labelsize':8,
    'figure.figsize': [8,8],
}

plt.rcParams.update(params)

''' Working code for one
plt.show()
axes = plt.gca()
axes.set_xlim(0,36000/60)
axes.set_ylim(0, 300)
line_mRNA, = axes.plot([], [], 'r-')
line_protein, = axes.plot([], [], 'b-')
### Ploting should be made prettier, some of the plt functions dont work with figures
#fig , axs = plt.subplots(2)
for i in range(simulation_len):

    line_mRNA.set_xdata(t[0:i]/60)
    line_mRNA.set_ydata(mRNA[0:i])
    #line_protein.set_ydata(Protein[0:i]/1000)
    plt.draw()
    plt.pause(1e-17)
    time.sleep(1e-12)
plt.show()
'''

step = 100

def update_plot(molecule, i, time, values, mole_type):
    if mole_type == 1:
        molecule.set_data(time[0:i]/60, values[0:i]/2000)
    else:
        molecule.set_data(time[0:i]/60, values[0:i])

data = [mRNA, Protein]

plt.figure()
plt.title("Variation of concentrations with time")
colors = ['orange', 'blue']
lines = [Line2D([0], [0], color=c, linewidth=1) for c in colors]
labels = ['protein', 'mRNA']
plt.legend(lines, labels)
plt.xlabel("time (mins)")
plt.ylabel("# of molecules")
plt.grid()
plt.xlabel('time(min)'); plt.ylabel('# molecules')
plt.axis([0, t1/60, 0, 1000])

molecules_plots = []
for i in range(2):
    mol_plot, = plt.plot([],[])
    molecules_plots.append(mol_plot)

for i in range (0, simulation_len, step):
    for j, mol_plot in enumerate(molecules_plots):
        update_plot(mol_plot, i, t, data[j],j)
    plt.draw()
    plt.pause(1e-20)
