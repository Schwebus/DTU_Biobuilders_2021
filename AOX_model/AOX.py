import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

####
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

def ODEs(variables, t):
    #variables = list of concentrations, so here, [mRNA , Protein]. t = time
    mRNA = variables[0] 
    Protein = variables[1]  #
    #hill_coefficient_MXR1_methanol = 1.539 
    hill_coeff_AOX_methanol = 3 # mxr1 and mit1, prm without competition according to Wang, X., Wang, Q., Wang, J., Bai, P., Shi, L., Shen, W., ... & Cai, M. (2016). Mit1 transcription factor mediates methanol signaling and regulates the alcohol oxidase 1 (AOX1) promoter in Pichia pastoris. Journal of Biological Chemistry, 291(12), 6245-6261.
    #Km_MXR! = 200 #nM
    methanol_concentration = 1 #mM
    # functionn to model to activity level of gene transcription depending on TF concentration
    #hill_eq_MXR1_vs_methanol = methanol_concentration**hill_coefficient_MXR1_methanol/(K**hill_coefficient_MXR1_methanol+methanol_concentration**hill_coefficient) #nM we can vary TF and so indirectly methanol here
    
    # Couderc, R., & Baratti, J. (1980). Oxidation of methanol by the yeast, Pichia pastoris. Purification and properties of the alcohol oxidase. Agricultural and biological chemistry, 44(10), 2279-2289.
    # this paper gives 1.4 mM as methanol Km for AOX activity on 0.19 mM O2 and 3.1 mM at 0.93 mM O2
    Km_AOX = 1.4 #mM, alternatively 3.1 at higher O2
    hill_eq_AOX_vs_methanol = methanol_concentration**hill_coeff_AOX_methanol/(Km_AOX**hill_coeff_AOX_methanol+methanol_concentration**hill_coeff_AOX_methanol)


# arbitrary numbers, try so that concentration is just a little above Kd (steep curve will make it big fast)

    coef_repr = 100
    K_repressor = 50
    conc_repr = 0
    repressor = K_repressor**coef_repr/(K_repressor**coef_repr+conc_repr**coef_repr)

    leakiness = 0.0000001

   
    dTF_dt = ktx - deg_mRNA*mRNA

    # RNA

    dmRNA_dt =       leakiness + (1-leakiness)*repressor*ktx*hill_eq_AOX_vs_methanol - deg_mRNA*mRNA
    
    # Protein

    max_Hemoglobin = 1000 # g/gDW, corresponding to 10% of cell protein
    factor = 100
    bound_term = np.heaviside(max_Hemoglobin-Protein,0)

    dProtein_dt =   (bound_term * ktl*mRNA) - deg_Protein*Protein

    return [dmRNA_dt, dProtein_dt] 


#####
#Solving the ODEs
#####
t0 = 0              #Initial time
t1 = 36000           #Final time
total =  1000000     #Number of time steps (larger the better)

initial_conditions = [0.0, 0.0]        #set the initial values for [mRNA] and [Protein]
t = sp.linspace(t0,t1,total)                       #set the array of time values to integrate over

solution = odeint(ODEs , initial_conditions , t) #Produces an 2d array of solutions
                                                 #for each variable wrt time

mRNA = solution[:,0]    #Index all values in first column
Protein = solution[:,1] #Index all values in second column


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
### Ploting should be made prettier, some of the plt functions dont work with figures
fig , axs = plt.subplots(2)
axs[0].plot(t/60 , mRNA, label = "mRNA # of molecules")
axs[1].plot(t/60 , Protein, label = "Protein # of molecules")
axs[0].set_title("mRNA # of molecules/time (min)")
axs[1].set_title("Protein # of molecules/time (min)")
fig.suptitle("Variation of concentrations with time")
plt.show()
plt.title("Variation of concentrations with time")
plt.xlabel("time (mins)")
plt.ylabel("# of molecules")
plt.grid()
plt.legend()
#plt.plot(t/60 , mRNA, label = "mRNA # of molecules")
#plt.plot(t/60 , Protein, label = "Protein # of molecules")
#plt.title("Variation of concentrations with time")
#plt.xlabel("time (mins)")
#plt.ylabel("# of molecules")
#plt.grid()
#plt.legend()
#plt.show()
#plt.save('test.png')
