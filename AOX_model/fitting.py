import numpy, scipy, matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings
import numpy as np


###
# Anggiani, M., Helianti, I., & Abinawanto, A. (2018, October). Optimization of methanol induction for expression of synthetic gene Thermomyces lanuginosus lipase in Pichia pastoris. In AIP Conference Proceedings (Vol. 2023, No. 1, p. 020157). AIP Publishing LLC.
###

conc = np.array([0.5,1,1.5,2,2.5,3])
acti = np.array([7,8.1,8.2,8.3,8.9,9])
plt.plot(conc,acti)
plt.show()
# times 9 as the max

def func(x, b, c):
    return  (numpy.power(x, b) / (numpy.power(c, b) + numpy.power(x, b)))*9 

#initialParameters = numpy.array([1.0, 1.0])

fittedParameters, _ = curve_fit(func, conc, acti)

modelPredictions = func([0,0.0001,0.001,0.01,0.1,0.2,0.5,1,1.5,2,2.5,3], *fittedParameters)

print(modelPredictions)
print(acti)
print(fittedParameters)
#print(modelPredictions - acti)
