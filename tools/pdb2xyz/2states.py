import sys
import numpy as np
import matplotlib.pyplot as plt 

from matplotlib import rc
#rc('text.latex', preamble=r'\usepackage[table]{xcolor}')
rc('font',**{'family':'arial','size':80})
#rc('text', usetex=True)
rc('figure', figsize=(20.0*2, 10.0*2))

e1=float(sys.argv[1])
s1=float(sys.argv[2])
e2=float(sys.argv[3])
s2=float(sys.argv[4])
fmod=float(sys.argv[5])

def twoStates(e1,s1,e2,s2,x):
    sx12 = np.power(s1/x,12)
    sx6  = np.power(s1/x,6)
    g = e2*np.exp(-((x-2.0*s1)**2)/s2)
    return 4.0*e1*(sx12-sx6)+g

def line(F,x):
    return -F*x

def lineInverse(F,x):
    return -0.5*F*F/x

x=np.linspace(-2,20,1000)


plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)

#plt.gca().axes.get_xaxis().set_ticklabels([])
#plt.gca().axes.get_yaxis().set_ticklabels([])

ax = plt.gca()
ax.set_xlim([2,15])
ax.set_ylim([-e1*1.2-fmod*10,e2*4])

#plt.title("Two states model")
ax.set_ylabel('Free energy')
ax.set_xlabel('Reaction coordinate')

if(fmod==0):
    plt.plot(x,twoStates(e1,s1,e2,s2,x),linewidth=10)
else:
    plt.plot(x,twoStates(e1,s1,e2,s2,x)+line(x,fmod),linewidth=10)
    plt.plot(x,line(x,fmod),'--',linewidth=10)

#plt.plot(x,twoStates(e1,s1,e2,s2,x)+lineInverse(x,fmod),linewidth=10)
#plt.plot(x,lineInverse(x,fmod),'--',linewidth=10)
plt.axhline(0, color='black')
plt.axvline(3, color='black')
plt.savefig('twoStates.svg',format='svg')
plt.show()
