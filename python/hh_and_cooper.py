import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

"""
    Values to alter / update
"""
v_impact=1.0    # impact speed between drop and ice particle
E=1.0           # collision efficiency, just use 1 for now - size independent
C=0.16          # C parameter (see hh and cooper)
T=-5.0          # temperature in deg C - determines HM coefficient
plot_flag=False # whether to generate figure of PSDs
plot_fits=False # whether to generate figure of PSDs
nPSD=1          # set to 1, 2, 3, or 4 - to do analysis of particular PSD
allRime=False   # if False it does the HH and Cooper analysis that depends on drop size
                # if True it calculates SIP based on all Rime accreted



"""
the data in figure 2
"""
dat1=np.array([11.531531531531531, 0.0007692307692307443,13.513513513513509, 0.0007692307692307443,14.41441441441441, 0.004615384615384577,14.95495495495495, 0.013846153846153841,15.135135135135137, 0.019999999999999962,15.31531531531531, 0.028461538461538455,15.85585585585585, 0.056153846153846165,16.036036036036037, 0.07076923076923075,16.21621621621621, 0.08384615384615382,16.396396396396398, 0.10153846153846152,16.57657657657657, 0.11692307692307691,16.75675675675675, 0.1276923076923077,16.936936936936938, 0.15230769230769228,17.11711711711711, 0.17307692307692307,17.2972972972973, 0.1892307692307692,17.65765765765766, 0.21307692307692305,17.83783783783784, 0.2238461538461538,18.01801801801801, 0.23307692307692307,18.378378378378372, 0.23769230769230767,18.73873873873874, 0.23076923076923073,18.918918918918912, 0.22153846153846152,19.279279279279272, 0.2015384615384615,19.45945945945946, 0.19076923076923075,19.819819819819813, 0.15923076923076923,20, 0.1392307692307692,20.36036036036036, 0.11384615384615382,20.72072072072072, 0.08615384615384614,21.26126126126126, 0.05307692307692308,21.62162162162162, 0.03615384615384615,22.882882882882875, 0.007692307692307665,23.783783783783782, 0.003076923076923088,24.504504504504503, 0.001538461538461544,25.405405405405403, 0.0007692307692307443])
dat2=np.array([9.549549549549546, 0.0007692307692307443,11.711711711711708, 0.0038461538461538325,12.972972972972968, 0.010769230769230753,14.594594594594593, 0.023846153846153822,15.31531531531531, 0.03307692307692309,16.21621621621621, 0.04384615384615384,16.936936936936938, 0.05307692307692308,17.47747747747747, 0.06,18.378378378378372, 0.06923076923076923,19.0990990990991, 0.07384615384615384,19.819819819819813, 0.07846153846153844,20.54054054054054, 0.07999999999999999,21.8018018018018, 0.07923076923076922,22.882882882882875, 0.07461538461538461,23.783783783783782, 0.06999999999999998,25.045045045045036, 0.061538461538461514,25.405405405405403, 0.057692307692307654,26.306306306306304, 0.05076923076923076,27.207207207207205, 0.04384615384615384,27.927927927927925, 0.038461538461538436,28.828828828828826, 0.03230769230769229,29.909909909909906, 0.026153846153846166,30.810810810810807, 0.020769230769230762,31.891891891891888, 0.01692307692307693,33.51351351351351, 0.011538461538461497,34.77477477477477, 0.008461538461538465,36.39639639639639, 0.005384615384615377,38.55855855855856, 0.003076923076923088,40.36036036036036, 0.001538461538461544,41.98198198198198, 0.001538461538461544,44.68468468468468, 0.0007692307692307443])
dat3=np.array([18.1981981981982, 0.0007692307692307443,19.45945945945946, 0.0007692307692307443,20.36036036036036, 0.0038461538461538325,20.72072072072072, 0.006923076923076921,21.081081081081074, 0.011538461538461497,21.44144144144144, 0.01615384615384613,21.62162162162162, 0.022307692307692306,21.8018018018018, 0.030769230769230743,22.34234234234234, 0.04153846153846155,22.52252252252252, 0.049230769230769245,22.7027027027027, 0.05461538461538459,22.882882882882875, 0.06076923076923074,23.063063063063062, 0.06846153846153846,23.423423423423422, 0.08153846153846153,23.603603603603602, 0.08923076923076922,23.783783783783782, 0.09692307692307692,24.144144144144143, 0.1069230769230769,24.504504504504503, 0.11538461538461536,24.684684684684683, 0.1192307692307692,25.045045045045036, 0.1246153846153846,25.765765765765764, 0.12384615384615383,26.306306306306304, 0.11692307692307691,26.486486486486484, 0.1123076923076923,26.846846846846844, 0.10461538461538461,27.207207207207205, 0.09692307692307692,27.387387387387385, 0.09153846153846151,27.747747747747745, 0.08307692307692305,28.108108108108105, 0.07461538461538461,28.648648648648646, 0.06538461538461535,29.009009009009006, 0.061538461538461514,29.369369369369366, 0.057692307692307654,29.729729729729726, 0.05461538461538459,30.630630630630627, 0.052307692307692305,31.891891891891888, 0.049230769230769245,32.252252252252255, 0.04769230769230767,33.153153153153156, 0.04461538461538461,33.87387387387387, 0.03999999999999998,34.77477477477477, 0.031538461538461515,35.49549549549549, 0.026153846153846166,36.21621621621621, 0.020769230769230762,36.93693693693693, 0.01615384615384613,37.83783783783783, 0.011538461538461497,38.73873873873873, 0.007692307692307665,39.63963963963963, 0.004615384615384577,40.72072072072071, 0.003076923076923088,42.16216216216216, 0.001538461538461544,43.78378378378378, 0.001538461538461544])
dat4=np.array([11.171171171171167, 0.0007692307692307443,14.054054054054053, 0.0007692307692307443,15.31531531531531, 0.0023076923076922884,17.2972972972973, 0.006923076923076921,18.73873873873874, 0.011538461538461497,20.180180180180173, 0.018461538461538474,21.62162162162162, 0.026923076923076883,23.603603603603602, 0.038461538461538436,25.765765765765764, 0.04769230769230767,27.567567567567565, 0.05307692307692308,29.909909909909906, 0.05384615384615382,31.351351351351354, 0.052307692307692305,32.79279279279279, 0.049230769230769245,34.41441441441441, 0.04461538461538461,35.85585585585586, 0.03923076923076921,37.29729729729729, 0.034615384615384576,38.73873873873873, 0.03,40.36036036036036, 0.024615384615384595,41.801801801801794, 0.019999999999999962,43.78378378378378, 0.01615384615384613,45.4054054054054, 0.012307692307692297,46.84684684684684, 0.010000000000000009,48.28828828828829, 0.008461538461538465,49.909909909909906, 0.006153846153846176,51.35135135135135, 0.005384615384615377,52.79279279279279, 0.0038461538461538325,55.85585585585584, 0.0023076923076922884,58.19819819819819, 0.001538461538461544,59.81981981981982, 0.0007692307692307443])

"""
 extract the data into arrays
"""
diams1  =dat1[0:-1:2]
dNdlogD1=dat1[1:2*len(diams1):2]
diams2  =dat2[0:-1:2]
dNdlogD2=dat2[1:2*len(diams2):2]
diams3  =dat3[0:-1:2]
dNdlogD3=dat3[1:2*len(diams3):2]
diams4  =dat4[0:-1:2]
dNdlogD4=dat4[1:2*len(diams4):2]



"""
    table 1, current paper
"""
N_DSD=[100.,100.,[100.,100.],100.] # generic values for N
D_DSD=[18.4e-6,22.1e-6,[25.3e-6,32.1e-6],31.0e-6] # values of D from table 1
D_DSD=[18.4e-6,20.8e-6,[25.3e-6,32.1e-6],30.0e-6] # values of D from my analysis of Fig 2
S_DSD=[1.6,5.5,[1.9,3.1],8.4] # values of sigma from table 1
S_DSD=[1.08,1.28,[1.1,1.1],1.3] # values of sigma from my analysis of Fig 2

"""
    plot the digitized data
"""
if(plot_flag):
    plt.ion()
    plt.plot(diams1,dNdlogD1,'r-')
    plt.plot(diams2,dNdlogD2,'g-')
    plt.plot(diams3,dNdlogD3,'b-')
    plt.plot(diams4,dNdlogD4,'y-')

d=np.logspace(-7,-4,1000)
n=0
f=1./(np.sqrt(2.0*np.pi)*np.log(S_DSD[n]))* \
    np.exp(-np.log(d/D_DSD[n])**2/(2.0*np.log(S_DSD[n])**2))/np.log(1e6)*0.7
if(plot_fits):
    plt.plot(d*1e6,f,'r--')
n=1
f=1./(np.sqrt(2.0*np.pi)*np.log(S_DSD[n]))* \
    np.exp(-np.log(d/D_DSD[n])**2/(2.0*np.log(S_DSD[n])**2))/np.log(1e6)*0.7
if(plot_fits):
    plt.plot(d*1e6,f,'g--')
n=2
f=1.1/(np.sqrt(2.0*np.pi)*np.log(S_DSD[n][0]))* \
    np.exp(-np.log(d/D_DSD[n][0])**2/(2.0*np.log(S_DSD[n][0])**2))/np.log(1e6)*0.7/2.
f=f+0.4/(np.sqrt(2.0*np.pi)*np.log(S_DSD[n][1]))* \
    np.exp(-np.log(d/D_DSD[n][1])**2/(2.0*np.log(S_DSD[n][1])**2))/np.log(1e6)*0.7/2.
if(plot_fits):
    plt.plot(d*1e6,f,'b--')
n=3
f=1./(np.sqrt(2.0*np.pi)*np.log(S_DSD[n]))* \
    np.exp(-np.log(d/D_DSD[n])**2/(2.0*np.log(S_DSD[n])**2))/np.log(1e6)*0.7*0.8
if(plot_fits):
    plt.plot(d*1e6,f,'y--')
    plt.legend(['DSD1','DSD2','DSD3','DSD4','DSD1-fit','DSD2-fit','DSD3-fit','DSD4-fit'])
    plt.xlim((0,60))
N_DSD=[100.,100.,[100.*1.1,100.*0.4],100.]

"""
analysis for particular PSD
"""
n=nPSD-1


"""
 HM temperature function
"""
def f_hm(T):
    if(T <= -5):
        f=np.minimum(np.maximum(0.,(T+8)/3.),1.0)
    else:
        f=np.minimum(np.maximum(0.,(-T-2)/3.),1.0)
    return f

"""
    this is just used to integrate the drop distribution (to check integration is working)
    the result should be equal to N_DSD
"""
def DropDist(x):
    if(n == 2):
        f=N_DSD[n][0]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n][0]))* \
            np.exp(-np.log(x/D_DSD[n][0])**2/(2.0*np.log(S_DSD[n][0])**2))
        f=f+N_DSD[n][1]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n][1]))* \
            np.exp(-np.log(x/D_DSD[n][1])**2/(2.0*np.log(S_DSD[n][1])**2))
    else:
        f=N_DSD[n]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n]))* \
            np.exp(-np.log(x/D_DSD[n])**2/(2.0*np.log(S_DSD[n])**2))
    
    return f

"""
    this is the size dependent integral in HH and cooper
"""
def g_func_hhc(x):
    if(n == 2):
        f=N_DSD[n][0]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n][0]))* \
            np.exp(-np.log(x/D_DSD[n][0])**2/(2.0*np.log(S_DSD[n][0])**2))
        f=f+N_DSD[n][1]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n][1]))* \
            np.exp(-np.log(x/D_DSD[n][1])**2/(2.0*np.log(S_DSD[n][1])**2))
    else:
        f=N_DSD[n]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n]))* \
            np.exp(-np.log(x/D_DSD[n])**2/(2.0*np.log(S_DSD[n])**2))
    
    return f*x**2/4.0*E

"""
    this is to do the integral for the mass-riming rate
"""
def riming_rate(x):
    if(n == 2):
        f=N_DSD[n][0]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n][0]))* \
            np.exp(-np.log(x/D_DSD[n][0])**2/(2.0*np.log(S_DSD[n][0])**2))
        f=f+N_DSD[n][1]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n][1]))* \
            np.exp(-np.log(x/D_DSD[n][1])**2/(2.0*np.log(S_DSD[n][1])**2))
    else:
        f=N_DSD[n]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n]))* \
            np.exp(-np.log(x/D_DSD[n])**2/(2.0*np.log(S_DSD[n])**2))
    
    return f*(1e-3+x)**2/4.0*E*np.pi*np.pi/6.0*x**3*1000.*v_impact

"""
    this is to do the integral in HH and cooper, to calculate the SIP due to RS
"""
def P_hhc(x,y,gR):
    if(n == 2):
        f1=N_DSD[n][0]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n][0]))* \
            np.exp(-np.log(x/D_DSD[n][0])**2/(2.0*np.log(S_DSD[n][0])**2))
        f1=f1+N_DSD[n][1]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n][1]))* \
            np.exp(-np.log(x/D_DSD[n][1])**2/(2.0*np.log(S_DSD[n][1])**2))

        f2=N_DSD[n][0]/(y*np.sqrt(2.0*np.pi)*np.log(S_DSD[n][0]))* \
            np.exp(-np.log(y/D_DSD[n][0])**2/(2.0*np.log(S_DSD[n][0])**2))
        f2=f2+N_DSD[n][1]/(y*np.sqrt(2.0*np.pi)*np.log(S_DSD[n][1]))* \
            np.exp(-np.log(y/D_DSD[n][1])**2/(2.0*np.log(S_DSD[n][1])**2))
    else:
        f1=N_DSD[n]/(x*np.sqrt(2.0*np.pi)*np.log(S_DSD[n]))* \
            np.exp(-np.log(x/D_DSD[n])**2/(2.0*np.log(S_DSD[n])**2))
    
        f2=N_DSD[n]/(y*np.sqrt(2.0*np.pi)*np.log(S_DSD[n]))* \
            np.exp(-np.log(y/D_DSD[n])**2/(2.0*np.log(S_DSD[n])**2))
    
    f2=1.0 # f2 needs to be 1, because this is for one ice particle in the experiment
    
    P=C*f_hm(T)*gR*np.pi*(x+1.e-3)**2*v_impact*f2*f1*E # note, there is a mistake here
                # I actually need to divide by 4.
    return P



if __name__=="__main__":
    # check the integration works
    I=integrate.quad(DropDist,0,1.e-2)
    print('The integral of the PSD for n=' + str(n) + ' is ' + str(I[0]))
    
    
    # g13
    G13=integrate.quad(g_func_hhc,0,13.e-6,epsabs=1.49e-20,epsrel=1.49e-20)
    print('The G13 integrated over the PSD for n=' + str(n) + ' is ' + str(G13[0]))
    
    # gall
    Gall=integrate.quad(g_func_hhc,0,1.e-3)
    print('The Gall integrated over the PSD for n=' + str(n) + ' is ' + str(Gall[0]))
    
    if allRime:
        gR=1.0
    else:
        gR=G13[0]/Gall[0]
    print('Fraction of rime accreted of sizes less than 13 microns ' + str(gR))
    
    # riming rate - 0 to 1mm - an effective max
    Rime=integrate.quad(riming_rate,0,1.e-3)
    print('The riming rate integrated over the PSD for n=' + str(n) + ' is ' + str(Rime[0]))
    
    # the integral for the production rate
    if allRime:
        P=integrate.dblquad(P_hhc,1.e-6,1.e-3,1.e-6, 1.e-3,args=(gR,))
    else:
        # integrate everything larger than 24 microns, up to 1mm - an effective maximum
        P=integrate.dblquad(P_hhc,24.e-6,1.e-3,1.e-6, 1.e-3,args=(gR,))
    print('Production rate per mg of rime ' + str(P[0]/Rime[0]/1.e6))
    
    
