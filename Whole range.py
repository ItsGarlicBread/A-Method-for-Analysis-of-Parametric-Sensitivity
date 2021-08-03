import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

#just to check run time
import time
start_time = time.time()


#setting variables
#Note: Variables ending in 'u' are for ramping up, 'd' are for ramping down, 'c' is for base equation, 'v' is for variable equation
k = 0.3
alpha = 1
beta = -1
omegauc = 1.25
omegadc = 1.25
omegauv = 1.25
omegadv = 1.25
rsv_up = []
rsv_down = []
rsc_up = []
rsc_down = []
distance_down = []
distance_up = []


#setting steps and iterations for gamma (variables end in g)
num_stepsg = 4000
stepg = 0.0001
intervalg = num_stepsg * stepg


cu, du = 1, 0
cd, dd = 1, 0
au, bu = 1, 0
ad, bd = 1, 0
nsg = np.linspace(0, num_stepsg, num_stepsg)


#setting steps and iterations for alpha (variables end in a)
num_stepsa = 8000
stepa = 0.0125
intervala = num_stepsa * stepa
nsa = np.linspace(0, num_stepsa, num_stepsa)


#ramping up for base equation
def duffinguc(xuc,tuc):
    return [xuc[1], xuc[0] - 0.3*xuc[1] - 1*xuc[0]**3 + gammauc * np.cos(omegauc*tuc)]

for ng in nsg:
    gammauc = stepg * ng
    tuc = np.linspace(0, (4*np.pi) / omegauc, 200)
    xsuc = odeint(duffinguc, [cu, du], tuc)
    for i in range(2):
        cu = xsuc[100, 0]
        du = xsuc[100, 1]
        ruc = np.sqrt(cu**2 + du**2)
        rsc_up.append([ng, ruc])

rsc_up = np.array(rsc_up)
rscupv = rsc_up[:,1]


#defining base functions
def duffingdc(xdc,tdc):
    return [xdc[1], xdc[0] - 0.3*xdc[1] - 1*xdc[0]**3 + gammadc * np.cos(omegadc*tdc)]

#ramping down for base equation
for ng in nsg:
    gammadc = intervalg - stepg * ng
    tdc = np.linspace(0, (4*np.pi) / omegadc, 200)
    xsdc = odeint(duffingdc, [cd, dd], tdc)
    for i in range(2):
        cd = xsdc[100, 0]
        dd = xsdc[100, 1]
        rdc = np.sqrt(cd**2 + dd**2)
        rsc_down.append([num_stepsg - ng, rdc])

rsc_down = np.array(rsc_down)
rscdownv = rsc_down[:,1]

fig, ax = plt.subplots()
xtick_labels = np.linspace(0, intervalg, 5)
ax.set_xticks([xc / intervalg * num_stepsg for xc in xtick_labels])
ax.set_xticklabels(['{:.1f}'.format(xtick) for xtick in xtick_labels])

plt.plot(rsc_up[:, 0], rsc_up[:,1], 'r.', markersize=0.1)
plt.plot(rsc_down[:, 0], rsc_down[:,1], 'b.', markersize=0.05)
plt.xlabel(r'$\gamma$', fontsize=15)
plt.ylabel('r', fontsize=15)
plt.tick_params(labelsize=15)
plt.title('Base (Alpha = 1.0)')

#LOOK HERE FOR CHANGE
alpha = 1.1

#defining variable functions
def duffinguv(xuv, tuv):
    return [xuv[1], -beta * xuv[0] - k*xuv[1] - alpha * xuv[0]**3 + gammauv * np.cos(omegauv*tuv)]

def duffingdv(xdv, tdv):
    return [xdv[1], -beta * xdv[0] - k*xdv[1] - alpha * xdv[0]**3 + gammadv * np.cos(omegadv*tdv)]

for na in nsa:
    alpha = 1 + na * stepa
    print(int(na))
    #ramping up for variable equation
    for ng in nsg:
        gammauv = stepg * ng
        tuv = np.linspace(0, (4*np.pi) / omegauv, 200)
        xsuv = odeint(duffinguv, [au, bu], tuv)
        for i in range(2):
            au = xsuv[100, 0]
            bu = xsuv[100, 1]
            ruv = np.sqrt(au**2 + bu**2)
            rsv_up.append([ng, ruv])

    rsv_up = np.array(rsv_up)

    #LOOK HERE FOR CHANGE 
    rsupv = rsv_up[:,1]
    #diff = rsupv - rscupv 
    Distanceu = np.dot(rsupv,rscupv)/np.sqrt(np.dot(rscupv,rscupv)*np.dot(rsupv,rsupv))
    print("--- %s seconds ---" % (time.time() - start_time))
    distance_up.append([alpha, Distanceu])


    #ramping down for variable equation
    for ng in nsg:
        gammadv = intervalg - stepg * ng
        tdv = np.linspace(0, (4*np.pi) / omegadv, 200)
        xsdv = odeint(duffingdv, [ad, bd], tdv)
        for i in range(2):
            ad = xsdv[100, 0]
            bd = xsdv[100, 1]
            rdv = np.sqrt(ad**2 + bd**2)
            rsv_down.append([num_stepsg - ng, rdv])

    rsv_down = np.array(rsv_down)

    #LOOK HERE FOR CHANGE 
    rsdownv = rsv_down[:,1]
    #diff = rsupd - rscupd 
    Distanced = np.dot(rsdownv,rscdownv)/np.sqrt(np.dot(rscdownv,rscdownv)*np.dot(rsdownv,rsdownv))
    print("--- %s seconds ---" % (time.time() - start_time))
    distance_down.append([alpha, Distanced])
    rsv_up = []
    rsv_down = []
distance_down = np.array(distance_down)
distance_up = np.array(distance_up)


fig, ax = plt.subplots()
xtick_labels = np.linspace(0, 100, 6)
ax.set_xticks([x + 1 for x in xtick_labels])
ax.set_xticklabels(['{:.1f}'.format(xtick) for xtick in xtick_labels])

plt.plot(distance_up[:, 0], distance_up[:,1], 'r.', markersize=0.1)
plt.plot(distance_down[:, 0], distance_down[:,1], 'b.', markersize=0.1)
plt.xlabel(r'$\alpha$', fontsize=15)
plt.ylabel('Similarity', fontsize=15)
plt.tick_params(labelsize=15)
plt.savefig('Distance_With_Varying_Alpha')








plt.show()
