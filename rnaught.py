from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

#Falta rho, beta y D

def plotseird(t, S, E, I, R, D, L=None, R0=None):
    fp, axp = plt.subplots(1,1,figsize=(10,4))
    axp.plot(t, S, 'b', linewidth=2, label='Suceptibles')
    axp.plot(t, E, 'y', linewidth=2, label='Expuestos')
    axp.plot(t, I, 'r', linewidth=2, label='Infectados')
    axp.plot(t, R, 'g', linewidth=2, label='Recuperados')
    axp.plot(t, D, 'k', linewidth=2, label='Muertos')
    axp.set_xlabel('Tiempo (dias)')
    legend = axp.legend(borderpad=2.0)
    legend.get_frame().set_alpha(0.5)
    if L is not None:
        plt.title("Encierro en el día {} ".format(L))
    plt.show()

    f, ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(t, E, 'y', linewidth=2, label='Expuestos')
    ax.plot(t, I, 'r', linewidth=2, label='Infectados')
    ax.plot(t, R, 'g', linewidth=2, label='Recuperados')
    ax.plot(t, D, 'k', alpha=0.7, linewidth=2, label='Muertos')
    ax.set_xlabel('Tiempo (dias)')
    legend = ax.legend(borderpad=2.0)
    legend.get_frame().set_alpha(0.5)
    plt.show()
  
    f = plt.figure(figsize=(12,4))
    ax0 = f.add_subplot()
    ax0.plot(t, S, 'b', linewidth=2, label='Susceptible')
    ax.set_xlabel('Tiempo (dias)')
    legend = ax0.legend(borderpad=2.0)
    legend.get_frame().set_alpha(0.5)
    plt.show()



def deriv(y, t, N, beta, gamma, delta, alpha, rho):
    Z = np.random.uniform(0,1)
    S, E, I, R, D = y
    dSdt = (-beta(t) * S * I / N ) 
    #print(dSdt)
    dEdt = (beta(t) * S * I / N - delta * E ) 
    #print(dEdt)
    dIdt = (delta * E - (1 - alpha) * gamma * I - alpha * rho * I )
    dRdt = (1 - alpha) * gamma * I
    dDdt = alpha * rho * I
    return dSdt, dEdt, dIdt, dRdt, dDdt

def R_0(t):
    return 2.0 if t < L else 0.9

def logistic_R_0(t):
    return R_0_start if t < L else (R_0_start-R_0_end) / (1 + np.exp(-k*(-t+x0))) + R_0_end

def beta(t):
    return logistic_R_0(t) * gamma

N = 6220145
#N = 62201
D = 10.0 
gamma = 1.0 / D
delta = 1.0 / 6.0
R_0_start, k, x0, R_0_end = 4.0, 0.1, 30, 0.9
alpha = 0.096733  # % death rate
rho = 1/10.5844  # days from infection until death

L = 30

S0, E0, I0, R0, D0 = N-8, 6, 2, 0, 0  # initial conditions
t = np.linspace(0, 200, 100) # Grid of time points (in days)
y0 = S0, E0, I0, R0, D0 # Initial conditions vector

# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta, alpha, rho))
#print(ret)
"""
for d in range(len(ret)-1):
    z = np.random.normal(0,1)
    if z > 0.9:
        ret[d+1] = ret[d]
        if z > 0.95:
            new = z * 10/2
            ret[d][0] -= new
            ret[d][2] += new
    if z < 0.05:
        new = z * 10/2
        ret[d][0] -= new
        ret[d][1] += new
    #print(ret[d])


    axp.plot(d, newRet[d][0], 'b', linewidth=2, label='Suceptibles')
    axp.plot(d, newRet[d][1], 'y', linewidth=2, label='Expuestos')
    axp.plot(d, newRet[d][2], 'r', linewidth=2, label='Infectados')
    axp.plot(d, newRet[d][3], 'g', linewidth=2, label='Recuperados')
    axp.plot(d, newRet[d][4], 'k', linewidth=2, label='Muertos')
    axp.set_xlabel('Tiempo (dias)')
    legend = axp.legend(borderpad=2.0)
    legend.get_frame().set_alpha(0.5)
    plt.show()
"""
"""
newRet = np.zeros((len(ret),5))
for d in range(len(ret)):
    s = sum(ret[d])
    prob = [ret[d][0]/s,ret[d][1]/s,ret[d][2]/s,ret[d][3]/s,ret[d][4]/s]
    #printProgressBar(ret[d][0]/s,1)
    for i in range (N):
        ch = np.random.choice([0,1,2,3,4], p= prob)
        newRet[d][ch] += 1*10
        #print(ch)
    
    print(newRet[d])
"""

#print(ret)
#S, E, I, R, D = newRet.T

S, E, I, R, D = ret.T
R0_over_time = [logistic_R_0(i) for i in range(len(t))]  # to plot R_0 over time: get function values

plotseird(t, S, E, I, R, D, R0=R0_over_time)


def data_time(t,data):
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    labels = 'Suceptibles','Expuestos', 'Infectados', 'Recuperados', 'Muertos'
    datos = data[t]
    print(datos)
    sizes = datos[1:]
    explode = (0,0, 0, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
    fig1, ax1 = plt.subplots()
    ax1.pie(datos, explode=explode, labels=labels, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.show()

    barrasd = []
    mE = []
    mI = []
    mR = []
    mD = []
    c = []
    for i in range(0,len(ret),10):
        barrasd.append(ret[i][1:])
        mE.append(ret[i][1])
        mI.append(ret[i][2])
        mR.append(ret[i][3])
        mD.append(ret[i][4])
        c.append(i)
    #print(barrasd)
    n_groups = len(barrasd)

    # create plot
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.20
    opacity = 0.8

    rects1 = plt.bar(index, mE, bar_width,
    alpha=opacity,
    color='b',
    label='Exposed')
    rects2 = plt.bar(index + bar_width, mI, bar_width,
    alpha=opacity,
    color='r',
    label='Infected')
    rects3 = plt.bar(index + bar_width*2, mR, bar_width,
    alpha=opacity,
    color='y',
    label='Recovered')
    rects4 = plt.bar(index + bar_width*3, mD, bar_width,
    alpha=opacity,
    color='g',
    label='Dead')

    plt.xlabel('Days')
    plt.ylabel('People')
    plt.title('Scores')
    plt.xticks(index + bar_width, c)
    plt.legend()

    plt.tight_layout()
    plt.show()

data_time(40,ret)




