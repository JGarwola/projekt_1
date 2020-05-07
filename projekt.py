# -*- coding: utf-8 -*-

import numpy as np
import random
import matplotlib.pyplot as plt

# 2.6-4  2.4-13  2.2-20 2-28 

"day","cum.cases","cum.deaths","inc.cases","inc.deaths","no.tests"

real = np.array([
    [0,0,0,0,0,559],
    [1,1,0,1,0,585],
    [2,1,0,0,0,667],
    [3,5,0,4,0,860],
    [4,6,0,1,0,861],
    [5,11,0,5,0,1157],
    [6,17,0,6,0,1385],
    [7,22,0,5,0,1652],
    [8,31,0,9,0,2030],
    [9,51,1,21,1,2238],
    [10,68,2,17,1,2899],
    [11,104,3,36,1,4425],
    [12,125,3,21,0,5507],
    [13,177,4,52,1,6699],
    [14,238,5,61,1,7932],
    [15,287,5,49,0,9556],
    [16,355,5,68,0,11246],
    [17,425,5,70,0,13119],
    [18,536,5,111,0,15168],
    [19,634,7,98,2,17678],
    [20,749,8,115,1,20227],
    [21,901,10,152,2,23025],
    [22,1051,13,150,3,26338],
    [23,1221,16,170,3,29700],
    ])

t = np.array([38386000, 2.6, 2, 1/4.6, 0.8, 1/0.5, 1/1.4, 2/3, 1/3.6, 1/8,
              1/8, 0.066, 0.3, 0.57])
[N, R_00, lambd, sigma, f, a, b, q, w_1, w_2, w_3, e_1, e_2, e_3] = t
gamma = a*b/(a+b)
beta = (gamma)/(N*(1+f*a/b)/(1+lambd*f*a/b))

class population:
        
    def __init__(self, status = np.array([N-1,0,1,0,0,0,0,0,0]), kappa = 1, R0 = 2.6):
        self.status = status
        self.S = self.status[0]
        self.E = self.status[1]
        self.IA = self.status[2]
        self.IS = self.status[3]
        self.Q = self.status[4]
        self.H = self.status[5]
        self.V = self.status[6]
        self.D = self.status[7]
        self.R = self.status[8]
        self.kappa = kappa
        self.R0 = R0
        
    def state(self):
        tab = np.array([[self.S, self.E, self.IA, self.IS, self.Q,
             self.H, self.V, self.D, self.R],])
        return tab

    def update_pop(self):
        zagrozenia = self.kappa*beta*self.R0*(self.IA + lambd*self.IS)*self.S
        zarazenia = sigma*self.E
        objawowi_1 = (1 + f*(a/b) )*gamma*self.IA
        objawowi_2 = f*a*self.IA
        objawowi_3 = objawowi_1 - objawowi_2 #(1-f)*gamma*self.IA
        bezobjaw_1 = b*self.IS
        bezobjaw_2 = q*bezobjaw_1
        bezobjaw_3 = bezobjaw_1 - bezobjaw_2
        kwarantanna_1 = w_1*self.Q
        kwarantanna_2 = e_1*kwarantanna_1
        kwarantanna_3 = kwarantanna_1 - kwarantanna_2
        hospital_1 = w_2*self.H
        hospital_2 = e_2*hospital_1
        hospital_3 = hospital_1 - hospital_2
        intense_1 = w_3*self.V
        intense_2 = e_3*intense_1
        intense_3 = intense_1 - intense_2
        
        self.S += -zagrozenia
        self.E += zagrozenia - zarazenia
        self.IA += zarazenia - objawowi_1
        self.IS += objawowi_2 - bezobjaw_1
        self.Q += bezobjaw_2 - kwarantanna_1
        self.H += kwarantanna_2 - hospital_1
        self.V += hospital_2 - intense_1
        self.D += intense_2
        self.R += ( objawowi_3 + bezobjaw_3 + kwarantanna_3
                   + hospital_3 + intense_3)
        
    def szum(self):
        x = random.uniform(0, np.sqrt(self.IA + self.IS))
        self.S += -x
        self.IS += x
        
    def run_sim(self, K, rr, czy_szum,T=0):
        zapis = np.array([ [N-1,0,1,0,0,0,0,0,0], ])
        self.R0 = rr
        for t in range(100):
            
            if t+T >= 10:
                self.kappa = K
            if t % 7 == 0 and czy_szum == 1:    
                self.szum()
                
            self.update_pop()
            zapis = np.append(zapis, self.state(), axis=0)
        return zapis
    
    def check(self):
        return ( self.S + self.E + self.IA + self.IS + self.Q + self.H +
                self.V + self.D + self.R )
    
real_cases = [ real[i][1] for i in range(23) ]
real_deaths = [ real[i][2] for i in range(23) ]
real_tests = [ real[i][5] for i in range(23) ]

t = range(23)

#PRZESUNIĘCIE ZWIĄZANE Z T_SEED
delay = np.array([4,13,20,28])
#PIERWSZA KOLUMNA
for r in [2.6, 2.4, 2.2, 2]:
    f2 = int(-5*r +13)
    def zarazeni(kappa,r):
        model = p.run_sim(kappa/10, r, 0, -delay[f2])
        return [ ( model[i][2] + model[i][3] ) for i in range(delay[f2], delay[f2] + 23) ]
 
    fig=plt.figure()
    ax=fig.add_subplot(111)
    
    for j in range(11):
        p = population()
        ax.plot(t, zarazeni(j,r), label=str(j/10))
    ax.plot(t, real_cases, label = "dane", marker = 'o', c='b')
    plt.legend(loc=2)
    plt.show()
    
#DRUGA KOLUMNA (nie pomnozylem przez 0.01)
for r in [2.6, 2.4, 2.2, 2]:
    f2 = int(-5*r +13)
    def wykrywalnosc(kappa,r):
        model = p.run_sim(kappa/10, r, 0, -delay[f2])
        return [ real[i][1]/(model[i + delay[f2]][2] + model[i + delay[f2]][3]) for i in range(23) ]
 
    fig=plt.figure()
    ax=fig.add_subplot(111)
    
    for j in range(11):
        p = population()
        ax.plot(t, wykrywalnosc(j,r), label=str(j/10))
    plt.legend(loc=2)
    plt.show()

#TRZECIA KOLUMNA
for r in [2.6, 2.4, 2.2, 2]:
    f2 = int(-5*r +13)
    def ilosc_testow(kappa,r):
        model = p.run_sim(kappa/10, r, 0, -delay[f2])
        return [ real[i][5]/(model[i + delay[f2]][2] + model[i + delay[f2]][3]) for i in range(23) ]
 
    fig=plt.figure()
    ax=fig.add_subplot(111)
    
    for j in range(11):
        p = population()
        ax.plot(t, ilosc_testow(j,r), label=str(j/10))
    plt.legend(loc=2)
    plt.show()
    
#DODATKOWA KOLUMNA

def podatni(kappa,r):
    model = p.run_sim(kappa/10, r, 0)
    return [ ( model[i][0] ) for i in range(100) ]
 
fig=plt.figure()
ax=fig.add_subplot(111)
    
for j in range(11):
    p = population()
    ax.plot(range(100), podatni(j,2.6), label=str(j/10))
plt.legend(loc=2)
plt.show()

def hospitalizowani(kappa,r):
    model = p.run_sim(kappa/10, r, 0)
    return [ ( model[i][5]) for i in range(100) ]
#CZY DAC TU JESZCZ LUDZI Z Q?
fig=plt.figure()
ax=fig.add_subplot(111)
    
for j in range(11):
    p = population()
    ax.plot(range(100), hospitalizowani(j,2.6), label=str(j/10))
plt.legend(loc=2)
plt.show()

def intensywna(kappa,r):
    model = p.run_sim(kappa/10, r, 0)
    return [ (model[i][6] ) for i in range(100) ]
#CZY DAC TU JESZCZ LUDZI Z Q?
fig=plt.figure()
ax=fig.add_subplot(111)
    
for j in range(11):
    p = population()
    ax.plot(range(100), intensywna(j,2.6), label=str(j/10))
plt.legend(loc=2)
plt.show()

def zmarli(kappa,r):
    model = p.run_sim(kappa/10, r, 0)
    return [ (model[i][8] ) for i in range(100) ]
#CZY DAC TU JESZCZ LUDZI Z Q?
fig=plt.figure()
ax=fig.add_subplot(111)
    
for j in range(11):
    p = population()
    ax.plot(range(100), zmarli(j,2.6), label=str(j/10))
plt.legend(loc=2)
plt.show()