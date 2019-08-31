import math
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

SP0=2614.45
RU0=1512.155
dt=1/252
T=1362/252
days=1362

def MC(r_usd,div_sp,div_ru,sigma_sp,sigma_ru,rho):
    SP=[SP0]*(days+1)
    RU=[RU0]*(days+1)
    cov=[[dt,rho*dt],[rho*dt,dt]]
    dw1,dw2=np.random.multivariate_normal([0,0],cov,days).T
    
    for i in range(1,days+1):
        SP[i]=SP[i-1]+(r_usd-div_sp)*SP[i-1]*dt+sigma_sp*SP[i-1]*dw1[i-1]
        RU[i]=RU[i-1]+(r_usd-div_ru)*RU[i-1]*dt+sigma_ru*RU[i-1]*dw2[i-1]
        
    mean_sp=np.average(SP[1298:1363])
    mean_ru=np.average(RU[1298:1363])
    
    V=0
    
    if (mean_sp>=1.21*SP0) and (mean_ru>=1.21*RU0) :
        V=min(1000+1000*(min(mean_sp/SP0,mean_ru/RU0)-1.21)*3.34+415,2116.4)*math.exp(-r_usd*1365/252)
    elif (mean_sp>=SP0) and (mean_ru>=RU0) and (min(mean_sp/SP0,mean_ru/RU0)<1.21):
        V=(1000+1000*(min(mean_sp/SP0,mean_ru/RU0)-1)*1.5+100)*math.exp(-r_usd*1365/252)
    elif (mean_sp>=0.95*SP0) and (mean_ru>=0.95*RU0) and (min(mean_sp/SP0,mean_ru/RU0)<1):
        V=(1000+1000*(min(mean_sp/SP0,mean_ru/RU0)-0.95)*2)*math.exp(-r_usd*1365/252)
    else:
        V=(1000*min(mean_sp/SP0,mean_ru/RU0)+50)*math.exp(-r_usd*1365/252)
        
    return(V)

def value(a,b,c,d,e,f,n_sim):
    value=[0]*n_sim
    for i in range(0,n_sim):
        value[i]=MC(0.026825,0.0192,0.013,0.17,0.2,0.906769)
    return(np.mean(value))

value0=value(0.026825,0.0192,0.013,0.17,0.2,0.906769,10000)

####Sensitivity Analysis1
fig = plt.figure()
ax = Axes3D(fig)
X = np.arange(0.15, 0.26, 0.01)
Y = np.arange(0.15, 0.26, 0.01)
X, Y = np.meshgrid(X, Y)
Z=np.ones((12,12))*0
for i in range(0,12,1):
    for j in range(0,12,1):
        Z[j-1,i-1]=value(0.026825,0.0192,0.013,0.15+0.01*i,0.15+0.01*j,0.906769,10000)
ax.plot_surface(X, Y, Z,cmap="rainbow",alpha=0.5)
ax.set_xlabel('sigma_SP', color='b')
ax.set_ylabel('sigma_Russell', color='b')
ax.set_zlabel('price', color='b')
ax.set_zlim(900, 1050)
plt.show()

####Sensitivity Analysis2
rho_list=np.arange(0.85, 0.96, 0.01)
v=[0]*11
for i in range(0,11):
    v[i]=value(0.026825,0.0192,0.013,0.17,0.2,rho_list[i],10000)