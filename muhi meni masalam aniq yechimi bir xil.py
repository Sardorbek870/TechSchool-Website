import math
import numpy
import matplotlib.pyplot as plt
import numpy as np


def T(n,x):
    if n==0:
        return 1
    else:
        T = list([1 for i in range(n+1)])
        T[0]=1
        T[1]=x
        for i in range(2,n+1):
            T[i]=2*x*T[i-1]-T[i-2]
        return T[n]

def delta(i,j):
    if i==j:
        return 1
    else:
        return 0

n=int(input('n='))
eps=0.01

def u(y):
    return ((eps-0.5)/(1-math.exp(-1/eps)))*(1- math.exp(-(y+1)/(2*eps)))-((eps*(y+1))/2)+(y+1)*(y+1)/8

y=list([0 for i in range(n+1)])

#Aniq yechimni hisoblash
ua=list([0 for i in range(n+1)])
for l in range(n+1):
    y[l]=-math.cos(math.pi*l/n)
    ua[l]=u(y[l])

c=list([1 for i in range(n+1)])
c[0]=c[n]=2

b=list([0 for i in range(n+1)])

#b vektorning qiymatlarini hisoblash
for i in range(0,n+1):
    s=0
    for l in range(0,n+1):
        s=s+(y[l]+1)*T(i,y[l])/c[l]
    b[i]=2*s/(n*8*c[i])

beta1=list([0 for i in range(n+1)])
beta0=list([0 for i in range(n+1)])

#beta1 va beta0 ning qiymatlarini hisoblash
for i in range(0,n+1):
    beta1[i]=c[i]/(2*i+2)
    beta0[i]=beta1[i]/(2*(i+2))

gamma0=list([0 for i in range(n+1)])
gamma1=list([0 for i in range(n+1)])

#gamma1 ning qiymatlarini hisoblash
for i in range(2,n+1):
    gamma1[i]=-1/(2*(i-1))

#gamma0 ning qiymatlarini hisoblash
for i in range(1,n+1):
    gamma0[i]=(-beta1[i]+gamma1[i])/(2*i)

delta0=list([0 for i in range(n+1)])

#delta0 ning qiymatlarini hisoblash
for i in range(3,n+1):
    delta0[i]=-gamma1[i]/(2*(i-2))

f1 = [[0 for i in range(n+1)] for j in range(n+1)]
f0 = [[0 for i in range(n+1)] for j in range(n+1)]

#f1 va f0 matritsalarni hosil qilish
for j in range(0,n+1):
    for i in range(0,n+1):
        f1[j][i]=delta(j,i+1)*beta1[i]+delta(j,i-1)*gamma1[i]
        f0[j][i]=delta(j,i+2)*beta0[i]+delta(j,i)*gamma0[i]+delta(j,i-2)*delta0[i]

sigma0=list([0 for i in range(n+1)])
sigma1=list([0 for i in range(n+1)])
deltachiz0=list([0 for i in range(n+1)])
deltachiz1=list([0 for i in range(n+1)])

#sigma0 va deltachiziqcha0 ning qiymatlarini hisoblash
for i in range(0,n+1):
    for j in range(0,n+1):
        sigma0[i]+=f0[j][i]
        deltachiz0[i]+=((-1)**j)*f0[j][i]

#sigma1 va deltachiziqcha1 ning qiymatlarini hisoblash
for i in range(0,n+1):
    for j in range(0,n+1):
        sigma1[i]+=f1[j][i]
        deltachiz1[i]+=((-1)**j)*f1[j][i]

#g0 va g1 matritsalarni hosil qilish
g0 = [[0 for i in range(n+1)] for j in range(n+1)]
g1 = [[0 for i in range(n+1)] for j in range(n+1)]
for j in range(0,n+1):
    for i in range(0,n+1):
        g1[j][i]=f1[j][i]+delta(j,0)*0.5*(deltachiz0[i]-sigma0[i])
        g0[j][i]=f0[j][i]+delta(j,1)*0.5*(deltachiz0[i]-sigma0[i])-delta(j,0)*0.5*(deltachiz0[i]+sigma0[i])

#g1 matritsani hosil qilish
for j in range(0,n+1):
    for i in range(0,n+1):
        g1[j][i]=0.5*g1[j][i]
        if i==j:
            g1[j][i]+=eps

#Tenglamalar sistemasini yechib a vektorni topish
a=numpy.linalg.solve(g1,b)
print('a vektor:')
print(a)


#Taqribiy yechimni hisoblash
u=list([0 for i in range(n+1)])
for l in range(0,n+1):
    for j in range(n+1):
        for i in range(n+1):
            u[l]+=g0[j][i]*a[i]*T(j,y[l])

c1,c2=0,0
for i in range(n+1):
    c1+=0.5*(deltachiz0[i]-sigma0[i])*a[i]
    c2+=-0.5*(sigma0[i]+deltachiz0[i])*a[i]


ua1=list([0 for i in range(n+1)])
ua2=list([0 for i in range(n+1)])
for l in range(n+1):
    ua1[l]=(eps-0.5)*math.exp(-(y[l]+1)/(2*eps))/(2*eps*(1-math.exp(-1/eps))) + (y[l]+1)/4 -eps/2
for l in range(n+1):
    ua2[l]=1/4-(eps-0.5)*math.exp(-(y[l]+1)/(2*eps))/(4*eps**2*(1-math.exp(-1/eps)))

u1=list([0 for i in range(n+1)])
u2=list([0 for i in range(n+1)])

for l in range(n+1):
    u1[l]=c1*T(0,y[l])
    for j in range(0,n+1):
        for i in range(n+1):
            u1[l]+=f1[j][i]*a[i]*T(j,y[l])

for l in range(n+1):
    for j in range(n+1):
        u2[l]+=a[j]*T(j,y[l])

print('Funksiya:')
for l in range(0,n+1):
    print('l=',l,' Aniq yechim=',ua[l], ' Taqribiy yechim=', u[l], ' Xatolik=', math.fabs(ua[l]-u[l]))

print('1-tartibli hosila:')
for l in range(0,n+1):
    print('l=',l,' Aniq yechim=',ua1[l], ' Taqribiy yechim=', u1[l], ' Xatolik=', math.fabs(ua1[l]-u1[l]))

print('2-tartibli hosila:')
for l in range(0,n+1):
    print('l=',l,' Aniq yechim=',ua2[l], ' Taqribiy yechim=', u2[l], ' Xatolik=', math.fabs(ua2[l]-u2[l]))


y_l_array = np.zeros(n+1)
for i in range(n+1):
    y_l_array[i] = -1 * np.cos(np.pi * i / n)

plt.plot(y_l_array,u,'o--r',label='Taqribiy yechim')
plt.plot(y_l_array,ua,':', label='Aniq yechim')
plt.title(r'$Yechimlar-grafigi$')
plt.ylabel('U')
plt.xlabel('l')
plt.legend(loc='best', fontsize=10)
plt.show()

plt.plot(u1,'o--r', label=r'$Taqribiy-yechim$')
plt.plot(ua1,':', label=r'$Aniq-yechim$')
plt.title(r'$1-tartibli hosila$')
plt.ylabel('U')
plt.xlabel('l')
plt.legend(loc='best', fontsize=10)
plt.show()

plt.plot(u2,'o--r', label=r'$Taqribiy-yechim$')
plt.plot(ua2, ':',  label=r'$Aniq-yechim$')
plt.title(r'$2-tartibli hosila$')
plt.ylabel('U')
plt.xlabel('l')
plt.legend(loc='best', fontsize=10)
plt.show()

