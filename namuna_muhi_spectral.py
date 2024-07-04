import math
import numpy
import matplotlib.pyplot as plt

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

eps=0.001
def f(y):
      return ( (eps**5-2*eps**2)*y**4 + (16*eps**4-16*eps)*y**3 + (72*eps**3-2*eps**5-24+4*eps**2)*y**2 + (96*eps**2-16*eps**4+16*eps)*y - 2*eps**2+8 )*math.exp(eps*y)

n=int(input('n='))


y=list([0 for i in range(n+1)])
for l in range(n+1):
    y[l]=math.cos(math.pi*l/n)

c=list([1 for i in range(n+1)])
c[0]=c[n]=2

b=list([0 for i in range(n+1)])

#b vektorning qiymatlarini hisoblash
for i in range(0,n+1):
    s=0
    for l in range(0,n+1):
        s=s+f(y[l])*T(i,y[l])/c[l]
    b[i]=2*s/(n*c[i])

beta3=list([0 for i in range(n+1)])
beta2=list([0 for i in range(n+1)])
beta1=list([0 for i in range(n+1)])
beta0=list([0 for i in range(n+1)])

#betalarning qiymatlarini hisoblash
for i in range(0,n+1):
    beta3[i]=c[i]/(2*i+2)
    beta2[i]=beta3[i]/(2*(i+2))
    beta1[i]=beta2[i]/(2*(i+3))
    beta0[i]=beta1[i]/(2*(i+4))

gamma3=list([0 for i in range(n+1)])
gamma2=list([0 for i in range(n+1)])
gamma1=list([0 for i in range(n+1)])
gamma0=list([0 for i in range(n+1)])

#gamma3 ning qiymatlarini hisoblash
for i in range(2,n+1):
    gamma3[i]=-1/(2*(i-1))

#gamma2 ning qiymatlarini hisoblash
for i in range(1,n+1):
    gamma2[i]=(-beta3[i]+gamma3[i])/(2*i)

#gamma1 ning qiymatlarini hisoblash
for i in range(0,n+1):
    gamma1[i]=(-beta2[i]+gamma2[i])/(2*(i+1))

#gamma0 ning qiymatlarini hisoblash
for i in range(0,n+1):
    gamma0[i]=(-beta1[i]+gamma1[i])/(2*(i+2))

delta2=list([0 for i in range(n+1)])
delta1=list([0 for i in range(n+1)])
delta0=list([0 for i in range(n+1)])

#delta2 ning qiymatlarini hisoblash
for i in range(3,n+1):
    delta2[i]=-gamma3[i]/(2*(i-2))

#delta1 ning qiymatlarini hisoblash
for i in range(2,n+1):
    delta1[i]=(-gamma2[i]+delta2[i])/(2*(i-1))

#delta0 ning qiymatlarini hisoblash
for i in range(1,n+1):
    delta0[i]=(-gamma1[i]+delta1[i])/(2*i)

p1=list([0 for i in range(n+1)])
p0=list([0 for i in range(n+1)])

for i in range(4,n+1):
    p1[i]=-delta2[i]/(2*(i-3))

for i in range(3,n+1):
    p0[i]=(-delta1[i]+p1[i])/(2*(i-2))

eta0=list([0 for i in range(n+1)])
for i in range(5,n+1):
    eta0[i]=-p1[i]/(2*(i-4))

f3 = [[0 for i in range(n+1)] for j in range(n+5)]
f2 = [[0 for i in range(n+1)] for j in range(n+5)]
f1 = [[0 for i in range(n+1)] for j in range(n+5)]
f0 = [[0 for i in range(n+1)] for j in range(n+5)]

#f0,f1, f2,f3 matritsalarni hosil qilish
for j in range(0,n+5):
    for i in range(0,n+1):
        f3[j][i]=delta(j,i+1)*beta3[i]+delta(j,i-1)*gamma3[i]
        f2[j][i]=delta(j,i+2)*beta2[i]+delta(j,i)*gamma2[i]+delta(j,i-2)*delta2[i]
        f1[j][i]=delta(j,i+3)*beta1[i]+delta(j,i+1)*gamma1[i]+delta(j,i-1)*delta1[i]+delta(j,i-3)*p1[i]
        f0[j][i]=delta(j,i+4)*beta0[i]+delta(j,i+2)*gamma0[i]+delta(j,i)*delta0[i]+delta(j,i-2)*p0[i]+delta(j,i-4)*eta0[i]

sigma0=list([0 for i in range(n+1)])
sigma1=list([0 for i in range(n+1)])
sigma2=list([0 for i in range(n+1)])
sigma3=list([0 for i in range(n+1)])

sigmachiz0=list([0 for i in range(n+1)])
sigmachiz1=list([0 for i in range(n+1)])
sigmachiz2=list([0 for i in range(n+1)])
sigmachiz3=list([0 for i in range(n+1)])

for i in range(n+1):
    for j in range(n+5):
        sigma0[i]+=f0[j][i]
        sigmachiz0[i]+=((-1)**j)*f0[j][i]

for i in range(n+1):
    for j in range(n+4):
        sigma1[i]+=f1[j][i]
        sigmachiz1[i]+=((-1)**j)*f1[j][i]

for i in range(n+1):
    for j in range(n+3):
        sigma2[i]+=f2[j][i]
        sigmachiz2[i]+=((-1)**j)*f2[j][i]

for i in range(n+1):
    for j in range(n+2):
        sigma3[i]+=f3[j][i]
        sigmachiz3[i]+=((-1)**j)*f3[j][i]

#g0 va g1 matritsalarni hosil qilish
g0 = [[0 for i in range(n+1)] for j in range(n+1)]
g1 = [[0 for i in range(n+1)] for j in range(n+1)]
g2 = [[0 for i in range(n+1)] for j in range(n+1)]
g3 = [[0 for i in range(n+1)] for j in range(n+1)]

for j in range(0,n+1):
    for i in range(0,n+1):
        g3[j][i]=f3[j][i] + delta(j,0)*1.5*((sigma0[i]-sigmachiz0[i])-(sigma1[i]+sigmachiz1[i]))
        g2[j][i]=f2[j][i] + delta(j,1)*(-0.5*(sigma1[i]-sigmachiz1[i])) + delta(j,0)*1.5*((sigma0[i]-sigmachiz0[i])-(sigma1[i]+sigmachiz1[i]))
        g1[j][i]=f1[j][i]+delta(j,2)*0.25*1.5*((sigma0[i]-sigmachiz0[i])-(sigma1[i]+sigmachiz1[i])) + delta(j,1)*(-0.5*(sigma1[i]-sigmachiz1[i])) + delta(j,0)*(-0.25*(3*(sigma0[i]-sigmachiz0[i])-(sigma1[i]+sigmachiz1[i])))
        g0[j][i]=f0[j][i]+delta(j,3)*1.5*((sigma0[i]-sigmachiz0[i])-(sigma1[i]+sigmachiz1[i]))/24 + delta(j,2)*(-0.5*(sigma1[i]-sigmachiz1[i]))/4 + delta(j,1)*( (-0.25*(3*(sigma0[i]-sigmachiz0[i])-(sigma1[i]+sigmachiz1[i]))) - 1.5*((sigma0[i]-sigmachiz0[i])-(sigma1[i]+sigmachiz1[i]))/8 )+delta(j,0)*0.5*(-(sigma0[i]+sigmachiz0[i])+0.25*(sigma1[i]-sigmachiz1[i]))

#g matritsani hosil qilish
g = [[0 for i in range(n+1)] for j in range(n+1)]

for i in range(0,n+1):
    for k in range(0,n+1):
        g[i][k]= eps*delta(i,k)-2*g2[i][k]+g0[i][k]

#Tenglamalar sistemasini yechib a vektorni topish
a=numpy.linalg.solve(g,b)
print('a vektor:')
print(a)

c1,c2,c3,c4=0,0,0,0
for i in range(0,n+1):
    c1+=1.5*( (sigma0[i]-sigmachiz0[i]) - (sigma1[i]+sigmachiz1[i]) )*a[i]
    c2+=-0.5*( sigma1[i] - sigmachiz1[i] )*a[i]
    c3+=-0.25*( 3*(sigma0[i]-sigmachiz0[i]) + (sigma1[i]+sigmachiz1[i]) )*a[i]
    c4+=0.5*( -(sigma0[i]+sigmachiz0[i])+0.25*(sigma1[i]-sigmachiz1[i]) )*a[i]

print('c1=',c1)
print('c2=',c2)
print('c3=',c3)
print('c4=',c4)

u=list([0 for i in range(n+1)])
u1=list([0 for i in range(n+1)])
u2=list([0 for i in range(n+1)])
u3=list([0 for i in range(n+1)])
u4=list([0 for i in range(n+1)])

for l in range(n+1):
    u[l]=c1*T(3,y[l])/24 + c2*T(2,y[l])/4 + (c3-c1/8)*T(1,y[l])+c4*T(0,y[l])
    for j in range(0,n+5):
        for i in range(n+1):
            u[l]+=f0[j][i]*a[i]*T(j,y[l])

for l in range(n+1):
    u1[l]=c1*T(2,y[l])/4 + c2*T(1,y[l])+c3*T(0,y[l])
    for j in range(0,n+4):
        for i in range(n+1):
            u1[l]+=f1[j][i]*a[i]*T(j,y[l])

for l in range(n+1):
    u2[l]=c1*T(1,y[l])+c2*T(0,y[l])
    for j in range(0,n+3):
        for i in range(n+1):
            u2[l]+=f2[j][i]*a[i]*T(j,y[l])

for l in range(n+1):
    u3[l]=c1*T(0,y[l])
    for j in range(0,n+2):
        for i in range(n+1):
            u3[l]+=f3[j][i]*a[i]*T(j,y[l])

for l in range(n+1):
    for j in range(n+1):
        u4[l]+=a[j]*T(j,y[l])


ua=list([0 for i in range(n+1)])
ua1=list([0 for i in range(n+1)])

for l in range(n+1):
    ua[l]=(1-y[l]**2)**2*math.exp(eps*y[l])

plt.plot(u, color='blue', label=r'$mmm$')
plt.plot(ua,":", label=r'$mmm$')
plt.show()

for l in range(0,n+1):
    print('l=',l,' Aniq yechim=',ua[l], ' Taqribiy yechim=', u[l], ' Xatolik=', math.fabs(ua[l]-u[l]))

for l in range(n+1):
    ua1[l]=(y[l]-1)*(eps*y[l]**2+4*y[l]-eps)*math.exp(eps*y[l])

plt.plot(u1, "r--", label=r'$mmm$')
plt.plot(ua1, ":", label=r'$mmm$')
plt.show()

plt.plot(u2, ":", label=r'$mmm$')


plt.plot(u3, "-", label=r'$mmm$')


plt.plot(u4, ":", label=r'$mmm$')
plt.show()

