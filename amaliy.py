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


eps=0.001


N=int(input('N='))

c=list([1 for i in range(N+1)])
c[0]=c[N]=2

b=list([0 for i in range(N+1)])

y=list([0 for i in range(N+1)])
for l in range(N+1):
    y[l]=math.cos(math.pi*l/N)

def f(y):
    d = (1 - y ** 2) ** 2 * math.exp(eps * y)
    d2 = (eps ** 2 * y ** 4 + 8 * eps * y ** 3 + (12 - 2 * eps ** 2) * y ** 2 - 8 * eps * y + eps ** 2 - 4) * math.exp(
        eps * y)
    d4 = (eps ** 4 * y ** 4 + 16 * eps ** 3 * y ** 3 + (72 * eps ** 2 - 2 * eps ** 4) * y ** 2 + (
            96 * eps - 16 * eps ** 3) * y + eps ** 4 - 24 * eps ** 2 + 24) * math.exp(eps * y)
    return eps * d4 - 2 * d2 + d
for i in y:
    print(f(i))
print("#" * 10)
#b vektorning qiymatlarini hisoblash
for i in range(0,N+1):
    s=0
    for l in range(0,N+1):
        s=s+2*f(y[l])*T(i,y[l])/c[l]
    b[i]=s/(N*c[i])

for i in range(N-3,N+1):
    b[i]=0
#a matritsani hisoblash
a= [[0 for i in range(N+1)] for j in range(N+1)]

for n in range(0,N-3):
    a[n][n]=c[n]
    for p in range(n+4,N+1):
        if (p-n)%2==0:
            a[n][p]=p*( p**2*(p**2-4)**2 - 3*n**2*p**4 + 3*n**4*p**2 -n**2*(n**2-4)**2 )*eps/(24*c[n])
    for p in range(n+2,N+1):
        if (p-n)%2==0:
            a[n][p]+=-2*p*(p**2-n**2)/c[n]

for i in range(0,N+1):
    a[N-3][i]=(-1)**i
    a[N - 2][i] = (-1) ** (i - 1) * i ** 2
    a[N-1][i] = 1
    a[N][i]=i**2


# for i in range(0,N+1):
#     for j in range(0,N+1):
#         print(a[i][j], end=' ')
#     print()

aa=numpy.linalg.solve(a,b)


#funksiya
ut=list([0 for i in range(N+1)])
u=list([0 for i in range(N+1)])

for l in range(0,N+1):
    for i in range(0,N+1):
        ut[l]+=aa[i]*T(i,y[l])

for l in range(N+1):
    u[l]=(1-y[l]**2)**2*math.exp(eps*y[l])

print('Funksiyaning ozi')
for l in range(0,N+1):

    print('l=',l,' Aniq yechim=',u[l], ' Taqribiy yechim=', ut[l], ' Xatolik=', math.fabs(u[l]-ut[l]))

plt.plot(ut,":", label=r'$Taqribiy-yechim$')
plt.plot(u, label=r'$Aniq-yechim$')
plt.title(r'$Funksiya-yechimlar-grafigi$')
plt.ylabel('U')
plt.xlabel('l')
plt.legend(loc='best', fontsize=10)
plt.show()
#funksiya_tugadi



#1-tartibli hosila
a1=list([0 for i in range(N+1)])
for n in range(0,N+1):
    for p in range(n+1,N+1):
        if (p+n-1)%2==0:
            a1[n]+=2*p*aa[p]/c[n]

ut1=list([0 for i in range(N+1)])
for l in range(0,N+1):
    for i in range(0,N+1):
        ut1[l]+=a1[i]*T(i,y[l])

u1=list([0 for i in range(N+1)])
for l in range(N+1):
    u1[l]=(y[l]**2-1)*(eps*y[l]**2 +4*y[l]-eps)*math.exp(eps*y[l])

print('1-tartibli hosila')
for l in range(0,N+1):
    print('l=',l,' Aniq yechim=',u1[l], ' Taqribiy yechim=', ut1[l], ' Xatolik=', math.fabs(u1[l]-ut1[l]))
plt.plot(ut1,":", label=r'$Taqribiy-yechim$')
plt.plot(u1, label=r'$Aniq-yechim$')
plt.title(r'$Birinchi-tartibli-hosila$')
plt.ylabel('U1')
plt.xlabel('l')
plt.legend(loc='best', fontsize=10)
plt.show()
#1-tartibli hosila tugadi


#2-tartibli hosila
a2=list([0 for i in range(N+1)])
for n in range(0,N+1):
    for m in range(n+1,N+1):
        if (m+n-1)%2==0:
            a2[n]+=2*m*a1[m]/c[n]

ut2=list([0 for i in range(N+1)])
for l in range(0,N+1):
    for i in range(0,N+1):
        ut2[l]+=a2[i]*T(i,y[l])

u2=list([0 for i in range(N+1)])
for l in range(N+1):
    u2[l]=(eps ** 2 * y[l] ** 4 + 8 * eps * y[l] ** 3 + (12 - 2 * eps ** 2) * y[l] ** 2 - 8 * eps * y[l] + eps ** 2 - 4) * math.exp(
        eps * y[l])
print('2-tartibli hosila')
for l in range(0,N+1):
    print('l=',l,' Aniq yechim=',u2[l], ' Taqribiy yechim=', ut2[l], ' Xatolik=', math.fabs(u2[l]-ut2[l]))
plt.plot(ut2,":", label=r'$Taqribiy-yechim$')
plt.plot(u2, label=r'$Aniq-yechim$')
plt.title(r'$Ikkinchi-tartibli-hosila$')
plt.ylabel('U2')
plt.xlabel('l')
plt.legend(loc='best', fontsize=10)
plt.show()
#2-tartibli hosila tugadi




#3-tartibli hosila
a3=list([0 for i in range(N+1)])
for n in range(0,N+1):
    for p in range(n+3,N+1):
        if (p+n-1)%2==0:
            a3[n]+=(1/(4*c[n]))*( p*( (p**2-1)**2-2*(p**2+1)*n**2+n**4)*aa[p]  )

ut3=list([0 for i in range(N+1)])
for l in range(0,N+1):
    for i in range(0,N+1):
        ut3[l]+=a3[i]*T(i,y[l])

u3=list([0 for i in range(N+1)])
for l in range(N+1):
    u3[l]=(eps**3*y[l]**4 + 12*eps**2*y[l]**3 + (36*eps-2*eps**2)*y[l]**2 + (24-12*eps**2)*y[l] + eps**3 -12*eps  )*math.exp(eps*y[l])
print('3-tartibli hosila')
for l in range(0,N+1):
    print('l=',l,' Aniq yechim=',u3[l], ' Taqribiy yechim=', ut3[l], ' Xatolik=', math.fabs(u3[l]-ut3[l]))
plt.plot(ut3,":", label=r'$Taqribiy-yechim$')
plt.plot(u3, label=r'$Aniq-yechim$')
plt.title(r'$Uchinchi-tartibli-hosila$')
plt.ylabel('U3')
plt.xlabel('l')
plt.legend(loc='best', fontsize=10)
plt.show()
#3-tartibli hosila tugadi


#4-tartibli hosila
a4=list([0 for i in range(N+1)])
for n in range(0,N+1):
    for p in range(n+3,N+1):
        if (p+n-1)%2==0:
            a4[n]+=(1/(24*c[n]))*( p*( p**2*(p**2-4)**2 - 3*n**4*p**4 + 3*n**4*p**2 - n**2*(n**2-4)**2 )*aa[p]  )

ut4=list([0 for i in range(N+1)])
for l in range(0,N+1):
    for i in range(0,N+1):
        ut4[l]+=a4[i]*T(i,y[l])

u4=list([0 for i in range(N+1)])
for l in range(N+1):
    u4[l]=(eps ** 4 * y[l] ** 4 + 16 * eps ** 3 * y[l] ** 3 + (72 * eps ** 2 - 2 * eps ** 4) * y[l] ** 2 + (
            96 * eps - 16 * eps ** 3) * y[l] + eps ** 4 - 24 * eps ** 2 + 24) * math.exp(eps * y[l])
print('4-tartibli hosila')
for l in range(0,N+1):
    print('l=',l,' Aniq yechim=',u4[l], ' Taqribiy yechim=', ut4[l], ' Xatolik=', math.fabs(u4[l]-ut4[l]))
plt.plot(ut4,":", label=r'$Taqribiy-yechim$')
plt.plot(u4, label=r'$Aniq-yechim$')
plt.title(r'$Tortinchi-tartibli-hosila$')
plt.ylabel('U4')
plt.xlabel('l')
plt.legend(loc='best', fontsize=10)
plt.show()
#4-tartibli hosila tugadi

