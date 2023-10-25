#!/usr/bin/env python
# coding: utf-8

# # ![image.png](attachment:image.png)

# In[1]:


from sympy.interactive import printing
printing.init_printing(use_latex=True)
from sympy import *
from sympy.solvers import solve
from time import time


import sympy as sp

start1=time()
a11,a12,a13,a21,a22,a23,a31,a32,a33=map(int,input("Enter the matrix coefficients : ").split())
r1,r2,r3=[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]
amat=sp.Matrix((r1,r2,r3))
print("Input matrix")
display(amat)
    

s1=a11+a22+a33
display('S1={0}+{1}+{2}'.format(a11,a22,a33))
display('s1={0}'.format(s1))
s2=(((a22*a33)-(a32*a23))+((a11*a33)-(a31*a13))+((a11*a22)-(a21*a12)))
display('s2=((({0}*{1})-({2}*{3}))+(({4}*{5})-({6}*{7}))+(({8}*{9})-({10}*{11}))'.format(a22,a33,a32,a23,a11,a33,a31,a13,a11,a22,a21,a12))
display('s2=({0}+{1}+{2})'.format(((a22*a33)-(a32*a23)),((a11*a33)-(a31*a13)),((a11*a22)-(a21*a12))))
s3=det(amat)
display('s3=det(amat)')
display("s3={0}".format(s3))
display('s1={0} s2={1} s3={2}'.format(s1,s2,s3))
x=sp.symbols('λ')
eq1=sp.Function("eq1")
eq1=Eq((x**3 - (s1)*x**2 +(s2)*x- s3))
display(eq1)
k=sp.factor(eq1)
display(k)
K=solve(k)
print("Eigen values are : ",K)

mat2=Matrix(([x,0,0],[0,x,0],[0,0,x]))

a,b,c,d,e,f,mi=sp.symbols('x1 x2 x3 * = -> -')
print(">>>>>>>>>>>>>>>>>> To find eigen vector <<<<<<<<<<<<<<<")
print("({0}-{1}{2}){3}={4}".format('A',x,'I','X',0))
cmat=sp.Matrix((0,0,0))
display(amat,mi,mat2,e,cmat,f)
display(amat-mat2)
xmat=sp.Matrix((a,b,c))
display(d,xmat,e,cmat)


for i in range(0,len(K)):
    print("\n\nIf λ={0}".format(K[i]))
    print("({0}-{1}{2}){3}={4}".format('A',K[i],'I','X',0))
    lmat1=Matrix(([K[i],0,0],[0,K[i],0],[0,0,K[i]]))
    display(amat,mi,lmat1,e,amat-lmat1)
    print('now')
    display(amat-lmat1,d,xmat,e,cmat,f)

    eq2=sp.Function("eq2")
    eq3=sp.Function("eq3")
    eq4=sp.Function("eq4")

    eq2=Eq((a11-K[i])*a+a12*b+a13*c)
    eq3=Eq(a21*a+(a22-K[i])*b+a23*c)
    eq4=Eq(a31*a+a32*b+(a33-K[i])*c)

    display(eq2,eq3,eq4)


    b1,b2,b3=[(a11-K[i]),a12,a13,0],[a21,(a22-K[i]),a23,0],[a31,a32,(a33-K[i]),0]

    system=Matrix((b1,b2,b3))
    print('\nBy Gaussian Elimination Method.....\n')
    n=sp.solve_linear_system(system,a,b,c)
    display(n)
    if(eq2!=eq3 and eq3!=eq4):
        print("\nBy Using Cofactor Method\n")
        xm=((a12*a23)-((a22-K[i])*a13))
        ym=(((a11-K[i])*a23)-(a21*a13))
        zm=(((a11-K[i])*(a22-K[i]))-(a21*a12))
        if(xm==0 and ym==0 and zm==0):
            subst=input("Which variable you want to assume the value : ")
            if(subst=='x1' or subst=='X1'):
                xm=1
            if(subst=='x2' or subst=='X2'):
                ym=-1
            if(subst=='x3' or subst=='X3'):
                zm=1
        print('({0}/{1})=(-({2}/{3}))=({4}/{5})'.format(a,xm,b,ym,c,zm))
        if(xm==(-ym)==zm):
            if(xm==(-ym)==zm==0):
                v1=Matrix([xm,-ym,zm])
            else:
                v1=Matrix([xm/xm,-ym/(-ym),zm/zm])
        elif((xm%2==0 and -ym%2==0 and zm%2==0)):
            v1=Matrix([xm/2,-ym/2,zm/2])
            if(xm!=1 and -ym!=1 and zm!=1):
                xm/=2
                ym/=2
                zm/=2
        elif((xm%3==0 and -ym%3==0 and zm%3==0)):
            v1=Matrix([xm/3,-ym/3,zm/3])
            if(xm!=1 and -ym!=1 and zm!=1):
                xm/=3
                ym/=3
                zm/=3
        elif((xm%5==0 and -ym%5==0 and zm%5==0)):
            v1=Matrix([xm/5,-ym/5,zm/5])
            if(xm!=1 and -ym!=1 and zm!=1):
                xm/=5
                ym/=5
                zm/=5
        elif((xm%7==0 and -ym%7==0 and zm%7==0)):
            v1=Matrix([xm/7,-ym/7,zm/7])
            if(xm!=1 and -ym!=1 and zm!=1):
                xm/=7
                ym/=7
                zm/=7
        else:
            v1=Matrix([xm,-ym,zm])
        print("\n\nEigen vector when {0}={1} : ".format(x,K[i]))
        display(v1)
    elif(eq2==eq3==eq4):
        print("\nBy Substitution Method.....\n")
        print('\nLet X2=0\n')
        ys=0
        xm=eq2.subs(b,ys)
        display(xm)
        display(solve(xm))
        xs=int(input('X1 substitution : '))
        display(xm.subs(a,xs))
        display(solve(xm.subs(a,xs)))
        zs=solve(xm.subs(a,xs))
        print('C=',*zs)
        v2=Matrix([xs,ys,*zs])
        print("\nEigen vectors when {0}={1} and when X2=0 \n".format(x,K[0]))
        display(v2)
        print('\n\nLet C=0\n')
        zs=0
        xm=eq2.subs(c,zs)
        display(xm)
        display(solve(xm))
        xs=int(input('X1 substitution : '))
        display(xm.subs(a,xs))
        display(solve(xm.subs(a,xs)))
        ys=solve(xm.subs(a,xs))
        print('B=',*ys)
        v3=Matrix([xs,*ys,zs])
        display(v3)
        print("\nEigen vectors when {0}={1} and when X3=0 \n".format(x,K[0]))
        display(v3)
end1=time()
print("Execution Time")
display(end1-start1)
print("Using Built in function")
start=time()
ev=amat.eigenvals()
j=amat.eigenvects()
print("Eigen value")
display(*ev)
print("Eigen Vector")
display(j)
end=time()
print("Execution time ")
display(end-start)




# 
# 

# In[4]:


y=int(input("How many sums going to check "))
for i in range(0,y):
    dim=int(input("Enter the dimension : "))
    if(dim==2):
        a11,a12,a21,a22=map(int,input("Enter the matrix coefficients : ").split())
        r1,r2=[a11,a12],[a21,a22]
        amat=sp.Matrix((r1,r2))
        print("Input matrix")
        display(amat)
        s1=a11+a22
        display('S1={0}+{1}'.format(a11,a22))
        display('s1={0}'.format(s1))
        display('s2=({0}*{1})-({2}*{3})'.format(a11,a22,a12,a21))
        s2=(a11*a22)-(a12*a21)
        display('s2={0}'.format(s2))
        x=sp.symbols('λ')
        eq1=sp.Function("eq1")
        eq1=Eq((x**2 +(s1)*x- s2))
        display(eq1)
        k=sp.factor(eq1)
        display(k)
        K=solve(k)
        print("Eigen values are : ",*K)
    if(dim==3):
        a11,a12,a13,a21,a22,a23,a31,a32,a33=map(int,input("Enter the matrix coefficients : ").split())
        r1,r2,r3=[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]
        amat=sp.Matrix((r1,r2,r3))
        print("Input matrix")
        display(amat)


        s1=a11+a22+a33
        display('S1={0}+{1}+{2}'.format(a11,a22,a33))
        display('s1={0}'.format(s1))
        s2=(((a22*a33)-(a32*a23))+((a11*a33)-(a31*a13))+((a11*a22)-(a21*a12)))
        display('s2=((({0}*{1})-({2}*{3}))+(({4}*{5})-({6}*{7}))+(({8}*{9})-({10}*{11}))'.format(a22,a33,a32,a23,a11,a33,a31,a13,a11,a22,a21,a12))
        display('s2=({0}+{1}+{2})'.format(((a22*a33)-(a32*a23)),((a11*a33)-(a31*a13)),((a11*a22)-(a21*a12))))
        s3=det(amat)
        display('s3=det(amat)')
        display("s3={0}".format(s3))
        display('s1={0} s2={1} s3={2}'.format(s1,s2,s3))
        x=sp.symbols('λ')
        eq1=sp.Function("eq1")
        eq1=Eq((x**3 - (s1)*x**2 +(s2)*x- s3))
        display(eq1)
        k=sp.factor(eq1)
        display(k)
        K=solve(k)
        print("Eigen values are : ",K)
    P,D=amat.diagonalize()
    print("_____Diagonalization_____")
    display(P,D,P**-1)


# In[ ]:




