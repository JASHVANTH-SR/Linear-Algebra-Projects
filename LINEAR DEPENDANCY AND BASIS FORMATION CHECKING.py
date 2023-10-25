#!/usr/bin/env python
# coding: utf-8

# In[25]:


from sympy.interactive import printing
printing.init_printing(use_latex=True)
from sympy import Eq,solve_linear_system,Matrix,det,transpose
from time import time

import sympy as sp

print("------------- Welcome! To Find linear combination dependancy and basis formation checking ---------------")
no=int(input("No of Times you going to execute the program : "))
for j in range(0,no):
    n=int(input('Number of vectors : '))
    m,a,e,g=sp.symbols('* + = -->')
    if n==2:
            x,y,m,a,e,g=sp.symbols("x1 x2 * + = -->")
            a11,a12=map(int,input("Enter 1st row : ").split())
            a21,a22=map(int,input("Enter 2nd row : ").split())
            vect1,vect2=[a11,a12],[a21,a22]
            display(x,m,vect1,a,y,m,vect2)
            amat=Matrix((vect1,vect2))
            print("Given Matrix")
            display(amat)
            print('now')
            display(amat,m,Matrix((x,y)),e,Matrix((0,0)),g)
            r1,r2=[a11,a21,0],[a12,a22,0]

            eq1,eq2,eqn,eqm=sp.Function("eq1"),sp.Function("eq2"),sp.Function("eqn"),sp.Function("eqm")
            eqm=(a11*x,a12*x)
            eqn=(a21*y,a22*y)
            display(eqm,a,eqn,e,(0,0))
            eq1,eq2=Eq(a11*x+a21*y,0),Eq(a12*x+a22*y,0)
            print("arranged matrix")
            display(transpose(amat))
            display(eq1,eq2)
            system=Matrix((r1,r2))
            print("augmented form")
            display(system)

            print("\n\n\nThe two vectors are ")
            if(det(amat)==0):
                print("linearly dependant")
                display(solve_linear_system(system,x,y))
                print("dim R2 != 2, The given vectors NOT form a basis of R^2")
            else:
                print("linearly independant")
                display(solve_linear_system(system,x,y))
                print("The given vectors does form a basis of R^2")
    if n==3:
            x,y,z=sp.symbols("x1 x2 x3")
            a11,a12,a13=map(int,input("Enter 1st row : ").split())
            a21,a22,a23=map(int,input("Enter 2nd row : ").split())
            a31,a32,a33=map(int,input("Enter 3rd row : ").split())
            vect1,vect2,vect3=[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]
            display(x,m,vect1,a,y,m,vect2,a,z,m,vect3,e,(0,0,0))
            amat=Matrix((vect1,vect2,vect3))
            print("Given Matrix")
            display(amat)
            print('now')
            display(amat,m,Matrix((x,y,z)),e,Matrix((0,0,0)),g)
            r1,r2,r3=[a11,a21,a31,0],[a12,a22,a32,0],[a13,a23,a33,0]
            eq1,eq2,eq3,eqn,eqm,eqo=sp.Function("eq1"),sp.Function("eq2"),sp.Function("eq3"),sp.Function("eqn"),sp.Function("eqm"),sp.Function("eqo")
            eqm=(a11*x,a12*x,a13*x)
            eqn=(a21*y,a22*y,a23*y)
            eqo=(a31*z,a32*z,a33*z)
            display(eqm,a,eqn,a,eqo,e,(0,0,0),g)
            eq1,eq2,eq3=Eq(a11*x+a21*y+a31*z,0),Eq(a12*x+a22*y+a32*z,0),Eq(a13*x+a23*y+a33*z,0)
            display(eq1,eq2,eq3)
            print("arranged matrix")
            display(transpose(amat))
            system=Matrix((r1,r2,r3))
            print("augmented form")
            display(system)

            print("\n\n\nThe two vectors are ")
            if(det(transpose(amat))==0):
                print("linearly dependant")
                display(solve_linear_system(system,x,y,z))
                print("dim R3 != 3, The given vectors NOT form a basis of R^3")
            else:
                print("linearly independant")
                display(solve_linear_system(system,x,y,z))
                print("The given vectors does form a basis of R^3")
    if n==4:
            x,y,z,w=sp.symbols("x1 x2 x3 x4")
            a11,a12,a13,a14=map(int,input("Enter 1st row : ").split())
            a21,a22,a23,a24=map(int,input("Enter 2nd row : ").split())
            a31,a32,a33,a34=map(int,input("Enter 3rd row : ").split())
            a41,a42,a43,a44=map(int,input("Enter 4th row : ").split())
            vect1,vect2,vect3,vect4=[a11,a12,a13,a14],[a21,a22,a23,a24],[a31,a32,a33,a34],[a41,a42,a43,a44]
            display(x,m,vect1,a,y,m,vect2,a,z,m,vect3,a,w,m,vect4)
            amat=Matrix((vect1,vect2,vect3,vect4))
            print("Given Matrix")
            display(amat)
            print('now')
            display(amat,m,Matrix((x,y,z,w)),e,Matrix((0,0,0,0)),g)
            r1,r2,r3,r4=[a11,a21,a31,a41,0],[a12,a22,a32,a42,0],[a13,a23,a33,a43,0],[a14,a24,a34,a44,0]
            eq1,eq2,eq3,eq4=sp.Function("eq1"),sp.Function("eq2"),sp.Function("eq3"),sp.Function("eq4")
            eqm=(a11*x,a12*x,a13*x,a14*x)
            eqn=(a21*y,a22*y,a23*y,a24*y)
            eqo=(a31*z,a32*z,a33*z,a34*z)
            eqp=(a41*z,a42*z,a43*z,a44*z)
            display(eqm,a,eqn,a,eqo,a,eqp,e,(0,0,0,0),g)
            eq1,eq2,eq3,eq4=Eq(a11*x+a21*y+a31*z+a41*w,0),Eq(a12*x+a22*y+a32*z+a42*w,0),Eq(a13*x+a23*y+a33*z+a43*w,0),Eq(a14*x+a24*y+a34*z+a44*w,0)
            display(eq1,eq2,eq3,eq4)
            system=Matrix((r1,r2,r3,r4))
            print('rearranged matrix')
            display(transpose(amat))
            print("augmented form")
            display(system)

            print("\n\n\nThe two vectors are ")
            if(det(transpose(amat))==0):
                print("linearly dependant")
                display(solve_linear_system(system,x,y,z,w))
                print("dim R4 != 4, The given vectors NOT form a basis of R^4")
            else:
                print("linearly independant")
                display(solve_linear_system(system,x,y,z,w))
                print("The given vectors does form a basis of R^4")


# In[ ]:




