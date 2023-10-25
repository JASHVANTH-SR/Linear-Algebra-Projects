#!/usr/bin/env python
# coding: utf-8

# ![image.png](attachment:image.png)

# In[95]:


from sympy.interactive import printing
printing.init_printing(use_latex=True)
from sympy import *
import sympy as sp
from sympy.physics.quantum import Bra,Ket
mv1,mv2,mv3,eq=sp.symbols('||v1|| ||v2|| ||v3|| =')
su1,su2,su3,sv1,sv2,sv3=sp.symbols('U1 U2 U3 V1 V2 V3')
bv1,bv2,bv3=Bra('V1'),Bra('V2'),Bra('V3')
ku1,ku2,ku3=Ket('U1'),Ket('U2'),Ket('U3')
kv1,kv2,kv3=Ket('V1'),Ket('V2'),Ket('V3')

u1=Matrix(list(map(int,input('Vector 1 :').split())))
u2=Matrix(list(map(int,input('Vector 2 :').split())))
u3=Matrix(list(map(int,input('Vector 3 :').split())))

print('Let the given vectors be {u1,u2,u3}')
display((u1,u2,u3))

eq1=Eq(su1,sv1)
display(eq1)
display(u1)
v1=u1
print('For our convinience we first find self inner products of v1')
v1s=v1.dot(v1)
ev1s=Bra(v1)*Ket(v1)
display((mv1**2,eq,bv1*kv1,eq,ev1s,eq,v1s))
eq2=Eq(((su2-((bv2*ku1)/(mu1**2))*su1)),sv2)
display(eq2)
u2v1=u2.dot(v1)
display((bv1*ku2,eq,Bra(u2)*Ket(v1),eq,u2v1))                    
div1=u2v1/v1s
v2=u2-(div1*v1)
display(v2)
v2s=v2.dot(v2)
display((bv2*kv2,eq,Bra(v2)*Ket(v2),eq,v2s))
eq3=((su3-(((bv3*ku1)/(mu1**2))*su1)-((bv3*ku2)/(mu2**2))*su2))
display((eq3,eq,sv3))
u3v1=u3.dot(v1)
u3v2=u3.dot(v2)
display((bv3*ku1,eq,Bra(u3)*Ket(v1),eq,u3v1))                    
display((bv3*ku2,eq,Bra(u3)*Ket(v2),eq,u3v2))   
v2s=v2.dot(v2)
ev2s=Bra(v2)*Ket(v2)
display((mv2**2,eq,bv2*kv2,eq,ev2s,eq,v2s))
div2=u3v1/v1s
div3=u2v1/v2s
v3=u3-(div2*v1)-(div3*v2)
display(v3)
v3s=v3.dot(v3)
ev3s=Bra(v3)*Ket(v3)
display((mv3**2,eq,bv3*kv3,eq,ev3s,eq,v3s))
print('Orthogonal Vectors We Got')
display((v1,v2,v3))
print('Orthonormal Basis')

w1,w2,w3=sp.symbols('w1 w2 w3')
display((w1,eq,(sv1/mv1)))
display((w2,eq,(sv2/mv2)))
display((w3,eq,(sv3/mv3)))
W1=v1/sqrt(v1s)
W2=v2/sqrt(v2s)
W3=v3/sqrt(v3s)
display((w1,eq,W1))
display((w2,eq,W2))
display((w3,eq,W3))
print('The obtained Orthonormal Basis is ')
display((W1,W2,W3))


# # 

# In[ ]:




