# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 21:44:43 2019

@author: Elvin
"""
import random
def MillerRabinTest(num, count):
    if num <= 3:
        raise ValueError("Too small number for Miller-Rabin Primality Test")
    print("Miller-Rabin Primality Test for", num, "number")

    r, s = 0, num - 1
    while s % 2 == 0:
        r, s = r + 1, s//2
    for _ in range(count):
        base = random.randint(2, num - 2)
        x = pow(base, s, num)
        if x == 1 or x == num - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, num)
            if x == 1:
                print("base =", base, "composite", "x == 1")
                return
            elif x == num - 1:
                break
        else:
            print("base =", base, "composite", "x == -1")
            return
    print("base =", base, "probably simple")
    return
def jacobi(a, n):
    s = 1
    while True:
        if n < 1: raise ValueError("Too small module for Jacobi symbol: " + str(n))
        if n & 1 == 0: raise ValueError("Jacobi is defined only for odd modules")
        if n == 1: return s
        a = a % n
        if a == 0: return 0
        if a == 1: return s

        if a & 1 == 0:
            if n % 8 in (3, 5):
                s = -s
            a >>= 1
            continue

        if a % 4 == 3 and n % 4 == 3:
            s = -s

        a, n = n, a
    return
def gcd(a, b):
    x=[]
    y=[]
    q=[]
    r=[]
    if(a>b):
        r=[a,b]
        x=[1,0]
        y=[0,1]
        q=[0,int(a/b)]
        
    else:
        r=[b,a]
        x=[1,0]
        y=[0,1]
        q=[0,int(b/a)]
    i=1
    #and i<10
    while (r[i]!=0) :
        r.append(r[i-1]-r[i]*q[i])
        x.append(x[i-1]-q[i]*x[i])
        y.append(y[i-1]-q[i]*y[i])
        if(r[i+1]!=0):
            q.append(int(r[i]/r[i+1]))
        i=i+1
    return r[i-1]
def lcf(a, b):
    x=[]
    y=[]
    q=[]
    r=[]
    if(a>b):
        r=[a,b]
        x=[1,0]
        y=[0,1]
        q=[0,int(a/b)]
        
    else:
        r=[b,a]
        x=[1,0]
        y=[0,1]
        q=[0,int(b/a)]
    i=1
    #and i<10
    while (r[i]!=0) :
        r.append(r[i-1]-r[i]*q[i])
        x.append(x[i-1]-q[i]*x[i])
        y.append(y[i-1]-q[i]*y[i])
        if(r[i+1]!=0):
            q.append(int(r[i]/r[i+1]))
        i=i+1
    return(r,x,y,q,i)
def SolovayStrassenTest(num, count):
    k=0
    counter=0
    if num <= 3:
         raise ValueError("Too small number for Solovay-Strassen Primality Test")  
    print("Solovay-Strassen Primality Test for", num, "number")
           
    for i in range(count):
        counter=counter+1
        if k==1:
            break
        base = random.randint(2, num - 2)
        d = gcd(num, base) 
        if d != 1:
            print("base =", base, "composite", "GCD condition")
        
        # calculate legendre symbol from euler criterion formula  
        y = pow(base, (num - 1) // 2, num)
        x = jacobi(base, num)
        if y != x % num:
            k=1
            print("base =", base, "composite", "Jacobi condition")
            print(base,"^",(num - 1) // 2,"mod",num,"=",pow(base, (num - 1) // 2, num))
        else:
            print("base =", base, "probably simple")
            print(base,"^",(num - 1) // 2,"mod",num,"=",pow(base, (num - 1) // 2, num))

    return
def FermatTest(num, count):
    k=0
    counter=0
    if num <= 3:
        raise ValueError("Too small number for Fermat Primality Test")
        print("Fermat Primality Test for", num, "number") 
	# check number for 10 different bases
    for i in range(count):
        counter=counter+1
        if k==1:
            break
        base = random.randint(2, num - 2)
        if pow(base, num - 1, num) != 1:
            k=1
            print("base =", base, "composite")
            print(base,"^",num-1,"mod",num,"=",pow(base, num - 1, num))
        else:
            if counter<=5:
                print("base =", base, "probably simple")
                print(base,"^",num-1,"mod",num,"=",pow(base, num - 1, num))
    return 

#FermatTest(39932034899759958979,20)
#testing if carlmichel is simple





def CarlMickleTest(num, count):#проверка на пседослючайность
    k=0
    counter=0
    if num <= 3:
        raise ValueError("Too small number for Fermat Primality Test")
        print("Fermat Primality Test for", num, "number") 
	# check number for 10 different bases
    for i in range(count):
        counter=counter+1
        if k==1:
            break
        base = random.randint(2, num - 2)
        #base=counter+1
        r,x,y,q,i=lcf(base,num)
        if r[i-1]==1:
        #base = random.randint(2, num - 2)
            if pow(base, num - 1, num) != 1:
                k=1 
                print("base =", base, "composite")
                print(base,"^",num-1,"mod",num,"=",pow(base, num - 1, num))
            else:
                if counter<5 or counter>num-6:
                    print("base =", base, "probably simple")
                    print(base,"^",num-1,"mod",num,"=",pow(base, num - 1, num))
    return 

def SolovayStrassenTest12(num, count):
    if num <= 3:
         raise ValueError("Too small number for Solovay-Strassen Primality Test")  
    print("Solovay-Strassen Primality Test for", num, "number")
           
    for i in range(count):
        base = random.randint(2, num - 2)
        
        d = gcd(num, base) 
        if d != 1:
            print("base =", base, "composite", "GCD condition failed")
        
        # calculate legendre symbol from euler criterion formula  
        y = pow(base, (num - 1) // 2, num)
        x = jacobi(base, num)
        if y != x % num:
            print("base =", base, "composite", "Jacobi condition failed")
        else:
            print("base =", base, "probably simple")
    return



#FermatTest(2810864562635368426005268142616001,20)
#FermatTest(7947876367161318130381303018250493665173,20)
#FermatTest(510241663229783503211079698102859446339,20)
#FermatTest(2819632050530705041831693710220059046232636245640560854597862247156896041042889,20)
#CarlMickleTest(379382381447399527322618466130154668512652910714224209601,200)
#SolovayStrassenTest(2819632050530705041831693710220059046232636245640560854597862247156896041042889,10)
MillerRabinTest(2810864562635368426005268142616001,100)



#r,x,y,q,i=lcf(26,30)
#print(r[i-1])

#def stepen(c,k):
#    b=1
#    while(k!=0):
#        if k%2==0:
#            k=k//2
#            c=c*c
#        else:
#            k=k-1
#            b=b*c
#    return b,c     
#    
#print(stepen(15,10258460621399486474))