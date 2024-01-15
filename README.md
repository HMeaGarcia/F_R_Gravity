# F_R_Gravity

This repository contains part of the implemented code for my master thesis titled *Non-linear analysis of solar system tests* in $f(R)$ gravity which can be found in  http://132.248.9.195/ptd2023/febrero/0835817/Index.html

The main objective of the code is to solve numerically the $f(R)$ field equations for a 
for a spherically symmetric and static (SSS) matter distribution.

$$
ds^2 =-n(r)dt^2+m(r)dr^2+  r^2 d\Omega^2
$$
<br/>

### F_R_Tools.py
   This file contains the main classes used for solving the RHS of the $f(R)$ field equations. 

   ```
   class STBS
   ```
Corresponds to the Staroninsky $f(R)$ model (Starobinsky, A. A. “Disappearing Cosmological Constant in f(R) Gravity”. JETP Letters, 86, (2007))

$$
    f(R)^\text{STBS}= R + \lambda R_S\left[\left(1+\frac{R^2}{R_S^2}\right)^{-q}-1\right]
$$
<br/>



### Additional Notes

The shooting method is not within this files. However, you can request the code by emailing me."


<br/>



