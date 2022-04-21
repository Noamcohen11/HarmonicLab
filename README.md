# HarmonicLab
Matlab code for all graphs needed on lab 2 - damped harmonic acceleration
The lab has an automatic distance measuring machine a string and weights. 
When the weights are added to the string a damped harmonic movment accures. 

## The graphes in this repo include:

* [ Resting position vs potential height energy](#-Resting-position-vs-potential-height-energy)
* [ Time vs distnace for each mass](#Time-vs-distnace-for-each-mass)

## Resting position vs potential height energy

The goal of this graph is to calculate the k constant of the string. 
Since the equation is  ```mg = k*(x-x0)```.
In the resting position x = 0 so  ```mg = k*(-x0)```.
As a result the k constant is calulated as the graph's slope.

##  Time vs distnace for each mass

This code generates graphs for each of the weights messured. 
It also fits using the damped harmonic formula:
```a*exp(b*x)*(cos(2*pi*x/c + 2*pi/d))+e```
