# Lab 11
### The Diffusion equation


In this lab we consider the diffusion equation
<p align="center">
<img src="stuffy_stuff/f1.png" width="100">
</p>
where *D* is a positive constant.

In the lecture we derived the implicit algorithm
<p align="center">
<img src="stuffy_stuff/f2.png" width="300">
</p>
which may be written as a system of linear equations
<p align="center">
<img src="stuffy_stuff/f4.png" width="130">
</p>
where the coefficient matrix *M* is given by

<p align="center">
<img src="stuffy_stuff/f3.png" width="350">
</p>

In the given source code implement the function **step()** that solves this system of linear equations in each time step.

Note that the matrix is tri-diagonal and that it is sufficient to store the entries on the three diagonals in vectors.
