# A-multivariate-numerical-integration-program-with-a-two-dimensional-Laplacian-operator-in-R
This repository stores the codes and data of our article (in submission) as well.

For many ecologists and environmental scientists, R is still the top priority data analysis tool. deSolve Package provides a relatively ideal integrator in R for a large group of R programmers who occasionally need to build some differential systems. To my knowledge, deSolve in R accompany with Intel Math Kernel Library (MKL) is shown to be able to meet most tasks that are not very computationally intensive, which help to save the time from learning other widely used but complex mathematical tools such as Matlab, Mathematica or general programming laguages such as Python and Julia, etc.

During researching in the present paper, I noticed that the profile of deSolve provided instruction in designing systems with two-dimensional Laplacian operators, which simulate the process of diffusion or dispersal of natural materials or populations. However, it becomes a problem when we need to arbitrarily change the number of variables. In other words, the algorithm given in the profile cannot deal with systems whose input and output are vectors of any dimension.

I expanded the usage of Laplacian operator to multivariate cases by manipulating arrays in R in these codes, and with them one can easily change the specific form of multivariate systems and inquire how the systems perform beautifully in a two-dimensional space.

It is also notable that this reprository is used to store the codes and data of our recent article in submission.
