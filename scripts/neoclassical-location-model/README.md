# The Neoclassical Model of Location, Albouy and Stuart (2020)

The Neoclassical Model of Location, as set up by Albouy and Stuart ([doi](https://doi.org/10.1111/iere.12419), [working paper](https://bryan-stuart.com/files/UrbanPopulation_Feb2019.pdf)), is a general equilibrium model with 2 output sectors-- one traded good and one home (housing) good-- as well as 3 factors of production-- land, labor, and capital. The full nonlinear model is shown in Supp. Appendix B and C. Because Julia has very powerful libraries for solving all kinds of systems of equations, I decided to sharpen my skills and write a function that solves the model for the 15 endogeneous variables, subject to 15 equilibrium conditions.

-   [This script](https://github.com/rcberg/macroeconomics/blob/main/scripts/neoclassical-location-model/large_city_small_city_simulator.jl) uses JuMP to set up the nonlinear system from Appendix B and C from the paper. Following the 'MPEC' method, the equilibrium conditions are expressed as constraints. For the large city baseline, the 'Ipopt' solver optimizes a loss function (the sum of the squared percentage difference from the target) of the CES parameters to match Table 1 targets for shares, subject to the equilibrium constraints. The important large city values are supplied for the small city, whose population and local prices are endogeneously determined. This project benefited from [Gemini 3 Pro](https://gemini.google.com/) to help better understanding aspects of the Julia language (and JuMP library).

## To-do list:

(Nothing as of now.)