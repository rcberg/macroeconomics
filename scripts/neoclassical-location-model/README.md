# The Neoclassical Model of Location, Albouy and Stuart (2020)

The Neoclassical Model of Location, as set up by Albouy and Stuart ([doi](https://doi.org/10.1111/iere.12419), [working paper](https://bryan-stuart.com/files/UrbanPopulation_Feb2019.pdf)), is a general equilibrium model with 2 output sectors-- one traded good and one home (housing) good-- as well as 3 factors of production-- land, labor, and capital. The full nonlinear model is shown in Supp. Appendix B and C. Because Julia has very powerful libraries for solving all kinds of systems of equations, I decided to sharpen my skills and write a function that solves the model for the 15 endogeneous variables, subject to 15 equilibrium conditions.

-   [This](https://github.com/rcberg/macroeconomics/blob/main/scripts/neoclassical-location-model/ge-symbolic-model-function.jl) is the workhorse function that establishes the nonlinear system of equations. What's really cool about it is how it just takes utility, expenditure, and cost functions and automatically differentiates them, so no need to manually write-in the first-order condition derivatives. To use, it needs to be passed to Julia's `NonlinearProblem` function from the NonlinearSolve package, and then to `solve`. Parameters of the model must be supplied, but due to the nature of the model, it only requires the population density, and not the total population and land. This is good because the workhorse function performs poorly when using large square mile and population numbers. To get the correct estimates for an area's population and land, after running the model, multiply outputs 7-15 by the land. (I have a couple scripts that prove this works.)

-   [Here](https://github.com/rcberg/macroeconomics/blob/main/scripts/neoclassical-location-model/economic-model-example.jl) is simple example, to show the workhorse function in action.

-   [...And more in-depth example](https://github.com/rcberg/macroeconomics/blob/main/scripts/neoclassical-location-model/economic-model-data-example.jl) that calibrates the model to population density across US Core Based Statistical Areas (CBSAs), from 2023 ACS data.

## To-do list:

~~1. Julia has particularly powerful native libraries for things like [symbolic math,](https://docs.sciml.ai/Symbolics/stable/) not to mention [auto-differentiation and forward differentiation](https://juliadiff.org/ForwardDiff.jl/stable/). As a result, I probably do not need to write-out the derivatives in the workhorse function; it would be cool to see if I can just pass the cost and expenditure functions through an algorithm.~~

[Completed!](https://github.com/rcberg/macroeconomics/blob/main/scripts/neoclassical-location-model/ge-symbolic-model-function.jl) ✅

1.  Right now, the workhorse function takes many important parameters normalized to one (quality of life; traded and home productivity). I want to see if I can tweak the code to solve for variables like quality of life and home/traded sector productivity, if I am able to supply wage and house price differentials from the data.