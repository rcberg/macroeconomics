# Macroeconomics
This repository was initiated to store models and code for macroeconomic modeling, data, and analysis.

## General equilibrium Julia model
Using Julia, I've implemented a static general equilibrium model with 1 household agent, 2 goods/firms, and 3 factors of production. Households optimize over consumption of the two goods, as well as labor supply. One of the factors ("**L**and") has a fixed supply; labor has a fixed cap with households and firms "negotiating" the final market supply; capital has a fixed rental price. I would like do a more detailed write-up of this model soon.

[Click here](https://github.com/rcberg/macroeconomics/blob/master/scripts/general_equilibrium_function_1hh2good3fct.jl) for the primary workhorse function, which takes as inputs 1 vector of model parameters and 1 vector with the initial "guess." The output is 1 vector of first-order condition residuals (should be very close to 0); 1 solution vector; and a 0-1 indicating if the solver converged.

[Click here](https://github.com/rcberg/macroeconomics/blob/master/scripts/general_equilibrium_function_1hh2good3fct_stress_test.jl) to check out a demonstration of the function. It has a for-loop which scales-up the fixed land $L$ and maximum labor supply $N_T$ from 1 to 100, and another for-loop which holds either $L$ or $N_T$ at 100, and scales-back the other down to 1. Note: Right now, as the loop scale approaches $\approx \frac{1}{10}^{th}$ of the scale max (100), the solutions start becoming numerically unstable. 

## Mosquito problem solution
A worked-through answer to an old macroeconomics preliminary exam question: [Click here](https://raw.githack.com/rcberg/macroeconomics/master/reports/mosquito-problem/mosquito_macro_problem.html)

## Overlapping Generations (OLG) economy simulation
Below is R code which produces simulations of a basic 2-period Overlapping Generations macroeconomic model with simple logarithmic utility (CES with $\sigma$ = 1). I am working on extending the code to a more general CES utility function.

The model incorporates a lump-sum wage tax during the agent's working period, which is rebated to the agent as a public good. The tax is set to the level which causes the agent's savings to become consistent with dynamically-efficient capital accumulation.

Code: [Click here](https://github.com/rcberg/macroeconomics/blob/master/scripts/efficient-olg-simulation-with-taxes.R)

## Ramsey-Cass-Koopmans (RCK) economy simulation
Below is R code which produces simulations of a basic Ramsey-Cass-Koopmans macroeconomic model. The model features a growing representative dynasty that lives to infinity, and plans as a single decision-maker. The code features a shooting algorithm that numerically solves for the dynamic saddle path to the steady state. (This algorithm is fairly sensitive to model parameters.)

Code: [Click here](https://github.com/rcberg/macroeconomics/blob/master/scripts/infinite-ramsey-model-sim.R)

For a walk-through of the Ramsey-Cass-Koopmans model that's oriented towards people with a decent mathematical background, but who might be less-familiar with macroeconomic models, make sure to check out my write-up below.

Write-up: [(Light mode](https://raw.githack.com/rcberg/macroeconomics/refs/heads/main/reports/rck-sim-writeup/rck_model_simulation.html)[| Dark mode)](https://raw.githack.com/rcberg/macroeconomics/refs/heads/main/reports/rck-sim-writeup/rck_model_simulation_dark.html)
