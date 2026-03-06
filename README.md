# Macroeconomics
This repository was initiated to store models and code for macroeconomic modeling, data, and analysis.

## General equilibrium Julia models
I've implemented a couple static general equilibrium models in Julia which have 2 goods/firms and 3 factors of production. One of the models has a single household, another has a 2-household setup. Households optimize over consumption of the two goods, and in the 2-household case, one of the households optimizes over labor supply as well. "**L**and" is a fixed factor; labor is constrained by the fixed "hours in the day," but otherwise households and firms negotiate the final market supply within that; capital has a fixed rental price.

[Click here](https://github.com/rcberg/macroeconomics/blob/main/reports/general-eqm-1x2x3-writeup/general-equilibrium-1x2x3-writeup.ipynb) to check out the 1 household-2 good-3 factor model. In this version, the representative household owns all factors *and* supplies all the labor.

[Click here](https://github.com/rcberg/macroeconomics/blob/main/reports/general-eqm-2x2x3-writeup/general-equilibrium-2x2x3-writeup.ipynb) to see the version extended with a second household who does not participate in the labor force, but gets all their income from being the sole owner of the non-labor factors.

To see the actual functions which do the heavy-lifting, [here is the 1-household version](https://github.com/rcberg/macroeconomics/blob/main/scripts/general_equilibrium_function_1hh2good3fct.jl) and [here is the 2-household version](https://github.com/rcberg/macroeconomics/blob/main/scripts/general_equilibrium_function_2hh2good3fct.jl).

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
