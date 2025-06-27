# Macroeconomics
This repository was initiated to store models and code for macroeconomic modeling, data, and analysis.

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
