# Statistical-Arbitrage

Run a backtest simulation using some classical alphas. The backtest is intended to be realistic due to the inclusion of transaction costs.

# Risk Model:
Every month, for the stocks in the universe where isactivenow is equal to one
on the first day of the month, compute the shrinkage estimator of the covariance matrix using the
past year of daily returns.

# Alphas: 
The alpha we will inject into the optimizer will be a weighted blend of alphas drawn from the four main categories of alphas:

1.  SHORT-TERM CONTRARIAN (TECHNICAL): the short-term mean-reversion alpha with
triangular decay. Assign to it a weight of 50%.
2.  SHORT-TERM PROCYCLICAL (FUNDAMENTAL):  Assign to it a weight of 25%.
3.  LONG-TERM CONTRARIAN (FUNDAMENTAL): the value alpha we will take the market-to-book ratio whose data is contained in the matrix. Assign
to it a weight of 15%.
4.  LONG-TERM PROCYCLICAL (TECHNICAL): the momentum alpha. We will use the last 12 months, kicking out the last month. Assign to it a weight of 10%.

Cross-sectionally demean, standardize and windsorize every day each of these four alphas individually before blending them up

# Optimizer:
 Matlab’s quadprog.m
 
