tspan = [0 5];
lags = [1 0.2];
sol = dde23(@ddefun, lags, @history, tspan);