# attitudedyn 0.1.0

* Estimate The Nepaul-Stern Model of Attitude Stability and Constraint (Nepaul and Stern, 2026) <doi:10.1017/pan.2026.10034> for several Likert-scale attitudes measured on the same people across survey waves. Estimation is by simulated maximum likelihood.
* `predict()` returns the probability of each Likert category for each person-wave (`type = "prob"`) or the most likely category (`type = "class"`).
* `simulate()` draws a new panel from a fitted model. `simulate_attitudes()` draws a panel from known parameters.
* Persistence, impulse responses, influence shares, a variance decomposition, Wald tests, and score tests.
* The R package and the Python package estimate the same model.
* The vignette *attitudedyn user guide* walks through the General Social Survey application. A PDF of that guide is `docs/attitudedyn_user_guide.pdf`. The project site is <https://nepaulstern.github.io/>.
