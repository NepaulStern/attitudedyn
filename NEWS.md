# attitudedyn 0.1.0

* Fit a dynamic ordered-probit model for several ordinal panel attitudes by simulated maximum likelihood.
* `predict()` returns Likert category probabilities for each person-wave (`type = "prob"`) or the most likely category (`type = "class"`).
* `simulate()` draws a new panel from a fitted model. `simulate_attitudes()` draws a panel from known parameters.
* Persistence, impulse responses, influence shares, a variance decomposition, Wald tests, and score tests.
* The numerical work is a dependency-free C++ core shared with the Python wrapper.
* The vignette *attitudedyn user guide* walks through the GSS application. A PDF of that guide is `docs/attitudedyn_user_guide.pdf`.
