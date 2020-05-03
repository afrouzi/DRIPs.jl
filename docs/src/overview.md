## The Problem
A LQG Dynamic Rational Inattention Problem (**Drip**) is defined as the following tracking problem, where at any point in time the agent chooses a vector of actions $\vec{a}_t\in\mathbb{R}^m$ to track a Gaussian stochastic process ``\vec{x}_t\in\mathbb{R}^n``:
```math
\begin{aligned}
    & \min_{\{\vec{a}_t\}_{t\geq 0}} \mathbb{E}\left[\sum_{t=0}^\infty \beta^t (\vec{a}_t - \vec{x}_t'\mathbf{H})'(\vec{a}_t - \vec{x}_t'\mathbf{H}) - \omega \mathbb{I}(\vec{a}^t;\vec{x}^t|\vec{a}^{t-1})\lvert \vec{a}^{-1}\right] \\
    s.t.\quad &
        \vec{x}_t=\mathbf{A}\vec{x}_{t-1}+\mathbf{Q}\vec{u}_t,\quad \vec{u}_t\sim \mathcal{N}(\mathbf{0},\mathbf{I}^{k\times k}) \\
        & \vec{a}^{-1} \text{ given.}
\end{aligned}
```

Here:
* ``\vec{a}_t\in\mathbb{R}^m`` is the vector of the agent's actions at time $t$ (a firms choosing a price, or a househld choosing consumption and labor etc.) We denote the number of actions by $m$.
* ``\vec{x}_t\in\mathbb{R}^n`` is the vector of the shocks that the agent faces at time $t$ that are exogenous to her decision, but could be endogenous to the GE model (marginal cost of a firm, real interest rates etc.) We denote the number of shocks by $n$.

### The Parameters
The LQG-DRI problem is characterized by the following parameters:
* ``\omega\in \mathbb{R}_+``: cost of 1 bit of information in units of the agent's payoff.
* ``\beta\in[0,1]``: rate of discounting information.
* ``\mathbf{A}\in \mathbb{R}^{n\times n}, \mathbf{Q}\in\mathbb{R}^{n\times k}``: Determine the state space representation of ``\vec{x}_t``.
* ``\mathbf{H}\in \mathbb{R}^{n\times m}``: interaction of payoffs with shocks. This comes from a second order approximation to the utility function and is such that under full information ``\vec{a}^*=\mathbf{H}'\vec{x}``.

### The Solution
The solution to the dynamic rational inattention problem is a joint stochastic process between the actions and the states: ``\{(\vec{a}_t,\vec{x}_t):t\geq 0\}``. Moreover, in some economic applications, we are also interested in the law of motion for the agent's belief about $\vec{x}_t$ under the optimal information structure ``\hat{x}_t=\mathbb{E}_t[\vec{x}_t]`` where the expectation is taken conditional on the agent's time $t$ information.

Theorem 2 and Proposition 3 in [Afrouzi and Yang (2019)](http://www.afrouzi.com/dynamic_inattention.pdf) characterize this joint distribution as a function of a tuple ``(\mathbf{K_t},\mathbf{Y_t},\mathbf{\Sigma}_{z,t})`` where
```math
\begin{aligned}
    \vec{a}_t &= \mathbf{H}'\hat{x}_t = \mathbf{H}'\mathbf{A}\hat{x}_{t-1}+\mathbf{Y_t}'(\vec{x}_t-\mathbf{A}\hat{x}_{t-1})+\vec{z}_t \\
    \hat{x}_t &= \mathbf{A}\hat{x}_{t-1}+\mathbf{K}_t\mathbf{Y}_t'(\vec{x}_t-\mathbf{A}\hat{x}_{t-1})+\mathbf{K}_t\vec{z}_t,\quad \vec{z}\sim\mathcal{N}(0,\Sigma_z)
\end{aligned}
```

Here,

* ``\mathbf{K}_t\in\mathbb{R}^{n\times m}`` is the Kalman-gain matrix of the agent in under optimal information acquisition at time $t$.
* ``\mathbf{Y}_t\in\mathbb{R}^{m\times m}`` is the loading of optimal signals on the state at time $t$.
* ``\mathbf{\Sigma}_{z,t}\in\mathbb{R}^{m\times m}`` is the variance-covariance matrix of the agent's rational inattention error at time $t$.

In addition to these, we might also be interested in the agent's prior and posterior subjective uncertainty, along with the continuation value that she assigns to information:
* ``\mathbf{\Sigma}_{p,t}=\mathbb{V}ar(\vec{x}_t|\vec{a}^{t})\in\mathbb{R}^{n\times n}``.
* ``\mathbf{\Sigma}_{-1,t}=\mathbb{V}ar(\vec{x}_t|\vec{a}^{t-1})\in\mathbb{R}^{n\times n}``,
* ``\bar{\Omega}_t\in\mathbb{R}^{n\times n}``.

where the matrix $\bar{\Omega}_t$ captures the value of information (see [Afrouzi and Yang (2019)](http://www.afrouzi.com/dynamic_inattention.pdf) for details)

## Steady State of DRIPs

The solver function is `solve_drip(ω,β,A,Q,H)`. It takes the primitives `(ω,β,A,Q,H)` as arguments and returns the solution of the model within a `Drip` structure that contains all the primitives and the solution objects of the model.

See the [syntax section for Drip methods](@ref DRIP_methods) for the definition of the `Drip` structure as well as for more information about `solve_drip`.

## Transition Dynamics of DRIPs

The Euler equation derived in [Afrouzi and Yang (2019)](http://www.afrouzi.com/dynamic_inattention.pdf) for the  also allows us to characterize the transition path of the information structure over time for an arbitrary initial prior.

The function `solve_trip(P::Drip,s0)` takes the steady state solution `P=solve_drip(.)` along with an initial condition `s0` as an input and returns a `Trip` structure that summarizes the transition path of the optimal information structure. The initial condition `s0` can be given either as an initial prior covariance matrix or alternatively as a one time signal about the state that perturbs the steady state prior.

See the [syntax section for Trip methods](@ref TRIP_methods) for the definition of the `Trip` structure as well as for more information about `solve_trip`.

## Impulse Response Functions

Once the model is solved, one can generate the impulse response functions of actions and beliefs using the laws of motion stated above.

We have also included a built-in function that generates these IRFs. The function `dripirfs(P::Drip)` takes a `Drip` structure as input and returns the irfs of the state, beliefs and actions to all structural shocks within a `Dripirfs` structure.

The function also returns IRFs for transition dynamics if an initial signal is specified. See the [syntax section Impulse Response Functions](@ref IRFs) for the definition of the `Dripirfs` structure as well as for more information about `dripirfs`.
