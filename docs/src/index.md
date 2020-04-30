# A Toolbox for Solving Dynamic Rational Inattention Problems

This notebook explains how to solve LQG Dynamic Rational Inattention models using the methods and the [solver](https://github.com/choongryulyang/dynamic_multivariate_RI) from [Afrouzi and Yang (2019)](http://www.afrouzi.com/dynamic_inattention.pdf) that was adapted for Julia by the authors and [Miguel Acosta](https://www.acostamiguel.com/home). 

To add the package, execute:

```
] DRIPs
```

To include the package, execute:

```julia
using DRIP;
```

## A Structure for LQG Dynamic Rational Inattention Problems (D.R.I.P.)

### The Problem
A LQG Dynamic Rational Inattention Problem (**drip**) is defined as the following tracking problem, where at any point in time the agent chooses a vector of actions $\vec{a}_t\in\mathbb{R}^m$ to track a Gaussian stochastic process $\vec{x}_t\in\mathbb{R}^n$:
\begin{align}
    & \min_{\{\vec{a}_t\}_{t\geq 0}} \mathbb{E}\left[\sum_{t=0}^\infty \beta^t (\vec{a}_t - \vec{x}_t'\mathbf{H})'(\vec{a}_t - \vec{x}_t'\mathbf{H}) - \omega \mathbb{I}(\vec{a}^t;\vec{x}^t|\vec{a}^{t-1})\lvert \vec{a}^{-1}\right] \\
    s.t.\quad &
        \vec{x}_t=\mathbf{A}\vec{x}_{t-1}+\mathbf{Q}\vec{u}_t,\quad \vec{u}_t\sim \mathcal{N}(\mathbf{0},\mathbf{I}^{k\times k})
\end{align}
Here:
* $\vec{a}_t\in\mathbb{R}^m$ is the vector of the agent's actions at time $t$ (a firms choosing a price, or a househld choosing consumption and labor etc.) We denote the number of actions by $m$.
* $\vec{x}_t\in\mathbb{R}^n$ is the vector of the shocks that the agent faces at time $t$ that are exogenous to her decision, but could be endogenous to the GE model (marginal cost of a firm, real interest rates etc.) We denote the number of shocks by $n$.

#### The Parameters
The LQG-DRI problem is characterized by the following parameters:
* $\omega\in \mathbb{R}_+$: cost of 1 bit of information in units of the agent's payoff. 
* $\beta\in[0,1]$: rate of discounting information.
* $\mathbf{A}\in \mathbb{R}^{n\times n}, \mathbf{Q}\in\mathbb{R}^{n\times k}$: Determine the state space representation of $\vec{x}_t$.
* $\mathbf{H}\in \mathbb{R}^{n\times m}$: interaction of payoffs with shocks. This comes from a second order approximation to the utility function and is such that under full information $\vec{a}^*=\mathbf{H}'\vec{x}$.

#### The Solution
The solution to the dynamic rational inattention problem is a joint stochastic process between the actions and the states: $\{(\vec{a}_t,\vec{x}_t):t\geq 0\}$. Moreover, in some economic applications, we are also interested in the law of motion for the agent's belief about $\vec{x}_t$ under the optimal information structure $\hat{x}_t=\mathbb{E}_t[\vec{x}_t]$ where the expectation is taken conditional on the agent's time $t$ information.

Theorem 2 and Proposition 3 in [Afrouzi and Yang (2019)](http://www.afrouzi.com/dynamic_inattention.pdf) characterize this joint distribution as a function of a tuple $(\mathbf{K},\mathbf{Y},\mathbf{\Sigma}_z)$ where
\begin{align}
    \vec{a}_t &= \mathbf{H}'\hat{x}_t = \mathbf{H}'\mathbf{A}\hat{x}_{t-1}+\mathbf{Y}'(\vec{x}_t-\mathbf{A}\hat{x}_{t-1})+\vec{z}_t \\
    \hat{x}_t &= \mathbf{A}\hat{x}_{t-1}+\mathbf{K}\mathbf{Y}'(\vec{x}_t-\mathbf{A}\hat{x}_{t-1})+\mathbf{K}\vec{z}_t,\quad \vec{z}\sim\mathcal{N}(0,\Sigma_z)
\end{align}

Here, 

* $\mathbf{K}\in\mathbb{R}^{n\times m}$ is the Kalman-gain matrix of the agent in under optimal information acquisition in the stationary solution.
* $\mathbf{Y}\in\mathbb{R}^{m\times m}$ maps the state to the agent's action under rational inattetion and govern the covariance of the two.
* $\mathbf{\Sigma}_z\in\mathbb{R}^{m\times m}$ is the variance-covariance matrix of the agent's rational inattention error.

In addition to these, we might also be interested in the agent's prior and posterior subjective undertainty, along with the continution value that she assigns to information:
* $\mathbf{\Sigma}_{p}=\mathbb{V}ar(\vec{x}_t|\vec{a}^{t})\in\mathbb{R}^{n\times n}$.
* $\mathbf{\Sigma}_{-1}=\mathbb{V}ar(\vec{x}_t|\vec{a}^{t-1})\in\mathbb{R}^{n\times n}$, 
* $\bar{\Omega}\in\mathbb{R}^{n\times n}$.

The matrix $\bar{\Omega}$, which is the continuation value of information, is also important for the number of iterations that the code needs to converge. When $n$ is large, accurate guesses for $\bar{\Omega}$ and $\mathbf{\Sigma}_{-1}$ speed up the convergence considerably.

### The Solver for the Steady State of D.R.I.P.

The solver function is `solve_drip(ω,β,A,Q,H;Σ0 = A*A'+Q*Q',Ω0 = H*H')`. It takes the primitives `(ω,β,A,Q,H;Σ0,Ω0)` as arguments, where the arguments `Σ0,Ω0` are the initial guesses for $\mathbf{\Sigma}_{-1}$ and $\bar{{\Omega}}$ and are optional for cases in which an accurate initial guess is at hand. If not specified, the default values are set as `Ω0=H*H'` and `Σ0=Q*Q'`. See the code for extra options regarding tolerance of convergence, maximum number of iterations, and updating weight for every iteration.

The function returns the solution of the model within a `drip` structure,  that contains all the primitives and the solution objects of the model:


```julia
mutable struct Drip 
    ω; β; A; Q; H;                    # primitives
    K; Y; Σ_z; Σ_p; Σ_1; Ω;           # solution
    Drip()          = new();
    Drip(ω,β,A,Q,H) = new(ω,β,A,Q,H);
end
```

### Impulse Resonse Functions using Steady State Information Structure

Once the model is solved, one can generate the impulse response functions of actions and beliefs using the laws of motion stated above. 

We have also included a built-in function that generates these IRFs. The function `dripirfs(P::drip,T::Int)` takes a `drip` structure as input along with the length of irfs `T` and returns the irfs of the state, beliefs and actions to all possible shocks in the following structure:


```julia
struct Dripirfs
    T     :: Int            # length of IRFs
    x     :: Array{Float64} # IRFs of the fundamental shocks
    x_hat :: Array{Float64} # IRFs of beliefs
    a     :: Array{Float64} # IRFs of actions
end
```

Here 
* `x` is of dimension $n\times k \times T$ where `x[i,j,:]` is the IRF of the $i$'th element of the state vector ($\vec{x}_t$) to the $j$'th shock in the vector $\vec{u}_t$.

* `x_hat` is of dimension $n\times k \times T$ where `x_hat[i,j,:]` is the IRF of the agent's belief about $i$'th element of the state vector  ($\hat{x}_t$) to the $j$'th shock in the vector $\vec{u}_t$.

* `x_hat` is of dimension $m\times k \times T$ where `a[i,j,:]` is the IRF of the $i$'th element of the action vector ($\vec{a}_t$) to the $j$'th shock in the vector $\vec{u}_t$.

### The Solver for the Transition dynamics of Rational Inattention Problems (T.R.I.P.)

The Euler equation derived in [Afrouzi and Yang (2019)](http://www.afrouzi.com/dynamic_inattention.pdf) for the  also allows us to characterize the transition path of the information structure over time for an arbitrary initial prior. The function `solve_trip(P::drip; T=100)` takes the solution `P=solve_drip(.)` as an input and returns a `trip` structure that summarizes the transition path of information benefit matrics, priors and a vector for the marginal values of information over time (the input `T` determines the number of periods that it takes to reach to the steady state and is set to 100 by default):


```julia
struct Trip 
    P::Drip;                # problem and solution in steady state
    T::Int;                 # length of T.R.I.P.
    Σ_1s; Σ_ps; Ωs;         # priors, posteriors and benefit matrices 
    Ds;                     # eigenvalues of Σ_t^(0.5)Ω_tΣ_t^(0.5)
end 
```
