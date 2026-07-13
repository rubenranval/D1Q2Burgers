# D1Q2Burgers

A Julia implementation of the D1Q2 lattice Boltzmann scheme for the inviscid Burgers equation, together with an exact reference solution and a convergence study in the relaxation parameter s.

$$\partial_t u + \partial_x\ \left(\frac{u^2}{2}\right) = 0, \qquad x \in [0,1],\ t \in [0,T].$$

The scheme is the one analysed by Graille as a lattice-Boltzmann reading of the Jin–Xin relaxation method [1]. The point of the repository is to show numerically the behaviour predicted by the equivalent-equation analysis: the scheme is first-order accurate for $s \neq 2$ and second-order accurate for $s = 2$.

## The scheme

Two populations $f_1, f_2$ live on a uniform lattice of step $\Delta x$, with velocities $\pm\lambda$ and time step $\Delta t = \Delta x/\lambda$ (so that free flight moves a population exactly one cell).

The conserved moment is the density, the second one carries the flux:

$$m_0 = f_1 + f_2 = u, \qquad m_1 = \lambda\,(f_2 - f_1), \qquad m_1^{\mathrm{eq}} = \frac{u^2}{2}.$$

Inverting the moment matrix gives the equilibria used in the code,

$$f_i^{\mathrm{eq}}(u) = \frac{u}{2} + c_i\,\frac{u^2}{4\lambda}, \qquad c_i = \pm 1 .$$

At eachh time step we have: a collision, then transport:

$$f_i \leftarrow f_i + s\left(f_i^{\mathrm{eq}}(u) - f_i\right), \qquad f_i(x, t+\Delta t) \leftarrow f_i(x - c_i\lambda\,\Delta t,\ t).$$

Since $\sum_i f_i^{\mathrm{eq}}(u) = u$, the density is conserved exactly by the collision; only $m_1$ relaxes, at rate $s = 1/\tau$.

Stability requires $0 < s \le 2$ together with the subcharacteristic condition $\lambda \ge \max|u|$. In the test case below $\max|u| = 1$, so the default $\lambda = 1$ is the marginal choice, larger $\lambda$ can be passed to d1q2.

**Order of accuracy.** The equivalent equation of the scheme is

$$\partial_t u + \partial_x\ \left(\frac{u^2}{2}\right) = \Delta x \left(\frac{1}{s} - \frac{1}{2}\right)\partial_x\ \Big[(\lambda^2 - u^2)\,\partial_x u\Big] + \mathcal{O}(\Delta x^2),$$

so the leading numerical viscosity is proportional to $1/s - 1/2$ and **vanishes exactly at $s = 2$**, hence order 2 there, order 1 elsewhere.

## Usage

```julia
using D1Q2Burgers

res  = d1q2(; k = 10, s = 2.0, T = 0.2)   # N = 2^k cells
u_ex = burgers_exact.(res.x, res.t)       # NB: res.t, not res.T 
```

`d1q2` returns a `NamedTuple` `(; x, u, k, s, λ, dx, dt, nt, T, t)`. The number of time steps is `nt = round(Int, T/dt)`, so the time actually reached, `t = nt*dt`, is in general slightly different from the requested `T`. **Errors must be measured against the exact solution at `t`**, otherwise the mismatch acts as an $\mathcal{O}(\Delta t)$ perturbation and pollutes the second-order rate at `s = 2`.

Convergence study:

```julia
res = convergence(; s_list = [2.0, 1.9, 1.75, 1.0, 0.75, 0.5], ks = 3:16, T = 0.2)
res.errors[2.0]   # relative L^2 errors
res.orders[2.0]   # estimated orders between successive meshes
```

## Results

Relative discrete $L^2$ error at $T = 0.2$, $\lambda = 1$:

| $\Delta x$ | $s = 2$ | order | $s = 1.75$ | order | $s = 1$ | order |
|---|---|---|---|---|---|---|
| $2^{-8}$  | 4.25e-04 | 1.92 | 1.17e-03 | 1.08 | 7.62e-03 | 0.98 |
| $2^{-9}$  | 1.08e-04 | 1.98 | 5.72e-04 | 1.03 | 3.85e-03 | 0.98 |
| $2^{-10}$ | 2.75e-05 | 1.97 | 2.84e-04 | 1.01 | 1.95e-03 | 0.99 |
| $2^{-11}$ | 6.87e-06 | 2.00 | 1.41e-04 | 1.01 | 9.77e-04 | 0.99 |
| $2^{-12}$ | 1.73e-06 | 1.99 | 7.03e-05 | 1.00 | 4.89e-04 | 1.00 |

Values of $s$ close to 2 (e.g. $s = 1.9$) show a preasymptotic second-order regime on coarse meshes, before the $\mathcal{O}(\Delta x)$ term, small but not zero, takes over. The results are similar to the ones obtained by B. Graille and al.


## References

[1] B. Graille, *Approximation of mono-dimensional hyperbolic systems: a lattice Boltzmann scheme as a relaxation method*, Journal of Computational Physics **266** (2014), 74–88. [doi:10.1016/j.jcp.2014.02.017](https://doi.org/10.1016/j.jcp.2014.02.017)


## License
