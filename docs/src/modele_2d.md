# Two-dimensional problems

In this document, we detail some aspects of the $2$-dimensional
semi-Lagrangian method and give some examples to validate the
implementation.

## General context

We consider $2D$ transport equation of the form 

```math
\tag{1}
\partial_t f + u_x \partial_x f + u_y \partial_y f = 0, f(t=0, x, y)= f_0(x, y), x, y\in \Omega\subset \mathbb{R}^2,
```

where the advection field $(u_x, u_y)(t, x, y)$ satisfies the
incompressibility condition $\partial_x u_x + \partial_y u_y=0$ which
implies (1) can be reformulated as

```math
\partial_t f + \partial_x( u_x  f )+ \partial_y(u_y  f) = 0, f(t=0, x, y)= f_0(x, y),
```

from which we deduce the mass conservation
$\int\!\!\int f(t, x, y) dxdy = \int\!\!\int f_0(x, y) dxdy$. To solve
numerically (1), we
will use a $2D$ semi-Lagrangian method which is based on the fact that
the solution $f$ is constant along the characteristics
$X(t)=(x(t), y(t))$ defined by

```math
\dot{X}(t) = U(t, X(t)), \;\; X(s) = X_g,
```

with
$U(t, X):=U(t, x, y)=(u_x(t, x, y), u_y(t, x, y))$, $s$ is a time and
$X_g$ is a prescribed condition (which will be a grid point). Hence, we
can write $f(x, X(s))=f(t, X(t))$ for all $t, s$. Considering a
discretization of the time $t^n=n\Delta t$ with $n\in \mathbb{N}$ and
$\Delta t>0$ the time step, we rewrite the latter equality with
$s=t^{n+1}$ and $t=t^n$ to get 

```math
f(t^{n+1}, X_g) = f(t^n, X(t^n)).
```

We want to compute the left hand side which corresponds to the numerical
solution at time $t^{n+1}$ and at the grid point $X_g$. To do so, we
assume (by induction) that the solution is known at time $t^n$ and at
the grid points $X_g$, thus, we have to interpolate (in $2D$) to compute
$f(t^n, X(t^n))$. The semi-Lagrangian methods can be split into two
steps

1.  compute $X(t^n)$

2.  compute $f(t^n, X(t^n))$: from the known values $f(t^n, X_g)$,
    interpolate at $X(t^n)$.

These two steps are details in the next section.

## Details on the $2D$ semi-Lagrangian method

The two steps of the semi-Lagrangian are detailed and some examples are
given.

## ODE solver

First, we need to compute $X(t^n)$ which the solution at time $t^n$ of

```math
\tag{2}
\dot{X}(t) = U(t, X(t)), \;\; X(t^{n+1}) = X_g.
```
When $U$ is simple
enough, $X(t^n)$ can be computed analytically but in general, we need a
solver of this differential equation. The main difficulty comes from the
fact that (2) has to
be solved backward in time and when the time dependency of $U$ is
nonlinearly coupled to the solution $f$ itself (see Examples 3 and 4
below), we do not know $U(t, \dot)$ for $t>t^n$ (and time extrapolation
has to be used [^filbet]).

### First order 

A simple scheme to compute $X(t^n)$ is the Euler scheme applied to (2)

```math
\frac{X(t^{n+1}) - X(t^n)}{\Delta t} = U(t^n, X(t^{n+1}),
```
and using
the condition $X(t^{n+1}) = X_g$, we get the following first order
approximation for $X(t^n)$ 

```math
X(t^n) = X_g - \Delta t U(t^n, X_g).
```

### Second order 

A mid-point scheme (which is second order accurate) can be used to solve (2):

```math
\frac{X(t^{n+1}) - X(t^n)}{\Delta t} = U\Big(t^{n+1/2}, \frac{X(t^{n+1}) + X(t^n)}{2}\Big),
```

which gives an implicit expression for $X(t^n)$ 

```math
\begin{equation}
\label{ode_2nd_imp}
X(t^n) = X(t^{n+1}) -\Delta t \, U\Big(t^{n+1/2}, \frac{X(t^{n+1}) + X(t^n)}{2}\Big).
\end{equation}
```

As mentioned above, we first need to extrapolate $U(t^{n+1/2}, \cdot)$.
To do so, we use $U(t^{n-1}, \cdot)$ and $U(t^{n}, \cdot)$ and we
construct a first order Lagrange polynomial
${\cal L}(t), t\in[t^{n-1}, t^n]$

```math
{\cal L}(t) = U(t^n, \cdot) \frac{t^{n-1}-t}{\Delta t} + U(t^{n-1}, \cdot) \frac{t - t^{n}}{\Delta t}.
```

We then approximate $U(t^{n+1/2}, \cdot)$ by ${\cal L}(t^{n+1/2})$. We then have to solve
\eqref{ode_2nd_imp} using a fixed point

```math
X^{n, k+1} = X_g -\Delta t \, {\cal L}\Big(t^{n+1/2}, \frac{X_g + X^{n,k}}{2}\Big), \;\; \mbox{ for } k\geq 0, X^{n, 0}=X_g,
```

up to convergence.

### Extension to higher order 

We can look at the schemes proposed in [^filbet] but we can also use the package `DifferentialEquations.jl`

## $2D$ Interpolation

Once $X(t^n)$ has been computed using the techniques detailed before,
one has to compute $f(t^n, X(t^n))$. Since the grid values $f(t^n, X_g)$
are known, we can reconstruct a piecewise $2D$ polynomial function
${\cal P}(x, y), x, y\in \Omega$ using Lagrange or splines such that
$f(t^n, X_g) = {\cal P}(X_g)$, and then we approximate
$$f(t^n, X(t^n)) \approx {\cal P}(X(t^n)).$$ We consider a cartesian
grid of $N_x\times N_y$ points $x_i = (i-1)\Delta x$ and
$y_j =(j-1)\Delta y$ for $i=1, \dots, N_x$ and $j=1, \dots, N_y$. We
assume periodicity in $x$ and $y$ which means
$f(x_{N_x+1}, y) = f(x_1, y)$ and $f(x, y_{N_y+1}) = f(x, y_1)$.

### Lagrange interpolation 

A tensor product can be done. If we denote $L_{x,i}$ and $L_{y,j}$ the
Lagrange polynomial of degree $2d+1$ in the $x$ and $y$ direction
$L_{x,i}(x) = \Pi_{k=i-d, k\neq i}^{i+d} \frac{(x-x_k)}{(x_i-x_k)}$ and
$L_{y,j}(x) = \Pi_{k=j-d, k\neq j}^{j+d} \frac{(y-y_k)}{(y_j-y_k)}$.
Then, we have

```math
f(x,y) \approx {\cal P}(x,y) =\sum_{i=1}^{N_x} \sum_{j=1}^{N_y} f(x_i,y_j) L_{x,i}(x)L_{y,j}(y).
```

### Splines interpolation 

For cubic splines, one can use

```julia
splinePP 
```

which I try to explain below.

We introduce the tridiagonal matrix $A\in {\cal M}_{Nx+1, N_x+1}$ with
$4$ on the diagonal and $1$ on the two extradiagonals. We consider
$f(x_i, y_j)$ for $i=1, \dots, N_x$ and $j=1, \dots, N_y$ and we denote
$F\in {\cal M}_{N_x+1, N_y}$ the matrix such that $F_{i,j} = f(x_i,y_j)$
for $i=1, \dots, N_x$ and $j=1, \dots, N_y$ and
$F_{N_x+1, j}= f(x_1, y_j)$ for periodicity. We compute the splines
coefficients as

- solve $A \eta_j = 6 f(:,y_j), \;\; \forall j=1, \dots, N_y$ with
  $\eta_j\in \mathbb{R}^{N_x+1}$

- gestion du bord (Pierre ?) pour obtenir $\eta(1:Nx+3, 1:N_y+1)$

- solve $A \,$coef$_i \, = 6 \eta(i,:), \;\; \forall i=1, \cdots, N_x+3 $ with coef$_i\in \mathbb{R}^{N_y+1}$

- gestion du bord (Pierre ?) pour obtenir coef$(1:Nx+3, 1:N_y+3)$

```math
f(x,y) \approx {\cal S}(x,y) = \sum_i\sum_j {\tt coef}_{i,j} B_i(x) B_j(y)
```

## Some examples

In this part, some examples of increasing difficulty are proposed to
validate our algorithms.

### Example 1 : rotation

```math
\partial_t f + y \partial_x f - x \partial_y f = 0, f(t=0, x, y)= f_0(x, y).
```

We already discussed this test and we have to check it again. In this
case, the characteristics can be solved exactly in time since
$X(t^n) = e^{-J \Delta t}X_g$ with $J$ the symplectic matrix. This test
will enable us to validate the $2D$ interpolation step. Here

```math
J=\Big(
\begin{matrix}
0 & 1 \\
-1 & 0
\end{matrix}
\Big) 
\mbox{ and } 
e^{-J \Delta t}=\Big(
\begin{matrix}
\cos \Delta t & -\sin\Delta t \\
\sin\Delta t & \cos\Delta t
\end{matrix}
\Big).
```

### Example 2 : swirling deformation flow

```math
\begin{aligned}
\partial_t f + \Big( \sin^2(\pi x) \sin(2\pi y) g(t)\Big) \partial_x f - \Big(\sin^2(\pi y)\sin(2\pi x) g(t)  \Big) \partial_y f = 0, \\
f(t=0, x, y)= f_0(x, y)
\end{aligned}
```

with $g(t)=\cos(\pi t/T)$. The solution slows down and reverses
direction in such a way that the initial condition should be recovered
at time $T$: $f(T, x, y)=f_0(x,y)$. This gives a very useful test to
validate our methods since we know the exact solution at time $T$. For
the initial condition, we consider 

```math
\begin{equation}
f_0(x, y) = 
\left\{ 
\begin{array}{cc}
1, & \mbox{ if } (x-1)^2+(y-1)^2 <0.8\\
0, & \mbox{ otherwise}.
\end{array}
\right.
\end{equation}
```

We can choose $T=1.5$ and the spatial domain is $[0, 1]^2$.
Some results are given in [^leveque] or [^qiu].

### Example 3: Vlasov-Poisson

```math
\partial_t f + y \partial_x f +E \partial_y f = 0, f(t=0, x, y)= f_0(x, y), x\in [0, 4\pi], y\in \mathbb{R},
```

where the electric field $E$ derives from a potential
$\phi(t, x)\in\mathbb{R}$ which satisfies a Poisson equation

```math
\partial_x^2 \phi = \int_{\mathbb{R}} f dy - 1.
```

The initial condition is

```math
f_0(x, y)= \frac{1}{\sqrt{2}\pi}e^{-y^2/2}(1 + 0.001\cos(x/2)).
```

For this problem, we define the electric energy ``{\cal E}_e = \int E^2dx ``
and the kinetic energy ``{\cal E}_k := \int\!\int y^2 f dx dy`` so that the
total energy ``{\cal E}_e+{\cal E}_k`` is preserved with time. We can also
consider the time evolution of ``{\cal E}_e`` for which we know the
behavior.

### Example 4: guiding-center

```math
\partial_t f + E_x \partial_x f +E_y \partial_y f = 0, f(t=0, x, y)= f_0(x, y)
```

where the electric field $E=(E_x, E_y)(t, x, y)$ derives from a
potential $\phi(t, x, y)\in\mathbb{R}$ which satisfies a Poisson
equation 

```math
\Delta \phi = f.
```

The spatial domain is $[0, 4\pi]\times [0, 2\pi]$ and the initial condition is

```math
f_0(x, y)= \sin(y) + 0.015\cos(x/2).
```

For this problem, the electric
energy ``{\cal E}_e :=\int\!\int (E_x^2+E_y^2)dxdy`` and the enstrophy
``{\cal E}_f:=\int\!\int f^2 dx dy`` are preserved with time. Some results
are given in [^qiu] or [^crouseilles].

## Etienne's models

The code for the Etienne's model is available at
<https://github.com/vressegu/sqgmu> (matlab). If I understood well, the
deterministic case is exactly the guiding-center model (with different
initial condition but with periodic boundary conditions). Stochastic
terms can be added that we perhaps can take into account in our
numerical method.

---

## References

[^filbet]: F. Filbet, C. Prouveur, *High order time discretization for backward semi-Lagrangian methods*, Journal of Computational and Applied Mathematics, vol 303, (2016) pp. 171-188.

[^leveque]: R. LeVeque, *High-resolution conservative algorithms for advection in incompressible flow*, SIAM Journal on Numerical Analysis, (1996), pp.  627-665.  <https://www.jstor.org/stable/2158391?seq=29#metadata_info_tab_contents>

[^qiu]: J. Qiu, C.-W. Shu, *Conservative high order semi-Lagrangian finite difference WENO methods for advection in incompressible flow*, Journal of Computational Physics, Volume 230, Issue 4, 20 (2011), pp. 863-889.

[^lauritzen]: P. Lauritzen, D. Ramachandran, P. Ullrich, *A conservative semi-Lagrangian multi-tracer transport scheme (CSLAM) on the cubed-sphere grid*, J. Comput. Phys. 229, (2010), pp. 1401-1424.

[^crouseilles]: N. Crouseilles, M. Mehrenberger, E. Sonnendrücker *Conservative semi-Lagrangian methods for the Vlasov equations*, J. Comput. Phys., 229, pp 1927-1953, (2010).

