# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Semi-Lagrangian inverse problem

import numpy as np
from numpy.linalg import solve
from numpy.linalg import qr
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from tqdm.notebook import trange

# +
problem = 'ts'
nx, nv = 2*163, 543
L = 20.0/3.0 * np.pi
LV = 12.0

num_steps = 300
deltat  = 0.1

xs, hx = np.linspace(0, L, nx, endpoint=False, retstep=True)
vs, hv = np.linspace(-LV, LV, nv, endpoint=False, retstep=True)

V, X = np.meshgrid(vs, xs)

f_eq = (0.9*np.exp(-V**2)/np.sqrt(np.pi)  + 0.1*np.exp(-2*(V-4.5)**2)/np.sqrt(np.pi/2))


print(hx, hv, deltat)
plt.contourf(V, X, f_eq)


# -

def plot_all(f, ee, H, hist_J):
    plt.figure()
    plt.subplot(2,2,1)
    plt.imshow(f.transpose(), extent=[xs[0], xs[-1], vs[0], vs[-1]])
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('v')

    plt.subplot(2,2,2)
    plt.plot(xs, compute_rho(f), label='rho')
    plt.plot(xs, compute_E(f), label='E')
    plt.plot(xs, H, label='H')
    plt.xlabel('x')
    plt.ylabel('E')
    plt.legend()

    plt.subplot(2,2,3)
    plt.semilogy(np.array(range(len(ee)))*deltat, ee)
    plt.xlabel('t')
    plt.ylabel('electric energy')
    
    plt.subplot(2,2,4)
    plt.semilogy(np.arange(len(hist_J)), hist_J)
    plt.xlabel("iterations")
    plt.ylabel("J")

    plt.show()


# +
def semilag_x(deltat, f):
    """
    1D Semi-Lagrangian solver
    """
    for i, v in enumerate(vs):
        f[:,i] = np.interp(xs - v*deltat, xs, f[:,i], period=L)

def semilag_v(deltat, f, E, s=None):
    if np.sum(E*deltat % hv == 0.0) != 0:
        raise ValueError("E*deltat must _not_ be a multiple of hv")
    for i in range(len(xs)):
        f[i,:] = np.interp(vs - E[i]*deltat, vs, f[i,:], period=2.0*LV)
        f[i,:] += deltat*s[i] if s is not None else 0.0


# +
def compute_rho(f):
    return 1.0 - hv*np.sum(f, axis=1)

def compute_E_from_rho(rho):
    rhohat = np.fft.fft(rho)
    Ehat = 1j*np.zeros(len(rhohat))
    # -\partial_xx \phi = 1 - \rho
    # k^2 \hat{\phi} = ...
    # i k
    Ehat[1:] = -1.0/(1j*2*np.pi/L*np.fft.fftfreq(len(rhohat), 1)[1:]*len(rhohat))*rhohat[1:]
    Ehat[0] = 0.0
    return np.real(np.fft.ifft(Ehat))

def compute_E(f):
    return compute_E_from_rho(compute_rho(f))

def electric_energy(E):
    return 0.5*np.sum(E**2)*hx


# +
def time_step_forward(dt, f, H):
    semilag_x(0.5*dt, f)
    f_star = f.copy()
    E_total = compute_E(f) + H
    semilag_v(dt, f, E_total)
    semilag_x(0.5*dt, f)
    return f_star, E_total

def run_forward(H):
    f = (1.0+0.03*np.cos(0.3*X))*f_eq

    t = 0.0
    ees = []
    fs = [f]
    E_stars = []
    f_stars = []
    for i in trange(num_steps):
        E = compute_E(f)
        ee = electric_energy(E)
        ees.append(ee)

        f_star, E_total_star = time_step_forward(deltat, f, H)
        fs.append(f)
        f_stars.append(f_star)
        E_stars.append(E_total_star)

        t += deltat
        
    return fs, ees, f_stars, E_stars


# -

def J(f):
    "energy"
    return np.sum((f[-1]-f_eq)**2)*hx*hv
    #return sum((f-f_eq)**2)*hx*hv*deltat


H = 0.0*xs
f0 = (1.0+0.03*np.cos(0.3*X))*f_eq
hist_J = []
# forward problem
fs, ee, f_stars, E_stars = run_forward(H)
hist_J.append(J(fs))
print(' goal: ', hist_J[-1])

plot_all(fs[-1], ee, H, hist_J)


# shift the second index by -n(EH)
def shift_by_n(dt, g, EH):
    out = np.zeros_like(g)
    for i in range(nx):
        n = int(np.floor(-dt*EH[i]/hv))
        out[i, :] = np.roll(g[i, :], -n)
    return out


def time_step_backward(dt, f, f_star, EH_star):
    semilag_x(-0.5*dt, f)
    g_starstar = f.copy()
    g_ss_s = shift_by_n(0*dt, f, EH_star)
    g_ss_d = (g_ss_s - np.roll(g_ss_s, 1, axis=1))/hv
    s = compute_E_from_rho( np.sum(g_ss_d*f_star, axis=1)*0*hv )
    #print(f's:',max(s))
    semilag_v(-dt, f, EH_star, s)
    semilag_x(-0.5*dt, f)
    return g_starstar


def run_backward(gT, f_stars, E_stars):
    f = gT.copy()
    t = num_steps*deltat # start from the final time
    gs = [gT]
    g_starstars = []
    for i in reversed(trange(num_steps)):
        g_starstar = time_step_backward(deltat, f, f_stars[i], E_stars[i])
        gs.append(f)
        g_starstars.append(g_starstar)
        t -= deltat
    return gs[::-1], g_starstars[::-1]


gT = 2*(fs[-1] - f_eq)/deltat
gs, g_starstars = run_backward(gT, f_stars, E_stars)

 plot_all(g_starstars[-1], ee, H, hist_J)

# gradient
f_ss_s = [shift_by_n(0*deltat, f, EH) for f, EH in zip(f_stars,E_stars)]
f_ss_d = [(np.roll(f, 1, axis=1) - np.roll(f, -1, axis=1))/2/hv for f in f_ss_s]
gradient = np.sum([np.sum(f, axis=1)*hv for f, g in zip(f_ss_s, g_starstars)], axis=0)*deltat**2
plt.plot(xs, gradient, label='gradient')
plt.legend()


