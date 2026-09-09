from typing import Protocol

import numpy as np
import numpy.typing as npt

Array = npt.NDArray[np.float64]


class OCPProtocol(Protocol):
    nx: int

    def eval_f(self, t: float, x: Array, u: Array) -> Array: ...


def forward_euler(
    ocp: OCPProtocol, t: float, dt: float, x: Array, u: Array
) -> Array:
    dx = np.zeros(ocp.nx)
    dx = ocp.eval_f(t, x, u)
    x1 = x + dt * dx
    return x1

def RK4(ocp: OCPProtocol, t: float, dt: float, x: Array, u: Array) -> Array:
    k1 = np.zeros(ocp.nx)
    k2 = np.zeros(ocp.nx)
    k3 = np.zeros(ocp.nx)
    k4 = np.zeros(ocp.nx)
    x1 = np.zeros(ocp.nx)
    k1 = ocp.eval_f(t, x, u)
    x1 = x + 0.5 * dt * k1
    k2 = ocp.eval_f(t+0.5*dt, x1, u)
    x1 = x + dt * 0.5 * (np.sqrt(2.0)-1.0) * k1 + dt*(1.0-(1.0/np.sqrt(2.0))) * k2
    k3 = ocp.eval_f(t+0.5*dt, x1, u)
    x1 = x - dt * 0.5 * np.sqrt(2.0) * k2 + dt * (1.0+(1.0/np.sqrt(2.0))) * k3
    k4 = ocp.eval_f(t+dt, x1, u)
    x1 = x + (dt/6.0) * (k1+(2.0-np.sqrt(2.0))*k2 + (2.0+np.sqrt(2.0))*k3+k4)
    return x1
