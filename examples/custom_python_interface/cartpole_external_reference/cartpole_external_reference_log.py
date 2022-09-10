import cgmres.cartpole_external_reference
import cgmres.common
import numpy as np


ocp = cgmres.cartpole_external_reference.OCP()

# this is actually std::shared_ptr<cgmres::OCP_cartpoleExternalReference::ExternalReference>
external_reference = cgmres.cartpole_external_reference.ExternalReference() 
external_reference.cart_position = 0.0
# set the external reference to ocp
ocp.external_reference = external_reference

horizon = cgmres.common.Horizon(Tf=2.0) # fixed length

settings = cgmres.common.SolverSettings()
settings.sampling_time = 0.001
settings.zeta = 1000
settings.finite_difference_epsilon = 1e-08
settings.max_iter = 50
settings.opterr_tol = 1e-06
settings.verbose_level = 1 # print opt error

t0 = 0.0
x0 = np.zeros(ocp.nx)

# Initialize solution using zero horizon OCP solution
initializer = cgmres.cartpole_external_reference.ZeroHorizonOCPSolver(ocp, settings)
uc0 = np.array([0.01])
initializer.set_uc(uc0)
initializer.solve(t0, x0)

# Create MPC solver and set the initial solution 
mpc = cgmres.cartpole_external_reference.MultipleShootingCGMRESSolver(ocp, horizon, settings)
mpc.set_uc(initializer.ucopt)
mpc.init_x_lmd(t0, x0)
mpc.init_dummy_mu()

from autogenu import Logger
logger = Logger(log_dir='logs/cartpole_external_reference', log_name='cartpole_external_reference')

# simple simulation with forward Euler 
tsim = 10.0
sampling_time = settings.sampling_time
t = t0
x = x0.copy()
for _ in range(int(tsim/sampling_time)):
    u = mpc.uopt[0]
    dx = ocp.eval_f(t, x, u)
    x1 = x + sampling_time * dx
    mpc.update(t, x)

    logger.save(t, x, u, mpc.opt_error())

    x = x1
    t = t + sampling_time
    print('t: ', t, ', x: ', x)

    # switch the reference cart position from 0.0 to 1.0
    if t > 5.0: 
      external_reference.cart_position = 1.0

logger.close()

print("\n======================= MPC used in this simulation: =======================")
print(mpc)

from autogenu import Plotter
plotter = Plotter(log_dir='logs/cartpole_external_reference', log_name='cartpole_external_reference')
plotter.show()
plotter.save()