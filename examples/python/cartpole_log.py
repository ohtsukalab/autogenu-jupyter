import cgmres.cartpole
import cgmres.common
import numpy as np
from autogenu import RK4, Logger, Plotter

ocp = cgmres.cartpole.OCP()

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
initializer = cgmres.cartpole.ZeroHorizonOCPSolver(ocp, settings)
uc0 = np.array([0.01])
initializer.set_uc(uc0)
initializer.solve(t0, x0)

# Create MPC solver and set the initial solution 
mpc = cgmres.cartpole.MultipleShootingCGMRESSolver(ocp, horizon, settings)
mpc.set_uc(initializer.ucopt)
mpc.init_x_lmd(t0, x0)
mpc.init_dummy_mu()

logger = Logger(log_dir='logs/cartpole', log_name='cartpole')

# simple simulation with forward Euler 
tsim = 10.0
sampling_time = settings.sampling_time
t = t0
x = x0.copy()
for _ in range(int(tsim/sampling_time)):
    u = mpc.uopt[0]
    x1 = RK4(ocp, t, sampling_time, x, u)
    mpc.update(t, x)

    logger.save(t, x, u, mpc.opt_error())

    x = x1
    t = t + sampling_time
    print('t: ', t, ', x: ', x)

logger.close()

print("\n======================= MPC used in this simulation: =======================")
print(mpc)

plotter = Plotter(log_dir='logs/cartpole', log_name='cartpole')
plotter.show()
plotter.save()