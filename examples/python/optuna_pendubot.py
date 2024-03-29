import cgmres.pendubot
import cgmres.common
import numpy as np
import optuna
from collections import namedtuple


SIMULATION_TIME = 10.0
SAMPLING_TIME = 0.001

TrainParams = namedtuple('TrainParams', ['q', 'r', 'q_terminal', 'Tf'])


def run_simulation(params, verbose=False):
    # set train params to OCP
    ocp = cgmres.pendubot.OCP()
    ocp.q = params.q
    ocp.r = params.r 
    ocp.q_terminal = params.q_terminal

    # set train params to MPC 
    Tf = params.Tf
    horizon = cgmres.common.Horizon(Tf) # fixed length
    settings = cgmres.common.SolverSettings()
    settings.sampling_time = SAMPLING_TIME
    settings.zeta = 1.0 / settings.sampling_time
    mpc = cgmres.pendubot.MultipleShootingCGMRESSolver(ocp, horizon, settings)

    # initial state of the simulation
    t0 = 0.0
    x0 = np.zeros(ocp.nx)

    # Initialize solution using zero horizon OCP solution
    initializer = cgmres.pendubot.ZeroHorizonOCPSolver(ocp, settings)
    uc0 = np.array([0.01])
    initializer.set_uc(uc0)
    initializer.solve(t0, x0)
    mpc.set_uc(initializer.ucopt)
    mpc.init_x_lmd(t0, x0)
    mpc.init_dummy_mu()

    # run simulation
    t = t0
    x = x0.copy()
    sampling_time = SAMPLING_TIME
    sim_steps = int(SIMULATION_TIME/SAMPLING_TIME)
    xs = []
    us = []
    opt_error = []
    for _ in range(sim_steps):
        u = mpc.uopt[0]
        xs.append(x)
        us.append(u)
        opt_error.append(mpc.opt_error())
        dx = ocp.eval_f(t, x, u)
        x1 = x + sampling_time * dx
        mpc.update(t, x)
        x = x1
        t = t + sampling_time
        if verbose:
            print('t: ', t, ', x: ', x)
    if verbose:
        print(mpc)
    return xs, us, opt_error


def eval_simulation(xs, us, opt_error):
    NAN_INF_PENALTY = 1.0e10
    goal_state = np.array([np.pi, 0., 0., 0.])
    q = np.array([1.0, 1.0, 0.1, 0.1])
    score = 0.
    for x in xs:
        diff = x - goal_state
        if not np.isnan(diff).any() and not np.isinf(diff).any():
            score += diff.T @ np.diag(q) @ diff
        else:
            score += NAN_INF_PENALTY
    for e in opt_error:
        if not np.isnan(e) and not np.isinf(e):
            score += e 
        else:
            score += NAN_INF_PENALTY
    return score


def objective(trial):
    q = np.array([trial.suggest_loguniform('q_0', 1.0e-04, 1.0e04), 
                    trial.suggest_loguniform('q_1', 1.0e-04, 1.0e04), 
                    trial.suggest_loguniform('q_2', 1.0e-04, 1.0e04), 
                    trial.suggest_loguniform('q_3', 1.0e-04, 1.0e04)])
    r = np.array([trial.suggest_loguniform('r', 1.0e-04, 1.0e04)])
    q_terminal = np.array([trial.suggest_loguniform('q_terminal_0', 1.0e-04, 1.0e04), 
                            trial.suggest_loguniform('q_terminal_1', 1.0e-04, 1.0e04), 
                            trial.suggest_loguniform('q_terminal_2', 1.0e-04, 1.0e04), 
                            trial.suggest_loguniform('q_terminal_3', 1.0e-04, 1.0e04)])
    Tf = trial.suggest_uniform('Tf', 0.1, 5.0)
    params = TrainParams(q, r, q_terminal, Tf)
    xs, us, opt_error = run_simulation(params)
    score = eval_simulation(xs, us, opt_error)
    return score


# parameter tuning via optuna
study = optuna.create_study(study_name='pendubot_param_tune',
                            storage='sqlite:///./pendubot_param_tune.db',
                            load_if_exists=True)
study.optimize(objective, n_trials=200)


# simulation with the best params
print('The best params: ', study.best_params)
print('Simulation with the best params')
q = np.array([study.best_params['q_0'], 
              study.best_params['q_1'], 
              study.best_params['q_2'], 
              study.best_params['q_3']])
r = np.array([study.best_params['r']])
q_terminal = np.array([study.best_params['q_terminal_0'], 
                       study.best_params['q_terminal_1'], 
                       study.best_params['q_terminal_2'], 
                       study.best_params['q_terminal_3']])
Tf = study.best_params['Tf']
params = TrainParams(q, r, q_terminal, Tf)
run_simulation(params, verbose=True)