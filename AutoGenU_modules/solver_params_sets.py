class SolverParams:
    def __init__(self, T_f, alpha, horizon_division_num, difference_increment, zeta, max_dim_krylov):
        self.T_f = T_f
        self.alpha = alpha
        self.horizon_division_num = horizon_division_num
        self.difference_increment = difference_increment
        self.zeta = zeta
        self.max_dim_krylov = max_dim_krylov


class InitializationParams:
    def __init__(self, initial_guess_solution, convergence_radius, maximum_itr_newton):
        self.initial_guess_solution = initial_guess_solution
        self.convergence_radius = convergence_radius
        self.maximum_itr_newton = maximum_itr_newton