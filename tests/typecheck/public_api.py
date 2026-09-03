"""Static-only examples that exercise the published type information."""

from pathlib import Path

import numpy as np
import numpy.typing as npt

from autogenu import RK4, AutoGenU, Logger, NLPType, forward_euler

Array = npt.NDArray[np.float64]


class LinearOCP:
    nx: int = 1

    def eval_f(self, t: float, x: Array, u: Array) -> Array:
        return -x + u


generator = AutoGenU("typed_problem", nx=1, nu=1)
x_symbols = generator.define_x()
u_symbols = generator.define_u()
generator.set_functions(
    f=[u_symbols[0]],
    C=[],
    h=[],
    L=(x_symbols[0] ** 2 + u_symbols[0] ** 2) / 2,
    phi=x_symbols[0] ** 2 / 2,
)
generator.set_nlp_type(NLPType.SingleShooting)
generator.set_horizon_params(Tf=1.0)
generator.set_solver_params(0.01, 10, 1.0e-8, 100.0, 3)
generator.set_initialization_params([0.0])
generator.set_simulation_params(0.0, [0.1], 1.0)

executable: Path = generator.get_executable_path()
installed_at: Path = generator.install_python_interface()

ocp = LinearOCP()
state = np.array([1.0], dtype=np.float64)
control = np.array([0.0], dtype=np.float64)
euler_state: Array = forward_euler(ocp, 0.0, 0.01, state, control)
rk4_state: Array = RK4(ocp, 0.0, 0.01, state, control)

logger = Logger("logs", "typed_problem")
logger.save(0.0, state, control, 0.0)
logger.close()
