## Python examples
Here we show how to use Python interfaces.
First, run `cartpole.ipynb`, `hexacopter.ipynb`, `mobilerobot.ipynb`, or `pendubot.ipynb` in the project root directory and install Python interfaces.
Second, set PYTHONPATH according to the generated messages in the notebooks.
Then you can run examples, e.g., via
```
python3 cartpole.py
```

## Documentation 
Python API is very similar to the C++ API.  
Please refere to the C++ API documentation https://mayataka.github.io/autogenu-jupyter/ as well as the following tips of the API conversions between C++ and Python.


## Tips of API conversions between C++ and Python

### Horizon
[C++ Horizon documentation](https://mayataka.github.io/autogenu-jupyter/classcgmres_1_1_horizon.html)
|                         |  C++  |  Python  |
| ----------------------- | ----- | -------- |
| Include                 | #include "cgmres/horizon.hpp"  |  import cgmres.common |
| Constructor             | const double Tf = ...; <br> const double alpha = ...; <br> cgmres::Horizon horizon(Tf, alpha);  |  horizon = cgmres.common.Horizon(Tf=..., alpha=...)  |
| Methods                 | const double t = ...; <br> const double T = horizon.T(t);  |  T = horizon.T(t=...) |
| Print out               | std::cout << horizon << std::endl;  |  print(horizon) |

### SolverSettings
[C++ SolverSettings documentation](https://mayataka.github.io/autogenu-jupyter/structcgmres_1_1_solver_settings.html)
|                         |  C++  |  Python  |
| ----------------------- | ----- | -------- |
| Include                 | #include "cgmres/solver_settings.hpp"  |  import cgmres.common |
| Constructor             | cgmres::SolverSettings settings;  |  settings = cgmres.common.SolverSettings()  |
| Scalar member variables | settings.max_iter = ...;     |  settings.max_iter = ... |
| Print out               | std::cout << settings << std::endl;  |  print(settings) |

### OCP (OCP_cartpole as an example)
[C++ OCP_cartpole documentation](https://mayataka.github.io/autogenu-jupyter/classcgmres_1_1_o_c_p__cartpole.html)
|                         |  C++  |  Python  |
| ----------------------- | ----- | -------- |
| Include                 | #include "ocp.hpp"  |  import cgmres.cartpole  |
| Constructor             | cgmres::OCP_cartpole ocp();  |  ocp = cgmres.cartpole.OCP()  |
| Scalar member variables | ocp.g = 0.80665;     |  ocp.g = 0.80665       |
| Array member variables  | ocp.q = std::array<double, 4>({2.5, 10, 0.1, 0.1});   |  ocp.q = np.array([2.5, 10, 0.1, 0.1])  |
| Static member variables | const int nx = cgmres::OCP_cartpole::nx;   |  nx = cgmres.cartpole.OCP.nx  |
| Member functions        | const double t = ...;<br> const cgmres::VectorX x = ...; <br> const cgmres::VectorX u = ...; <br> cgmres::VectorX dx = ...; <br> ocp.eval_f(t, x, u, dx);  | t = ... <br> x = np.array([...]) <br> u = np.array([...]) <br> dx = ocp.eval_f(t, x, u) |
| Print out               | std::cout << ocp << std::endl;  |  print(ocp) |

### ZeroHorizonOCPSolver (OCP_cartpole as an example)
[C++ ZeroHorizonOCPSolver documentation](https://mayataka.github.io/autogenu-jupyter/classcgmres_1_1_zero_horizon_o_c_p_solver.html)
|                         |  C++  |  Python  |
| ----------------------- | ----- | -------- |
| Include                 | #include "cgmres/zero_horizon_ocp_solver.hpp"<br> #include "ocp.hpp"  |  import cgmres.cartpole  |
| Constructor             | cgmres::OCP_cartpole ocp(); <br> constexpr int kmax = ...; <br> cgmres::SolverSettings settings; <br> settings.max_iter = ...; <br> cgmres::ZeroHorizonOCPSolver<cgmres::OCP_cartpole, kmax> solver(ocp, settings); |  ocp = cgmres.cartpole.OCP() <br> settins = cgmres.common.SolverSettings() <br> settings.max_iter = ... <br> solver = cgmres.cartpole.ZeroHorizonOCPSolver(ocp, settings) |
| Member functions        | const double t = ...;<br> const cgmres::VectorX x = ...; <br> solver.solve(t, x);  | t = ... <br> x = np.array([...]) <br> solver.solve(t, x) |
| Setter functions        | const cgmres::VectorX u = ...;  <br> solver.set_u(u); | u = np.array([...]) <br> solver.set_u(u) |
| Getter functions        | const auto& uopt = solver.uopt(); | uopt = solver.uopt |
| Print out               | std::cout << solver << std::endl;  |  print(solver) |

### MultipleShootingCGMRESSolver (OCP_cartpole as an example)
[C++ MultipleShootingCGMRESSolver documentation](https://mayataka.github.io/autogenu-jupyter/classcgmres_1_1_multiple_shooting_c_g_m_r_e_s_solver.html)
|                         |  C++  |  Python  |
| ----------------------- | ----- | -------- |
| Include                 | #include "cgmres/multiple_shooting_cgmres_solver.hpp"<br> #include "ocp.hpp"  |  import cgmres.cartpole  |
| Constructor             | cgmres::OCP_ocpname ocp(); <br> constexpr int N = ...; <br> constexpr int kmax = ...; <br> cgmres::SolverSettings settings; <br> settings.sampling_time = ...; <br> cgmres::Horizon horizon(...); <br> cgmres::MultipleShootingCGMRESSolver<cgmres::OCP_ocpname, N, kmax> mpc(ocp, horizon, settings); |  ocp = cgmres.cartpole.OCP() <br> settings = cgmres.common.SolverSettings() <br> settings.sampling_time = ... <br> horizon = Horizon(Tf=..., alpha=...) <br> mpc = cgmres.cartpole.MultipleShootingCGMRESSolver(ocp, horizon, settings) |
| Member functions        | const double t = ...;<br> const cgmres::VectorX x = ...; <br> mpc.update(t, x);  | t = ... <br> x = np.array([...]) <br> mpc.update(t, x) |
| Setter functions        | const cgmres::VectorX u = ...;  <br> mpc.set_u(u); | u = np.array([...]) <br> mpc.set_u(u) |
| Getter functions        | const auto& uopt0 = solver.uopt()[0]; | uopt0 = solver.uopt[0] |
| Print out               | std::cout << mpc << std::endl;  |  print(mpc) |