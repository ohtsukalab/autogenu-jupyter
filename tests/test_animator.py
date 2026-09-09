import subprocess
import sys

import numpy as np


def test_mobile_robot_animation_accepts_modern_matplotlib(tmp_path):
    """The reference marker must pass sequences to Line2D.set_data()."""
    log_name = "mobile_robot"
    np.savetxt(tmp_path / f"{log_name}_t.log", [0.0, 0.1])
    np.savetxt(
        tmp_path / f"{log_name}_x.log",
        [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]],
    )
    script = """
import sys

import matplotlib

matplotlib.use("Agg")

from autogenu.animator import FuncAnimation, MobileRobot

FuncAnimation.save = lambda self, *args, **kwargs: None
animator = MobileRobot(
    log_dir=sys.argv[1],
    log_name="mobile_robot",
    vx_ref=1.0,
    X1=1.0,
    Y1=0.25,
    R1=0.5,
    X2=2.0,
    Y2=-0.25,
    R2=0.5,
)
animator.generate_animation()
"""

    subprocess.run([sys.executable, "-c", script, str(tmp_path)], check=True)
