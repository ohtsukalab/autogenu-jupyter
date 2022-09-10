import numpy as np
import os

class Logger(object):
    def __init__(self, log_dir, log_name: str):
        os.makedirs(log_dir, exist_ok=True)
        self.t_log = open(os.path.join(log_dir, log_name+"_t.log"), mode='w')
        self.x_log = open(os.path.join(log_dir, log_name+"_x.log"), mode='w')
        self.u_log = open(os.path.join(log_dir, log_name+"_u.log"), mode='w')
        self.opterr_log = open(os.path.join(log_dir, log_name+"_opterr.log"), mode='w')
        self.open = True

    def __del__(self):
        self.close()

    def save(self, t: float, x: np.ndarray, u: np.ndarray, opterr: float):
        if self.open:
            np.savetxt(self.t_log, np.array([t]))
            np.savetxt(self.x_log, [x])
            np.savetxt(self.u_log, [u])
            np.savetxt(self.opterr_log, np.array([opterr]))

    def close(self):
        self.t_log.close()
        self.x_log.close()
        self.u_log.close()
        self.opterr_log.close()
        self.open = False