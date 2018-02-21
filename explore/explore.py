from pypet import Environment
from pypet.utils.explore import cartesian_product
from pypet import pypetconstants
import subprocess
import numpy as np

def cross_tune(traj, w=(100, 1)):
    """ dX_upstream
        dY_upstream
        dZ_upstream
        dTheta_upstream
        dPhi_upstream
        dX_downstream
        dY_downstream
        dZ_downstream
        dTheta_downstream
        dPhi_downstream
    """
    global step
    print(traj.dx)
    exe = ["./LaserCrossTune"]
    pars = [traj.dx, traj.dy, traj.dz, traj.dThetha, traj.dPhi,
            traj.ux, traj.uy, traj.uz, traj.uThetha, traj.uPhi]

    pars = [str(item) for item in pars]

    output = subprocess.check_output(exe + pars + [str(step)])
    mean, rms = output.decode("utf-8").strip().split(" ")[-2:]

    traj.f_add_result('mean', mean)
    traj.f_add_result('rms', rms)

step = 0
filename = 'exploration.hdf5'
env = Environment(trajectory='Multiplication',
                  filename=filename,
                  overwrite_file=True,
                  file_title='Example_01_First_Steps',
                  comment='The first example!',
                  large_overview_tables=True,
                  wrap_mode=pypetconstants.WRAP_MODE_QUEUE,
                  # To see a nice overview of all
                  # computed `z` values in the resulting HDF5 file.
                  # Per default disabled for more compact HDF5 files.
                  )

# The environment has created a trajectory container for us
traj = env.trajectory

# Add both parameters
traj.f_add_parameter('dx', np.float64(1))
traj.f_add_parameter('dy', np.float64(1))
traj.f_add_parameter('dz', np.float64(1))
traj.f_add_parameter('dThetha', np.float64(0))
traj.f_add_parameter('dPhi', np.float64(0))
traj.f_add_parameter('ux', np.float64(1))
traj.f_add_parameter('uy', np.float64(1))
traj.f_add_parameter('uz', np.float64(1))
traj.f_add_parameter('uThetha', np.float64(0))
traj.f_add_parameter('uPhi', np.float64(0))

# Explore the parameters with a cartesian product
space_range = np.arange(-2,2,0.5)
dx_range = 2.5 + space_range
dy_range = 0.8 + space_range
dz_range = [0]

trajectory = {'dx': dx_range,
              'dy': dx_range,
              'dz': dx_range,
              'dThetha': [np.float64(0.006)],
              'dPhi': [np.float64(0.1)],
              'ux': [np.float64(0)],
              'uy': [np.float64(0)],
              'uz': [np.float64(0)],
              'uThetha': [np.float64(-0.001)],
              'uPhi': [np.float64(-0.003 )]
              }

traj.f_explore(cartesian_product(trajectory))

# Run the simulation
env.run(cross_tune)