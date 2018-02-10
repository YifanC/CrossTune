import subprocess
from scipy.optimize import minimize
from functools import partial
import datetime
import csv

def cross_tune(params, w=[100, 1]):
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
    global resul
    print(params)

    exe = ["./LaserCrossTune"]
    pars = [str(p) for p in params]

    output = subprocess.check_output(exe + pars + [str(step)])
    mean, rms = output.decode("utf-8").strip().split(" ")[-2:]
    resul = [mean, rms]
    print(mean, rms)
    step += 1
    return abs(w[0]*float(mean) + w[1]*float(rms))


step = 0

x0 = 10 * [0]

w = [1000, 1]
cross_tune_part = partial(cross_tune, w=w)



bounds = [(-4, 4), (-4, 4), (-4, 4), (-0.1, 0.1), (-0.1, 0.1), \
          (-4, 4), (-4, 4), (-4, 4), (-0.1, 0.1), (-0.1, 0.1)]

print(bounds)

res = minimize(cross_tune_part, x0, bounds=bounds)


with open("results.csv", "a") as result:
    now = datetime.datetime.now().strftime("%d-%m-%y-%H:%M")
    writer = csv.writer(result)
    row = w + res.x.tolist() + resul
    writer.writerow(row)