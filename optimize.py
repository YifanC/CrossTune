import subprocess
from scipy.optimize import minimize
from functools import partial
import datetime
import csv

import matplotlib.pyplot as plt

def cross_tune(params, w=(100, 1)):
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
    global para
    global part_res

    exe = ["./../LaserCrossTune"]
    pars = [str(p) for p in params]
    print(params)
    para.append(params)

    output = subprocess.check_output(exe + pars + [str(step)])
    mean, rms = output.decode("utf-8").strip().split(" ")[-2:]
    resul = [mean, rms]
    part_res.append(resul)
    print(mean, rms)
    step += 1
    res = abs(w[0]*float(mean)+w[1]*float(rms))
    print(res)
    return res

para = []
part_res = []
step = 0

x0 = [4, 0, 0, 0, 0, 4, 0, 0, 0, 0]

w = [1, 100]
cross_tune_part = partial(cross_tune, w=w)

space_limits = (-5., 5.)
rad_limit = (-0.5, 0.5)

bounds = [space_limits, space_limits, space_limits, rad_limit, rad_limit,
          space_limits, space_limits, space_limits, rad_limit, rad_limit]

options = {'eps': 2, 'maxiter': 100}

method = 'BFGS'

res = minimize(cross_tune_part, x0, method=method, bounds=bounds, options=options)

with open("results.csv", "a") as result:
    now = datetime.datetime.now().strftime("%d-%m-%y-%H:%M")
    writer = csv.writer(result)
    row = res.x.tolist() + resul
    row.insert(0, now)
    row.insert(1, method)
    writer.writerow(row)

f, axarr = plt.subplots(3, sharex=True)

for idx, (px, part) in enumerate(zip(para, part_res)):
    axarr[0].plot(idx, px[0], "x", color="r")
    axarr[0].plot(idx, px[1], "x", color="g")
    axarr[0].plot(idx, px[2], "x", color="b")
    axarr[0].plot(idx, px[5], "o", color="r")
    axarr[0].plot(idx, px[6], "o", color="g")
    axarr[0].plot(idx, px[7], "o", color="b")
    axarr[1].plot(idx, px[3], "x", color="r")
    axarr[1].plot(idx, px[4], "x", color="g")
    axarr[1].plot(idx, px[8], "o", color="r")
    axarr[1].plot(idx, px[9], "o", color="g")
    axarr[2].plot(idx, part[0], "x", color="r")
    axarr[2].plot(idx, part[1], "o", color="g")

axarr[0].set_ylabel("cm")
axarr[1].set_ylabel("rad")
axarr[1].set_xlabel("step")

plt.savefig("res-{}-{}.png".format(w[0], w[1]))