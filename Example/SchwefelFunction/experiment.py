import subprocess
import os
import matplotlib.pyplot as plt
import numpy as np

DIMENSIONS = 6
# 2d optimal
# x_true = [420.966, 420.982]
# f_true = -837.966
# 6d optimal
x_true = [420.973, 420.996, 420.967, 420.976, 420.95, 420.97] 
f_true = -2513.9
repeat = 10000
l2_limit = 10


failed_runs = 0

bins = np.linspace(-2600, -1400, 100)
parameter_paths = [
    ["SA", "Example/SchwefelFunction/parameters.json", "Release/SA_run"],
    ["GA", "GAparameters.json", "./GA_run"],
]

for desc, parameter_path, program_path in parameter_paths:
    f_curr_list = []
    f_best_list = []
    runtime_list = []
    numSuccess = 0
    for i in range(repeat):
        print(str(i), "/", repeat, end='\r', flush=True)
        # run the SA algor
        res = subprocess.run([program_path, parameter_path], cwd=os.getcwd(), stdout=subprocess.PIPE)
        s = str(res.stdout)

        try:
            # runtime
            start = s.find("Optimisation took ") + 18
            end = s.find("ms", start)
            runtime_list.append(float(s[start:end]))

            # # final solution
            # start = s.find("current solution: x: [") + 22
            # end = s.find("]", start)
            # x = [float(xi) for xi in s[start:end].split(", ")]
            # f = float(s[s.find("f: ", start)+3:s.find("best", start)-2])
            # f_curr_list.append(f)

            # best solution found
            start = s.find("best solution: x: [") + 19
            end = s.find("]", start)
            x = [float(xi) for xi in s[start:end].split(", ")]
            f = float(s[s.find("f: ", start)+3:s.find("\n", start)-2])
            f_best_list.append(f)

            # see if it is optimal
            l2_norm = 0
            for i in range(DIMENSIONS): l2_norm += (x[i]-x_true[i])**2
            if l2_norm**0.5 < l2_limit: numSuccess += 1
        except:
            failed_runs += 1
        
    # print(x_true, f_true)
    print("prob success:", numSuccess/repeat)
    print("average runtime: ", np.mean(runtime_list), "ms")
    plt.hist(f_best_list, bins, label=desc, alpha=0.5)

print("failed number of runs: ", failed_runs)
plt.xlabel("optimisation result")
plt.ylabel(str("number of instances out of " + str(repeat) + " repeats"))
plt.legend()
plt.savefig("outcome_dist.png")
