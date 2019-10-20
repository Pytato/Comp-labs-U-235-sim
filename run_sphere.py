import uranium_sim
import pylab as pl
import seaborn
import multiprocessing
import time
import numpy
import tqdm


def magnitude(x):
    return numpy.floor(numpy.log10(x))


def step_to_true(initial, old_initial, ctx_data, prev_out=0):
    ctx_data["run_count"] += 1

    with uranium_sim.Uranium(initial, shape="sphere") as uranium_c:
        crit_prop_a_l = 0
        sec_dec_a_l = 0

        for i in range(int(ctx_data["average_length"])):
            sec_prop = uranium_c.fission(ctx_data["neutron_i"])[1]
            crit_prop_a_l += sec_prop

        sec_prop_a = crit_prop_a_l / ctx_data["average_length"]

        ctx_data["critical_mass_prop"].append(sec_prop_a)
        ctx_data["lengths"].append(initial)

    new_out = sec_prop_a

    if new_out > ctx_data["ideal_out"] > prev_out or new_out < ctx_data["ideal_out"] < prev_out:
        ctx_data["step"] *= -0.25
    else:
        ctx_data["step"] *= 1.25

    if round(new_out, ctx_data["num_dp"]) == 1.0 or ctx_data["step"] == 0:
        return initial, ctx_data

    else:
        # print(new_out, initial)

        if ctx_data["run_count"] > 25:
            ctx_data["num_dp"] -= 2
            ctx_data["run_count"] = 0

        return step_to_true(initial + ctx_data["step"], initial, ctx_data, prev_out=new_out)


def gen_sphere_rad(current_step):
    import numpy

    context_dict = {
        "step": 0.015,
        "l_i": (0.09 - 0.060) * numpy.random.random() + 0.060,
        "neutron_i": 100,
        "average_length": 50.0,
        "ideal_out": 1.0,
        "critical_mass_prop": [],
        "lengths": [],
        "run_count": 0,
        "num_dp": 6,
    }

    # final_rad, ctx_data_out = step_to_true(context_dict["l_i"], 0, context_dict)

    # pl.plot(numpy.arange(0.0, len(ctx_data_out["lengths"])), ctx_data_out["lengths"],
    #        label=("Starting at L = {}m".format(ctx_data_out["l_i"])))
    out = list(step_to_true(context_dict["l_i"], 0, context_dict))
    return out


def handle_outputs(results_array):
    final_rad, ctx_data_out = results_array
    pl.plot(numpy.arange(0.0, len(ctx_data_out["lengths"])), ctx_data_out["lengths"],
            label=("Starting at L = {}m".format(ctx_data_out["l_i"])), linewidth=0.5)
    return final_rad


if __name__ == "__main__":
    seaborn.set()
    plot_colours = ["red", "green", "blue"]
    critical_length_a = []
    num_iter = 20
    print("Please note, thread windup-time is around 20 seconds on eight core 3.9GHz "
          "CPU for the default system, little to no activity will be visible before this time.")

    worker_pool = multiprocessing.Pool(multiprocessing.cpu_count()-2)
    print(f"Processing pool of {multiprocessing.cpu_count()-2} workers established.")
    sim_results = []

    async_sim_result = worker_pool.imap_unordered(gen_sphere_rad, list(range(num_iter)))
    time.sleep(0.05)
    main_prog_bar = tqdm.tqdm(async_sim_result, desc="Simulation Progress", ncols=100,
                              unit="sim", smoothing=0, total=num_iter)

    for simulation_result in main_prog_bar:
        critical_rad = handle_outputs(simulation_result)
        critical_length_a.append(critical_rad)

    time.sleep(0.05)

    print("Mean:", numpy.mean(critical_length_a))
    print("Range:", (max(critical_length_a) - min(critical_length_a))/2)
    print("Stdev:", numpy.std(critical_length_a))

    pl.xlabel("Iteration Step")
    pl.ylabel("Radius of U-235 Sphere (m)")
    pl.title("Plot of critical value calculation\n for sphere of uranium with radius\n"
             " L, given in meters.")
    pl.savefig("./length_oscillation_multi_plot_sphere_2.png", dpi=300)
