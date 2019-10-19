import uranium_sim
import pylab as pl
import seaborn
import multiprocessing
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


def gen_sphere_rad(run_point):
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
    global critical_length_a
    final_rad, ctx_data_out = results_array
    critical_length_a.append(final_rad)
    pl.plot(numpy.arange(0.0, len(ctx_data_out["lengths"])), ctx_data_out["lengths"],
            label=("Starting at L = {}m".format(ctx_data_out["l_i"])), linewidth=0.5)


if __name__ == "__main__":
    seaborn.set()
    plot_colours = ["red", "green", "blue"]
    critical_length_a = []
    num_iter = 40

    worker_pool = multiprocessing.Pool(multiprocessing.cpu_count()-2)
    sim_results = []

    for iterable in range(num_iter):
        async_sim_result = worker_pool.apply_async(gen_sphere_rad, args=[iterable], callback=handle_outputs)
        sim_results.append(async_sim_result)

    main_prog_bar = tqdm.tqdm(sim_results, desc="Simulation Progress", ncols=100, unit="sim", smoothing=0.5)
    for result in main_prog_bar:
        result.wait()

    print("Mean:", numpy.mean(critical_length_a))
    print("Range:", (max(critical_length_a) - min(critical_length_a))/2)
    print("Stdev:", numpy.std(critical_length_a))

    pl.xlabel("Iteration Step")
    pl.ylabel("Radius of U-235 Sphere (m)")
    pl.title("Plot of critical value calculation\n for sphere of uranium with radius\n"
             " L, given in meters.")
    pl.savefig("./length_oscillation_multi_plot_sphere_2.png", dpi=300)
