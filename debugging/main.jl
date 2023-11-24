#%%
using TulipaEnergyModel
using BenchmarkTools
using Profile

# input_dir = "test/inputs/Norse"
# input_dir = "debugging/LargeNorse/"
input_dir = "debugging/input_NLBE/"

#%%

@time graph, representative_periods =
    create_graph_and_representative_periods_from_csv_folder(input_dir);

#%%

@time constraints_partitions = compute_constraints_partitions(graph, representative_periods);

#%%

@time model = create_model(graph, representative_periods, constraints_partitions);
# @profview create_model(graph, representative_periods, constraints_partitions);

#%%

@time solution = solve_model(model);
