#%%
using Chain, DataFrames, JuMP, MetaGraphsNext, TulipaEnergyModel

println("""
Instructions:

Change the `input_dir` to the desired problem and run again
- Tiny is very small
- Norse is still small wrt `time_block`s
- EU is very large

Search for `#%%` to find the sections in this file.

""")

#%% Strategies

function current_best(tbl_cons, tbl_flows)
    incoming = Vector{AffExpr}(undef, length(tbl_cons.asset))
    for row in Tables.rows(tbl_cons)
        incoming[row.index] = AffExpr(0.0)
        idx = findall(
            (tbl_flows.rp .== row.rp) .&&
            (tbl_flows.to .== row.asset) .&&
            (last.(tbl_flows.time_block) .≥ row.time_block[1]) .&&
            (first.(tbl_flows.time_block) .≤ row.time_block[end]),
        )
        if length(idx) > 0
            incoming[row.index] = sum(
                tbl_flows.flow[idx] .* length.(Ref(row.time_block) .∩ tbl_flows.time_block[idx]) * 3.14row.rp,
            )
        end
    end
end

function also_decent(tbl_cons, df_flows)
    [
        begin
            df = @view df_flows[
                (df_flows.rp.==row.rp).&&(df_flows.to.==row.asset).&&(last.(
                    df_flows.time_block
                ).≥row.time_block[1]).&&(first.(df_flows.time_block).≤row.time_block[end]),
                :,
            ]
            if size(df, 1) > 0
                sum(df.flow .* length.(Ref(row.time_block) .∩ df.time_block) * 3.14row.rp)
            else
                AffExpr(0.0)
            end
        end for row in Tables.rows(tbl_cons)
    ]
end

function older_strategy(df_cons, df_flows)
    transform!(
        df_cons,
        [:rp, :asset, :time_block] =>
            ByRow(
                (rp, a, time_block) -> begin
                    df = @view df_flows[
                        (df_flows.rp.==rp).&&(df_flows.to.==a).&&(last.(
                            df_flows.time_block
                        ).≥first.(
                            Ref(time_block)
                        )).&&(first.(df_flows.time_block).≤last.(Ref(time_block))),
                        :,
                    ]
                    coalesce(
                        sum(df.flow .* length.(Ref(time_block) .∩ df.time_block) * 3.14rp),
                        AffExpr(0.0),
                    )
                end,
            ) => :incoming_term,
    )
end

function leftjoin_strategy(df_cons, df_flows)
    @chain df_cons begin
        leftjoin(df_flows; on = [:rp, :asset => :to], makeunique = true)
        groupby([:asset, :rp, :time_block, :index])
        combine(
            [:time_block, :time_block_1, :flow, :rp] =>
                (
                    (t, t1, fl, rp) -> sum(
                        3.14 * length(t[i] ∩ t1[i]) * rp[i] * fl[i] for
                        i = 1:length(fl) if !ismissing(fl[i]);
                        init = AffExpr(0.0),
                    )
                ) => :incoming_term,
        )
    end
end

#%% Reading and preparing data

input_dir = joinpath(@__DIR__, "test/inputs/Tiny")
# input_dir = joinpath(@__DIR__, "test/inputs/Norse")
# input_dir = joinpath(@__DIR__, "benchmark/EU")

graph, representative_periods = create_graph_and_representative_periods_from_csv_folder(input_dir)
constraints_partitions = compute_constraints_partitions(graph, representative_periods)

df_flows = DataFrame(
    (
        ((from = u, to = v, rp = rp, time_block = B) for B ∈ graph[u, v].partitions[rp]) for
        (u, v) ∈ edge_labels(graph), rp = 1:length(representative_periods)
    ) |> Iterators.flatten,
)
df_flows.index = 1:size(df_flows, 1)

df_cons = DataFrame(
    (
        (
            (asset = a, rp = rp, time_block = B) for
            B ∈ constraints_partitions[:lowest_resolution][(a, rp)]
        ) for a ∈ labels(graph), rp = 1:length(representative_periods)
    ) |> Iterators.flatten,
)
df_cons.index = 1:size(df_cons, 1)

model = Model()
flow =
    model[:flow] =
        df_flows.flow = [
            begin
                u, v, rp, time_block = row.from, row.to, row.rp, row.time_block
                @variable(model, base_name = "flow[($u, $v), $rp, $time_block]")
            end for row in eachrow(df_flows)
        ]

tbl_flows = Tables.columntable(df_flows)
tbl_cons = Tables.columntable(df_cons)

#%% Timing

println("Strategies timing on the data from $input_dir")

if input_dir |> endswith("EU")
    println("Skipping other strategies")
else
    println("Current best:")
    current_best(tbl_cons, tbl_flows)
    @time current_best(tbl_cons, tbl_flows)

    println("Also decent:")
    also_decent(tbl_cons, df_flows)
    @time also_decent(tbl_cons, df_flows)

    println("Older strategy")
    older_strategy(df_cons, df_flows)
    @time older_strategy(df_cons, df_flows)
end

println("Leftjoin strategy")
leftjoin_strategy(df_cons, df_flows)
@time leftjoin_strategy(df_cons, df_flows)

println("Done")
