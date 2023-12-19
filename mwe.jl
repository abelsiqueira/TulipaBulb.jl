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
    outgoing = Vector{AffExpr}(undef, length(tbl_cons.asset))
    for row in Tables.rows(tbl_cons)
        incoming[row.index] = AffExpr(0.0)
        outgoing[row.index] = AffExpr(0.0)
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
        idx = findall(
            (tbl_flows.rp .== row.rp) .&&
            (tbl_flows.from .== row.asset) .&&
            (last.(tbl_flows.time_block) .≥ row.time_block[1]) .&&
            (first.(tbl_flows.time_block) .≤ row.time_block[end]),
        )
        if length(idx) > 0
            outgoing[row.index] = sum(
                tbl_flows.flow[idx] .* length.(Ref(row.time_block) .∩ tbl_flows.time_block[idx]) * 3.14row.rp,
            )
        end
    end
    return incoming, outgoing
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

function dan_strategy(df_cons, df_flows)
    Tmin = min(minimum(first.(df_flows.time_block)), minimum(first.(df_cons.time_block)))
    Tmax = max(maximum(last.(df_flows.time_block)), maximum(last.(df_cons.time_block)))
    Tspan = Tmin:Tmax
    nT = length(Tspan)
    from_list = unique(df_flows.from)
    from_lookup = Dict(a => i for (i, a) in enumerate(from_list))
    nA = length(from_list)

    pre_incoming = Matrix{AffExpr}(undef, nT, nA)
    df_cons.incoming_term_dan .= AffExpr(0.0)

    g_flows = groupby(df_flows, [:rp, :to])
    g_cons = groupby(df_cons, [:rp, :asset])

    for ((rp, to), sdf) in pairs(g_cons)
        haskey(g_flows, (; rp, to)) || continue
        pre_incoming .= AffExpr(0.0)
        for row in eachrow(g_flows[(; rp, to)]), t in row.time_block
            j = from_lookup[row.from]
            pre_incoming[t-Tmin+1, j] = row.flow * 3.14 * rp
        end
        for row in eachrow(sdf)
            row.incoming_term_dan =
                sum(sum(pre_incoming[t-Tmin+1, :]) for t in row.time_block; init = AffExpr(0.0))
        end
    end
end

function rocco_strategy(df_cons, df_flows)
    df_cons.incoming_term_rocco .= AffExpr(0)
    # sort!(df_cons, [:rp, :asset])
    # sort!(df_flows, [:rp, :to])

    function cons_inc_term!(frp, fto, crp, cass, ftb, ctb, flow, it)
        rp_to = collect(zip(frp, fto))
        idr = [searchsorted(rp_to, (a, r)) for (a, r) in zip(crp, cass)]
        id = findall(!isempty, idr)
        for i in id
            len = [length(∩(ctb[i], fl)) * 3.14 * crp[i] for fl in ftb[idr[i]]]
            it[i] = len' * flow[idr[i]]
        end
    end

    cons_inc_term!(
        df_flows.rp,
        df_flows.to,
        df_cons.rp,
        df_cons.asset,
        df_flows.time_block,
        df_cons.time_block,
        df_flows.flow,
        df_cons.incoming_term_rocco,
    )
end

#%% Reading and preparing data

# input_dir = joinpath(@__DIR__, "test/inputs/Tiny")
# input_dir = joinpath(@__DIR__, "test/inputs/Norse")
# input_dir = joinpath(@__DIR__, "test/inputs/Variable Resolution/")
input_dir = joinpath(@__DIR__, "benchmark/EU")
# input_dir = mktempdir()

if input_dir |> startswith("/tmp")
    NORSE_PATH = joinpath(@__DIR__, "test/inputs/Norse")
    new_rp_length = 100
    for file in readdir(NORSE_PATH; join = false)
        cp(joinpath(NORSE_PATH, file), joinpath(input_dir, file))
    end
    # Add another line to rep-periods-data.csv
    open(joinpath(input_dir, "rep-periods-data.csv"), "a") do io
        println(io, "3,$new_rp_length,0.1")
    end
    open(joinpath(input_dir, "rep-periods-mapping.csv"), "a") do io
        println(io, "216,3,1.0")
    end
    # Add profiles to flow and asset
    open(joinpath(input_dir, "flows-profiles.csv"), "a") do io
        for (u, v) in [("Asgard_E_demand", "Valhalla_E_balance")]
            for i = 1:new_rp_length
                println(io, "$u,$v,3,$i,0.95")
            end
        end
    end
    open(joinpath(input_dir, "assets-profiles.csv"), "a") do io
        for a in ["Asgard_E_demand"]
            for i = 1:new_rp_length
                println(io, "$a,3,$i,0.95")
            end
        end
    end
end

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
sort!(df_cons, [:rp, :asset])
sort!(df_flows, [:rp, :to])

# df_flows = DataFrame(;
#     from = ["p1", "p1", "p2", "p2", "p1", "p2"],
#     to = ["d", "d", "d", "d", "d", "d"],
#     rp = [1, 1, 1, 1, 2, 2],
#     time_block = [1:3, 4:6, 1:2, 3:6, 1:12, 1:12],
#     index = 1:6,
# )
# df_cons = DataFrame(;
#     asset = ["p1", "p1", "p1", "p2", "p2", "d", "d", "d"],
#     rp = [1, 1, 2, 1, 2, 1, 1, 2],
#     time_block = [1:3, 4:6, 1:12, 1:6, 1:12, 1:4, 5:6, 1:12],
#     index = 1:8,
# )

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

    # println("Also decent:")
    # also_decent(tbl_cons, df_flows)
    # @time also_decent(tbl_cons, df_flows)

    # println("Older strategy")
    # older_strategy(df_cons, df_flows)
    # @time older_strategy(df_cons, df_flows)
end

# println("Leftjoin strategy")
# leftjoin_strategy(df_cons, df_flows)
# @time leftjoin_strategy(df_cons, df_flows)

println("Dan strategy")
dan_strategy(df_cons, df_flows)
@time dan_strategy(df_cons, df_flows)

println("Rocco strategy")
rocco_strategy(df_cons, df_flows)
@time rocco_strategy(df_cons, df_flows)

expression_is_zero(e::AffExpr) = e == 0 || all(collect(values(e.terms)) .== 0)
expressions_are_equal(e1, e2) = all(expression_is_zero.(e1 - e2))
@info expressions_are_equal(df_cons.incoming_term_dan, df_cons.incoming_term_rocco)

println("Done")
