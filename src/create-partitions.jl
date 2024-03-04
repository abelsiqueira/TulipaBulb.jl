export create_partitions

function create_partitions(graph, representative_periods)
    create_partitions(
        MetaGraphsNext.labels(graph) |> collect,
        MetaGraphsNext.edge_labels(graph) |> collect,
        [rp.time_steps for rp in representative_periods],
    )
end

function create_partitions(
    assets::Vector{String},
    flows::Vector{Tuple{String,String}},
    representative_periods::Vector{UnitRange{Int}},
)
    num_assets = length(assets)
    num_rps = length(representative_periods)
    df_assets = DataFrame(;
        asset = repeat(assets; outer = num_rps),
        rps = repeat(1:length(representative_periods); inner = num_assets),
        specification = :uniform,
        partition = "1",
    )

    num_flows = length(flows)
    flows = repeat(flows; outer = num_rps)
    df_flows = DataFrame(;
        from = getfield.(flows, 1),
        to = getfield.(flows, 2),
        rps = repeat(1:length(representative_periods); inner = num_flows),
        specification = :uniform,
        partition = "1",
    )

    return df_assets, df_flows
end
