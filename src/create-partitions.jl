export create_partitions

"""
    df_assets, df_flows = create_partitions(graph, representative_periods; callback=...)

Creates the dataframes storing the partition specification of the representative periods per asset/flow.

The `callback` function will be called as `callback(df_assets, df_flows)`.
Use [`callback_partition_list`](@ref) to see some examples of implemented callback functions.
"""
function create_partitions(graph, representative_periods; kwargs...)
    _create_partitions(
        MetaGraphsNext.labels(graph) |> collect,
        MetaGraphsNext.edge_labels(graph) |> collect,
        [rp.time_steps for rp in representative_periods];
        kwargs...,
    )
end

function _create_partitions(
    assets::Vector{String},
    flows::Vector{Tuple{String,String}},
    representative_periods::Vector{UnitRange{Int}};
    callback = (args...) -> nothing,
)
    num_assets = length(assets)
    num_rps = length(representative_periods)
    df_assets = DataFrame(;
        asset = repeat(assets; outer = num_rps),
        rep_period = repeat(1:length(representative_periods); inner = num_assets),
        specification = :uniform,
        partition = "1",
    )

    num_flows = length(flows)
    flows = repeat(flows; outer = num_rps)
    df_flows = DataFrame(;
        from = getfield.(flows, 1),
        to = getfield.(flows, 2),
        rep_period = repeat(1:length(representative_periods); inner = num_flows),
        specification = :uniform,
        partition = "1",
    )

    callback(df_assets, df_flows)

    return df_assets, df_flows
end

"""
    callback_partition_list()

returns a list of callback functions. See their help to understand how to use them.
"""
function callback_partition_list()
    return [callback_partition_read_from_csv]
end

function callback_partition_read_from_csv(filepath)
    return (df_assets, df_flows) -> begin
        df_new_assets = rename(
            read_csv_with_schema(
                joinpath(filepath, "assets-rep-periods-partitions.csv"),
                AssetRepPeriodPartitionData,
            ),
        )
        df_new_flows = rename(
            read_csv_with_schema(
                joinpath(filepath, "flows-rep-periods-partitions.csv"),
                FlowRepPeriodPartitionData,
            ),
            :from_asset => :from,
            :to_asset => :to,
        )

        # TODO: Is there a better way?
        for row in eachrow(df_new_assets)
            i = findfirst(df_assets.asset .== row.asset .&& df_assets.rep_period .== row.rep_period)
            @assert !isnothing(i)
            df_assets.specification[i] = row.specification
            df_assets.partition[i] = row.partition
        end

        for row in eachrow(df_new_flows)
            i = findfirst(
                df_flows.from .== row.from .&&
                df_flows.to .== row.to .&&
                df_flows.rep_period .== row.rep_period,
            )
            @assert !isnothing(i)
            df_flows.specification[i] = row.specification
            df_flows.partition[i] = row.partition
        end

        return nothing
    end
end

function callback_partition_filter_and_set(asset_filter, flow_filter, specification, partition)
    return (df_assets, df_flows) -> begin
        df_filtered = filter(asset_filter, df_assets; view = true)
        df_filtered.specification .= specification
        df_filtered.partition .= partition
        df_filtered = filter(flow_filter, df_flows; view = true)
        df_filtered.specification .= specification
        df_filtered.partition .= partition
    end
end
