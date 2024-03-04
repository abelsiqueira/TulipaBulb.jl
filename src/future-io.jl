function read_input_data_from_csv_folder(input_folder)
    fillpath(filename) = joinpath(input_folder, filename)
    time_invariant_input = Dict(
        :assets => DataFrames.rename(
            read_csv_with_schema(fillpath("assets-data.csv"), AssetData),
            :name => :asset,
        ),
        :flows => DataFrames.rename(
            read_csv_with_schema(fillpath("flows-data.csv"), FlowData),
            :from_asset => :from,
            :to_asset => :to,
        ),
    )
    clustered_profiles = Dict(
        :assets => read_csv_with_schema(fillpath("assets-profiles.csv"), AssetProfiles),
        :flows => read_csv_with_schema(fillpath("flows-profiles.csv"), FlowProfiles),
    )
    clustering_info = Dict(
        :rep_periods => read_csv_with_schema(fillpath("rep-periods-data.csv"), RepPeriodData),
        :rp_mapping_df =>
            read_csv_with_schema(fillpath("rep-periods-mapping.csv"), RepPeriodMapping),
    )

    time_invariant_input, clustered_profiles, clustering_info
end

function create_graph(time_invariant_input)
    df_assets, df_flows = time_invariant_input[:assets], time_invariant_input[:flows]

    num_assets = size(df_assets, 1)
    asset_lookup = Dict(row.asset => i for (i, row) in enumerate(eachrow(df_assets)))
    _graph = Graphs.DiGraph(num_assets)
    for row in eachrow(df_flows)
        Graphs.add_edge!(_graph, asset_lookup[row.from], asset_lookup[row.to])
    end

    graph = MetaGraphsNext.MetaGraph(
        _graph,
        [row.asset => (index = i) for (i, row) in enumerate(eachrow(df_assets))],
        [(row.from, row.to) => (index = i) for (i, row) in enumerate(eachrow(df_flows))],
        nothing,
        nothing,
        nothing,
    )

    return graph
end

function create_partition(df_assets, df_flows, df_rep_periods; callback = identity)
    df_assets_partitions = crossjoin(
        df_rep_periods[!, [:id]],
        df_assets[!, [:asset]];
        renamecols = (_ -> :rep_period) => identity,
    )
    df_assets_partitions.specification .= :uniform
    df_assets_partitions.partition .= "1"

    df_flows_partitions = crossjoin(
        df_rep_periods[!, [:id]],
        df_flows[!, [:from, :to]];
        renamecols = (_ -> :rep_period) => identity,
    )
    df_flows_partitions.specification .= :uniform
    df_flows_partitions.partition .= "1"

    return Dict(:assets => df_assets_partitions, :flows => df_flows_partitions)
end

function callback_partition_read_from_csv!(dfs_partitions, filepath)
    df_assets, df_flows = dfs_partitions[:assets], dfs_partitions[:flows]
    df_new_assets = rename(
        read_csv_with_schema(joinpath(filepath, "assets-partitions.csv"), AssetPartitionData),
        :rep_period_id => :rep_period,
    )
    df_new_flows = rename(
        read_csv_with_schema(joinpath(filepath, "flows-partitions.csv"), FlowPartitionData),
        :from_asset => :from,
        :to_asset => :to,
        :rep_period_id => :rep_period,
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

function create_variable_and_constraints(
    dfs_time_independent,
    dfs_profiles,
    dfs_partitions,
    df_rep_periods,
)
    # Generate assets and flows investment
    #  these are the easy ones, since they just depend on a single set
    df_assets_investment = filter(row -> row.investment, dfs_time_independent[:assets])[!, [:asset]]
    df_flows_investment =
        filter(row -> row.investment, dfs_time_independent[:flows])[!, [:from, :to]]

    # Flows, storage intra, storage inter

    df_assets_unrolled_partition = flatten(
        transform(
            dfs_partitions[:assets],
            [:specification, :partition, :rep_period] =>
                ByRow(
                    (specificatoin, partition, rep_period) -> _parse_rp_partition(
                        Val(specificatoin),
                        partition,
                        1:df_rep_periods.num_time_steps[rep_period],
                    ),
                ) => :time_block,
        ),
        :time_block,
    )
    #
    # lowest, highest_in_out, highest_in and highest_out
    #
end
