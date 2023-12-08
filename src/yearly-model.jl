export create_model_yearly

"""
    model = create_model_yearly(graph, representative_periods)

Create the energy model given the graph and representative_periods.
"""
function create_model_yearly(
    graph,
    representative_periods,
    constraints_partitions;
    verbose = false,
    write_lp_file = false,
)

    ## Helper functions
    # Computes the duration of the `block` that is within the `period`, and
    # multiply by the resolution of the representative period `rp`.
    # It is equivalent to finding the indexes of these values in the matrix.
    # function duration(B1, B2, rp)
    #     return length(B1 ∩ B2) * representative_periods[rp].resolution
    # end

    # function duration(B, rp)
    #     return length(B) * representative_periods[rp].resolution
    # end

    # # Sums the profile of representative period rp over the time block B
    # # Uses the default_value when that profile does not exist.
    # function profile_sum(profiles, rp, B, default_value)
    #     if haskey(profiles, rp)
    #         return sum(profiles[rp][B])
    #     else
    #         return length(B) * default_value
    #     end
    # end

    # function assets_profile_sum(a, t, default_value)
    #     return profile_sum(graph[a].profiles, rp, B, default_value)
    # end

    # # Same as above but for flow
    # function flows_profile_sum(u, v, rp, B, default_value)
    #     return profile_sum(graph[u, v].profiles, rp, B, default_value)
    # end

    ## Sets unpacking
    A = labels(graph) |> collect
    F = edge_labels(graph) |> collect
    filter_assets(key, value) =
        Iterators.filter(a -> getfield(graph[a], key) == value, A) |> collect
    filter_flows(key, value) =
        Iterators.filter(f -> getfield(graph[f...], key) == value, F) |> collect

    Ac = filter_assets(:type, "consumer")
    Ap = filter_assets(:type, "producer")
    Ai = filter_assets(:investable, true)
    As = filter_assets(:type, "storage")
    Ah = filter_assets(:type, "hub")
    Acv = filter_assets(:type, "conversion")
    Fi = filter_flows(:investable, true)
    Ft = filter_flows(:is_transport, true)
    # RP = 1:length(representative_periods)
    P = 1:(365*24)

    ## Model
    model = Model(HiGHS.Optimizer)
    set_attribute(model, "output_flag", verbose)

    ## Variables
    @variable(model, flow[(u, v) ∈ F, P])
    @variable(model, 0 ≤ assets_investment[Ai], Int)  #number of installed asset units [N]
    @variable(model, 0 ≤ flows_investment[Fi], Int)
    @variable(model, 0 ≤ storage_level[a ∈ As, P])

    # TODO: Fix storage_level[As, RP, 0] = 0

    ## Expressions for the objective function
    assets_investment_cost = @expression(
        model,
        sum(graph[a].investment_cost * graph[a].capacity * assets_investment[a] for a ∈ Ai)
    )

    flows_investment_cost = @expression(
        model,
        sum(
            graph[u, v].investment_cost * graph[u, v].unit_capacity * flows_investment[(u, v)]
            for (u, v) ∈ Fi
        )
    )

    flows_variable_cost =
        @expression(model, sum(graph[u, v].variable_cost * flow[(u, v), t] for (u, v) ∈ F, t ∈ P))

    ## Objective function
    @objective(model, Min, assets_investment_cost + flows_investment_cost + flows_variable_cost)

    ## Expressions for balance constraints
    @expression(
        model,
        incoming_flow[a ∈ A, t ∈ P],
        sum(flow[(u, a), t] for u in inneighbor_labels(graph, a))
    )

    @expression(
        model,
        outgoing_flow[a ∈ A, t ∈ P],
        sum(flow[(a, v), t] for v in outneighbor_labels(graph, a))
    )

    @expression(
        model,
        incoming_flow_w_efficiency[a ∈ A, t ∈ P],
        sum(flow[(u, a), t] * graph[u, a].efficiency for u in inneighbor_labels(graph, a))
    )

    @expression(
        model,
        outgoing_flow_w_efficiency[a ∈ A, t ∈ P],
        sum(flow[(a, v), t] / graph[a, v].efficiency for v in outneighbor_labels(graph, a))
    )

    @expression(
        model,
        energy_limit[a ∈ As∩Ai],
        graph[a].energy_to_power_ratio * graph[a].capacity * assets_investment[a]
    )

    @expression(
        model,
        storage_inflows[a ∈ As, t ∈ P],
        3.14 * (graph[a].initial_storage_capacity + (a ∈ Ai ? energy_limit[a] : 0.0))
    )

    ## Balance constraints (using the lowest resolution)
    # - consumer balance equation
    @constraint(
        model,
        consumer_balance[a ∈ Ac, t ∈ P],
        incoming_flow[a, t] - outgoing_flow[a, t] == 3.14 * graph[a].peak_demand
    )

    # - storage balance equation
    @constraint(
        model,
        storage_balance[a ∈ As, t ∈ P],
        storage_level[a, t] ==
        (t > 1 ? storage_level[a, t-1] : graph[a].initial_storage_level) +
        storage_inflows[a, t] +
        incoming_flow_w_efficiency[a, t] - outgoing_flow_w_efficiency[a, t]
    )

    # - hub balance equation
    @constraint(model, hub_balance[a ∈ Ah, t ∈ P], incoming_flow[a, t] == outgoing_flow[a, t])

    # - conversion balance equation
    @constraint(
        model,
        conversion_balance[a ∈ Acv, t ∈ P],
        incoming_flow_w_efficiency[a, t] == outgoing_flow_w_efficiency[a, t]
    )

    ## Expression for capacity limit constraints
    @expression(
        model,
        assets_profile_times_capacity[a ∈ A, t ∈ P],
        3.14 *
        (graph[a].initial_capacity + (a ∈ Ai ? (graph[a].capacity * assets_investment[a]) : 0.0))
    )

    ## Capacity limit constraints (using the highest resolution)
    # - maximum output flows limit
    @constraint(
        model,
        max_output_flows_limit[a ∈ Acv∪As∪Ap, t ∈ P],
        outgoing_flow[a, t] ≤ assets_profile_times_capacity[a, t]
    )

    # - maximum input flows limit
    @constraint(
        model,
        max_input_flows_limit[a ∈ As, t ∈ P],
        incoming_flow[a, t] ≤ assets_profile_times_capacity[a, t]
    )

    # - define lower bounds for flows that are not transport assets
    for f ∈ F, t ∈ P
        if f ∉ Ft
            set_lower_bound(flow[f, t], 0.0)
        end
    end

    ## Expressions for transport flow constraints
    @expression(
        model,
        upper_bound_transport_flow[(u, v) ∈ F, t ∈ P],
        3.14 * (
            graph[u, v].initial_capacity +
            (graph[u, v].investable ? graph[u, v].export_capacity * flows_investment[(u, v)] : 0.0)
        )
    )

    @expression(
        model,
        lower_bound_transport_flow[(u, v) ∈ F, t ∈ P],
        3.14 * (
            graph[u, v].initial_capacity +
            (graph[u, v].investable ? graph[u, v].import_capacity * flows_investment[(u, v)] : 0.0)
        )
    )

    ## Constraints that define bounds for a transport flow Ft
    @constraint(
        model,
        max_transport_flow_limit[f ∈ Ft, t ∈ P],
        flow[f, t] ≤ upper_bound_transport_flow[f, t]
    )

    @constraint(
        model,
        min_transport_flow_limit[f ∈ Ft, t ∈ P],
        flow[f, t] ≥ -lower_bound_transport_flow[f, t]
    )

    ## Extra constraints for storage assets
    # - maximum storage level limit
    @constraint(
        model,
        max_storage_level_limit[a ∈ As, t ∈ P],
        storage_level[a, t] ≤ graph[a].initial_storage_capacity + (a ∈ Ai ? energy_limit[a] : 0.0)
    )

    # - cycling condition for storage level
    for a ∈ As
        set_lower_bound(storage_level[a, end], graph[a].initial_storage_level)
    end

    ## Expressions for the extra constraints of investment limits
    @expression(
        model,
        flow_max_capacity[(u, v) ∈ Fi],
        max(graph[u, v].export_capacity, graph[u, v].import_capacity)
    )

    ## Extra constraints for investment limits
    # - maximum (i.e., potential) investment limit for assets
    for a ∈ Ai
        if graph[a].capacity > 0 && !ismissing(graph[a].investment_limit)
            set_upper_bound(assets_investment[a], graph[a].investment_limit / graph[a].capacity)
        end
    end

    # - maximum (i.e., potential) investment limit for flows
    for (u, v) ∈ Fi
        if flow_max_capacity[(u, v)] > 0 && !ismissing(graph[u, v].investment_limit)
            set_upper_bound(
                flows_investment[(u, v)],
                graph[u, v].investment_limit / flow_max_capacity[(u, v)],
            )
        end
    end

    if write_lp_file
        write_to_file(model, "model.lp")
    end

    return model
end
