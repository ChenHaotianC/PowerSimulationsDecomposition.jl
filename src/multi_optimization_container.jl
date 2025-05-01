Base.@kwdef mutable struct MultiOptimizationContainer{T <: DecompositionAlgorithm} <:
                           ISOPT.AbstractOptimizationContainer
    main_problem::PSI.OptimizationContainer
    subproblems::Dict{String, PSI.OptimizationContainer}
    time_steps::UnitRange{Int}
    resolution::Dates.TimePeriod
    settings::PSI.Settings
    settings_copy::PSI.Settings
    variables::Dict{ISOPT.VariableKey, AbstractArray}
    aux_variables::Dict{ISOPT.AuxVarKey, AbstractArray}
    duals::Dict{ISOPT.ConstraintKey, AbstractArray}
    constraints::Dict{ISOPT.ConstraintKey, AbstractArray}
    objective_function::PSI.ObjectiveFunction
    expressions::Dict{ISOPT.ExpressionKey, AbstractArray}
    parameters::Dict{ISOPT.ParameterKey, PSI.ParameterContainer}
    primal_values_cache::PSI.PrimalValuesCache
    initial_conditions::Dict{ISOPT.InitialConditionKey, Vector{<:PSI.InitialCondition}}
    initial_conditions_data::PSI.InitialConditionsData
    base_power::Float64
    optimizer_stats::ISOPT.OptimizerStats  # TODO: needs custom struct
    built_for_recurrent_solves::Bool
    metadata::ISOPT.OptimizationContainerMetadata
    default_time_series_type::Type{<:PSY.TimeSeriesData}  # Maybe isn't needed here
    mpi_info::Union{Nothing, MpiInfo}
    if_coordination::Int
end

function MultiOptimizationContainer(
    ::Type{T},
    sys::PSY.System,
    settings::PSI.Settings,
    ::Type{U},
    subproblem_keys::Vector{String},
) where {T <: DecompositionAlgorithm, U <: PSY.TimeSeriesData}
    resolutions = PSY.get_time_series_resolutions(sys)

    if length(resolutions) > 1
        error(
            "Multiple time series resolutions not supported for MultiOptimizationContainer",
        )
    else
        resolution = resolutions[1]
    end

    if isabstracttype(U)
        error("Default Time Series Type $U can't be abstract")
    end

    # define dictionary containing the optimization container for the subregion
    subproblems = Dict(
        k => PSI.OptimizationContainer(sys, settings, nothing, U) for k in subproblem_keys
    )

    return MultiOptimizationContainer{T}(;
        main_problem=PSI.OptimizationContainer(sys, settings, nothing, U),
        subproblems=subproblems,
        time_steps=1:1,
        resolution=IS.time_period_conversion(resolution),
        settings=settings,
        settings_copy=PSI.copy_for_serialization(settings),
        variables=Dict{ISOPT.VariableKey, AbstractArray}(),
        aux_variables=Dict{ISOPT.AuxVarKey, AbstractArray}(),
        duals=Dict{ISOPT.ConstraintKey, AbstractArray}(),
        constraints=Dict{ISOPT.ConstraintKey, AbstractArray}(),
        objective_function=PSI.ObjectiveFunction(),
        expressions=Dict{ISOPT.ExpressionKey, AbstractArray}(),
        parameters=Dict{ISOPT.ParameterKey, PSI.ParameterContainer}(),
        primal_values_cache=PSI.PrimalValuesCache(),
        initial_conditions=Dict{ISOPT.InitialConditionKey, Vector{PSI.InitialCondition}}(),
        initial_conditions_data=PSI.InitialConditionsData(),
        base_power=PSY.get_base_power(sys),
        optimizer_stats=ISOPT.OptimizerStats(),
        built_for_recurrent_solves=false,
        metadata=ISOPT.OptimizationContainerMetadata(),
        default_time_series_type=U,
        mpi_info=nothing,
        if_coordination=2,
    )
end

function get_container_keys(container::MultiOptimizationContainer)
    return Iterators.flatten(keys(getfield(container, f)) for f in STORE_CONTAINERS)
end

PSI.get_default_time_series_type(container::MultiOptimizationContainer) =
    container.default_time_series_type
PSI.get_duals(container::MultiOptimizationContainer) = container.duals
PSI.get_expressions(container::MultiOptimizationContainer) = container.expressions
PSI.get_initial_conditions(container::MultiOptimizationContainer) =
    container.initial_conditions
PSI.get_initial_conditions_data(container::MultiOptimizationContainer) =
    container.initial_conditions_data
PSI.get_initial_time(container::MultiOptimizationContainer) =
    PSI.get_initial_time(container.settings)
PSI.get_jump_model(container::MultiOptimizationContainer) =
    PSI.get_jump_model(container.main_problem)
PSI.get_metadata(container::MultiOptimizationContainer) = container.metadata
PSI.get_optimizer_stats(container::MultiOptimizationContainer) = container.optimizer_stats
PSI.get_parameters(container::MultiOptimizationContainer) = container.parameters
PSI.get_resolution(container::MultiOptimizationContainer) = container.resolution
PSI.get_settings(container::MultiOptimizationContainer) = container.settings
PSI.get_time_steps(container::MultiOptimizationContainer) = container.time_steps
PSI.get_variables(container::MultiOptimizationContainer) = container.variables

PSI.set_initial_conditions_data!(container::MultiOptimizationContainer, data) =
    container.initial_conditions_data = data
PSI.get_objective_expression(container::MultiOptimizationContainer) =
    container.objective_function
PSI.is_synchronized(container::MultiOptimizationContainer) =
    container.objective_function.synchronized
PSI.set_time_steps!(container::MultiOptimizationContainer, time_steps::UnitRange{Int64}) =
    container.time_steps = time_steps

PSI.get_aux_variables(container::MultiOptimizationContainer) = container.aux_variables
PSI.get_base_power(container::MultiOptimizationContainer) = container.base_power
PSI.get_constraints(container::MultiOptimizationContainer) = container.constraints

get_if_coordination(container::MultiOptimizationContainer) = container.if_coordination

function get_subproblem(container::MultiOptimizationContainer, id::String)
    return container.subproblems[id]
end

function check_optimization_container(container::MultiOptimizationContainer)
    for subproblem in values(container.subproblems)
        PSI.check_optimization_container(subproblem)
    end
    PSI.check_optimization_container(container.main_problem)
    return
end

function _finalize_jump_model!(
    container::MultiOptimizationContainer,
    settings::PSI.Settings,
)
    @debug "Instantiating the JuMP model" _group = PSI.LOG_GROUP_OPTIMIZATION_CONTAINER
    PSI._finalize_jump_model!(container.main_problem, settings)
    return
end

function init_optimization_container!(
    container::MultiOptimizationContainer,
    network_model::PSI.NetworkModel{<:PM.AbstractPowerModel},
    sys::PSY.System,
)
    PSY.set_units_base_system!(sys, "SYSTEM_BASE")
    # The order of operations matter
    settings = PSI.get_settings(container)

    if PSI.get_initial_time(settings) == PSI.UNSET_INI_TIME
        if PSI.get_default_time_series_type(container) <: PSY.AbstractDeterministic
            PSI.set_initial_time!(settings, PSY.get_forecast_initial_timestamp(sys))
        elseif PSI.get_default_time_series_type(container) <: PSY.SingleTimeSeries
            ini_time, _ = PSY.check_time_series_consistency(sys, PSY.SingleTimeSeries)
            PSI.set_initial_time!(settings, ini_time)
        else
            error("Bug: unhandled $(PSI.get_default_time_series_type(container))")
        end
    end

    # TODO: what if the time series type is SingleTimeSeries?
    if PSI.get_horizon(settings) == PSI.UNSET_HORIZON
        PSI.set_horizon!(settings, PSY.get_forecast_horizon(sys))
    end
    horizon_count = (PSI.get_horizon(settings) ÷ PSI.get_resolution(settings))
    @assert horizon_count > 0
    container.time_steps = 1:horizon_count

    stats = PSI.get_optimizer_stats(container)
    stats.detailed_stats = PSI.get_detailed_optimizer_stats(settings)

    # need a special method for the main problem to initialize the optimization container
    # without actually caring about the subnetworks
    # PSI.init_optimization_container!(subproblem, network_model, sys)

    for (index, subproblem) in container.subproblems
        @debug "Initializing Container Subproblem $index" _group =
            PSI.LOG_GROUP_OPTIMIZATION_CONTAINER
        subproblem.settings = deepcopy(settings)
        PSI.init_optimization_container!(subproblem, network_model, sys)
        subproblem.built_for_recurrent_solves = true
    end
    _finalize_jump_model!(container, settings)
    return
end

function PSI.serialize_optimization_model(
    container::MultiOptimizationContainer,
    save_path::String,
) end
