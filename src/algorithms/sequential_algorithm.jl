function build_impl!(
    container::MultiOptimizationContainer{SequentialAlgorithm},
    template::MultiProblemTemplate,
    sys::PSY.System,
)
    for (index, sub_template) in get_sub_templates(template)
        @info "Building Subproblem $index" _group = PSI.LOG_GROUP_OPTIMIZATION_CONTAINER
        PSI.build_impl!(get_subproblem(container, index), sub_template, sys)
        
        if container.if_coordination > 0
            ## ch add the Relief variables in the model
            subproblem_container = get_subproblem(container, index)
            time_steps = PSI.get_time_steps(subproblem_container)
            # println("Built Subproblem $index", time_steps)
            # println(typeof(time_steps))
            if length(time_steps) == 1
            # uc_sub_network = PSI.get_network_model(sub_template)
                if container.if_coordination == 1
                    segments = [1]
                end
                if container.if_coordination == 2
                    segments = collect(1:15)
                end
                
                monitor_upper_key = PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.RateLimitConstraint, PSY.MonitoredLine}("ub")
                monitor_lower_key = PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.RateLimitConstraint, PSY.MonitoredLine}("lb")
                ####revise the tieline constraint of the monitor_line###
                if haskey(subproblem_container.constraints, monitor_lower_key) && haskey(subproblem_container.constraints, monitor_upper_key)
                    monitored_lines = PSI.get_constraint(subproblem_container,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb").axes[1]
                    num_monitored_lines = length(monitored_lines)
                    monitored_lines_index = Dict(line => index for (index, line) in enumerate(monitored_lines))
    
                    monitored_lines_segments = ["$(l)_$(s)_$(i)" for l in monitored_lines for s in segments for i in [1,2]]
                    
                    region_relief = PSI.add_variable_container!(
                        subproblem_container,
                        PowerSimulationsDecomposition.Relief(),
                        PSY.MonitoredLine,
                        monitored_lines_segments,
                        time_steps,
                        meta = "",
                    )
                    
                    for l in monitored_lines
                        for s in segments
                            for i in [1,2]
                            
                                index_r = "$(l)_$(s)_$(i)"
                                region_relief[index_r,1] = JuMP.@variable(
                                    PSI.get_jump_model(subproblem_container),
                                    base_name = "Relief_{$(s)_$(i), $(l)}",
                                    lower_bound = 0.0
                                )
                                JuMP.set_upper_bound(region_relief[index_r,1], 0)
                                JuMP.set_lower_bound(region_relief[index_r,1], 0)
                            end
                        end
                    end
                    

                    line_capacity = 0.1
                    for l in monitored_lines
                        monitor_lower_con = PSI.get_constraint(subproblem_container,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb")[l,1]
                        monitor_upper_con = PSI.get_constraint(subproblem_container,PSI.RateLimitConstraint(),PSY.MonitoredLine,"ub")[l,1]
                        JuMP.set_normalized_rhs(monitor_lower_con, -Inf)
                        JuMP.set_normalized_rhs(monitor_upper_con, -line_capacity)
                        for s in segments
                            for i in [1,2]
                                index_r = "$(l)_$(s)_$(i)"
                                JuMP.set_normalized_coefficient(monitor_lower_con, region_relief[index_r,1], -1)
                                JuMP.set_normalized_coefficient(monitor_upper_con, region_relief[index_r,1], -1)
                                JuMP.set_objective_coefficient(subproblem_container.JuMPmodel, region_relief[index_r,1], 10000000)
                            end
                        end
                    end

                end
            end
        end      
    end

    build_main_problem!(container, template, sys)

    check_optimization_container(container)

    return
end

function build_main_problem!(
    container::MultiOptimizationContainer{SequentialAlgorithm},
    template::MultiProblemTemplate,
    sys::PSY.System,
) end

# The drawback of this approach is that it will loop over the results twice
# once to write into the main container and a second time when writing into the
# store. The upside of this approach is that doesn't require overloading write_model_XXX_results!
# methods from PowerSimulations.
function write_results_to_main_container(container::MultiOptimizationContainer)
    # TODO: This process needs to work in parallel almost right away
    # TODO: This doesn't handle the case where subproblems have an overlap in axis names.
    for (k, subproblem) in container.subproblems
        for field in CONTAINER_FIELDS
            subproblem_data_field = getproperty(subproblem, field)
            main_container_data_field = getproperty(container, field)
            for (key, src) in subproblem_data_field
                if src isa JuMP.Containers.SparseAxisArray
                    @debug "Skip SparseAxisArray" field key
                    continue
                end
                num_dims = ndims(src)
                num_dims > 2 && error("ndims = $(num_dims) is not supported yet")
                data = nothing
                data = PSI.jump_value.(src)
                dst = main_container_data_field[key]
                if num_dims == 1
                    dst[1:length(axes(src)[1])] = data
                elseif num_dims == 2
                    #ychen fix horizontal passing ACbusinjection issue
                    if field == :expressions
                        field1 = :parameters
                        subproblem_data_field1 = getproperty(subproblem, field1)
                        src1 =
                            subproblem_data_field1[InfrastructureSystems.Optimization.ParameterKey{
                                PowerSimulationsDecomposition.StateEstimationInjections,
                                PSY.ACBus,
                            }(
                                "",
                            )]
                        B = parse.(Int, axes(src1.parameter_array)[1])
                        A = axes(src)[1]
                        C = filter(x -> !(x in B), A)
                        columns = C
                    else
                        columns = axes(src)[1]
                    end
                    #ychen end                    
                    len = length(axes(src)[2])
                    dst[columns, 1:len] = PSI.jump_value.(src[:, :])
                    #try 
                    #   println("======111  dst,k,",dst[203, :],",subproblem,",k)
                    #catch e end   
                elseif num_dims == 3
                    # TODO: untested
                    axis1 = axes(src)[1]
                    axis2 = axes(src)[2]
                    len = length(axes(src)[3])
                    dst[axis1, axis2, 1:len] = PSI.jump_value.(src[:, :, :])
                end
            end
        end
        _write_parameter_results_to_main_container(container, subproblem)
    end
    return
end

function _write_parameter_results_to_main_container(
    container::MultiOptimizationContainer,
    subproblem,
)
    for (key, parameter_container) in subproblem.parameters
        num_dims = ndims(parameter_container.parameter_array)
        num_dims > 2 && error("ndims = $(num_dims) is not supported yet")
        src_param_data = PSI.jump_value.(parameter_container.parameter_array)
        src_mult_data = PSI.jump_value.(parameter_container.multiplier_array)
        dst_param_data = container.parameters[key].parameter_array
        dst_mult_data = container.parameters[key].multiplier_array
        #println("*****11,subproblem,",subproblem)
        #println("*****12,key,",key)
        #println("*****13,num_dims,",num_dims)
        if num_dims == 1
            dst_param_data[1:length(axes(src_param_data)[1])] = src_param_data
            dst_mult_data[1:length(axes(src_mult_data)[1])] = src_mult_data

        elseif num_dims == 2
            param_columns = axes(src_param_data)[1]
            mult_columns = axes(src_mult_data)[1]
            len = length(axes(src_param_data)[2])
            @assert_op len == length(axes(src_mult_data)[2])
            dst_param_data[param_columns, 1:len] = PSI.jump_value.(src_param_data[:, :])
            dst_mult_data[mult_columns, 1:len] = PSI.jump_value.(src_mult_data[:, :])
        else
            error("Bug")
        end
    end
end

function solve_impl!(
    container::MultiOptimizationContainer{SequentialAlgorithm},
    sys::PSY.System,
)
    # Solve main problem
    status = ISSIM.RunStatus.RUNNING

    cost_curve = Dict()
    cost_upper = Dict()
    cost_lower = Dict()
    relief_price = Dict()
    adjust_rate = 0.1
    monitor_upper_key = PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.RateLimitConstraint, PSY.MonitoredLine}("ub")
    monitor_lower_key = PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.RateLimitConstraint, PSY.MonitoredLine}("lb")
    for (index, subproblem) in container.subproblems
        @debug "Solving problem $index"
        status = PSI.solve_impl!(subproblem, sys)
        println("yc -- solving problem,", index," ",status)
        # if status != ISSIM.RunStatus.SUCCESSFULLY_FINALIZED
        #     return status
        # end
        # if status != ISSIM.RunStatus.SUCCESSFULLY_FINALIZED
        #     println("not fail yet")
        #     line_capacity = 0.1
        #     monitored_lines = subproblem.constraints[monitor_lower_key].axes[1]
        #     mon_line_name = subproblem.constraints[monitor_lower_key].axes[1][1]
        #     for l in monitored_lines
        #         monitor_upper_con = subproblem.constraints[monitor_upper_key][l,1]
        #         # monitor_upper_bound = JuMP.normalized_rhs(monitor_upper_con)  
        #         monitor_lower_con = subproblem.constraints[monitor_lower_key][l,1]
        #         # monitor_lower_bound = JuMP.normalized_rhs(monitor_lower_con)  
        #         # JuMP.set_normalized_rhs(monitor_upper_con,line_capacity)
        #         # JuMP.set_normalized_rhs(monitor_lower_con,-line_capacity) 
        #         # if JuMP.normalized_rhs(monitor_lower_con) == line_capacity
        #         #     JuMP.set_normalized_rhs(monitor_upper_con,-line_capacity)
        #         #     JuMP.set_normalized_rhs(monitor_lower_con,-Inf)  
        #         # else
        #         #     JuMP.set_normalized_rhs(monitor_upper_con,Inf)
        #         #     JuMP.set_normalized_rhs(monitor_lower_con,line_capacity) 
        #         # end
        #     end
        # end
        status = PSI.solve_impl!(subproblem, sys)
        if status != ISSIM.RunStatus.SUCCESSFULLY_FINALIZED
            println("indeed failure")
            return status
        end
        
        ########ch compute the dual price
        if container.if_coordination > 0 
            if haskey(subproblem.constraints, monitor_lower_key) && haskey(subproblem.constraints, monitor_upper_key) && length(PSI.get_time_steps(subproblem)) == 1
                monitored_lines = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb").axes[1]
                if length(monitored_lines) == 1
                    mon_line_name = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb").axes[1][1]

                    monitor_upper_con = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"ub")[mon_line_name,1]
                    monitor_lower_con = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb")[mon_line_name,1]
                    
                    if status == ISSIM.RunStatus.SUCCESSFULLY_FINALIZED
                        
                        relief_price[index] = JuMP.dual(PSI.get_constraint(subproblem,PSI.NetworkFlowConstraint(),PSY.MonitoredLine,"")[mon_line_name,1])
                        cost_upper[index] = -JuMP.dual(monitor_upper_con)
                        cost_lower[index] = -JuMP.dual(monitor_lower_con)
                        # cost_upper[index] = -JuMP.shadow_price(monitor_upper_con)
                    else
                        relief_price[index] = 10000000
                        cost_upper[index] = 10000000
                        cost_lower[index] = 10000000
                    end
                else
                    # monitored_lines_index = Dict(line => index for (index, line) in enumerate(monitored_lines))
                    for l in monitored_lines

                        monitor_upper_con = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"ub")[l,1]
                        monitor_lower_con = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb")[l,1]
                        

                        if status == ISSIM.RunStatus.SUCCESSFULLY_FINALIZED
                            relief_price[(index, l)] = JuMP.dual(PSI.get_constraint(subproblem,PSI.NetworkFlowConstraint(),PSY.MonitoredLine,"")[l,1])
                            cost_upper[(index, l)] = -JuMP.dual(monitor_upper_con)
                            cost_lower[(index, l)] = -JuMP.dual(monitor_lower_con)
                            # cost_upper[index] = -JuMP.shadow_price(monitor_upper_con)
                        else
                            relief_price[(index, l)] = 10000000
                            cost_upper[(index, l)] = 10000000
                            cost_lower[(index, l)] = 10000000
                        end    
                    end
                end     
            end  
        end
    end
    write_results_to_main_container(container)
    ########ch set relief
    if container.if_coordination == 1 
        for (index, subproblem) in container.subproblems
            if haskey(subproblem.constraints, monitor_lower_key) && haskey(subproblem.constraints, monitor_upper_key) && length(PSI.get_time_steps(subproblem)) == 1
                
                ##### for multiple lines######
                monitored_lines = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb").axes[1]
                monitored_lines_index = Dict(line => index for (index, line) in enumerate(monitored_lines))
                for l in monitored_lines
                    monitor_upper_con = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"ub")[l,1]
                    monitor_lower_con = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb")[l,1]
                    monitor_upper_bound = JuMP.normalized_rhs(monitor_upper_con)  
                    monitor_lower_bound = JuMP.normalized_rhs(monitor_lower_con)  
                    current_flow = JuMP.value(subproblem.variables[PSY.InfrastructureSystems.Optimization.VariableKey{PSI.FlowActivePowerVariable, PSY.MonitoredLine}("")][l,1])
                    relief_name = "$(l)_$(1)"
                    relief = subproblem.variables[PSY.InfrastructureSystems.Optimization.VariableKey{PowerSimulationsDecomposition.Relief, PSY.MonitoredLine}("")][relief_name,1]
                    current_flow = current_flow - JuMP.value(relief)
                    other_region_relief_cost_upper = 0
                    relief_cost_upper = cost_upper[(index, l)]
                    for (key,value) in cost_upper
                        if key[1] != index
                            other_region_relief_cost_upper = value
                        end
                    end
                    if current_flow <= monitor_upper_bound + 0.0001 && current_flow >= monitor_upper_bound - 0.0001
                        if relief_cost_upper > other_region_relief_cost_upper
                            JuMP.set_upper_bound(relief, 0)
                            JuMP.set_lower_bound(relief, -adjust_rate*monitor_upper_bound)
                        end
                        if relief_cost_upper < other_region_relief_cost_upper
                            JuMP.set_upper_bound(relief, adjust_rate*monitor_upper_bound)
                            JuMP.set_lower_bound(relief, 0)
                        end
                    end

                    if current_flow >= monitor_upper_bound + 0.0001
                        if relief_cost_upper > other_region_relief_cost_upper
                            JuMP.set_upper_bound(relief, 0.000001)
                            JuMP.set_lower_bound(relief, 0)
                        end
                        if relief_cost_upper < other_region_relief_cost_upper
                            JuMP.set_upper_bound(relief, current_flow - monitor_upper_bound)
                            JuMP.set_lower_bound(relief, 0)
                        end
                    end
    
                    if current_flow <= monitor_upper_bound - 0.0001 && current_flow >= monitor_lower_bound + 0.0001
                        if relief_cost_upper < other_region_relief_cost_upper
                            JuMP.set_upper_bound(relief, 0)
                            JuMP.set_lower_bound(relief, current_flow - monitor_upper_bound)              
                        end
                        if relief_cost_upper > other_region_relief_cost_upper
                            JuMP.set_upper_bound(relief, 0)
                            JuMP.set_lower_bound(relief, -0.000001)
                        end
                        if relief_cost_upper == other_region_relief_cost_upper
                            JuMP.set_upper_bound(relief, 0)
                            JuMP.set_lower_bound(relief, 0.5*(current_flow - monitor_upper_bound))
                        end
                    end

                    #######set the relief based on the lower bound cost###
                    other_region_relief_cost_lower = 0
                    relief_cost_lower = cost_lower[(index, l)]
                    for (key,value) in cost_lower
                        if key[1] != index
                            other_region_relief_cost_lower = value
                        end
                    end

                    if current_flow <= monitor_lower_bound + 0.0001 && current_flow >= monitor_lower_bound - 0.0001
                        if -relief_cost_lower < -other_region_relief_cost_lower
                            JuMP.set_upper_bound(relief, -adjust_rate*monitor_lower_bound)
                            JuMP.set_lower_bound(relief, 0)
                        end
                        if -relief_cost_lower > -other_region_relief_cost_lower
                            JuMP.set_upper_bound(relief, 0)
                            JuMP.set_lower_bound(relief, adjust_rate*monitor_lower_bound)
                        end
                    end
    
                    if current_flow <= monitor_lower_bound - 0.0001
                        if -relief_cost_lower < -other_region_relief_cost_lower
                            JuMP.set_upper_bound(relief, 0.000001)
                            JuMP.set_lower_bound(relief, 0)
                        end
                        if -relief_cost_lower > -other_region_relief_cost_lower
                            JuMP.set_upper_bound(relief, 0)
                            JuMP.set_lower_bound(relief, current_flow - monitor_lower_bound)
                        end
                    end
    
                    if current_flow <= monitor_upper_bound - 0.0001 && current_flow >= monitor_lower_bound + 0.0001
                        if relief_cost_lower < other_region_relief_cost_lower
                            JuMP.set_upper_bound(relief, 0)
                            JuMP.set_lower_bound(relief, current_flow - monitor_upper_bound)              
                        end
                        if relief_cost_lower > other_region_relief_cost_lower
                            JuMP.set_upper_bound(relief, 0)
                            JuMP.set_lower_bound(relief, -0.0001)
                        end
                        if relief_cost_lower == other_region_relief_cost_lower
                            JuMP.set_upper_bound(relief, 0)
                            JuMP.set_lower_bound(relief, 0.5*(current_flow - monitor_upper_bound))
                        end
                    end

                    if current_flow > monitor_upper_bound-0.0001
                        JuMP.set_objective_coefficient(subproblem.JuMPmodel, relief, other_region_relief_cost_upper)
                    elseif current_flow < monitor_lower_bound+0.0001
                        JuMP.set_objective_coefficient(subproblem.JuMPmodel, relief, other_region_relief_cost_lower)
                    else
                        JuMP.set_objective_coefficient(subproblem.JuMPmodel, relief, 0)
                    end

                end
                ##############################
            end       
        end
    end

    if container.if_coordination == 2 
        #############################ch get the cost curve of relief####   
        println("generate curve")  
        current_flow = Dict()   
        for (index, subproblem) in container.subproblems
            if (haskey(subproblem.constraints, monitor_lower_key) || haskey(subproblem.constraints, monitor_upper_key)) && length(PSI.get_time_steps(subproblem)) == 1
                println("computing ",index)
                
                cost_points = Dict()
                line_capacity = 0.1
                # marginal_cost_range = line_capacity
                trans_limit_change = collect(-line_capacity:line_capacity/7:line_capacity)
                segments = collect(1:15)
                monitored_lines = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb").axes[1]
                for l in monitored_lines
                    monitor_upper_con = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"ub")[l,1]
                    monitor_lower_con = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb")[l,1]
                    current_flow = JuMP.value(subproblem.variables[PSY.InfrastructureSystems.Optimization.VariableKey{PSI.FlowActivePowerVariable, PSY.MonitoredLine}("")][l,1])
                    println("flow ", current_flow)
                    
                    relief = subproblem.variables[PSY.InfrastructureSystems.Optimization.VariableKey{PowerSimulationsDecomposition.Relief, PSY.MonitoredLine}("")]
                    current_total_flow = current_flow - sum(JuMP.value(relief["$(l)_$(s)_$(i)", 1]) for s in segments for i in [1, 2])

                    
                    for s in segments
                        for i in [1,2]
                            relief_name = "$(l)_$(s)_$(i)"
                            var = relief[relief_name,1]

                            println("relief $l $s $i ", JuMP.value(var)," cost ", JuMP.coefficient(JuMP.objective_function(PSI.get_jump_model(subproblem)), var))
                            JuMP.set_normalized_coefficient(monitor_lower_con, var, 0)
                            JuMP.set_normalized_coefficient(monitor_upper_con, var, 0)
                        end
                    end
                    JuMP.set_normalized_rhs(monitor_upper_con,line_capacity)
                    for diff in trans_limit_change
                        JuMP.set_normalized_rhs(monitor_upper_con, current_flow - diff) 
                        status = PSI.solve_impl!(subproblem, sys)
                        if status == ISSIM.RunStatus.SUCCESSFULLY_FINALIZED
                            cost_points[diff] = -JuMP.dual(monitor_upper_con)
                        else
                            cost_points[diff] = 1000000
                        end        
                        # JuMP.set_normalized_rhs(monitor_upper_con,line_capacity)  
                    end
                    # if current_flow>0
                    #     # JuMP.set_normalized_rhs(subproblem.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.NetworkFlowConstraint, PSY.MonitoredLine}("")][mon_line_name,1], line_capacity)
                    #     for diff in trans_limit_change
                    #         JuMP.set_normalized_rhs(monitor_upper_con, current_flow - diff) 
                    #         JuMP.set_normalized_rhs(monitor_lower_con, -Inf)  
                    #         status = PSI.solve_impl!(subproblem, sys)
                    #         if status == ISSIM.RunStatus.SUCCESSFULLY_FINALIZED
                    #             cost_points[diff] = -JuMP.dual(monitor_upper_con)
                    #         else
                    #             cost_points[diff] = 1000000
                    #         end        
                    #         # JuMP.set_normalized_rhs(monitor_upper_con,line_capacity)  
                    #     end
                    # else
                    #     # JuMP.set_normalized_rhs(subproblem.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.NetworkFlowConstraint, PSY.MonitoredLine}("")][mon_line_name,1], -line_capacity)
                    #     for diff in trans_limit_change
                    #         JuMP.set_normalized_rhs(monitor_upper_con, Inf)
                    #         JuMP.set_normalized_rhs(monitor_lower_con, current_flow + diff)  
                    #         status = PSI.solve_impl!(subproblem, sys)
                    #         if status == ISSIM.RunStatus.SUCCESSFULLY_FINALIZED
                    #             cost_points[diff] = -JuMP.dual(monitor_lower_con)
                    #         else
                    #             cost_points[diff] = -1000000
                    #         end
                    #         # JuMP.set_normalized_rhs(monitor_lower_con,-line_capacity)          
                    #     end
                    # end
                    unique_values = Dict()
                    epsilon = 1e-5
                    normalize_value = x -> abs(x) < epsilon ? 0.0 : x
                    sorted_cost_points = sort(collect(cost_points), by = first)
                    println("sorted_cost_points ", sorted_cost_points)
                    sub_cost_curve = []
                    for (key, value) in sorted_cost_points
                        norm_value = normalize_value(value)
                        if !haskey(unique_values, norm_value) && length(sub_cost_curve) < 15
                            unique_values[norm_value] = true  
                            push!(sub_cost_curve, (key, value))  
                        end
                    end
                    last_key, last_value = last(sorted_cost_points)
                    if !any(p -> p[1] == last_key && p[2] == last_value, sub_cost_curve)
                        push!(sub_cost_curve, (last_key, last_value))
                    end
                    cost_curve[(l,index)] = sub_cost_curve

                    for s in segments
                        for i in [1,2]
                            relief_name = "$(l)_$(s)_$(i)"
                            var = relief[relief_name,1]
                            JuMP.set_normalized_coefficient(monitor_lower_con, var, -1)
                            JuMP.set_normalized_coefficient(monitor_upper_con, var, -1)
                        end
                    end
                    JuMP.set_normalized_rhs(monitor_upper_con,-line_capacity)

                    # if current_flow>0                 
                    #     # JuMP.set_normalized_rhs(monitor_upper_con,0)
                    #     JuMP.set_normalized_rhs(monitor_upper_con,-line_capacity)
                    #     # JuMP.set_normalized_rhs(subproblem.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.NetworkFlowConstraint, PSY.MonitoredLine}("")][mon_line_name,1], line_capacity)
                    # else
                    #     # JuMP.set_normalized_rhs(monitor_lower_con,0)
                    #     JuMP.set_normalized_rhs(monitor_lower_con,line_capacity)
                    #     # JuMP.set_normalized_rhs(subproblem.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.NetworkFlowConstraint, PSY.MonitoredLine}("")][mon_line_name,1], -line_capacity)
                    # end
                    # if current_total_flow + 2*line_capacity> 0 
                    # # if current_flow>0                
                    #     # JuMP.set_normalized_rhs(monitor_upper_con,0)
                    #     JuMP.set_normalized_rhs(monitor_upper_con,-line_capacity)
                    #     JuMP.set_normalized_rhs(monitor_lower_con,-3*line_capacity)
                    #     # JuMP.set_normalized_rhs(subproblem.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.NetworkFlowConstraint, PSY.MonitoredLine}("")][mon_line_name,1], line_capacity)
                    # else
                    #     # JuMP.set_normalized_rhs(monitor_lower_con,0)
                    #     JuMP.set_normalized_rhs(monitor_upper_con,3*line_capacity)
                    #     JuMP.set_normalized_rhs(monitor_lower_con,line_capacity)
                    #     # JuMP.set_normalized_rhs(subproblem.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.NetworkFlowConstraint, PSY.MonitoredLine}("")][mon_line_name,1], -line_capacity)
                    # end

                end
                println("cost curve ",cost_curve)
            end
        end

        #### set the bound
        for (index, subproblem) in container.subproblems
            if (haskey(subproblem.constraints, monitor_lower_key) || haskey(subproblem.constraints, monitor_upper_key)) && length(PSI.get_time_steps(subproblem)) == 1
                println("setting bound ",index)
                segments = collect(1:15)
                line_capacity = 0.1
                monitored_lines = PSI.get_constraint(subproblem,PSI.RateLimitConstraint(),PSY.MonitoredLine,"lb").axes[1]
                relief = subproblem.variables[PSY.InfrastructureSystems.Optimization.VariableKey{PowerSimulationsDecomposition.Relief, PSY.MonitoredLine}("")]
                for l in monitored_lines

                    #consider the case with multiple regions and lines
                    sub_cost_curve = []
                    for (key,value) in cost_curve
                        if key[2] != index && key[1] == l
                            push!(sub_cost_curve, value)
                        end
                    end
                    println("sub_cost_curve ",sub_cost_curve)

                    
                    current_flow = JuMP.value(subproblem.variables[PSY.InfrastructureSystems.Optimization.VariableKey{PSI.FlowActivePowerVariable, PSY.MonitoredLine}("")][l,1])
                    i = 1
                    for cost_curve_value in sub_cost_curve
                        println("cost_curve_value", cost_curve_value)
                        num_pieces = length(cost_curve_value)
                        for s in segments      
                            relief_name = "$(l)_$(s)_$(i)"
                            var = relief[relief_name,1]
                            if s < num_pieces
                                JuMP.set_upper_bound(var, cost_curve_value[s+1][1] - cost_curve_value[s][1])
                            else
                                JuMP.set_upper_bound(var, 0)
                            end
                            if s == 1 && num_pieces == 1
                                JuMP.set_upper_bound(var, line_capacity)
                            end
                            coeff = s <= num_pieces ? cost_curve_value[s][2] : 1000000  

                            JuMP.set_objective_coefficient(subproblem.JuMPmodel, var, coeff)
                        end
                        i += 1
                    end
                end
            end

        end
 
    end
    return status
end
