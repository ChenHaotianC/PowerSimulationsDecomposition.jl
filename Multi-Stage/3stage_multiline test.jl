using Pkg
#Pkg.activate(".")
Pkg.activate(@__DIR__)
#run(`bash -c "module load gurobi"`)

#ENV["GUROBI_HOME"] = "/nopt/nrel/apps/software/gurobi/gurobi1100/linux64"
#ENV["PATH"] *= ":$ENV[\"GUROBI_HOME\"]/bin"
#ENV["GRB_LICENSE_FILE"] = "/nopt/nrel/apps/software/gurobi/tlicense/gurobi.lic"  

#ENV["XPRESSDIR"] = "C:\\xpressmp"
#ENV["XPAUTH_PATH"] = "C:\\xpressmp\\bin"

#ENV["XPRESSDIR"] = "/nopt/nrel/apps/software/xpressmp/"
#ENV["XPAUTH_PATH"] ="/nopt/nrel/apps/software/xpressmp/9.2.5/bin"
using PowerSystems
using PowerSimulations
using PowerSystemCaseBuilder
using PowerNetworkMatrices
using HydroPowerSimulations
using StorageSystemsSimulations
using JuMP
using Dates
#using Gurobi
using PowerSimulationsDecomposition
using Revise
using Logging
using CSV
using DataFrames
using Xpress
using HiGHS
using InfrastructureSystems

const HPS = HydroPowerSimulations
const PSI = PowerSimulations
const SSS = StorageSystemsSimulations
const PSY = PowerSystems

function write_cons_coeff(uc2a,con2a,fnm)
    open(fnm,"w") do io
        redirect_stdout(io) do
            for k in all_variables(uc2a.JuMPmodel)  
                if normalized_coefficient(con2a, k)!=0         
                    println("coef;",name(k),";",normalized_coefficient(con2a, k), ";is_fixed;",is_fixed(k),";value;",value(k))
                end  
            end 
        end
    end
end

function write_model(uc0,fnm)
    open(fnm,"w") do io
        redirect_stdout(io) do
            println(objective_function(uc0.JuMPmodel))
            for k in all_constraints(uc0.JuMPmodel,; include_variable_in_set_constraints = true)           
                println(name(k),",",k) 
            end    
        end    
    end
end

function objterms(uc0,ff)
    objd0=Dict(); sol0=Dict(); objv0=0
    obj0=objective_function(uc0.JuMPmodel); objd0=Dict(zip(name.(keys(obj0.terms)),values(obj0.terms))); 
    av0 = JuMP.all_variables(uc0.JuMPmodel); sol0=Dict(zip(name.(av0), value.(av0)))

    open(ff,"w") do io
        redirect_stdout(io) do
            for (k,v) in objd0
                objv0=objv0+sol0[k]*objd0[k]
                println("name;value;objcoef;obj_contribution;",k,";",sol0[k],";",objd0[k],";", sol0[k]*objd0[k])      
            end
        end
    end  
    return(objd0,sol0, objv0)  
end    

function read_bus_df(results_rt,se=0,standardload=0)
    ac = read_realized_variable(results_rt, "ActivePowerVariable__ThermalStandard")
    ac1 = read_realized_variable(results_rt, "ActivePowerVariable__RenewableDispatch")
    ac2 = read_realized_variable(results_rt, "ActivePowerVariable__HydroDispatch")

    ab = read_realized_variable(results_rt, "ActivePowerBalance__ACBus")
    if se==1  al = read_realized_variable(results_rt, "StateEstimationInjections__ACBus") end
    if standardload==1
        ald = read_realized_variable(results_rt, "ActivePowerTimeSeriesParameter__StandardLoad")
    else
        ald = read_realized_variable(results_rt, "ActivePowerTimeSeriesParameter__PowerLoad")
    end        
    #if length(line)>0 am = read_realized_variable(results_rt, "FlowActivePowerVariable__MonitoredLine") end
    ahvdc=DataFrame()
    try 
        ahvdc=read_realized_variable(results_rt, "FlowActivePowerVariable__TwoTerminalHVDCLine")
    catch e end

    bus_ActivePowerVariable__ThermalStandard=Dict()
    bus_ActivePowerVariable__RenewableDispatch=Dict()
    bus_ActivePowerVariable__HydroDispatch=Dict()
    bus_ActivePowerTimeSeriesParameter__PowerLoad=Dict()
    bus_ActivePowerBalance__ACBus=Dict()
    bus_StateEstimationInjections__ACBus=Dict()
    bus_ActivePowerVariable__HVDC=Dict()

    bus_df=DataFrames.DataFrame(;
    t =Int64[], bus=[], area = [], ThermalStandard=[],Renewable=[],Hydro=[],
    PowerLoad=[], HVDC=[],ActivePowerBalance__ACBus=[], StateEstimationInjections__ACBus=[]
    #,ptdf=[], flowcontribution=[], loopflowcontribution=[]
    )
    NT=size(ac,1)
    for t in 1:NT
        for b in PSY.get_components(PSY.Bus,sys)
            bn=get_number(b); bus_ActivePowerVariable__ThermalStandard[t,bn]=0
            bus_ActivePowerVariable__RenewableDispatch[t,bn]=0
            bus_ActivePowerVariable__HydroDispatch[t,bn]=0
            bus_ActivePowerVariable__HVDC[t,bn]=0
            bus_ActivePowerTimeSeriesParameter__PowerLoad[t,bn]=0; bus_ActivePowerBalance__ACBus[t,bn]=0
            bus_StateEstimationInjections__ACBus[t,bn]=0
        end
    end

    for i in names(ab)
        try b=parse(Int,i)
            for t=1:NT bus_ActivePowerBalance__ACBus[t,b]=ab[t,i] end
        catch e println("error= ",e,",i=",i)  end        
    end
    if se==1
        for i in names(al)
            try b=parse(Int,i)
                for t=1:NT bus_StateEstimationInjections__ACBus[t,b]=al[t,i] end
            catch e println("error= ",e,",i=",i) end     
        end
    end    

    for i in names(ald)
        if standardload==1
            try b=get_number(get_bus(get_component(StandardLoad, sys, i)))
                for t=1:NT  bus_ActivePowerTimeSeriesParameter__PowerLoad[t,b]=bus_ActivePowerTimeSeriesParameter__PowerLoad[t,b]+ald[t,i] end  
            catch e println("error= ",e,",i=",i) end 
        else
            try b=get_number(get_bus(get_component(PowerLoad, sys, i)))
                for t=1:NT  bus_ActivePowerTimeSeriesParameter__PowerLoad[t,b]=bus_ActivePowerTimeSeriesParameter__PowerLoad[t,b]+ald[t,i] end  
            catch e println("error= ",e,",i=",i) end 
        end                    
    end

    if size(ahvdc,1)>0
        for i in names(ahvdc)
            try b=get_number(get_from(get_arc(get_component(TwoTerminalHVDCLine, sys, i))))
            for t=1:NT 
                #bus_ActivePowerVariable__ThermalStandard[t,b]=bus_ActivePowerVariable__ThermalStandard[t,b]-ahvdc[t,i] 
                bus_ActivePowerVariable__HVDC[t,b]=bus_ActivePowerVariable__HVDC[t,b]-ahvdc[t,i] 
            end  
            catch e println("error ",e,",i=",i) end 
            try b=get_number(get_to(get_arc(get_component(TwoTerminalHVDCLine, sys, i))))
                for t=1:NT 
                    #bus_ActivePowerVariable__ThermalStandard[t,b]=bus_ActivePowerVariable__ThermalStandard[t,b]+ahvdc[t,i] 
                    bus_ActivePowerVariable__HVDC[t,b]=bus_ActivePowerVariable__HVDC[t,b]+ahvdc[t,i] 
                end  
            catch e println("error ",e,",i=",i) end 
        end
    end    
    
    for i in names(ac)
        try b=get_number(get_bus(get_component(ThermalStandard, sys, i)))
        for t=1:NT bus_ActivePowerVariable__ThermalStandard[t,b]=bus_ActivePowerVariable__ThermalStandard[t,b]+ac[t,i] end  
        catch e println("error ",e,",i=",i) end 
    end
    for i in names(ac1)
        try b=get_number(get_bus(get_component(RenewableDispatch, sys, i)))
        for t=1:NT bus_ActivePowerVariable__RenewableDispatch[t,b]=bus_ActivePowerVariable__RenewableDispatch[t,b]+ac1[t,i] end  
        catch e println("error ",e,",i=",i) end 
    end
    for i in names(ac2)
        try b=get_number(get_bus(get_component(HydroDispatch, sys, i)))
        for t=1:NT bus_ActivePowerVariable__HydroDispatch[t,b]=bus_ActivePowerVariable__HydroDispatch[t,b]+ac2[t,i] end  
        catch e println("error ",e,",i=",i) end 
    end

    for t in 1:NT
        for b in PSY.get_components(PSY.Bus,sys)
            bn=get_number(b); area=get_name(get_area(b)); 
            #gsf=Dict(); facbus=Dict(); facbus_se=Dict()
            #for ln in line
            #    gsf[ln]=ptdf[ln,bn]
            #    facbus[ln]=bus_ActivePowerBalance__ACBus[t,bn]*gsf[ln]
            #    facbus_se[ln]=bus_StateEstimationInjections__ACBus[t,bn]*gsf[ln]
            #end 
            push!(bus_df,(t,bn,area, bus_ActivePowerVariable__ThermalStandard[t,bn], 
            bus_ActivePowerVariable__RenewableDispatch[t,bn],bus_ActivePowerVariable__HydroDispatch[t,bn],
            bus_ActivePowerTimeSeriesParameter__PowerLoad[t,bn],bus_ActivePowerVariable__HVDC[t,bn],
            bus_ActivePowerBalance__ACBus[t,bn], bus_StateEstimationInjections__ACBus[t,bn] #, gsf,
            ##bus_ActivePowerBalance__ACBus[t,bn]*gsf,bus_StateEstimationInjections__ACBus[t,bn]*gsf
            #facbus,facbus_se
            ))
        end
    end

    bus_df[!,"buscheck"]=bus_df[!,"ActivePowerBalance__ACBus"]-bus_df[!,"ThermalStandard"]/100-bus_df[!,"Renewable"]/100-bus_df[!,"Hydro"]/100-bus_df[!,"PowerLoad"]/100
    println("bus check ActivePowerBalance__ACBus<>ThermalStandard/100-PowerLoad/100 ,",filter([:t,:buscheck] => (t,buscheck)-> (t>=1) && (abs(buscheck)>0.000001),  bus_df))
    #println(filter([:t,:bus] => (t,bus)-> (t>=1) && (bus==203),  bus_df))
    return(bus_df)
end    

function buildsubsystem1(sys, area_region_dict, selected_line, monitoredline_subsys=[])
    subsys=unique(values(area_region_dict))

    for b in PSY.get_components(PSY.AreaInterchange, sys)
        println(b); 
        for v in subsys
            add_component_to_subsystem!(sys, v, b) 
        end
    end

    for b in PSY.get_components(PSY.Area, sys)
        region=area_region_dict[get_name(b)]
        PSY.set_ext!(b, Dict("subregion" => Set([region])))
        add_component_to_subsystem!(sys, region, b)
    end
    for b in PSY.get_components(PSY.StaticInjection, sys)
        region=area_region_dict[get_name(get_area(get_bus(b)))]
        PSY.set_ext!(b, Dict("subregion" => Set([region])))
        add_component_to_subsystem!(sys, region, b)
    end
    for b in PSY.get_components(PSY.Bus, sys)
        region=area_region_dict[get_name(get_area(b))]
        PSY.set_ext!(b, Dict("subregion" => Set([region])))
        add_component_to_subsystem!(sys, region, b)
    end
    
    if length(selected_line)>0
        for b in selected_line 
            l=PSY.get_component(ACBranch, sys,b)
            for v in subsys
                PSY.set_ext!(l, Dict("subregion" => Set([v])))
                add_component_to_subsystem!(sys, v, l) 
            end            
        end  
    end    
    
    if length(monitoredline_subsys) >0
        for (k,s) in monitoredline_subsys 
            l=PSY.get_component(ACBranch, sys,k)
            for v in s
                PSY.set_ext!(l, Dict("subregion" => Set([v])))
                add_component_to_subsystem!(sys, v, l) 
            end            
        end  
    end

    for dc in get_components(TwoTerminalHVDCLine, sys)
        tbus=get_to(get_arc(dc)); fbus=get_from(get_arc(dc))
        tregion=area_region_dict[get_name(get_area(tbus))]
        fregion=area_region_dict[get_name(get_area(fbus))]
        PSY.set_ext!(dc, Dict("subregion" => Set([tregion])))
        add_component_to_subsystem!(sys, tregion, dc) 
        
        #Option 1, add HVDC also to fregion. Can cause double counting of fixed and dispatchable HVDC MW
        #if fregion!=tregion
        #    PSY.set_ext!(dc, Dict("subregion" => Set([fregion])))
        #    add_component_to_subsystem!(sys, fregion, dc)
        #end                

        #Option 2, add fbus to tregion
        #if fregion!=tregion
        #    PSY.set_ext!(fbus, Dict("subregion" => Set([tregion])))
        #    remove_component_from_subsystem!(sys, fregion, fbus)
        #    add_component_to_subsystem!(sys, tregion, fbus)
        #end           
    end    
end


# EI
# include("final_EI_model.jl")

#RTS
include("final_RTS.jl")

mipgap = 0.001
#=
optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "Threads" => (length(Sys.cpu_info()) ÷ 2) - 1,
    "MIPGap" => mipgap,
    "TimeLimit" => 3000,
)
=#

optimizer = optimizer_with_attributes(
    Xpress.Optimizer,
    #"Threads" => (length(Sys.cpu_info()) ÷ 2) - 1,
    #"MIPGap" => mipgap,
    #"TimeLimit" => 3000,
)


#optimizer = optimizer_with_attributes(
#    HiGHS.Optimizer,
#)

ptdf = VirtualPTDF(sys; tol = 1e-4, max_cache_size = 10000)        
#monitoredlined_line=[ln]

for b in monitoredlined_line
    l = PSY.get_component(ACBranch, sys, b)
    set_rating!(l,limit[b])
    #PSY.get_component(ACBranch, sys, b)
    line = PSY.get_component(Line, sys, b)
    PSY.convert_component!(sys, line, MonitoredLine)
end

r = run_df[3, :]
CF,NS,netmodels,copperplates,NTS,area_region_dicts,monitoredline_subsys=r[["CF","NS","netmodels","copperplates","NTS","area_region_dicts","monitoredline_subsys"]]
syss=[sys]
templates=[]
NS
for i in 1:NS
    NT=NTS[i]
    if i>1                               
        sys2=deepcopy(sys)
        push!(syss,sys2)
    end
    transform_single_time_series!(syss[i], Hour(NT), Hour(NT))

    area_region_dict=area_region_dicts[i]  
    if netmodels[i]=="SplitAreaPTDFPowerModel"     
        subsys=unique(values(area_region_dict))
        for v in subsys
            println("add sys",i," subsystem,",v)
            add_subsystem!(sys2, v)
        end
        template = 
        MultiProblemTemplate(NetworkModel(SplitAreaPTDFPowerModel; use_slacks=true,PTDF_matrix = ptdf,), subsys)
    else          
        template = ProblemTemplate(NetworkModel(AreaPTDFPowerModel; use_slacks = true, PTDF_matrix = ptdf,))
    end

    push!(templates,template)    

    if standardload==1
        PSI.set_device_model!(template, StandardLoad , StaticPowerLoad)
    else
        PSI.set_device_model!(template, PowerLoad, StaticPowerLoad)
    end

    PSI.set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
    PSI.set_device_model!(template, DeviceModel(HydroDispatch, HPS.HydroDispatchReservoirBudget,
                                    time_series_names = Dict{Any, String}(
                                        PSI.ActivePowerTimeSeriesParameter => "max_active_power",
                                        HPS.EnergyBudgetTimeSeriesParameter => "hydro_budget",)
                                            )
                    )
    if hvdc==1
        PSI.set_device_model!(template, TwoTerminalHVDCLine, HVDCTwoTerminalLossless)
    end

    set_device_model!(template, DeviceModel(MonitoredLine, StaticBranchUnbounded, use_slacks = true))
    
    if NT>1              
        PSI.set_device_model!(template, ThermalStandard, ThermalBasicUnitCommitment)
        PSI.set_device_model!(template, ThermalMultiStart, ThermalBasicUnitCommitment)
        #link_nsi_to_tielines. Need to set selected_line with tie lines
        #PSI.set_device_model!(template, AreaInterchange, StaticBranch)
        # not link_nsi_to_tielines. Need to activate HVDC        
        #PSI.set_device_model!(template, AreaInterchange, StaticBranch)
        PSI.set_device_model!(template, AreaInterchange, StaticBranchUnbounded)
    else
        PSI.set_device_model!(template, ThermalStandard, ThermalBasicDispatch)
        PSI.set_device_model!(template, AreaInterchange, StaticBranchUnbounded)
    end                  
    
    if copperplates[i]==0
        set_device_model!(template, DeviceModel(MonitoredLine, StaticBranch, use_slacks = true))
    end
                
    set_device_model!(template,
        DeviceModel(Line, StaticBranchUnbounded; attributes=Dict("filter_function" => x -> get_name(x) in selected_line),))             

    subsys=unique(values(area_region_dict))

    if length(subsys)>1
    #buildsubsystem1(sys2, area_region_dict, union(selected_line,monitoredlined_line))
    buildsubsystem1(sys2, area_region_dict, selected_line,monitoredline_subsys[i])
    end
end

dm=[]
for i in 1:NS
    if netmodels[i]=="SplitAreaPTDFPowerModel"
        push!(dm, DecisionModel(MultiRegionProblem, templates[i], syss[i], name="UC"*string(i),optimizer=optimizer,
            store_variable_names=true,initialize_model=false,optimizer_solve_log_print=false, 
            direct_mode_optimizer=true,check_numerical_bounds=false, 
            calculate_conflict=true,rebuild_model=false,system_to_file = true))
    else    
        push!(dm, DecisionModel(templates[i],syss[i], name="UC"*string(i), optimizer=optimizer,
        store_variable_names=true,initialize_model=false,optimizer_solve_log_print=false, 
        direct_mode_optimizer=true,check_numerical_bounds=false, 
        calculate_conflict=true,rebuild_model=false, system_to_file = true))
    end
end    


if NS==2
    models = PSI.SimulationModels(decision_models=[dm[1],dm[2]])
else    
    models = PSI.SimulationModels(decision_models=[dm[1],dm[2],dm[3]])
end    
uc_simulation_ffs=[]
feedforwards_dict=Dict()
for i in 1:NS-1
    uc_simulation_ff = Vector{PowerSimulations.AbstractAffectFeedforward}()
#feed forward area interchange
    FVFF_area_interchange = FixValueFeedforward(;component_type=AreaInterchange,source=FlowActivePowerVariable,affected_values=[FlowActivePowerVariable],)
    push!(uc_simulation_ff, FVFF_area_interchange)
    #feed forward HVDC (doesn't seem to work)
    #FVFF_hvdc = FixValueFeedforward(;component_type=TwoTerminalHVDCLine,source=FlowActivePowerVariable,affected_values=[FlowActivePowerVariable],)
    #push!(uc_simulation_ff, FVFF_hvdc)

    if NTS[i+1]==1
        SCFF = SemiContinuousFeedforward(;
        component_type=ThermalStandard, source=OnVariable, affected_values=[ActivePowerVariable],)
        push!(uc_simulation_ff, SCFF) 
        #LBFF = FixValueFeedforward(;component_type=ThermalStandard,source=OnVariable, affected_values=[OnVariable],)
        #push!(uc_simulation_ff, LBFF)
    end
    push!(uc_simulation_ffs,uc_simulation_ff)
    feedforwards_dict["UC"*string(i+1)]=uc_simulation_ff
end

# models.decision_models[3].internal.container.if_coordination = true

sequence = SimulationSequence(;
    models=models,
    feedforwards=feedforwards_dict, #Dict("UC2" => uc_simulation_ff[1], "UC3" => uc_simulation_ff[2],),
    ini_cond_chronology=InterProblemChronology(),
);

# Specify the simulation setup
# Here we specify the simulation name, the initial/start time, number of steps/days to execute, and the simulation folder.
sim = PSI.Simulation(
    name="ntps_3stageDA_DASplit_RTSplit",
    steps=1,
    models=models,
    sequence=sequence,
    simulation_folder=mktempdir(), #".", #
)


# println("before PSI.build")
PSI.build!(sim, serialize=false)
# println("after PSI.build")
t1 = time()
PSI.execute!(sim, enable_progress_bar=false)
t2 = time()
results = SimulationResults(sim,ignore_status = true)

# check results and calculate flows in different optimization models and the actual flow called calculated_seflow
#include("final_check_results.jl")  
results_uc=[]; uc_areainterchange=[];uc_monitoredline=[];
bus_df_uc=[];area_df=[];
tm=Array(1:NTS[1])
areaflow_df=Dict();af=Dict()
for ln in monitoredlined_line
    areaflow_df[ln]=[]; 
    af[ln]=DataFrames.DataFrame(;tm)
end    

for i in 1:NS
    push!(results_uc,get_decision_problem_results(results, "UC"*string(i)))
    push!(uc_areainterchange,read_realized_variable(results_uc[i], "FlowActivePowerVariable__AreaInterchange"))
    #push!(uc_monitoredline,read_realized_variable(results_uc[i], "FlowActivePowerVariable__MonitoredLine"))

    if netmodels[i]=="SplitAreaPTDFPowerModel"
        push!(bus_df_uc,read_bus_df(results_uc[i],1,standardload)) 
    else    
        push!(bus_df_uc,read_bus_df(results_uc[i],0,standardload))
    end    

    bus_df_uc[i][!,"buscheck"]=bus_df_uc[i][!,"ActivePowerBalance__ACBus"]-bus_df_uc[i][!,"HVDC"]/100-bus_df_uc[i][!,"ThermalStandard"]/100-bus_df_uc[i][!,"Renewable"]/100-bus_df_uc[i][!,"Hydro"]/100-bus_df_uc[i][!,"PowerLoad"]/100
    # println("UC0 bus check ActivePowerBalance__ACBus<>Gen/100-PowerLoad/100 ,",filter([:t,:buscheck] => (t,buscheck)-> (t>=1) && (abs(buscheck)>0.000001),  bus_df_uc[i]))
    # if some bus has non-zero buscheck value, then check the bus. Some MW may not be added to power balance correctly.
    #println(filter([:t,:bus] => (t,bus)-> (t==1) && (bus==113 || bus==316),  bus_df_uc0))

###### area MW check #########
#area level net Gen, load and ActivePowerBalance__ACBus_sum check
    push!(area_df,combine(groupby(bus_df_uc[i], [:t, :area]), [:ThermalStandard, :Renewable, :Hydro, :PowerLoad, :ActivePowerBalance__ACBus, :StateEstimationInjections__ACBus] .=> sum))
    area_df[i][!,"areacheck"]=area_df[i][!,"ActivePowerBalance__ACBus_sum"]-area_df[i][!,"ThermalStandard_sum"]/100-area_df[i][!,"Renewable_sum"]/100-area_df[i][!,"Hydro_sum"]/100-area_df[i][!,"PowerLoad_sum"]/100
    # println("i,area_df check,",sum(area_df[i][!,"areacheck"]))  #should be close to 0

    area_df[i][!,"TotalGen"]=area_df[i][!,"ThermalStandard_sum"]+area_df[i][!,"Renewable_sum"]+area_df[i][!,"Hydro_sum"]
    area_df[i][!,"NetInterchange"]=area_df[i][!,"TotalGen"]+area_df[i][!,"PowerLoad_sum"]
    area_df[i][!,"InterchangeCheck"]=area_df[i][!,"NetInterchange"]-100*area_df[i][!,"ActivePowerBalance__ACBus_sum"]
    # println(area_df[i][:,["t","area","ActivePowerBalance__ACBus_sum","NetInterchange","InterchangeCheck"]])

####### area flow check #########
#area flow contribution by area and interval
    bus_df_uc[i][!,:gsf].=0.0; bus_df_uc[i][!,:flowcontribution].=0.0; bus_df_uc[i][!,:loopflowcontribution].=0.0
    for ln in monitoredlined_line
        for b in eachrow(bus_df_uc[i])
           gsf=ptdf[ln,b[:bus]]; 
           b[:gsf]=gsf; b[:flowcontribution]=gsf*b[:ActivePowerBalance__ACBus]; b[:loopflowcontribution]=gsf*b[:StateEstimationInjections__ACBus]
        end

        push!(areaflow_df[ln], combine(groupby(bus_df_uc[i], [:t, :area]), [:flowcontribution, :loopflowcontribution] .=> sum))

        if length(area_region_dicts[i])==0
            area_region_dict=area_dict
        else   
            area_region_dict=area_region_dicts[i]
        end      
        
        subsys=unique(values(area_region_dict))

        calculated_flow_subsys=Dict()
        calculated_loopflow_subsys=Dict()
        NT=NTS[1]
        for s in subsys
            #calculated_flow_subsys[s]=zeros(NT)
            #calculated_loopflow_subsys[s]=zeros(NT)
            af[ln][!,"uc"*string(i)*"_flow_subsys_"*s]=zeros(NT)
            af[ln][!,"uc"*string(i)*"_loopflow_subsys_"*s]=zeros(NT)
        end    

        if netmodels[i]!="SplitAreaPTDFPowerModel"
            for r in keys(area_region_dict)
                s=area_region_dict[r]
                name="uc"*string(i)*"_flow_subsys_"*s
                af[ln][!,name]=af[ln][!,name]+filter([:t,:area]=>(t,area)->(area==r),areaflow_df[ln][i])[!,:flowcontribution_sum]
                #calculated_flow_subsys[s]=calculated_flow_subsys[s]+ filter([:t,:area]=>(t,area)->(area==r),areaflow_df[i])[!,:flowcontribution_sum]
            end
            #calculated_flow_uc0=zeros(NT)
            #for s in subsys
            #    calculated_flow_uc0=calculated_flow_uc0+calculated_flow_uc0_subsys[s]
            #end  
            af[ln][!,"uc"*string(i)*"_flow"]= zeros(NT) #calculated_flow_uc0

            #for (k,v) in calculated_flow_uc0_subsys
            #    name="uc"*string(i)*"_flow_subsys"*k
            #    af[!,name]=calculated_flow_uc0_subsys[k]
            for s in subsys
                af[ln][!,"uc"*string(i)*"_flow"]=af[ln][!,"uc"*string(i)*"_flow"]+af[ln][!,"uc"*string(i)*"_flow_subsys_"*s]#calculated_flow_uc0_subsys[k]
            end 
        else
            for r in keys(area_region_dict)
                s=area_region_dict[r]
                name="uc"*string(i)*"_flow_subsys_"*s
                af[ln][!,name]=af[ln][!,name]+filter([:t,:area]=>(t,area)->(area==r),areaflow_df[ln][i])[!,:flowcontribution_sum]
                name="uc"*string(i)*"_loopflow_subsys_"*s
                af[ln][!,name]=af[ln][!,name]+filter([:t,:area]=>(t,area)->(area==r),areaflow_df[ln][i])[!,:loopflowcontribution_sum]
            end

            calculated_flow_uc2=Dict()
            af[ln][!,"uc"*string(i)*"_seflow"]=zeros(NT)
            for s in subsys
                af[ln][!,"uc"*string(i)*"_flow_"*s]=zeros(NT)
                af[ln][!,"uc"*string(i)*"_seflow"]=af[ln][!,"uc"*string(i)*"_seflow"]+af[ln][!,"uc"*string(i)*"_flow_subsys_"*s]#calculated_flow_subsys[s]
                for s1 in subsys
                    if s==s1
                        #calculated_flow_uc2[s]=calculated_flow_uc2[s]+calculated_flow_subsys[s1]
                        name="uc"*string(i)*"_flow_"*s
                        af[ln][!,name]=af[ln][!,name]+af[ln][!,"uc"*string(i)*"_flow_subsys_"*s1]#calculated_flow_subsys[s1]
                    else
                        #calculated_flow_uc2[s]=calculated_flow_uc2[s]+calculated_loopflow_subsys[s1]  
                        name="uc"*string(i)*"_flow_"*s
                        af[ln][!,name]=af[ln][!,name]+af[ln][!,"uc"*string(i)*"_loopflow_subsys_"*s1]#calculated_loopflow_subsys[s1]  
                    end    
                end  
            end
        end       
        CSV.write("af_"*case*"_"*CF*"_"*ln*".csv",af[ln])
    end
end 
println("time is ", t2-t1)

m1 = models.decision_models[1].internal.container
m2a = models.decision_models[2].internal.container.subproblems["a"]
m2b = models.decision_models[2].internal.container.subproblems["b"]
m2c = models.decision_models[2].internal.container.subproblems["c"]
m3a = models.decision_models[3].internal.container.subproblems["a"]
m3b = models.decision_models[3].internal.container.subproblems["b"]
m3c = models.decision_models[3].internal.container.subproblems["c"]

println(objective_value(m1.JuMPmodel))
println(objective_value(m2a.JuMPmodel))
println(objective_value(m2b.JuMPmodel))
println(objective_value(m2c.JuMPmodel))
println(objective_value(m3a.JuMPmodel))
println(objective_value(m3b.JuMPmodel))
println(objective_value(m3c.JuMPmodel))

results_uc[3]
typeof(results_uc[3])
CSV.write("uc2_local_o.csv",read_realized_variable(results_uc[2], "OnVariable__ThermalStandard"))
CSV.write("start.csv",read_realized_variable(results_uc[1], "StartVariable__ThermalStandard"))
CSV.write("ACBus.csv",read_realized_variable(results_uc[1], "ActivePowerBalance__ACBus"))

read_realized_variable(results_uc[1], "StartVariable__ThermalStandard")

g=read_realized_expression(get_decision_problem_results(results, "UC"*string(3)), "ProductionCostExpression__ThermalStandard")
sum(sum(eachcol(g[!, Not("DateTime")])))

g=read_realized_expression(get_decision_problem_results(results, "UC"*string(3)), "ProductionCostExpression__RenewableDispatch")
sum(sum(eachcol(g[!, Not("DateTime")])))

g=read_realized_expression(get_decision_problem_results(results, "UC"*string(3)), "ProductionCostExpression__HydroDispatch")
sum(sum(eachcol(g[!, Not("DateTime")])))


println("####")
results_uc
results_uc[1]
results_uc[2]
results_uc[3]

read_realized_variable(results_uc[2], "OnVariable__ThermalStandard")


area_df[1]."PowerLoad_sum"
bus_df_uc[2]
area_3 = groupby(bus_df_uc[2], :area)[2]


models.decision_models[1]
res = OptimizationProblemResults(models.decision_models[1])
typeof(results_uc)
typeof(results_uc[2])
results_uc[3]
# models
# typeof(models)
# typeof(models.decision_models[3])
# models.decision_models[1].constraints
# models.decision_models[2]
# models.decision_models[3]
models.decision_models[1]
m1 = models.decision_models[1].internal.container
m2a = models.decision_models[2].internal.container.subproblems["a"]
m2b = models.decision_models[2].internal.container.subproblems["b"]
m2c = models.decision_models[2].internal.container.subproblems["c"]
m3a = models.decision_models[3].internal.container.subproblems["a"]
m3b = models.decision_models[3].internal.container.subproblems["b"]
m3c = models.decision_models[3].internal.container.subproblems["c"]

JuMP.write_to_file(m3a.JuMPmodel, "m3a.lp")
JuMP.write_to_file(m3b.JuMPmodel, "m3b.lp")
JuMP.write_to_file(m3c.JuMPmodel, "m3c.lp")

JuMP.objective_coefficient


JuMP.coefficient(JuMP.objective_function(PSI.get_jump_model(m2a)) , relief["A18_1_1",1])


# models.decision_models[3].internal.container.subproblems
monitor_upper_key = PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.RateLimitConstraint, PSY.MonitoredLine}("ub")
monitor_lower_key = PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.RateLimitConstraint, PSY.MonitoredLine}("lb")
mon_line_name = m3b.constraints[monitor_lower_key].axes[1][1]
JuMP.normalized_rhs(m3b.constraints[monitor_upper_key][mon_line_name,1]) 
JuMP.set_normalized_rhs(m3b.constraints[monitor_upper_key][mon_line_name,1], 1 - 0.1) 
JuMP.normalized_constant(m3b.constraints[monitor_upper_key][mon_line_name,1]) 
typeof(PSI.get_time_steps(m2c))
length(PSI.get_time_steps(m3c))

monitor_lower_key = PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.RateLimitConstraint, PSY.MonitoredLine}("lb")
m2b.constraints[monitor_lower_key]["A28",1]
m3b.constraints[monitor_lower_key]["A28",1]
m2b.constraints[monitor_lower_key].axes[1][1]
m3c.constraints[monitor_lower_key].axes[1][1]
m3b.constraints[monitor_lower_key].axes[1]
PSI.get_constraint_keys(m3b)
PSI.get_variable_keys(m3b)
relief = m3a.variables[PSY.InfrastructureSystems.Optimization.VariableKey{PowerSimulationsDecomposition.Relief, PSY.MonitoredLine}("")]
# m3a.constraints[InfrastructureSystems.Optimization.ConstraintKey{NetworkFlowConstraint, MonitoredLine}("")]

flow = m3a.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerVariable, MonitoredLine}("")][mon_line_name,1]
JuMP.value(flow)
relief
JuMP.value.(relief)


m1.constraints
typeof(m1)
PSI.get_time_steps(m3b)

for (k,v) in m3a.variables
    println(k)
end

for (k,v) in m3b.constraints
    println(k)
end

m3a.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{NetworkFlowConstraint, MonitoredLine}("")]["A28",1]

flow = m3b.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.NetworkFlowConstraint, PSY.MonitoredLine}("")]["A28",1]
JuMP.dual(flow)

m3b.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.RateLimitConstraint, PSY.MonitoredLine}("ub")]["A28",1]
m3b.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.RateLimitConstraint, PSY.MonitoredLine}("lb")]["A28",1]
flow = m3b.constraints[PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.NetworkFlowConstraint, PSY.MonitoredLine}("")]["A28",1]
r = m3b.variables[PSY.InfrastructureSystems.Optimization.VariableKey{Relief, MonitoredLine}("")][1,1]
m3b.variables[PSY.InfrastructureSystems.Optimization.VariableKey{FlowActivePowerSlackLowerBound, MonitoredLine}("")]["A28",1]
m3b.variables[PSY.InfrastructureSystems.Optimization.VariableKey{FlowActivePowerSlackUpperBound, MonitoredLine}("")]["A28",1]

length(PSI.get_time_steps(m2a))

monitor_lower_key = PSY.InfrastructureSystems.Optimization.ConstraintKey{PSI.RateLimitConstraint, PSY.MonitoredLine}("lb")
haskey(m2a.constraints, monitor_lower_key)
m2a.constraints[monitor_lower_key].axes[1]
keys(m2b.constraints[monitor_lower_key].axes[1])

mon_line_name = m2a.constraints[monitor_lower_key].axes[1][1]
"A28" in m2a.constraints[monitor_lower_key].axes[1]

if "A28" in keys(m2a.constraints[monitor_lower_key].axes[1]) && index2 in keys(subproblem_container.variables[variable_key].axes[2])
haskey(m2a.constraints[monitor_lower_key], ["A28",1])
m2a.constraints[monitor_lower_key]

if haskey(subproblem_container.constraints, monitor_lower_key) && haskey(subproblem_container.constraints[monitor_lower_key], "A28")
    monitor_lower_con = subproblem_container.constraints[monitor_lower_key]["A28",1]
    monitor_lower_bound = JuMP.normalized_rhs(monitor_lower_con)  
    flow = subproblem_container.variables[PSY.InfrastructureSystems.Optimization.VariableKey{PSI.FlowActivePowerVariable, PSY.MonitoredLine}("")]["A28",1]
    JuMP.delete(subproblem_container.JuMPmodel, monitor_lower_con)
    monitor_lower_con = JuMP.@constraint(subproblem_container.JuMPmodel, flow - sum(region_relief[:, 1]) >= monitor_lower_bound)   
    subproblem_container.constraints[monitor_lower_key]["A28", 1] = monitor_lower_con
end











for r in eachrow(run_df)
    CF,NS,netmodels,copperplates,NTS,area_region_dicts,monitoredline_subsys=r[["CF","NS","netmodels","copperplates","NTS","area_region_dicts","monitoredline_subsys"]]
#    CF,NS,netmodels,copperplates,NTS,area_region_dicts,monitoredline_subsys=run_df[1,:][["CF","NS","netmodels","copperplates","NTS","area_region_dicts","monitoredline_subsys"]]
    syss=[sys]
    templates=[]

    for i in 1:NS
        NT=NTS[i]
        if i>1                               
            sys2=deepcopy(sys)
            push!(syss,sys2)
        end
        transform_single_time_series!(syss[i], Hour(NT), Hour(NT))

        area_region_dict=area_region_dicts[i]  
        if netmodels[i]=="SplitAreaPTDFPowerModel"     
            subsys=unique(values(area_region_dict))
            for v in subsys
                println("add sys",i," subsystem,",v)
                add_subsystem!(sys2, v)
            end
            template = 
            MultiProblemTemplate(NetworkModel(SplitAreaPTDFPowerModel; use_slacks=true,PTDF_matrix = ptdf,), subsys)
        else          
            template = ProblemTemplate(NetworkModel(AreaPTDFPowerModel; use_slacks = true, PTDF_matrix = ptdf,))
        end

        push!(templates,template)    

        if standardload==1
            PSI.set_device_model!(template, StandardLoad , StaticPowerLoad)
        else
            PSI.set_device_model!(template, PowerLoad, StaticPowerLoad)
        end

        PSI.set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
        PSI.set_device_model!(template, DeviceModel(HydroDispatch, HPS.HydroDispatchReservoirBudget,
                                        time_series_names = Dict{Any, String}(
                                            PSI.ActivePowerTimeSeriesParameter => "max_active_power",
                                            HPS.EnergyBudgetTimeSeriesParameter => "hydro_budget",)
                                                )
                        )
        if hvdc==1
            PSI.set_device_model!(template, TwoTerminalHVDCLine, HVDCTwoTerminalLossless)
        end

        set_device_model!(template, DeviceModel(MonitoredLine, StaticBranchUnbounded, use_slacks = true))
        
        if NT>1              
            PSI.set_device_model!(template, ThermalStandard, ThermalBasicUnitCommitment)
            PSI.set_device_model!(template, ThermalMultiStart, ThermalBasicUnitCommitment)
            #link_nsi_to_tielines. Need to set selected_line with tie lines
            #PSI.set_device_model!(template, AreaInterchange, StaticBranch)
            # not link_nsi_to_tielines. Need to activate HVDC        
            #PSI.set_device_model!(template, AreaInterchange, StaticBranch)
            PSI.set_device_model!(template, AreaInterchange, StaticBranchUnbounded)
        else
            PSI.set_device_model!(template, ThermalStandard, ThermalBasicDispatch)
            PSI.set_device_model!(template, AreaInterchange, StaticBranchUnbounded)
        end                  
        
        if copperplates[i]==0
            set_device_model!(template, DeviceModel(MonitoredLine, StaticBranch, use_slacks = true))
        end
                    
        set_device_model!(template,
            DeviceModel(Line, StaticBranchUnbounded; attributes=Dict("filter_function" => x -> get_name(x) in selected_line),))             

        subsys=unique(values(area_region_dict))

        if length(subsys)>1
        #buildsubsystem1(sys2, area_region_dict, union(selected_line,monitoredlined_line))
        buildsubsystem1(sys2, area_region_dict, selected_line,monitoredline_subsys[i])
        end
    end

    #include("build_model_sim_run.jl")

    dm=[]
    for i in 1:NS
        if netmodels[i]=="SplitAreaPTDFPowerModel"
            push!(dm, DecisionModel(MultiRegionProblem, templates[i], syss[i], name="UC"*string(i),optimizer=optimizer,
                store_variable_names=true,initialize_model=false,optimizer_solve_log_print=false, 
                direct_mode_optimizer=true,check_numerical_bounds=false, 
                calculate_conflict=true,rebuild_model=false,system_to_file = true))
        else    
            push!(dm, DecisionModel(templates[i],syss[i], name="UC"*string(i), optimizer=optimizer,
            store_variable_names=true,initialize_model=false,optimizer_solve_log_print=false, 
            direct_mode_optimizer=true,check_numerical_bounds=false, 
            calculate_conflict=true,rebuild_model=false, system_to_file = true))
        end
    end    

    if NS==2
        models = PSI.SimulationModels(decision_models=[dm[1],dm[2]])
    else    
        models = PSI.SimulationModels(decision_models=[dm[1],dm[2],dm[3]])
    end    

    uc_simulation_ffs=[]
    feedforwards_dict=Dict()
    for i in 1:NS-1
        uc_simulation_ff = Vector{PowerSimulations.AbstractAffectFeedforward}()
    #feed forward area interchange
        FVFF_area_interchange = FixValueFeedforward(;component_type=AreaInterchange,source=FlowActivePowerVariable,affected_values=[FlowActivePowerVariable],)
        push!(uc_simulation_ff, FVFF_area_interchange)
        #feed forward HVDC (doesn't seem to work)
        #FVFF_hvdc = FixValueFeedforward(;component_type=TwoTerminalHVDCLine,source=FlowActivePowerVariable,affected_values=[FlowActivePowerVariable],)
        #push!(uc_simulation_ff, FVFF_hvdc)

        if NTS[i+1]==1
            SCFF = SemiContinuousFeedforward(;
            component_type=ThermalStandard, source=OnVariable, affected_values=[ActivePowerVariable],)
            push!(uc_simulation_ff, SCFF) 
            #LBFF = FixValueFeedforward(;component_type=ThermalStandard,source=OnVariable, affected_values=[OnVariable],)
            #push!(uc_simulation_ff, LBFF)
        end
        push!(uc_simulation_ffs,uc_simulation_ff)
        feedforwards_dict["UC"*string(i+1)]=uc_simulation_ff
    end

    sequence = SimulationSequence(;
        models=models,
        feedforwards=feedforwards_dict, #Dict("UC2" => uc_simulation_ff[1], "UC3" => uc_simulation_ff[2],),
        ini_cond_chronology=InterProblemChronology(),
    );

    # Specify the simulation setup
    # Here we specify the simulation name, the initial/start time, number of steps/days to execute, and the simulation folder.
    sim = PSI.Simulation(
        name="ntps_3stageDA_DASplit_RTSplit",
        steps=1,
        models=models,
        sequence=sequence,
        simulation_folder=mktempdir(), #".", #
    )

    println("before PSI.build")
    PSI.build!(sim, serialize=false)
    println("after PSI.build")
    PSI.execute!(sim, enable_progress_bar=false)

    results = SimulationResults(sim,ignore_status = true)

    # check results and calculate flows in different optimization models and the actual flow called calculated_seflow
    #include("final_check_results.jl")  
    results_uc=[]; uc_areainterchange=[];uc_monitoredline=[];
    bus_df_uc=[];area_df=[];
    tm=Array(1:NTS[1])
    areaflow_df=Dict();af=Dict()
    for ln in monitoredlined_line
        areaflow_df[ln]=[]; 
        af[ln]=DataFrames.DataFrame(;tm)
    end    

    for i in 1:NS
        push!(results_uc,get_decision_problem_results(results, "UC"*string(i)))
        push!(uc_areainterchange,read_realized_variable(results_uc[i], "FlowActivePowerVariable__AreaInterchange"))
        #push!(uc_monitoredline,read_realized_variable(results_uc[i], "FlowActivePowerVariable__MonitoredLine"))

        if netmodels[i]=="SplitAreaPTDFPowerModel"
            push!(bus_df_uc,read_bus_df(results_uc[i],1,standardload)) 
        else    
            push!(bus_df_uc,read_bus_df(results_uc[i],0,standardload))
        end    

        bus_df_uc[i][!,"buscheck"]=bus_df_uc[i][!,"ActivePowerBalance__ACBus"]-bus_df_uc[i][!,"HVDC"]/100-bus_df_uc[i][!,"ThermalStandard"]/100-bus_df_uc[i][!,"Renewable"]/100-bus_df_uc[i][!,"Hydro"]/100-bus_df_uc[i][!,"PowerLoad"]/100
        println("UC0 bus check ActivePowerBalance__ACBus<>Gen/100-PowerLoad/100 ,",filter([:t,:buscheck] => (t,buscheck)-> (t>=1) && (abs(buscheck)>0.000001),  bus_df_uc[i]))
        # if some bus has non-zero buscheck value, then check the bus. Some MW may not be added to power balance correctly.
        #println(filter([:t,:bus] => (t,bus)-> (t==1) && (bus==113 || bus==316),  bus_df_uc0))

    ###### area MW check #########
    #area level net Gen, load and ActivePowerBalance__ACBus_sum check
        push!(area_df,combine(groupby(bus_df_uc[i], [:t, :area]), [:ThermalStandard, :Renewable, :Hydro, :PowerLoad, :ActivePowerBalance__ACBus, :StateEstimationInjections__ACBus] .=> sum))
        area_df[i][!,"areacheck"]=area_df[i][!,"ActivePowerBalance__ACBus_sum"]-area_df[i][!,"ThermalStandard_sum"]/100-area_df[i][!,"Renewable_sum"]/100-area_df[i][!,"Hydro_sum"]/100-area_df[i][!,"PowerLoad_sum"]/100
        println("i,area_df check,",sum(area_df[i][!,"areacheck"]))  #should be close to 0

        area_df[i][!,"TotalGen"]=area_df[i][!,"ThermalStandard_sum"]+area_df[i][!,"Renewable_sum"]+area_df[i][!,"Hydro_sum"]
        area_df[i][!,"NetInterchange"]=area_df[i][!,"TotalGen"]+area_df[i][!,"PowerLoad_sum"]
        area_df[i][!,"InterchangeCheck"]=area_df[i][!,"NetInterchange"]-100*area_df[i][!,"ActivePowerBalance__ACBus_sum"]
        println(area_df[i][:,["t","area","ActivePowerBalance__ACBus_sum","NetInterchange","InterchangeCheck"]])

    ####### area flow check #########
    #area flow contribution by area and interval
        bus_df_uc[i][!,:gsf].=0.0; bus_df_uc[i][!,:flowcontribution].=0.0; bus_df_uc[i][!,:loopflowcontribution].=0.0
        for ln in monitoredlined_line
            for b in eachrow(bus_df_uc[i])
               gsf=ptdf[ln,b[:bus]]; 
               b[:gsf]=gsf; b[:flowcontribution]=gsf*b[:ActivePowerBalance__ACBus]; b[:loopflowcontribution]=gsf*b[:StateEstimationInjections__ACBus]
            end

            push!(areaflow_df[ln], combine(groupby(bus_df_uc[i], [:t, :area]), [:flowcontribution, :loopflowcontribution] .=> sum))

            if length(area_region_dicts[i])==0
                area_region_dict=area_dict
            else   
                area_region_dict=area_region_dicts[i]
            end      
            
            subsys=unique(values(area_region_dict))

            calculated_flow_subsys=Dict()
            calculated_loopflow_subsys=Dict()
            NT=NTS[1]
            for s in subsys
                #calculated_flow_subsys[s]=zeros(NT)
                #calculated_loopflow_subsys[s]=zeros(NT)
                af[ln][!,"uc"*string(i)*"_flow_subsys_"*s]=zeros(NT)
                af[ln][!,"uc"*string(i)*"_loopflow_subsys_"*s]=zeros(NT)
            end    

            if netmodels[i]!="SplitAreaPTDFPowerModel"
                for r in keys(area_region_dict)
                    s=area_region_dict[r]
                    name="uc"*string(i)*"_flow_subsys_"*s
                    af[ln][!,name]=af[ln][!,name]+filter([:t,:area]=>(t,area)->(area==r),areaflow_df[ln][i])[!,:flowcontribution_sum]
                    #calculated_flow_subsys[s]=calculated_flow_subsys[s]+ filter([:t,:area]=>(t,area)->(area==r),areaflow_df[i])[!,:flowcontribution_sum]
                end
                #calculated_flow_uc0=zeros(NT)
                #for s in subsys
                #    calculated_flow_uc0=calculated_flow_uc0+calculated_flow_uc0_subsys[s]
                #end  
                af[ln][!,"uc"*string(i)*"_flow"]= zeros(NT) #calculated_flow_uc0

                #for (k,v) in calculated_flow_uc0_subsys
                #    name="uc"*string(i)*"_flow_subsys"*k
                #    af[!,name]=calculated_flow_uc0_subsys[k]
                for s in subsys
                    af[ln][!,"uc"*string(i)*"_flow"]=af[ln][!,"uc"*string(i)*"_flow"]+af[ln][!,"uc"*string(i)*"_flow_subsys_"*s]#calculated_flow_uc0_subsys[k]
                end 
            else
                for r in keys(area_region_dict)
                    s=area_region_dict[r]
                    name="uc"*string(i)*"_flow_subsys_"*s
                    af[ln][!,name]=af[ln][!,name]+filter([:t,:area]=>(t,area)->(area==r),areaflow_df[ln][i])[!,:flowcontribution_sum]
                    name="uc"*string(i)*"_loopflow_subsys_"*s
                    af[ln][!,name]=af[ln][!,name]+filter([:t,:area]=>(t,area)->(area==r),areaflow_df[ln][i])[!,:loopflowcontribution_sum]
                end

                calculated_flow_uc2=Dict()
                af[ln][!,"uc"*string(i)*"_seflow"]=zeros(NT)
                for s in subsys
                    af[ln][!,"uc"*string(i)*"_flow_"*s]=zeros(NT)
                    af[ln][!,"uc"*string(i)*"_seflow"]=af[ln][!,"uc"*string(i)*"_seflow"]+af[ln][!,"uc"*string(i)*"_flow_subsys_"*s]#calculated_flow_subsys[s]
                    for s1 in subsys
                        if s==s1
                            #calculated_flow_uc2[s]=calculated_flow_uc2[s]+calculated_flow_subsys[s1]
                            name="uc"*string(i)*"_flow_"*s
                            af[ln][!,name]=af[ln][!,name]+af[ln][!,"uc"*string(i)*"_flow_subsys_"*s1]#calculated_flow_subsys[s1]
                        else
                            #calculated_flow_uc2[s]=calculated_flow_uc2[s]+calculated_loopflow_subsys[s1]  
                            name="uc"*string(i)*"_flow_"*s
                            af[ln][!,name]=af[ln][!,name]+af[ln][!,"uc"*string(i)*"_loopflow_subsys_"*s1]#calculated_loopflow_subsys[s1]  
                        end    
                    end  
                end
            end       
            CSV.write("af_"*case*"_"*CF*"_"*ln*".csv",af[ln])
        end
    end        
end

using Plots
ln="A28"
af0=DataFrame(CSV.File("af_RTS_NS3-0_"*ln*".csv"))
af1=DataFrame(CSV.File("af_RTS_NS3-1_"*ln*".csv"))
af2=DataFrame(CSV.File("af_RTS_NS3-2_"*ln*".csv"))
af3=DataFrame(CSV.File("af_RTS_NS3-3_"*ln*".csv"))

p0=plot(af0.tm, [af0.uc1_flow  af0.uc2_flow  af0.uc3_flow], label=["NS3-0_uc1_flow" "NS3-0_uc2_flow" "NS3-0_uc3_flow"], line=(3,[:solid :solid :dash]))
p1=plot(af1.tm, [af1.uc1_flow  af1.uc2_seflow  af1.uc3_seflow], label=["NS3-1_uc1_flow" "NS3-1_uc2_seflow" "NS3-1_uc3_seflow"], line=(3,[:solid :solid :dash]))
p2=plot(af2.tm, [af2.uc1_flow  af2.uc2_seflow  af2.uc3_seflow], label=["NS3-2_uc1_flow" "NS3-2_uc2_seflow" "NS3-2_uc3_seflow"], line=(3,[:solid :solid :dash]))
p3=plot(af3.tm, [af3.uc1_flow  af3.uc2_flow  af3.uc3_seflow], label=["NS3-3_uc1_flow" "NS3-3_uc2_flow" "NS3-3_uc3_seflow"], line=(3,[:solid :solid :dash]))
p4=plot(af0.tm, [af0.uc3_flow  af1.uc3_seflow  af2.uc3_seflow af3.uc3_seflow], label=["NS3-0_uc3_flow" "NS3-1_uc3_seflow" "NS3-2_uc3_seflow" "NS3-3_uc3_seflow"], line=(3,[:solid :solid :dash :dash]))

#=
af0=DataFrame(CSV.File("af_RTS_NS3-0.csv"))
af1=DataFrame(CSV.File("af_RTS_NS3-1.csv"))
af2=DataFrame(CSV.File("af_RTS_NS3-2.csv"))
af3=DataFrame(CSV.File("af_RTS_NS3-3.csv"))

p0=plot(af0.tm, [af0.uc1_flow  af0.uc2_flow  af0.uc3_flow], label=["NS3-0_uc1_flow" "NS3-0_uc2_flow" "NS3-0_uc3_flow"], line=(3,[:solid :solid :dash]))
p1=plot(af1.tm, [af1.uc1_flow  af1.uc2_seflow  af1.uc3_seflow], label=["NS3-1_uc1_flow" "NS3-1_uc2_seflow" "NS3-1_uc3_seflow"], line=(3,[:solid :solid :dash]))
p2=plot(af2.tm, [af2.uc1_flow  af2.uc2_seflow  af2.uc3_seflow], label=["NS3-2_uc1_flow" "NS3-2_uc2_seflow" "NS3-2_uc3_seflow"], line=(3,[:solid :solid :dash]))
#p3=plot(af3.tm, [af3.uc1_flow  af3.uc2_seflow  af3.uc3_seflow], label=["NS3-3_uc1_flow" "NS3-3_uc2_seflow" "NS3-3_uc3_seflow"], line=(3,[:solid :solid :dash]))
p3=plot(af3.tm, [af3.uc1_flow  af3.uc2_flow  af3.uc3_seflow], label=["NS3-3_uc1_flow" "NS3-3_uc2_seflow" "NS3-3_uc3_seflow"], line=(3,[:solid :solid :dash]))
p4=plot(af0.tm, [af0.uc3_flow  af1.uc3_seflow  af2.uc3_seflow af3.uc3_seflow], label=["NS3-0_uc3_flow" "NS3-1_uc3_seflow" "NS3-2_uc3_seflow" "NS3-3_uc3_seflow"], line=(3,[:solid :solid :dash :dash]))
=#
plot(p0,p1,p2,p3,p4, layout = (3, 2))
plot!(size=(1200,1200))
#=
uc0=sim.models.decision_models[1].internal.container
uc2=sim.models.decision_models[2].internal.container
uc3=sim.models.decision_models[3].internal.container
sum(value.(uc0.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerVariable, MonitoredLine}("")]["A28",:])-af[!,"uc1_flow"])
sum(value.(uc2.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerVariable, MonitoredLine}("")]["A28",:])-af[!,"uc2_flow"])

uc2a=sim.models.decision_models[2].internal.container.subproblems["a"]
uc2b=sim.models.decision_models[2].internal.container.subproblems["b"]
uc2c=sim.models.decision_models[2].internal.container.subproblems["c"]
sum(value.(uc0.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerVariable, MonitoredLine}("")]["B28",:])-af[!,"uc1_flow"])
sum(value.(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerVariable, MonitoredLine}("")]["B28",:])-af[!,"uc2_flow_a"])
sum(value.(uc2b.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerVariable, MonitoredLine}("")]["A28",:])-af[!,"uc2_flow_b"])
sum(value.(uc2c.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerVariable, MonitoredLine}("")]["A28",:])-af[!,"uc2_flow_c"])

uc3a=sim.models.decision_models[3].internal.container.subproblems["a"]
uc3b=sim.models.decision_models[3].internal.container.subproblems["b"]
sum(value.(uc3a.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerVariable, MonitoredLine}("")]["A28",:])-af[!,"uc3_flow_a"])
sum(value.(uc3b.variables[InfrastructureSystems.Optimization.VariableKey{FlowActivePowerVariable, MonitoredLine}("")]["A28",:])-af[!,"uc3_flow_b"])

uc2a.constraints[InfrastructureSystems.Optimization.ConstraintKey{NetworkFlowConstraint, MonitoredLine}("")] 
uc2a.constraints[InfrastructureSystems.Optimization.ConstraintKey{RateLimitConstraint, MonitoredLine}("ub")]

uc2a.constraints[InfrastructureSystems.Optimization.ConstraintKey{CopperPlateBalanceConstraint, Area}("")] 

objective_value(uc2a.JuMPmodel)+objective_value(uc2b.JuMPmodel)+objective_value(uc2c.JuMPmodel)-objective_value(uc0.JuMPmodel)

objterms(uc2,"obj_uc2a_"*CF)
write_coeff(uc2,"A28",1,"A28_uc2a_"*CF,1)

objterms(uc2a,"obj_uc2a_fix_"*CF)
write_coeff(uc2a,"A28",1,"A28_uc2a_fix_"*CF,1)


fix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,RenewableDispatch}("")]["122_WIND_1",1],1.8596521707224; force=true)
fix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["107_CC_1",1],1.71962875337936; force=true)
fix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["123_STEAM_3",1],2.14337921320769; force=true)
fix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["123_CT_5",1],0.22; force=true)
fix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["121_NUCLEAR_1",1],4; force=true)

fix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["115_STEAM_3",1],0; force=true)
fix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["116_STEAM_1",1],0; force=true)
fix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["118_CC_1",1],0; force=true)

fix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["123_CT_1",1],0; force=true)
fix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["123_CT_4",1],0; force=true)

unfix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["115_STEAM_3",1])
unfix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["116_STEAM_1",1])
unfix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["118_CC_1",1])
unfix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["123_CT_1",1])
unfix(uc2a.variables[InfrastructureSystems.Optimization.VariableKey{ActivePowerVariable,ThermalStandard}("")]["123_CT_4",1])


con2a=uc2a.constraints[InfrastructureSystems.Optimization.ConstraintKey{CopperPlateBalanceConstraint, Area}("")]["1",1]
write_cons_coeff(uc2a,con2a,"pb_1_1_fix_"*CF)

PSI.compute_conflict!(uc2a)

con2=uc2.constraints[InfrastructureSystems.Optimization.ConstraintKey{CopperPlateBalanceConstraint, Area}("")]["1",1]
write_cons_coeff(uc2,con2,"pb_1_1_"*CF)


=#