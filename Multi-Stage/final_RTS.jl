sys = build_system(PSISystems, "modified_RTS_GMLC_DA_sys")

selected_line=[]#["CA-1", "CB-1", "AB1"] 
line_name = "A18"#,"B28"
monitoredlined_line=[line_name]
# monitoredlined_line=["CA-1", "C35"]
#monitoredlined_line=["B28"]
# limit=Dict("A3" => 0.5)
limit=Dict(line_name => 1.0)
# limit=Dict("CA-1" => 1.0, "C35" => 1.0)
hvdc=1

standardload=0
exchange_1_2 = AreaInterchange(;
    name="1_2", available=true, active_power_flow=0.0, from_area=get_component(Area, sys, "1"), to_area=get_component(Area, sys, "2"),
    flow_limits=(from_to=99999, to_from=99999),
)
add_component!(sys, exchange_1_2)

exchange_1_3 = AreaInterchange(;
    name="1_3", available=true, active_power_flow=0.0, from_area=get_component(Area, sys, "1"), to_area=get_component(Area, sys, "3"),
    flow_limits=(from_to=99999, to_from=99999),
)
add_component!(sys, exchange_1_3)

exchange_2_3 = AreaInterchange(;
    name="2_3", available=true, active_power_flow=0.0, from_area=get_component(Area, sys, "2"), to_area=get_component(Area, sys, "3"),
    flow_limits=(from_to=99999, to_from=99999),
)
add_component!(sys, exchange_2_3)

area_dict=Dict("1" => "1", "2" => "2", "3" =>"3")

l=PSY.get_component(TwoTerminalHVDCLine, sys, "DC1")
set_active_power_limits_from!(l, (min = 0.0, max = 0.0))

run_df=DataFrames.DataFrame(;
CF=[], NS=[], netmodels = [], copperplates=[],NTS=[],area_region_dicts=[],monitoredline_subsys=[])

CF="NS3-0"
NS=3
netmodels=["AreaPTDFPowerModel","AreaPTDFPowerModel","AreaPTDFPowerModel"]
copperplates=[1,0,0]
NTS=[24,24,1]
area_region_dicts=[Dict(), Dict(), Dict()]
monitoredline_subsys=[Dict(),Dict(),Dict()]
push!(run_df,(CF,NS,netmodels,copperplates,NTS,area_region_dicts,monitoredline_subsys))

CF="NS3-1"
NS=3
netmodels=["AreaPTDFPowerModel","SplitAreaPTDFPowerModel","SplitAreaPTDFPowerModel"]
copperplates=[1,0,0]
NTS=[24,24,1]
area_region_dicts=[Dict(),
                  Dict("1" => "a", "2" => "b", "3" =>"c"),
                  Dict("1" => "a", "2" => "b", "3" =>"c")]

# monitoredline_subsys=[Dict(),Dict("A28" => ["a"],"B28" => ["b"]),Dict("A28" => ["a"],"B28" => ["b"])]
monitoredline_subsys=[Dict(),Dict(line_name => ["a"],line_name => ["b"]),Dict(line_name => ["a"],line_name => ["b"])]
push!(run_df,(CF,NS,netmodels,copperplates,NTS,area_region_dicts,monitoredline_subsys))

CF="NS3-2"
NS=3
netmodels=["AreaPTDFPowerModel","SplitAreaPTDFPowerModel","SplitAreaPTDFPowerModel"]
copperplates=[1,0,0]
NTS=[24,24,1]
area_region_dicts=[Dict(),
                  Dict("1" => "a", "2" => "b", "3" =>"c"),
                  Dict("1" => "a", "2" => "b", "3" =>"c")]
# monitoredline_subsys=[Dict(),Dict("A28" => ["a","b","c"],"B28" => ["a","b","c"]),Dict("A28" => ["a","b","c"],"B28" => ["a","b","c"])]
# monitoredline_subsys=[Dict(),Dict("A28" => ["a","b"],"B28" => ["a","b"]),Dict("A28" => ["a"],"B28" => ["b"])]
# monitoredline_subsys=[Dict(),Dict("CA-1" => ["a","b"],"C35" => ["a","b"]),Dict("CA-1" => ["a","b"],"C35" => ["a","b"])]
monitoredline_subsys=[Dict(),Dict(line_name => ["a","b","c"]),Dict(line_name => ["a","b","c"])]
push!(run_df,(CF,NS,netmodels,copperplates,NTS,area_region_dicts,monitoredline_subsys))

CF="NS3-3"
NS=3
netmodels=["AreaPTDFPowerModel","AreaPTDFPowerModel","SplitAreaPTDFPowerModel"]
copperplates=[1,0,0]
NTS=[24,24,1]
area_region_dicts=[Dict(),
                  Dict(),
                  Dict("1" => "ab", "2" => "ab", "3" =>"c")]
# monitoredline_subsys=[Dict(),Dict(),Dict("A18" => ["ab"],"A20" => ["ab"])]
monitoredline_subsys=[Dict(),Dict(),Dict(line_name => ["ab"],line_name => ["ab"])]
push!(run_df,(CF,NS,netmodels,copperplates,NTS,area_region_dicts,monitoredline_subsys))

case="RTS"