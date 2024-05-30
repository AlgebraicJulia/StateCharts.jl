module TestStateChartsABMsInterfaces

using Test
using Makie, CairoMakie
using StateCharts
using AlgebraicABMs
using Catlab


ENV["JULIA_DEBUG"] = "AlgebraicABMs"; 

### example code starts here ###

totalPopulation=5.0
β = 0.1
p_becomingInfective = 5.0

# create the pertussis state chart
pertussis = StateChartF((:S,:E,:I,:R₄,:V₄), # states
(:newExposure=>(:S=>:E,:Pattern=>1),:becomingInfective=>(:E=>:I,:TimeOut=>p_becomingInfective),:recovery=>(:I=>:R₄,:Rate=>1.0/21.0),:vaccinated=>(:S=>:V₄,:Rate=>0.01)), #non-start transitions
(), # alternatives for non-start transitions
(:PertussisStateChart=>:S), # start transitions
(:PertussisStateChart=>((:initialInfective,1.0/totalPopulation)=>:I)) # alternatives for start transitions
)

StateCharts.Graph(pertussis)

########### use ABM model schema of one object P #####################
########### the states of state charts are attributes ################

# define the new_infectious rule
transitions_rules_SingleObject=[(ContinuousHazard(1/β),make_infectious_rule_SingleObject(pertussis,[:S,:I],[:I,:I],[:I]))]

# Initial state for schema of single object
init_SingleObject = StateChartCset_SingleObject(pertussis)
add_parts!(init_SingleObject,:P,Int(totalPopulation-1),PPertussisStateChart=:S)
add_part!(init_SingleObject,:P,PPertussisStateChart=:I)

# need to update after AlgebraicABMs bugs fixed
#res = run!(make_ABM(pertussis,transitions_rules_SingleObject,is_schema_singObject=true), init_SingleObject; maxtime=20);

########### use ABM model schema of multiple objects#####################
########### each state of state charts is an object, they are all #######
########### map to a collection object, e.g., total population ##########

transitions_rules_MultipleObjects=[(ContinuousHazard(1/β),make_infectious_rule_MultipleObjects(pertussis,[:S,:I],[:I,:I],[:I]))]

# Initial state for schema of multiple objects
init_MultipleObjects = StateChartCset_MultipleObjects(pertussis)
add_parts!(init_MultipleObjects,:P,Int(totalPopulation))
add_parts!(init_MultipleObjects,:S,Int(totalPopulation-1),SP=1)
add_part!(init_MultipleObjects,:I,IP=1)

res = run!(make_ABM(pertussis,transitions_rules_MultipleObjects,is_schema_singObject=false), init_MultipleObjects; maxtime=10);

Makie.plot(res; Dict(o=>X->nparts(X,o) for o in [:S,:E,:I,:R₄,:V₄])...)

end