using StateCharts
using AlgebraicABMs, AlgebraicRewriting, Catlab
using Fleck
using Random
import AlgebraicRewriting.Incremental.IncrementalSum: IncSumHomSet
ENV["JULIA_DEBUG"] = "AlgebraicABMs"; 

############################################## the interface functions ####################################################
## Functions in this section are all included in the file "StateChartsInterface.jl", but this file is not #################
## usable yet .After the extension dependency of "StateCharts.jl" registered, just delete these functions #################
## in this section.                                                                                       #################
###########################################################################################################################
const DST = :TimeOut # the transtition type for discrete hazard
const CTT = :Rate # the transtition type for continueous hazard
const RLT = :Pattern # the transtition type for user defined rewrite rule

""" automatically generate the ABM model schema from the state chart. This schema contains one object, and the name of this object is given by the input 
argument: obn::Symbol. And each unit state chart (with one entrance transition) contributes one attributes of the object "obn" in the schema

"""

function StateChartCset_SingleObject(ss::AbstractStateChart,obn::Symbol=:P)
    attrsType = enames(ss) # the attributes' type is the name of each unit state chart
    obns = repeat([obn], length(attrsType))
    attrs = map((x,y)->Symbol(String(x)*String(y)),obns,attrsType) # define the name of the map from the object "obn" to each attribute type
    attrsTuple = map((x,y,z)->(x,y,z),attrs,obns,attrsType)
    schema = BasicSchema([obn], [], attrsType, attrsTuple)
    tp_asgn = Dict(zip(attrsType,repeat([Symbol], length(attrsType))))
    acset = AnonACSet(schema;type_assignment=tp_asgn)
    acset
end

function StateChartCset_MultipleObject(ss::AbstractStateChart)   
end

function make_infectious_rule(ss::AbstractStateChart;make_acset::Function=StateChartCset_SingleObject,r::Symbol=:I,obn::Symbol=:P)
    Lᵢ,Iᵢ,Rᵢ = [make_acset(ss,obn) for _ in 1:3]
    # return the attribute names 
    ## TODO: need to make "attrn" support multiple unit state chart, this needs to extend the 
    ## schema of the state charts
    attrn = attrs(acset_schema(Lᵢ);just_names=true)[1]
    add_parts!(Lᵢ,obn,2;attrn=>[:S,:I])
    add_part!(Iᵢ,obn;attrn=>:I)
    add_parts!(Rᵢ,obn,2;attrn=>[r,:I])
    Rule(homomorphism(Iᵢ,Lᵢ),homomorphism(Iᵢ,Rᵢ)) 
end

# currently, suppose we only have one unit state chart
# TODO: make it support multiple unit state charts
# generate the rewrite rule
# t is the index of a transition
# transitions_rules is an array with pair of (timer, rule) as its elements
function get_rule(ss::AbstractStateChart,t,transitions_rules,obn::Symbol=:P)
    tt = ttype(ss,t)
    if tt == RLT
        tidx=Int(texpr(ss,t))
        rule = transitions_rules[tidx][2]
    elseif tt in [DST, CTT]
        L,I,R = LIR = [StateChartCset_SingleObject(ss,obn) for _ in 1:3]
        # return the attribute names 
        ## TODO: need to make "attrn" and "attrt" support multiple unit state chart, this needs to extend the 
        ## schema of the state charts
        attrn = attrs(acset_schema(L);just_names=true)[1]
        attrt = attrtypes(acset_schema(L))[1]
        add_part!(L,obn;attrn=>sname(ss,tsource(ss,t)))
        add_part!(I,obn;attrn=>AttrVar(add_part!(I,attrt)))
        add_part!(R,obn;attrn=>sname(ss,ttarget(ss,t)))
        rule = Rule(homomorphism(I,L),homomorphism(I,R))  
    else
        tn = tname(ss,t)
        throw("transition:$tn does not have valid transition type!")
    end
end

function get_timer(ss::AbstractStateChart,t,obn::Symbol=:P)
    tt=ttype(ss,t)
    if tt == DST
        timeout = texpr(ss,t)
        timer = DiscreteHazard(timeout)
    elseif tt == CTT
        r = texpr(ss,t)
        timer = ContinuousHazard(1/r)
    elseif tt == RLT
        tidx=Int(texpr(ss,t))
        timer = transitions_rules[tidx][1]
    else
        tn = tname(ss,t)
        throw("transition:$tn does not have valid transition type!")
    end
    timer
end

function make_ABM(ss::AbstractStateChart) 
    return ABM(map(parts(ss, :T)) do t
        tn = tname(ss, t)
        ABMRule(tn, get_rule(ss, t,transitions_rules), get_timer(ss,t))
    end, [])
end

#############################################################################################
######################### end of the dependency function part ###############################


### example code starts here ###

totalPopulation=10.0
β = 0.001

# create the pertussis state chart
pertussis = StateChartF((:S,:E,:I,:R₄,:V₄), # states
(:newExposure=>(:S=>:E,:Pattern=>1),:becomingInfective=>(:E=>:I,:TimeOut=>5.0),:recovery=>(:I=>:R₄,:Rate=>1.0/21.0),:vaccinated=>(:S=>:V₄,:Rate=>0.01)), #non-start transitions
(), # alternatives for non-start transitions
(:PertussisStateChart=>:S), # start transitions
(:PertussisStateChart=>((:initialInfective,1.0/totalPopulation)=>:I)) # alternatives for start transitions
)

StateCharts.Graph(pertussis)

# define the new_infectious rule
transitions_rules=[(ContinuousHazard(1/β),make_infectious_rule(pertussis,r=:E))]

# Initial state
init = StateChartCset_SingleObject(pertussis)
add_parts!(init,:P,Int(totalPopulation-1),PPertussisStateChart=:S)
add_part!(init,:P,PPertussisStateChart=:I)

res = run!(make_ABM(pertussis), init; maxtime=20);

#Makie.plot(res; Dict(o=>X->nparts(X,o) for o in [:S,:I])...)