using StateCharts
using AlgebraicABMs, AlgebraicRewriting, Catlab
using Fleck
using Random
import AlgebraicRewriting.Incremental.IncrementalSum: IncSumHomSet


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

homs_name(cod::Symbol,dom::Symbol)=Symbol(String(cod)*String(dom))

vectorify(n::Vector) = collect(n)
vectorify(n::Tuple) = length(n) == 1 ? [n] : collect(n)
vectorify(n::SubArray) = collect(n)
vectorify(n) = [n]


# function return the ABM model Schema based on the input state charts
# the schema's structure is:
# 1. one object: usually is Person
# 2. each unit state chart contributes to one attribute

function StateChartABMSchema_SingleObject(ss::AbstractStateChart,obn::Symbol=:P)
    attrsType = enames(ss) # the attributes' type is the name of each unit state chart
    obns = repeat([obn], length(attrsType))
    attrs = map((x,y)->homs_name(x,y),obns,attrsType) # define the name of the map from the object "obn" to each attribute type
    attrsTuple = map((x,y,z)->(x,y,z),attrs,obns,attrsType)
    schema = BasicSchema([obn], [], attrsType, attrsTuple)
    schema
end

# function return ABM model schema based on the input state charts
# the schema's structure is:
# 1. Objects: total number of states in state charts + 1 (the coproduct of the objects of states), the name of the coproduct
#             object is given by input argument "obn"
# 2. Morphisms: each object of state has a morphism to the "1" -- the coproduct of states

function StateChartABMSchema_MultipleObjects(ss::AbstractStateChart,obn::Symbol=:P)
    # generate the list of objects
    obs = vcat(snames(ss),[obn])
    # generate the list of morphisms
    morphisms_names = map((x,y)-> homs_name(x,y), snames(ss), repeat([obn],length(snames(ss))))
    morphisms = map((x,y,z)->(x,y,z),morphisms_names,snames(ss),repeat([obn],length(snames(ss))))
    # create the schema
    schema = BasicSchema(obs, morphisms)
    # return the cset generated
    schema       
end

function StateChartCset_SingleObject(ss::AbstractStateChart,obn::Symbol=:P)
    attrsType = enames(ss)
    schema = StateChartABMSchema_SingleObject(ss,obn)
    tp_asgn = Dict(zip(attrsType,repeat([Symbol], length(attrsType))))
    acset = AnonACSet(schema;type_assignment=tp_asgn)
    acset
end

StateChartCset_MultipleObjects(ss::AbstractStateChart,obn::Symbol=:P) = AnonACSet(StateChartABMSchema_MultipleObjects(ss,obn))        

# simulate a single transition, e.g., :I->:E
function make_rule_SingleObject(ss::AbstractStateChart, l, r, i=[]; obn::Symbol=:P)
    l = vectorify(l)
    i = vectorify(i)
    r = vectorify(r)
    Lᵢ,Iᵢ,Rᵢ = [StateChartCset_SingleObject(ss,obn) for _ in 1:3]
    # return the attribute names 
    ## TODO: need to make "attrn" support multiple unit state chart, this needs to extend the 
    ## schema of the state charts
    attrn = attrs(acset_schema(Lᵢ);just_names=true)[1]
    attrt = attrtypes(acset_schema(Lᵢ))[1]
    add_parts!(Lᵢ,obn,length(l);attrn=>l)
    add_part!(Iᵢ,obn;attrn=>AttrVar(add_part!(Iᵢ,attrt)))
    add_parts!(Rᵢ,obn,length(r);attrn=>r)
    Rule(homomorphism(Iᵢ,Lᵢ),homomorphism(Iᵢ,Rᵢ)) 
end

function make_infectious_rule_SingleObject(ss::AbstractStateChart,l,r,i,obn::Symbol=:P)
    l = vectorify(l)
    i = vectorify(i)
    r = vectorify(r)
    Lᵢ,Iᵢ,Rᵢ = [StateChartCset_SingleObject(ss,obn) for _ in 1:3]
    # return the attribute names 
    ## TODO: need to make "attrn" support multiple unit state chart, this needs to extend the 
    ## schema of the state charts
    attrn = attrs(acset_schema(Lᵢ);just_names=true)[1]
    add_parts!(Lᵢ,obn,2;attrn=>l)
    ic = findall(in(l), i)
    @assert length(ic)==1 "There should only be 1 infective individual in the infectious rule!"
    add_part!(Iᵢ,obn;attrn=>i[ic][1])
    attrt = attrtypes(acset_schema(Lᵢ))[1]
    add_part!(Iᵢ,obn;attrn=>AttrVar(add_part!(Iᵢ,attrt)))
    add_parts!(Rᵢ,obn,2;attrn=>r)
    Rule(homomorphism(Iᵢ,Lᵢ),homomorphism(Iᵢ,Rᵢ))
end

# l, include the objects of l, e.g., l=[:S,:I]
# i, include the objects of i, e.g., i=[:S]
# r, include the objects of r, e.g., r=[:I,:I]
function make_rule_MultipleObjects(ss::AbstractStateChart, l, r, i=[]; obn::Symbol=:P)
    l = vectorify(l)
    i = vectorify(i)
    r = vectorify(r)
    Lᵢ,Iᵢ,Rᵢ = LIR = [StateChartCset_MultipleObjects(ss,obn) for _ in 1:3]
    # add the object obn to L, I and R
    map(x->add_part!(x,obn),LIR)

    # add the (state) objects to L, I and R
    map(x->add_part!(Lᵢ, x; homs_name(x,obn)=>1),l)
    map(x->add_part!(Iᵢ, x; homs_name(x,obn)=>1),i)
    map(x->add_part!(Rᵢ, x; homs_name(x,obn)=>1),r)

    Rule(homomorphism(Iᵢ,Lᵢ),homomorphism(Iᵢ,Rᵢ)) 
end

function make_infectious_rule_MultipleObjects(ss::AbstractStateChart, l, r, i; obn::Symbol=:P) 
    make_rule_MultipleObjects(ss, l, r, i, obn=obn)
end

# currently, suppose we only have one unit state chart
# TODO: make it support multiple unit state charts
# generate the rewrite rule
# t is the index of a transition
# transitions_rules is an array with pair of (timer, rule) as its elements
function get_rule(ss::AbstractStateChart,t,transitions_rules;is_schema_singObject::Bool=true,obn::Symbol=:P)
    tt = ttype(ss,t)
    if tt == RLT
        tidx=Int(texpr(ss,t))
        return transitions_rules[tidx][2]
    elseif tt in [DST, CTT]
        # return the attribute names 
        ## TODO: need to make "attrn" and "attrt" support multiple unit state chart, this needs to extend the 
        ## schema of the state charts
        l = sname(ss,tsource(ss,t))
        r = sname(ss,ttarget(ss,t))
        return is_schema_singObject ? make_rule_SingleObject(ss,l,r,obn=obn) : make_rule_MultipleObjects(ss,l,r,obn=obn)
    else
        tn = tname(ss,t)
        throw("transition:$tn does not have valid transition type!")
    end
end

function get_timer(ss::AbstractStateChart, transitions_rules,t;obn::Symbol=:P)
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

function make_ABM(ss::AbstractStateChart, transitions_rules;is_schema_singObject::Bool=true) 
    return ABM(map(parts(ss, :T)) do t
        tn = tname(ss, t)
        ABMRule(tn, get_rule(ss, t,transitions_rules,is_schema_singObject=is_schema_singObject), get_timer(ss,transitions_rules,t))
    end, [])
end

#############################################################################################
######################### end of the dependency function part ###############################

