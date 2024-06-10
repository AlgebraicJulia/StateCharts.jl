module StateChartsABMsInterfaces

export StateChartABMSchema_SingleObject, StateChartABMSchema_MultipleObjects, StateChartCset_SingleObject, StateChartCset_MultipleObjects,
make_rule_SingleObject, make_infectious_rule_SingleObject, make_rule_MultipleObjects, make_infectious_rule_MultipleObjects, get_rule,
get_timer, make_ABM

using AlgebraicABMs, AlgebraicRewriting, Catlab
using Random
using ..StateChartsSchema

#import AlgebraicRewriting.Incremental.IncrementalSum: IncSumHomSet

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

function StateChartABMSchema_SingleObject(ss::AbstractUnitStateChart,obn::Symbol=:P)
    attrsType = dnames(ss) # the attributes' type is the name of each unit state chart
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

function StateChartABMSchema_MultipleObjects(ss::AbstractUnitStateChart,obn::Symbol=:P)
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

function StateChartCset_SingleObject(ss::AbstractUnitStateChart,obn::Symbol=:P)
    attrsType = dnames(ss)
    schema = StateChartABMSchema_SingleObject(ss,obn)
    tp_asgn = Dict(zip(attrsType,repeat([Symbol], length(attrsType))))
    acset = AnonACSet(schema;type_assignment=tp_asgn)
    acset
end

StateChartCset_MultipleObjects(ss::AbstractUnitStateChart,obn::Symbol=:P) = AnonACSet(StateChartABMSchema_MultipleObjects(ss,obn))  

## create a representable of single object schema
## this is an instance with one person, and the attributes values given by a vecoter of pairs e.g.,[:Inf=>:R,:Healthy=>:H] (Attr=>state). Then the rest attributes
## are not included in this vector are all given variables.
function representable_SingleObject(ss::AbstractUnitStateChart,as=[],obn::Symbol=:P)
    rep = StateChartCset_SingleObject(ss,obn)
    as = vectorify(as)
    println(as)
    # add the attributes have vaules
    [add_part!(rep,obn;x[1]=>x[2]) for x in as]
    # add the attributes without values
    asall = attrs(acset_schema(rep);just_names=true)
    asvar = asall[asall .∉ Ref(map(first,as))] # ∉ is the unicode of \notin
    # based on the schema design logic, the Attr name is the cocatenate of the codomain (obn) and domain (AttrType), so we can get the AttrType simple chop the obn from the Attr
    ats = map(x->Symbol(chopprefix(String(x),String(obn))),asvar)
    map((attrn,attrt)->add_part!(rep,obn;attrn=>AttrVar(add_part!(rep,attrt))),asvar,ats)
    return rep
end

# this function defines rewrite rules generate by single transition
# here i should be []. and it generate only one single person without state
function singleTransitionLIR_SingleObject(ss::AbstractUnitStateChart, l, i, r, obn::Symbol=:P)
    l = vectorify(l)
    i = vectorify(i)
    r = vectorify(r)
    L,R = [representable_SingleObject(ss,x,obn) for x in [l,r]]
    I = representable_SingleObject(ss,i,obn)
    return [L,I,R]
end


function mk_rule(LIR; use_DataMigration::Bool = false, acset=nothing, migration_rule=nothing)
    L,I,R = LIR
    if !use_DataMigration
        return Rule(homomorphism(I,L;monic=true),homomorphism(I,R;monic=true)) 
    else
        Lm,Im,Rm = [migrate(acset, lir, migration_rule) for lir in LIR]
        println(Lm)
        println(Im)
        println(Rm)
        return Rule(homomorphism(Im,Lm;monic=true),homomorphism(Im,Rm;monic=true)) 
    end  
end

# this function defines rewrite rules generate by single transition
# here i should be []. and it generate only one single person without state
function make_rule_SingleObject(ss::AbstractUnitStateChart, l, i, r, obn::Symbol=:P; use_DataMigration::Bool = false, acset=nothing, migration_rule=nothing)
    LIR = singleTransitionLIR_SingleObject(ss,l,i,r,obn)
    mk_rule(LIR;use_DataMigration = use_DataMigration, acset = acset, migration_rule=migration_rule)   
end

function infectiousruleLIR_SingleObject(ss::AbstractUnitStateChart, l, i, r, obn::Symbol=:P)
    l = vectorify(l)
    i = vectorify(i)
    r = vectorify(r)
    L = coproduct([representable_SingleObject(ss,p,obn) for p in l]) |> apex
    I = coproduct([representable_SingleObject(ss,[],obn),representable_SingleObject(ss,i,obn)]) |> apex
    R = coproduct([representable_SingleObject(ss,p,obn) for p in r]) |> apex
    return [L,I,R] 
end

function make_infectious_rule_SingleObject(ss::AbstractUnitStateChart, l, i, r, obn::Symbol=:P; use_DataMigration::Bool = false, acset=nothing, migration_rule=nothing)
    LIR = infectiousruleLIR_SingleObject(ss, l, i, r, obn)
    mk_rule(LIR;use_DataMigration = use_DataMigration, acset = acset, migration_rule=migration_rule)  
end

## create a representable of multiple object schema
## this is an instance with one person, and each state (may accross multiple state charts) is an object point to the person
## as is a vector of the state the person is in. e.g., as =[:S,:H], which indicates this person is Susceptible and Healthy
## NOTE: if input [], then create a single person without any states
function representable_MultipleObjects(ss::AbstractUnitStateChart,as=[],obn::Symbol=:P)
    rep = StateChartCset_MultipleObjects(ss,obn)
    as = vectorify(as)

    # add the object of person
    add_part!(rep,obn)
    # add the (state) objects to L, I and R
    println(as)
    [add_part!(rep, s; homs_name(s,obn)=>i) for (i,s) in enumerate(as)]

    return rep
end

# return the L, I and R object for single transition rewrite rules
function singleTransitionLIR_MultipleObjects(ss::AbstractUnitStateChart, l, i, r, obn::Symbol=:P)
    l = vectorify(l)
    i = vectorify(i)
    r = vectorify(r)
    return [representable_MultipleObjects(ss,x,obn) for x in [l,i,r]] # array of LIR   
end

# this function defines rewrite rules generate by single transition
# here i should be []. and it generate only one single person without state
function make_rule_MultipleObjects(ss::AbstractUnitStateChart, l, i, r, obn::Symbol=:P; use_DataMigration::Bool = false, acset=nothing, migration_rule=nothing)
    LIR = singleTransitionLIR_MultipleObjects(ss,l,i,r,obn)
    mk_rule(LIR;use_DataMigration = use_DataMigration, acset = acset, migration_rule=migration_rule)     
end

function infectiousruleLIR_MultipleObjects(ss::AbstractUnitStateChart, l, i, r, obn::Symbol=:P)
    l = vectorify(l)
    i = vectorify(i)
    r = vectorify(r)
    L = coproduct([representable_MultipleObjects(ss,p,obn) for p in l]) |> apex
    I = coproduct([representable_MultipleObjects(ss,[],obn),representable_MultipleObjects(ss,i,obn)]) |> apex
    R = coproduct([representable_MultipleObjects(ss,p,obn) for p in r]) |> apex
    return [L,I,R]
end

function make_infectious_rule_MultipleObjects(ss::AbstractUnitStateChart, l, i, r, obn::Symbol=:P; use_DataMigration::Bool = false, acset=nothing, migration_rule=nothing)
    LIR = infectiousruleLIR_MultipleObjects(ss,l,i,r,obn)
    mk_rule(LIR;use_DataMigration = use_DataMigration, acset = acset, migration_rule=migration_rule)  
end

# generate the rewrite rule
# t is the index of a transition
# transitions_rules is an array with pair of (timer, rule) as its elements
function get_rule(ss::AbstractUnitStateChart,t,transitions_rules,obn::Symbol=:P;is_schema_singObject::Bool=true, use_DataMigration::Bool = false, acset=nothing, migration_rule=nothing)
    tt = ttype(ss,t)
    if tt == RLT
        tidx=Int(texpr(ss,t))
        return transitions_rules[tidx][2]
    elseif tt in [DST, CTT]
        # we assume all the states of other state charts are all unkown (set as variables)
        l = sname(ss,tsource(ss,t))
        r = sname(ss,ttarget(ss,t))
        statechart = dname(ss,l)
        attrtype = homs_name(obn,statechart)
        if is_schema_singObject
            if !use_DataMigration
                return make_rule_SingleObject(ss,[attrtype=>l],[],[attrtype=>r],obn)
            else 
                return make_rule_SingleObject(ss,[attrtype=>l],[],[attrtype=>r],obn; use_DataMigration = true, acset = acset, migration_rule=migration_rule)
            end
        else # multiple object schema
            if !use_DataMigration
                return make_rule_MultipleObjects(ss,[l],[],[r],obn)
            else 
                return make_rule_MultipleObjects(ss,[l],[],[r],obn; use_DataMigration = true, acset = acset, migration_rule=migration_rule)
            end  
        end          
    else
        tn = tname(ss,t)
        throw("transition:$tn does not have valid transition type!")
    end
end

function get_timer(ss::AbstractUnitStateChart, transitions_rules,t)
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

function make_ABM(ss::AbstractUnitStateChart, transitions_rules, obn=:P;is_schema_singObject::Bool=true, use_DataMigration::Bool = false, acset=nothing, migration_rule=nothing) 
    return ABM(map(parts(ss, :T)) do t
        tn = tname(ss, t)
        ABMRule(tn, get_rule(ss, t,transitions_rules,obn,is_schema_singObject=is_schema_singObject, use_DataMigration=use_DataMigration, acset = acset, migration_rule=migration_rule), get_timer(ss,transitions_rules,t))
    end, [])
end
   
end

#############################################################################################
######################### end of the dependency function part ###############################
