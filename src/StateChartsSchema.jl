""" Some description of ths package
"""
module StateChartsSchema

export UnitStateChartF, TheoryUnitStateChart, AbstractUnitStateChart, UnitStateChart,
StateChartsF, TheoryStateCharts, AbstractStateCharts, StateCharts,
ns, nt, nal, sname, snames, tname, tnames, aname, anames, texpr, dname, dnames,
ttype, tsource, ttarget, etarget, alsource, alsources, altarget, altargets, spartition

using Catlab

vectorify(n::Vector) = collect(n)
vectorify(n::Tuple) = length(n) == 1 ? [n] : collect(n)
vectorify(n::SubArray) = collect(n)
vectorify(n) = [n]

state_dict(n) = Dict(s=>i for (i, s) in enumerate(n))

# define the schema of unit state charts
@present TheoryUnitStateChart(FreeSchema) begin
    S::Ob #state
    T::Ob #transitions
    ALT::Ob #alternative branches

    ts::Hom(T,S)
    tt::Hom(T,S)
    abs::Hom(ALT,T)
    abt::Hom(ALT,S)
    
    Name::AttrType
    Expr::AttrType
    TT::AttrType # the transition's type
    
    texpr::Attr(T, Expr)
    aexpr::Attr(ALT, Expr)
    aname::Attr(ALT, Name)
    sname::Attr(S, Name)
    tname::Attr(T, Name)
    ttype::Attr(T, TT)
end

@abstract_acset_type AbstractUnitStateChart
@acset_type UnitStateChartUntype(TheoryUnitStateChart) <: AbstractUnitStateChart
const UnitStateChart = UnitStateChartUntype{Symbol,Float16,Symbol}

# define functions of adding components of state charts
add_transition!(ss::AbstractUnitStateChart,s,t;kw...) = add_part!(ss,:T;ts=s,tt=t,kw...)
add_transitions!(ss::AbstractUnitStateChart,s,t,n;kw...) = add_parts!(ss,:T,n;ts=s,tt=t,kw...)

add_state!(ss::AbstractUnitStateChart;kw...) = add_part!(ss,:S;kw...) 
add_states!(ss::AbstractUnitStateChart,n;kw...) = add_parts!(ss,:S,n;kw...)

add_alternative!(ss::AbstractUnitStateChart,s,t;kw...) = add_part!(ss,:ALT;abs=s,abt=t,kw...)
add_alternatives!(ss::AbstractUnitStateChart,s,t,n;kw...) = add_parts!(ss,:ALT,n;abs=s,abt=t,kw...)

""" UnitStateChartF(s::Tuple, t::Tuple, alt::Tuple)

Returns an instance of unit State Chart
"""
UnitStateChartF(s, t, alt) = begin
    
    ss = UnitStateChart()

    # add states
    states=vectorify(s)
    add_states!(ss,length(states),sname=states)

    # add transitions
    transitions=vectorify(t)
    tns=map(first,transitions)
    tst=map(first,map(last,transitions))
    tnr=map(last,map(last,transitions))
    s_idx=state_dict(states)
    t_idx=state_dict(tns)
    ts=map(x->s_idx[x], map(first,tst))
    tt=map(x->s_idx[x], map(last,tst))
    add_transitions!(ss,ts,tt,length(tns),tname=tns,ttype=map(first,tnr),texpr=map(last,tnr))

    # add alternative branches of transitions
    if length(alt)>0 # if have branches
        for (i,(tn,rps)) in enumerate(vectorify(alt))
            rps = vectorify(rps)
            t = t_idx[tn]
            ps = map(first,rps)
            an = map(first,ps)
            p = map(last,ps)
            sns = map(last,rps)
            s = map(x->s_idx[x],sns)
            add_alternatives!(ss,repeat([t], length(s)),s,length(s),aexpr=p,aname=an)
        end
    end

   ss
end

ns(ss::AbstractUnitStateChart) = nparts(ss,:S)
nt(ss::AbstractUnitStateChart) = nparts(ss,:T)
nal(ss::AbstractUnitStateChart) = nparts(ss,:ALT) # return the numbers of branch alternatives -- one branch counts to one alternative

sname(ss::AbstractUnitStateChart,s) = subpart(ss,s,:sname) 
snames(ss::AbstractUnitStateChart) = [sname(ss, s) for s in 1:ns(ss)]

tname(ss::AbstractUnitStateChart,t) = subpart(ss,t,:tname) 
tnames(ss::AbstractUnitStateChart) = [tname(ss, t) for t in 1:nt(ss)]

aname(ss::AbstractUnitStateChart,al) = subpart(ss,al,:aname) 
anames(ss::AbstractUnitStateChart) = [aname(ss, al) for al in 1:nal(ss)]

ttype(ss::AbstractUnitStateChart,t) = subpart(ss,t,:ttype) 

# return a transition's source state
tsource(ss::AbstractUnitStateChart,t) = subpart(ss,t,:ts) 
# return a transition's target state
ttarget(ss::AbstractUnitStateChart,t) = subpart(ss,t,:tt) 
# return an alternative's source transition
alsource(ss::AbstractUnitStateChart,al) = subpart(ss,al,:abs) 
alsources(ss::AbstractUnitStateChart) = [alsource(ss, al) for al in 1:nal(ss)]
# return an alternative's target state
altarget(ss::AbstractUnitStateChart,al) = subpart(ss,al,:abt) 
altargets(ss::AbstractUnitStateChart) = [altarget(ss, al) for al in 1:nal(ss)]

texpr(ss::AbstractUnitStateChart,t) = subpart(ss,t,:texpr) 


# define the schema of multiple state charts
@present TheoryStateCharts <: TheoryUnitStateChart begin

    D::Ob # Partitions of each unit state chart
    sd::Hom(S,D) # map states to their own unit state chart    
    dname::Attr(D, Name)

end

@abstract_acset_type AbstractStateCharts <: AbstractUnitStateChart
@acset_type StateChartsUntype(TheoryStateCharts) <: AbstractStateCharts
const StateCharts = StateChartsUntype{Symbol,Float16,Symbol}

add_state!(ss::AbstractStateCharts,d;kw...) = add_part!(ss,:S;sd=d,kw...) 
add_states!(ss::AbstractStateCharts,n,d;kw...) = add_parts!(ss,:S,n;sd=d,kw...)

add_partition!(ss::AbstractStateCharts;kw...) = add_part!(ss,:D;kw...) 
add_partitions!(ss::AbstractStateCharts,n;kw...) = add_parts!(ss,:D,n;kw...) 

# s = (partitionname1=>(set of state names of this partition),
#      partitionname2=>(set of state names of this partition))
StateChartsF(s, t, alt) = begin
    
    ss = StateCharts()

    s = vectorify(s)
    # add states
    for (i,(d,states)) in enumerate(s)
        states = vectorify(states)
        add_partition!(ss,dname=d)
        add_states!(ss,length(states),repeat([i],length(states)),sname=states)
    end

    # add transitions
    transitions=vectorify(t)
    tns=map(first,transitions)
    tst=map(first,map(last,transitions))
    tnr=map(last,map(last,transitions))
    s_idx=state_dict(snames(ss))
    ts=map(x->s_idx[x], map(first,tst))
    tt=map(x->s_idx[x], map(last,tst))
    add_transitions!(ss,ts,tt,length(tns),tname=tns,ttype=map(first,tnr),texpr=map(last,tnr))

    # add alternative branches of transitions
    if length(alt)>0 # if have branches
        for (i,(tn,rps)) in enumerate(vectorify(alt))
            rps = vectorify(rps)
            t = t_idx[tn]
            ps = map(first,rps)
            an = map(first,ps)
            p = map(last,ps)
            sns = map(last,rps)
            s = map(x->s_idx[x],sns)
            add_alternatives!(ss,repeat([t], length(s)),s,length(s),aexpr=p,aname=an)
        end
    end

   ss
end

nd(ss::AbstractStateCharts) = nparts(ss,:D)

spartition(ss::AbstractStateCharts,s) = subpart(ss,s,:sd)

dname(ss::AbstractStateCharts,d) = subpart(ss,d,:dname) 
dnames(ss::AbstractStateCharts) = [dname(ss, d) for d in 1:nd(ss)]

# hard coded a state chart name for unit state chart, is to adapt to convieniently use
# in the model running
dname(ss::AbstractUnitStateChart,d) = :UnitStateChart
dnames(ss::AbstractUnitStateChart) = [:UnitStateChart]

# schema having the hazard rates of transitions can depend on states. And the logic of the expression of the hazard rates of transitions
# are captured by GatExpr. Object representing Links from states to transitions is also included in this schema.
#@present TheoryLinkedStateCharts <: TheoryStateCharts begin

#    L::Ob # links from objects to transitions
    
#    lsrc::Hom(L,S)    
#    ltgt::Attr(L, T)
#end

#@abstract_acset_type AbstractLinkedStateCharts <: AbstractStateCharts
#@acset_type LinkedStateChartsUntype(TheoryLinkedStateCharts) <: AbstractLinkedStateCharts
#const LinkedStateCharts = LinkedStateChartsUntype{Symbol,GATExpr{T},Symbol} # TODO: replace {T} to concrete expressions





end



