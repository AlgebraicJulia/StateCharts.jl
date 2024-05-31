""" Some description of ths package
"""
module StateChartsSchema

export StateChartF, TheoryStateChart, AbstractStateChart, StateChart,
ns, nt, ne, neal, nal, sname, snames, tname, tnames, ename, enames, eaname, eanames, aname, anames, texpr,
ttype, ealsource, ealsources, ealtarget, ealtargets, tsource, ttarget, etarget, alsource, alsources, altarget, altargets

using Catlab

vectorify(n::Vector) = collect(n)
vectorify(n::Tuple) = length(n) == 1 ? [n] : collect(n)
vectorify(n::SubArray) = collect(n)
vectorify(n) = [n]

state_dict(n) = Dict(s=>i for (i, s) in enumerate(n))

# define the schema of state charts
@present TheoryStateChart(FreeSchema) begin
    S::Ob #state
    T::Ob #transitions
    ALT::Ob #alternative branches
    E::Ob # start transition (in this statechart schema, it can have multiple start transitions -- with multiple state chart in parallel)
    EALT::Ob # alternative branches for the start transition
    
    ts::Hom(T,S)
    tt::Hom(T,S)
    abs::Hom(ALT,T)
    abt::Hom(ALT,S)
    es::Hom(E,S)
    eabs::Hom(EALT,E)
    eabt::Hom(EALT,S)
    
    Name::AttrType
    Expr::AttrType
    TT::AttrType # type of the non-start transitions
    
    ttype::Attr(T, TT)
    texpr::Attr(T, Expr)
    aexpr::Attr(ALT, Expr)
    eaexpr::Attr(EALT, Expr)
    aname::Attr(ALT, Name)
    eaname::Attr(EALT, Name)
    sname::Attr(S, Name)
    tname::Attr(T, Name)
    ename::Attr(E, Name)
end

@abstract_acset_type AbstractStateChart
@acset_type StateChartUntype(TheoryStateChart) <: AbstractStateChart
#const StateChart = StateChartUntype{Symbol,Union{Float16,Sampleable},Symbol}
const StateChart = StateChartUntype{Symbol,Float16,Symbol}

# define functions of adding components of state charts
add_transition!(ss::AbstractStateChart,s,t;kw...) = add_part!(ss,:T;ts=s,tt=t,kw...)
add_transitions!(ss::AbstractStateChart,s,t,n;kw...) = add_parts!(ss,:T,n;ts=s,tt=t,kw...)

add_state!(ss::AbstractStateChart;kw...) = add_part!(ss,:S;kw...) 
add_states!(ss::AbstractStateChart,n;kw...) = add_parts!(ss,:S,n;kw...)

add_alternative!(ss::AbstractStateChart,s,t;kw...) = add_part!(ss,:ALT;abs=s,abt=t,kw...)
add_alternatives!(ss::AbstractStateChart,s,t,n;kw...) = add_parts!(ss,:ALT,n;abs=s,abt=t,kw...)

add_startTransition!(ss::AbstractStateChart,t;kw...) = add_part!(ss,:E;es=t,kw...)
add_startTransitions!(ss::AbstractStateChart,t,n;kw...) = add_parts!(ss,:E,n;es=t,kw...)

add_startAlternative!(ss::AbstractStateChart,s,t;kw...) = add_part!(ss,:EALT;eabs=s,eabt=t,kw...)
add_startAlternatives!(ss::AbstractStateChart,s,t,n;kw...) = add_parts!(ss,:EALT,n;eabs=s,eabt=t,kw...)

""" StateCharts(name::String)
StateCharts(s::Tuple, t::Tuple, alt::Tuple, e::Tuple, ealt::Tuple)

Returns an instance of State Charts
"""
StateChartF(s, t, alt, e, ealt) = begin
    
    ss = StateChart()

    # add states
    states=vectorify(s)
    add_states!(ss,length(states),sname=states)

    # add non-start transitions
    transitions=vectorify(t)
    tns=map(first,transitions)
    tst=map(first,map(last,transitions))
    tnr=map(last,map(last,transitions))
    s_idx=state_dict(states)
    t_idx=state_dict(tns)
    ts=map(x->s_idx[x], map(first,tst))
    tt=map(x->s_idx[x], map(last,tst))
    add_transitions!(ss,ts,tt,length(tns),tname=tns,ttype=map(first,tnr),texpr=map(last,tnr))

    # add start transitions
    es=vectorify(e)
    esname=map(first,es)
    esstates=map(last,es)
    e_idx=state_dict(esname)
    tes=map(x->s_idx[x], map(last,es))
    add_startTransitions!(ss,tes,length(esname),ename=esname)

    # add branches of start transitions
    if length(ealt)>0 # if have branches
        for (i,(tn,rps)) in enumerate(vectorify(ealt))
            rps = vectorify(rps)
            t = e_idx[tn]
            ps = map(first,rps)
            an = map(first,ps)
            p = map(last,ps)
            sns = map(last,rps)
            s = map(x->s_idx[x],sns)
            add_startAlternatives!(ss,repeat([t], length(s)),s,length(s),eaexpr=p,eaname=an)
        end
    end

    # add branches of non-start transitions
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

ns(ss::AbstractStateChart) = nparts(ss,:S)
nt(ss::AbstractStateChart) = nparts(ss,:T)
ne(ss::AbstractStateChart) = nparts(ss,:E)
neal(ss::AbstractStateChart) = nparts(ss,:EALT) # return the numbers of start branch alternatives -- one branch counts to one alternative
nal(ss::AbstractStateChart) = nparts(ss,:ALT) # return the numbers of non-start branch alternatives -- one branch counts to one alternative

sname(ss::AbstractStateChart,s) = subpart(ss,s,:sname) 
snames(ss::AbstractStateChart) = [sname(ss, s) for s in 1:ns(ss)]

tname(ss::AbstractStateChart,t) = subpart(ss,t,:tname) 
tnames(ss::AbstractStateChart) = [tname(ss, t) for t in 1:nt(ss)]

ename(ss::AbstractStateChart,e) = subpart(ss,e,:ename) 
enames(ss::AbstractStateChart) = [ename(ss, e) for e in 1:ne(ss)]

eaname(ss::AbstractStateChart,eal) = subpart(ss,eal,:eaname) 
eanames(ss::AbstractStateChart) = [eaname(ss, eal) for eal in 1:neal(ss)]

aname(ss::AbstractStateChart,al) = subpart(ss,al,:aname) 
anames(ss::AbstractStateChart) = [aname(ss, al) for al in 1:nal(ss)]

ttype(ss::AbstractStateChart,t) = subpart(ss,t,:ttype) 

ealsource(ss::AbstractStateChart,eal) = subpart(ss,eal,:eabs)
ealsources(ss::AbstractStateChart) = [ealsource(ss, eal) for eal in 1:neal(ss)]

ealtarget(ss::AbstractStateChart,eal) = subpart(ss,eal,:eabt)
ealtargets(ss::AbstractStateChart) = [ealtarget(ss, eal) for eal in 1:neal(ss)]

# return a transition's source state
tsource(ss::AbstractStateChart,t) = subpart(ss,t,:ts) 
# return a transition's target state
ttarget(ss::AbstractStateChart,t) = subpart(ss,t,:tt) 
# return a start transition's target state
etarget(ss::AbstractStateChart,e) = subpart(ss,e,:es) 
# return an alternative's source transition
alsource(ss::AbstractStateChart,al) = subpart(ss,al,:abs) 
alsources(ss::AbstractStateChart) = [alsource(ss, al) for al in 1:nal(ss)]
# return an alternative's target state
altarget(ss::AbstractStateChart,al) = subpart(ss,al,:abt) 
altargets(ss::AbstractStateChart) = [altarget(ss, al) for al in 1:nal(ss)]

texpr(ss::AbstractStateChart,t) = subpart(ss,t,:texpr) 

end


