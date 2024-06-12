# this moduel defines the basic networks
module Networks

export SchDirectedNetwork, SchUndirectedNetwork, SchUndirectedReflectiveNetwork, SchPVArrow, SchPVSpan,
       smallworldNetWork, Graph, res_groups

#using Catlab.Graphs.BasicGraphs
using Catlab
using Catlab: codom, right
import Graphs as NormalGraphs
using GraphPlot
import ..Visualization: Graph

# using basic graph schemas as the schema of basic networks
SchDirectedNetwork = SchGraph
# # The Category of Graphs
#
# The Theory of Graphs is given by the following Schema:
# ```julia
# @present SchGraph(FreeSchema) begin
#   V::Ob
#   E::Ob
#   src::Hom(E,V)
#   tgt::Hom(E,V)
# end

SchUndirectedNetwork = SchSymmetricGraph 


# # Symmetric graphs
#
# @present SchSymmetricGraph <: SchGraph begin
#  inv::Hom(E,E)
#  compose(inv,inv) == id(E)
#  compose(inv,src) == tgt
#  compose(inv,tgt) == src
#end

SchUndirectedReflectiveNetwork = SchSymmetricReflexiveGraph

# Symmetric reflexive graphs
#
# @present SchSymmetricReflexiveGraph <: SchSymmetricGraph begin
#  refl::Hom(V,E)
#  compose(refl, src) == id(V)
#  compose(refl, tgt) == id(V)
#  compose(refl, inv) == refl # Reflexive loop fixed by involution.
#end

# some interface schemas 
# schema assign each person a vertice: P -> V
@present SchPVArrow(FreeSchema) begin
    (P,V)::Ob
    PV::Hom(P,V)    
end

# schema represents the relationships between P and V: P monic<- Conn -> V
@present SchPVSpan(FreeSchema) begin
    (P,V,Conn)::Ob
    connP::Hom(Conn,P)
    connV::Hom(Conn,V)       
end

# Create network instances
# n is the number of the total vertices (usually the population), one person is in one vertex
# d is the average degree of vertices
# p is the probablity of random connection
# return a small-world network: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.60.7332
smallworldNetWork(n, d, p) = begin
    g = NormalGraphs.newman_watts_strogatz(n, d, p)
    ges = collect(NormalGraphs.edges(g))

    nw = SymmetricReflexiveGraph(NormalGraphs.nv(g))
    [Catlab.add_edge!(nw, NormalGraphs.src(ges[i]), NormalGraphs.dst(ges[i])) for i in 1:length(ges)]
    return g=>nw
end

##### visualizations #####
# return each time's model result's instance
function res_state(res,t)
    t = Int(t)
    if t == 0
        state = res.init
    else
        state = codom.(right.(res.hist))[t]
    end
    return state
end

# t: time step
# s: symbol of state
# obn: symbol, the object name of persons

# return the list of persons of each state (s)
function res_state_agents(res,t,s,obn)
    states = res_state(res,t)
    vs = [subpart(states,n,s*obn) for n in 1:nparts(states,s)]
    return vs
end

# return the list of each person's ID
res_state_agents_id(res,t,s,obn) = [subpart(states,p,obn*:ID) for p in res_state_agents(res,t,s,obn)]

function statecolor(colorgroups, colors, s)
    idx_colors = Dict(c=>i for (i, c) in enumerate(colors))
    return idx_colors[colorgroups[s]]
end

# return the input argument of "membership" for Graph function
# res: the results of run!()
# t: time step 
# states: the states of the model
# pop: total population
# colors: the set of colors
# colorgroups: the Dict of state=>color
# obn: the symbol of the object of persons
function res_groups(res,t,states,pop,colors,colorgroups,obn)
    groups = zeros(Int64,Int(pop))
    for s in states
        for v in res_state_agents(res,t,s,obn)
            groups[v] = Int(statecolor(colorgroups, colors, s))
        end
    end
    groups
end

## functions to plot out the graphs
Graph(g::NormalGraphs.SimpleGraph, res, t, membership, nodecolor, nodesizescale) = begin    
    nodefillc = nodecolor[membership]
    states = res_state(res,t)
    nodelabel = [subpart(states,p,obn*:ID) for p in collect(1:NormalGraphs.nv(g))] # labelled by the ID
    gplot(g, nodefillc=nodefillc, nodelabel=nodelabel, nodelabelsize=ones(NormalGraphs.nv(g))*nodesizescale)
end

Graph(g::HasGraph) = to_graphviz(g)

end