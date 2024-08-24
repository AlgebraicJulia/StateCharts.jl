# this moduel defines the basic networks
module Networks

export SchDirectedNetwork, SchUndirectedNetwork, SchUndirectedReflectiveNetwork, SchPVArrow, SchPVSpan,
       smallworldNetWork, Graph, Graph_init, res_state

using Catlab
using Catlab: codom, right
using Catlab.Graphics.Graphviz
import Graphs as NormalGraphs
import ..Visualization: Graph
import Base: *

*(x::Symbol, y::Symbol) = Symbol(String(x)*String(y))

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
    return nw
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

get_state(g, states, v, obn=:V) = begin
    for s in states
        if length(incident(g, v, s*obn)) > 0
            return s
        end
    end
end

get_state_color(g, statesColors, v, obn=:V) = begin
    states = keys(statesColors) |> collect
    return statesColors[get_state(g, states, v, obn)]
end

def_node(g, statesColors, v, obn=:V) = ("n$v", Attributes(:label=>"$(subpart(g,v,obn*:ID))",
                                     :shape=>"circle", 
                                     :color=>"black", 
                                     :style=>"filled", 
                                     :fillcolor=>get_state_color(g, statesColors, v, obn),
                                     :width=>"0.25",
                                     :fontsize=>"8pt",
                                     :fixedsize => "true"))
def_edge(g, s, t) =  ([s, t],Attributes(:arrowhead => "none",
                                        :weight => "3.0"))

subGraph(g, n, statesColors::Dict{Symbol, String}, obn::Symbol=:V) = begin
    
    stmts = Graphviz.Statement[]

    for v in 1:nparts(g,obn)
        push!(stmts, Node(def_node(g,statesColors,v,obn)...))
    end

    for k in 1:nparts(g,:E)
        if !(k ∈ refl(g)) && (k <= inv(g, k))
            s = "n$(src(g,k))"
            t = "n$(tgt(g,k))"
            push!(stmts,Edge(def_edge(g,s,t)...))
        end
    end

    graph_attrs = Attributes(:rankdir=>"LR", :rank => "same", 
                             :label=>"time step = $(string(n))", :labelloc=>"t")
    
    gr = Graphviz.Digraph("G", stmts; graph_attrs=graph_attrs, prog="neato")

    return gr
end

Graph_init(g, statesColors::Dict{Symbol, String}, obn::Symbol=:V) = subGraph(g, 0, statesColors, obn)

Graph(res, n::Int,  statesColors::Dict{Symbol, String}, obn::Symbol=:V) = subGraph(res_state(res,n), n, statesColors, obn)

# algorithms generating graph node positions
# the algorithm plans to be cited from the paper of "Efficient ahd High Quality Force-Directed Graph Drawing": http://yifanhu.net/PUB/graph_draw_small.pdf
# in the future, probably we can use the julia package "NetworkLayout": "https://github.com/JuliaGraphs/NetworkLayout.jl/tree/master" to generate the positions
# of nodes for different layouts needs.



end