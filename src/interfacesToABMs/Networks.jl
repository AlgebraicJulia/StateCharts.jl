# this moduel defines the basic networks
module Networks

export SchDirectedNetwork, SchUndirectedNetwork, SchUndirectedReflectiveNetwork, SchPVArrow, SchPVSpan,
       smallworldNetWork, Graph

#using Catlab.Graphs.BasicGraphs
using Catlab
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

## functions to plot out the graphs
Graph(g::NormalGraphs.SimpleGraph, membership, nodecolor) = begin
    #membership = [1,1,1,1,1,1,1,1,2,1,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2]
    #nodecolor = [colorant"lightseagreen", colorant"red"]
    # membership color
    nodefillc = nodecolor[membership]
    nodelabel = collect(1:NormalGraphs.nv(g))
    gplot(g, nodefillc=nodefillc, nodelabel=nodelabel)
end

Graph(g::HasGraph) = to_graphviz(g)

end