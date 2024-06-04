# this moduel defines the basic networks
module Networks

export SchDirectedNetwork, SchUndirectedNetwork, SchUndirectedReflectiveNetwork, SchPVArrow, SchPVSpan

using Catlab

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

end