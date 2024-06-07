module Visualization

export Graph

using Catlab.Graphics.Graphviz
import Catlab.Graphics.Graphviz: Graph, Subgraph
using Catlab.Graphics
import Base.Iterators: flatten

using ..StateChartsSchema

# visualization
TTYPE=Dict(:Rate=>"📉", :TimeOut=>"⏰", :Conditional=>"❓", :Pattern=>"🚩")
# for Pattern, it indicates pattern match. User will provide the pattern match rule, and the state chart schema structure only
# capture the index of the array of the user-defined pattern match rule.

def_state(ss, s) = ("s$s", Attributes(:label=>"$(sname(ss, s))",
                                     :shape=>"square", 
                                     :color=>"black", 
                                     :style=>"filled", 
                                     :fillcolor=>"gold2"))

def_start(ss, se) = ("se$se", Attributes(:label=>"$(ename(ss, se))",
                                     :shape=>"point"))

def_branchNode(ss, b, isStartTransition=true) = (isStartTransition ? "eb$b" : "b$b", Attributes(:label=>"",
                                                      :shape=>"diamond",
                                                      :height=>"0.2",
                                                      :width=>"0.2",
                                                      :color=>"black",
                                                      :fillcolor=>"white"))

def_transition(ss, s, t, f, ttype) = ([s, t],Attributes(:taillabel=>ttype, 
                                                        :label=>"$(tname(ss,f))",
                                                        :taillabelfontsize=>"4",
                                                        :labeldistance=>"1.7",
                                                        :labelangle=>"0"))

# def_start_transition(ss, s, t) =  ([s, t],Attributes(:label=>""))

def_default_branch(ss, s, t) =  ([s, t],Attributes(:label=>"",
                                                   :style=>"dashed"))

def_alternative(ss, s, t, e, isStartTransition=true) =  ([s, t], isStartTransition ? Attributes(:label=>"$(eaname(ss,e))") : Attributes(:label=>"$(aname(ss,e))"))

function Graph(ss::AbstractUnitStateChart)

    stateNodes = [Node(def_state(ss,s)...) for s in 1:ns(ss)]
    branchNodes = [Node(def_branchNode(ss,b,false)...) for b in unique(alsources(ss))]

    smts_Noes=vcat(stateNodes,branchNodes)

    edges_T=map(1:nt(ss)) do k
        state_idx_outfrom = tsource(ss,k)
        state_idx_into = ttarget(ss,k)
        t_type=ttype(ss,k)
        if k in alsources(ss)
            branch_idx_into = k
            [Edge(def_transition(ss,"s$state_idx_outfrom","b$branch_idx_into", k, TTYPE[t_type])...),
             Edge(def_default_branch(ss,"b$branch_idx_into","s$state_idx_into")...)]
        else
            [Edge(def_transition(ss,"s$state_idx_outfrom", "s$state_idx_into", k, TTYPE[t_type])...)]
        end
        end |> flatten |> collect

    edges_alternatives = map(1:nal(ss)) do k
        start_branch = alsource(ss,k)
        state_idx_into = altarget(ss,k)
        [Edge(def_alternative(ss,"b$start_branch","s$state_idx_into", k, false)...)]
        end |> flatten |> collect

    smts_edges=vcat(edges_T,edges_alternatives)

    smts=vcat(smts_Noes,smts_edges)

    graph_attrs = Attributes(:rankdir=>"TB")
    
    g = Graphviz.Digraph("G", smts; graph_attrs=graph_attrs)

    return g
end

end