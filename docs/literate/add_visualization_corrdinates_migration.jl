using Makie, CairoMakie
using StateCharts
using AlgebraicABMs
using Catlab

using NetworkLayout

using Catlab.Graphics.Graphviz
using DataMigrations

using AlgebraicRewriting

# TODO: add positions of nodes and have the positions stable
# https://stackoverflow.com/questions/13938770/how-to-get-the-coordinates-from-layout-from-graphviz
# solution of get the positions of nodes of graphviz

# TODO:
# https://stackoverflow.com/questions/5343899/how-to-force-node-position-x-and-y-in-graphviz
# solution of plot graphviz with the node positions

ENV["JULIA_DEBUG"] = "AlgebraicABMs"; 

# SIR Agent-based model
# This model is the same  as the AnyLogic example model "SIR Agent-Based Networks.". 
# This model currently support three types of networks:
#   1. Random (if the parameter "p_random_connect" = 1.0)
#   2. Ring lattice (if the parameter "p_random_connect" = 0.0)
#   3. Small world (if the parameter "p_random_connect":  0.0 < p_random_connect < 1.0 )
#       it indicates that there are "p_random_connect" percentage connections are radom connections, and 
#       ( 1 - "p_random_connect" ) percentage connections are connections are Ring lattice connections.
# This example model in Anylogic can be found: https://cloud.anylogic.com/model/7088e817-1dab-42b1-89fe-b19bc0a823e1?mode=SETTINGS

# Step 1: define parameters
# parameters of state transition
total_population = 5 # Unit: persons
frac_initial_infected = 0.2
contact_rate = 5 # Unit: contacts per day
infectivity = 0.01
average_illness_duration = 15 # days
# parameters define the network
average_connections = 2
p_random_connect = 0.9

# Step 2: create the StateChart

# function UnitStateChartF takes in 3 input argummapents:
# states: tuple/array of all states
# transitions: tuple/array of all transitions. The syntax of each transition:
#              transition_name => (source_state => target_state, tranisiton_type => values)
# alternative transitions: tuple/array of all alternative transitions. The syntax of each altanative transition:
#              :source_transition_name=>((alternative_transition_name, probability_value) => target_state)
states = [:S, :I, :R] # define the states
transitions = [ :Infection => (:S => :I, :Pattern => 1), :Recovery => (:I => :R, :Rate => 1.0 / average_illness_duration)]
alternatives = [] # this model does not include alternative transitions
# 2.1 create the Statechart
SIRStatechart = UnitStateChartF(states, transitions, alternatives)
# 2.2 Visualization of the Statechart
stateColors = Dict(:S => "green",:I => "red",:R => "gray"); # argument define colors of each state in StateChart
StateCharts.Graph(SIRStatechart, stateColors = stateColors)

# Step 3: create the model schema
### generate the model schema by composing auto state chart schema and the network schema

# 3.1 define the persons object named :V, because we plan to compose the schema_statechart with graph schema by identifying ":V"
schema_statechart = StateChartABMSchema_MultipleObjects(SIRStatechart,:V) |> schemaACSet |> schemaPresent
# 3.2 compose the state chart schema with the network schema (symmetric reflective graph) as the ABM model schema (without equations)
schema_model_without_equations = compose(Open([:S],schemaACSet(schema_statechart),[:V]), Open([:V],schemaACSet(SchUndirectedReflectiveNetwork),[:E])) |> apex |> schemaPresent 
# 3.3 add the composition equations of the model schema, since those equations disapper after composition. This would be fixed in the future using GATlab
@present schema_model <: schema_model_without_equations begin
    compose(inv,inv) == id(E)
    compose(inv,src) == tgt
    compose(inv,tgt) == src
    compose(refl, src) == id(V)
    compose(refl, tgt) == id(V)
    compose(refl, inv) == refl 
end
@acset_type SIRModelStaticNet(schema_model, index=[:src,:tgt]){Symbol} <: AbstractSymmetricReflexiveGraph
# show the model schema
to_graphviz(schema_model; prog="dot")
# 3.4 create the model schema with nodes' positions. This schema is only used for model results visualization
@present schema_model_view <: schema_model begin
    Pos::AttrType
    VPos::Attr(V, Pos)
end
@acset_type SIRModelStaticNetView(schema_model_view, index=[:src,:tgt]){Symbol,Tuple{Float32,Float32}}
# show this schema
to_graphviz(schema_model_view; prog="dot")

addPos = Migrate(schema_model, SIRModelStaticNet, schema_model_view, SIRModelStaticNetView; delta=false)

# Step 4: Data Migration from the ACSet with schema -- schema_statechart to the ACSet with schema -- schema_model
# 4.1 define the data migration rule to automatically migrate ACSets from pure state chart schema to the composed model shema    
const migrate_rule = @migration schema_model schema_statechart begin
    V => V; 
    ID=>ID; VID => VID
    S => S; SV => SV
    I => I; IV => IV
    R => R; RV => RV

    E => @product begin
        p1::V
        p2::V
    end
    src => p1
    tgt => p2
    inv => begin
        p1 => p2
        p2 => p1
    end
end

# 4.2 define the user_defined rewrite rules: rewrite rules for contacts of Infectives
# Note that each transitions_rule has a pair: timer=>rule
transition_rules = [ ContinuousHazard( 1.0 / (contact_rate * infectivity)) => make_infectious_rule_MultipleObjects(SIRStatechart, [[:S],[:I]],[:I],[[:I],[:I]], :V; use_DataMigration = true, acset = SIRModelStaticNet, migration_rule = migrate_rule)]

# Step 5: define the network
network = smallworldNetWork(Int(total_population), average_connections, p_random_connect);

# Step 6: Initialization
# In the initialization, randomly assign "total_population * frac_initial_infected" persons as Infectives (state I), and the rest persons are all Susceptibles (state S)
init = radomlyAssignInitialInfectives(SIRModelStaticNet(), Int(total_population), Int(frac_initial_infected * total_population); obn=:V, network = network)
# show the networks with the initial model state
init_graphviz = StateCharts.Graph_init(init,stateColors)

# add positions
init_positions=addPos(init)

# Step 7: run the ABM model
abm = make_ABM(SIRStatechart, transition_rules, :V; is_schema_singObject=false, use_DataMigration=true, acset = SIRModelStaticNet, migration_rule = migrate_rule)
abm_positions = addPos(abm)

########################################################################################
#### stable layout of network visualization

using SparseArrays


# Note: the reflective edges are not included, since we will not plot out the refl edges
function ajacent_matrix(g::AbstractSymmetricReflexiveGraph)
    n = nv(g)

    amtx = spzeros(n, n)

    for k in 1:nparts(g,:E)
        if !(k ∈ refl(g)) && (k <= inv(g, k))
            s = src(g,k)
            t = tgt(g,k)
            amtx[s,t]=1
            amtx[t,s]=1
        end
    end
    return amtx

end

adj_matrix=ajacent_matrix(init)
positions = @time SFDP(; dim=2, Ptype=Float32, tol=0.1, K=1)(adj_matrix)
n = length(positions)

using GeometryBasics.GeoInterface

pos=[(GeoInterface.coordinates(positions[i])[1],GeoInterface.coordinates(positions[i])[2]) for i in 1:n]

set_subpart!(init_positions, 1:Int(total_population), :VPos, pos)
###############################################################################

res = run!(abm_positions, init_positions; maxtime=10.0); # solutions with stable positions
res = run!(abm, init; maxtime=10.0); # solutions withiout stable positions

# Step 8: Visualization of results
# 8.1 plot out the time series of each state
Makie.plot(res; Dict(o=>X->nparts(X,o) for o in states)...)
# 8.2 show the networks at time t
t = 15
StateCharts.Graph(res,t,stateColors)
