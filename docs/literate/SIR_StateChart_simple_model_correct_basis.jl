using Makie, CairoMakie
using StateCharts
using AlgebraicABMs
using Catlab

ENV["JULIA_DEBUG"] = "AlgebraicABMs"; 

# SIR Agent-based model
# This model is a simplified version of the AnyLogic example model "SIR Agent-Based Networks." It only includes the statechart of the state transitions among S, I, and R, and does not include the networks.
# This example model in Anylogic can be found: https://cloud.anylogic.com/model/7088e817-1dab-42b1-89fe-b19bc0a823e1?mode=SETTINGS

# Step 1: define parameters

total_population = 4 # Unit: persons
frac_initial_infected = 1.0/4.0
contact_rate = 5 # Unit: contacts per day
infectivity = 0.01
average_illness_duration = 15 # days

# Step 2: create the StateChart

# function UnitStateChartF takes in 3 input arguments:
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
# 2.3 define the user_defined rewrite rules: rewrite rules for contacts of Infectives
# Note that each transitions_rule has a pair: timer=>rule (if no basis)
#                                         or: timer=> ( rule => basis ) (if has basis: L' -> L)

# l = i = [[],:I], indicates a pair of persons, one generic person, and one infective
# r = [:I,:I], indicates two infectives
# lp = [:I], indicates an infective person
# ac = [:S, :I], indicates a positive condition of a susceptible person and an infective person
transition_rules = [ ContinuousHazard( 1.0 / (contact_rate * infectivity)) => make_infectious_rule_MultipleObjects(SIRStatechart, [[],:I],[[],:I],[:I,:I],[:I],[:S,:I],:P)]

# Plot out the model schema
# Note: this step is not nessesary. Here it is only used to show the model schema
schema_statechart = StateChartABMSchema_MultipleObjects(SIRStatechart,:P)
to_graphviz(schema_statechart |> schemaACSet |> schemaPresent; prog="dot")

# Step 3: Initialization
# In the initialization, randomly assign "total_population * frac_initial_infected" persons as Infectives (state I), and the rest persons are all Susceptibles (state S)
init = radomlyAssignInitialInfectives(StateChartCset_MultipleObjects(SIRStatechart), Int(total_population), Int(frac_initial_infected * total_population))

# Step 4: run the ABM model
res = run!(make_ABM(SIRStatechart,transition_rules,is_schema_singObject=false), init; maxtime=8);

# Step 5: plot out the results of each state
Makie.plot(res; Dict(o=>X->nparts(X,o) for o in states)...)



