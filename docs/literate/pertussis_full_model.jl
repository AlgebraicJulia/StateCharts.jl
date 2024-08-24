using Makie, CairoMakie
using StateCharts
using AlgebraicABMs
using Catlab
import Base: *

ENV["JULIA_DEBUG"] = "AlgebraicABMs"; 

*(x::Symbol, y::Symbol) = Symbol(String(x)*String(y))

# define the parameter values
totalPopulation=200.0
c = 0.1411 * 10 # calculated by the average value among the 32 age groups in Hethocote model

# the transmission probability seems very high to me. So I temporarily times 0.2 of it. 
β = [0.8,0.4,0.2] * 0.2 # 0.8 for full Infectives I
                  # 0.4 for mild Infectives Im
                  # 0.2 for weak Infectives Iw
p_vaccinationRate = 0.76 / 365.0 # per day, auxiliary value
# The values of "p_incubationPeriodInDay_Im" and "p_incubationPeriodInDay_Iw" refers to the parameter settings in the published paper
# (in Table 1: https://peerj.com/articles/2337/#table-1), instead of the specfic Anylogic model.
p_incubationPeriodInDay_I = 10 # unit: days
p_incubationPeriodInDay_Im = 14 # unit: days
p_incubationPeriodInDay_Iw = 21 # unit: days
p_gam = 1 / 21.0
p_alp = 1 / (5 * 365.0) # waning rate of R
p_tau = 1 / (2 * 365.0) # waning rate of V
p_eps = 1 / (100 * 365.0) # waning rate from R1 or V1 to S

# define helper arrays for subgroups
Incubations = [:Incubation_I, :Incubation_Im, :Incubation_Iw]
Infectives = [:Infective_I, :Infective_Im, :Infective_Iw]
# define the index of transitions due to contacting of infectives
IContacts_idx = Dict(:StoI=>1:3,
                        :V1toI=>4:6,
                        :R1toIm=>7:9,
                        :V2toIm=>10:12,
                        :R2toIw=>13:15,
                        :R3toR4=>16:18,
                        :V3toR4=>19:21)

# define the states
states = [:Susceptible,
# 3 incubation states, with decreasing levels of coming virulence
:Incubation_I,
:Incubation_Im, 
:Incubation_Iw,
# 3 infective states, with decreasing levels of virulence
:Infective_I,
:Infective_Im, 
:Infective_Iw,
# 3 recovered states, with INCREASING levels of protection
:Recovered_1,
:Recovered_2,
:Recovered_3,
:Recovered_4,
# 3 Vaccinated states, with INCREASING levels of protection
:Vaccinated_1,
:Vaccinated_2,
:Vaccinated_3,
:Vaccinated_4]

# create the pertussis state chart

# UnitStateChartF takes in 3 input arguments:
# states: tuple/array of all states
# transitions: tuple/array of all transitions. The syntax of each transition:
#              transition_name => (source_state => target_state, tranisiton_type => values)
# alternative transitions: tuple/array of all alternative transitions. The syntax of each altanative transition:
#              :source_transition_name=>((alternative_transition_name, probability_value) => target_state)
pertussisStatechart = UnitStateChartF(states, # states
    (   
        # we deal here with transitions from Susceptible.
        # first, we deal with the single vaccination transition for a susceptible person
        :t_transition42_vaccinationFromSusceptible=>(:Susceptible=>:Vaccinated_1,:Rate=>p_vaccinationRate),
        # next, we deal successively with the infection transitions to incubation states 
        # for each of full-strength, mild and weak infections in turn
        map((x,y)->:t_S_I_newExposure_FullStrength_ * x =>(:Susceptible=>:Incubation_I,:Pattern=>y), Infectives, IContacts_idx[:StoI])...,
        map((x,y)->:t_V1_I_Vaccinated_1_I_newExposure_FullStrength_ * x =>(:Vaccinated_1=>:Incubation_I,:Pattern=>y), Infectives, IContacts_idx[:V1toI])...,
        map((x,y)->:t_R1_mI_newExposure_mild_ * x =>(:Recovered_1=>:Incubation_Im,:Pattern=>y), Infectives, IContacts_idx[:R1toIm])...,
        map((x,y)->:t_V2_mI_Vaccinated_2_I_newExposure_mild_ * x =>(:Vaccinated_2=>:Incubation_Im,:Pattern=>y), Infectives, IContacts_idx[:V2toIm])...,
        map((x,y)->:t_R2_wI_newExposure_weak_ * x =>(:Recovered_2=>:Incubation_Iw,:Pattern=>y), Infectives, IContacts_idx[:R2toIw])...,       
        # next, we deal successively with the post-incubation transitions to infectivity for each of 
        # full-strength, mild and weak infections in turn
        :t_transition41_becomingInfective_FullStrength=>(:Incubation_I=>:Infective_I,:TimeOut=>p_incubationPeriodInDay_I),        
        :becomingInfective_mild=>(:Incubation_Im=>:Infective_Im,:TimeOut=>p_incubationPeriodInDay_Im),       
        :becomingInfective_weak=>(:Incubation_Iw=>:Infective_Iw,:TimeOut=>p_incubationPeriodInDay_Im),
        # now we deal with infection transitions that are too weak to 
        # trigger incubation & infectivity.  These occur from Recovered_2 and
        # and Vaccinated_3; note that lower levels of both natural & vaccine-induced
        # immunity lead to infectivity, and pathogen exposure induces no change at higher 
        # levels of immunity.
        # First, boosting of immunity from natural recovery
        map((x,y)->:t_R3_R4_exposureBoostingRecovery_3_ * x =>(:Recovered_3=>:Recovered_4,:Pattern=>y), Infectives, IContacts_idx[:R3toR4])...,
        # NB: The model posits that exposure to pathogen in Vaccinated_3 indeed 
        # confers immunity comparable to from natural infection (i.e., to Recovery_4)
        map((x,y)->:t_V3_R4_exposureBoostingWanedVaccinated_3_ * x =>(:Vaccinated_3=>:Recovered_4,:Pattern=>y), Infectives, IContacts_idx[:V3toR4])...,
        # having dealt with infection transitions, 
        # we next deal with recovery (departure from states of infectivity) for
        # each virulence level; all are governed by hazard rate gam and go to Recovered_4 
        :t_gam_I_recoveryFromFullStrengthInfection=>(:Infective_I=>:Recovered_4,:Rate=>p_gam),       
        :t_gam_ml_recoveryFromMildInfection=>(:Infective_Im=>:Recovered_4,:Rate=>p_gam),
        :t_gam_wI_recoveryFromWeakInfection=>(:Infective_Iw=>:Recovered_4,:Rate=>p_gam),
        # transitions involving waning of immunity from recovery states  
        :t_alp_R4_waningImmunityFromRecovered_4=>(:Recovered_4=>:Recovered_3,:Rate=>p_alp),
        :t_alp_R3_waningImmunityFromRecovered_3=>(:Recovered_3=>:Recovered_2,:Rate=>p_alp),        
        :t_alp_R2_waningImmunityFromRecovered_2=>(:Recovered_2=>:Recovered_1,:Rate=>p_alp),
        :t_alp_R1_waningImmunityFromRecovered_1=>(:Recovered_1=>:Susceptible,:Rate=>p_eps),
        # now we deal with the transitions involving receipt of vaccines whilst in a recovered state
        # NB: There is no opportunity to boost immunity through vaccination (or natural exposure) 
        # in the highest state of natural induced immunity (Recovered_4)
        :t_transition6_vaccinationFromRecovered_1=>(:Recovered_1=>:Recovered_2,:Rate=>p_vaccinationRate),
        :t_transition8_vaccinationFromRecovered_2=>(:Recovered_2=>:Recovered_3,:Rate=>p_vaccinationRate),
        :t_R3_R4_vaccinationFromRecovered_3=>(:Recovered_3=>:Recovered_4,:Rate=>p_vaccinationRate),
        # transitions involving waning of immunity from vaccinated states  
        :t_tau_V4_waningImmunityFromVaccinated_4=>(:Vaccinated_4=>:Vaccinated_3,:Rate=>p_tau),
        :t_tau_V3_waningImmunityFromVaccinated_3=>(:Vaccinated_3=>:Vaccinated_2,:Rate=>p_tau),
        :t_tau_V2_waningImmunityFromVaccinated_2=>(:Vaccinated_2=>:Vaccinated_1,:Rate=>p_tau),
        :t_eps_S_waningImmunityFromVaccinated_1=>(:Vaccinated_1=>:Susceptible,:Rate=>p_eps),        
        # finally, we have transitions involving receipt of vaccines whilst in a vaccinated state
        # NB: There is no opportunity to boost immunity through vaccination (or natural exposure) 
        # in the highest state of vaccine-induced immunity (Vaccinated_4)
        :t_transition10_vaccinationFromVaccinated_1=>(:Vaccinated_1=>:Vaccinated_2,:Rate=>p_vaccinationRate),
        :t_transition11_vaccinationFromVaccinated_2=>(:Vaccinated_2=>:Vaccinated_3,:TimeOut=>p_vaccinationRate),
        :t_transition12_vaccinationFromVaccinated_3=>(:Vaccinated_3=>:Vaccinated_4,:TimeOut=>p_vaccinationRate),    
    ), #transitions
()# alternatives for transitions
)
# plot out the state chart
StateCharts.Graph(pertussisStatechart)

######  define the user_defined rewrite rules: rewrite rules for contacts of Infectives
######  Note that each transitions_rule has a pair timer=>rule

# helper function takes in two arguments:
# S: state of being infected: :Susceptible, :Vaccinated_1, :Vaccinated_2, :Vaccinated_3, :Recovered_1, :Recovered_2 or :Recovered_3
# I: state of infectives: :I, :Im or :Iw
# returns three arrays: states including in the ACSet of L, I and R, where L <- I -> R is the rewrite rule
f_infective_collect(S::Symbol,I::Symbol) = map(x->[[[S],[x]],[x],[[I],[x]]],Infectives) 

# define the user_defined rewrite rules
transitions_rules = [map((x,y)->ContinuousHazard(1/(c*x))=>make_infectious_rule_MultipleObjects(pertussisStatechart,y...),β,f_infective_collect(:Susceptible,:Incubation_I))...,
                     map((x,y)->ContinuousHazard(1/(c*x))=>make_infectious_rule_MultipleObjects(pertussisStatechart,y...),β,f_infective_collect(:Vaccinated_1,:Incubation_I))...,
                     map((x,y)->ContinuousHazard(1/(c*x))=>make_infectious_rule_MultipleObjects(pertussisStatechart,y...),β,f_infective_collect(:Recovered_1,:Incubation_Im))...,
                     map((x,y)->ContinuousHazard(1/(c*x))=>make_infectious_rule_MultipleObjects(pertussisStatechart,y...),β,f_infective_collect(:Vaccinated_2,:Incubation_Im))...,
                     map((x,y)->ContinuousHazard(1/(c*x))=>make_infectious_rule_MultipleObjects(pertussisStatechart,y...),β,f_infective_collect(:Recovered_2,:Incubation_Iw))...,
                     map((x,y)->ContinuousHazard(1/(c*x))=>make_infectious_rule_MultipleObjects(pertussisStatechart,y...),β,f_infective_collect(:Recovered_3,:Recovered_4))...,
                     map((x,y)->ContinuousHazard(1/(c*x))=>make_infectious_rule_MultipleObjects(pertussisStatechart,y...),β,f_infective_collect(:Vaccinated_3,:Recovered_4))...
];

# Plot out the model schema
# Note: this step is not nessesary. Here it is only used to show the model schema
schema_statechart = StateChartABMSchema_MultipleObjects(pertussisStatechart,:P)
to_graphviz(schema_statechart |> schemaACSet |> schemaPresent; prog="dot")

# define the initial state
init = StateChartCset_MultipleObjects(pertussisStatechart)
add_parts!(init,:P,Int(totalPopulation),PID=1:totalPopulation)
add_parts!(init,:Susceptible,Int(totalPopulation-1),SusceptibleP=2:totalPopulation)
add_part!(init,:Infective_I,Infective_IP=1)
init

# run the ABM model
res = run!(make_ABM(pertussisStatechart,transitions_rules,is_schema_singObject=false), init; maxtime=50);

# plot out the results of each state
Makie.plot(res; Dict(o=>X->nparts(X,o) for o in states)...)