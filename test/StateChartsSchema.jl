module TestStateChartsSchema

using Test
using StateCharts

# Build a state chart of smoking, including three states: 
# NS: Never Smoking
# CS: Current Smoking
# PS: Previous Smoking

smoking = UnitStateChartF((:NS,:CS,:PS), # states
(:initiate=>(:NS=>:CS,:Rate=>0.01),:quit=>(:CS=>:PS,:Rate=>0.1),:relapse=>(:PS=>:CS,:Rate=>0.05)), #non-start transitions
() # alternatives for non-start transitions
)

@test snames(smoking) == [:NS,:CS,:PS]
@test tnames(smoking) == [:initiate, :quit, :relapse]

Graph(smoking)

# Build a state chart of breastfeeding, including two states:
# NonBreastfeeding
# Breastfeeding

breastfeeding = UnitStateChartF((:NonBreastfeeding,:Breastfeeding), # states
(:babyBorn=>(:NonBreastfeeding=>:NonBreastfeeding,:TimeOut=>30),:returnToWork=>(:Breastfeeding=>:NonBreastfeeding,:Rate=>0.05),:lactationProblem=>(:Breastfeeding=>:NonBreastfeeding,:Rate=>0.1)), #non-start transitions
(:babyBorn=>((:optToBreastfeeding,0.6)=>:Breastfeeding)) # alternatives for non-start transitions
)

@test snames(breastfeeding) == [:NonBreastfeeding,:Breastfeeding]
@test tnames(breastfeeding) == [:babyBorn, :returnToWork, :lactationProblem]

Graph(breastfeeding)

########### create two state charts ##################
smoke_healthy = StateChartsF((:SmokeStateChart=>(:NS,:CS,:PS),(:HealthyStateChart=>(:UH,:H))), # partition and states
(:initiate=>(:NS=>:CS,:Rate=>0.01),:quit=>(:CS=>:PS,:Rate=>0.1),:relapse=>(:PS=>:CS,:Rate=>0.05),:sick=>(:H=>:UH,:Rate=>0.001),:recover=>(:UH=>:H,:TimeOut=>7)), #non-start transitions
() # alternatives for non-start transitions
)

@test snames(smoke_healthy) == [:NS,:CS,:PS,:UH,:H]
@test tnames(smoke_healthy) == [:initiate,:quit,:relapse,:sick,:recover]


Graph(smoke_healthy)

end