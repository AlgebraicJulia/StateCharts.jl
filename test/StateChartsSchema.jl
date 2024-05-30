module TestStateChartsSchema

using Test
using StateCharts

# Build a state chart of smoking, including three states: 
# NS: Never Smoking
# CS: Current Smoking
# PS: Previous Smoking

smoking = StateChartF((:NS,:CS,:PS), # states
(:initiate=>(:NS=>:CS,:Rate=>0.01),:quit=>(:CS=>:PS,:Rate=>0.1),:relapse=>(:PS=>:CS,:Rate=>0.05)), #non-start transitions
(), # alternatives for non-start transitions
(:SmokingStateChart=>:NS), # start transitions
(:SmokingStateChart=>((:initializeToCS,0.2)=>:CS,(:initializeToPS,0.1)=>:PS)) # alternatives for start transitions
)

@test snames(smoking) == [:NS,:CS,:PS]
@test tnames(smoking) == [:initiate, :quit, :relapse]

Graph(smoking)

# Build a state chart of breastfeeding, including two states:
# NonBreastfeeding
# Breastfeeding

breastfeeding = StateChartF((:NonBreastfeeding,:Breastfeeding), # states
(:babyBorn=>(:NonBreastfeeding=>:NonBreastfeeding,:TimeOut=>30),:returnToWork=>(:Breastfeeding=>:NonBreastfeeding,:Rate=>0.05),:lactationProblem=>(:Breastfeeding=>:NonBreastfeeding,:Rate=>0.1)), #non-start transitions
(:babyBorn=>((:optToBreastfeeding,0.6)=>:Breastfeeding)), # alternatives for non-start transitions
(:BreastFeedingStateChart=>:NonBreastfeeding), # start transitions
() # alternatives for start transitions
)

@test snames(breastfeeding) == [:NonBreastfeeding,:Breastfeeding]
@test tnames(breastfeeding) == [:babyBorn, :returnToWork, :lactationProblem]

Graph(breastfeeding)

end