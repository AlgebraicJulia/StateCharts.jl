using Test

@testset "Code Quality (Aqua.jl)" begin
  #include("aqua.jl")
end

@testset "StateChartsSchema" begin
  include("StateChartsSchema.jl")
end

@testset "StateChartsABMsInterfaces" begin
  include("interfacesToAMBs/StateChartsABMsInterfaces.jl")
end

