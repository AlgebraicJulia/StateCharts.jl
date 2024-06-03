using Test

@testset "Code Quality (Aqua.jl)" begin
  #include("aqua.jl")
end

@testset "StateChartsSchema" begin
  include("StateChartsSchema.jl")
end

@testset "StateChartsABMsInterfaces" begin
  include("interfacesToABMs/StateChartsABMsInterfaces.jl")
end

@testset "NetworkSchemaInterfaces" begin
  include("interfacesToABMs/NetworkSchemaInterfaces.jl")
end