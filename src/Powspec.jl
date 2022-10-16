module Powspec

using Parameters

export power_spectrum

include("structs.jl")
include("box_periodic.jl")
include("box_randoms.jl")
include("survey.jl")


end
