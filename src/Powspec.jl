module Powspec

using Parameters

export power_spectrum
println("Trying to locate library in ", ENV["LIBPOWSPEC_PATH"])
include("structs.jl")
include("box_periodic.jl")
include("box_randoms.jl")
include("survey.jl")


end
