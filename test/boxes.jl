using Revise
using Powspec
using StaticArrays
using Parameters
using NPZ
using Plots
using CSV
using Statistics
using FFTW
using DataFrames

FFTW.set_num_threads(128)

data_cat_fn = "/global/cfs/projectdirs/desi/mocks/UNIT/HOD_Shadab/HOD_boxes/redshift0.9873/UNIT_DESI_Shadab_HOD_snap97_ELG_v0.txt"
data_cat_pos = DataFrame(CSV.File(data_cat_fn, delim = ' ', ignorerepeated = true, header=[:x, :y, :d, :z], types = [Float64 for _ in 1:4]))
data_cat_pos = map(values, (data_cat_pos[!,:x], data_cat_pos[!,:y], data_cat_pos[!,:z]))

data = data_cat_pos
data_w = zero(data[1]) .+ 1
pk = power_spectrum(data, data_w, "test/powspec_auto.conf"; output_file = "test.dat", precision = :single)

p1 = plot(pk[:k], pk[:k] .* view(pk[:multipoles], 1,:))
p2 = plot(pk[:k], pk[:k] .* view(pk[:multipoles], 2,:))
p3 = plot(pk[:k], pk[:k] .* view(pk[:multipoles], 3,:))
plot(p1, p2, p3, layout=(3,1), xscale = :log10)
savefig("test/test.png")


random = Tuple(1000rand(Float64, 10size(data[1], 1)) for _ in 1:3)
random_w = zero(random[1]) .+ 1
pk = power_spectrum(data, data_w, random, random_w, "test/powspec_auto.conf"; output_file = "test.dat", precision = :single)

plot!(p1, pk[:k], pk[:k] .* view(pk[:multipoles], 1,:))
plot!(p2, pk[:k], pk[:k] .* view(pk[:multipoles], 2,:))
plot!(p3, pk[:k], pk[:k] .* view(pk[:multipoles], 3,:))
plot(p1, p2, p3, layout=(3,1), xscale = :log10)
savefig("test/test.png")


pk = power_spectrum((data, data), (data_w, data_w), "test/powspec_auto.conf", precision = :single)

plot!(p1, pk[:k], pk[:k] .* view(pk[:auto_multipoles][1], 1,:))
plot!(p2, pk[:k], pk[:k] .* view(pk[:auto_multipoles][1], 2,:))
plot!(p3, pk[:k], pk[:k] .* view(pk[:auto_multipoles][1], 3,:))

plot!(p1, pk[:k], pk[:k] .* view(pk[:auto_multipoles][2], 1,:))
plot!(p2, pk[:k], pk[:k] .* view(pk[:auto_multipoles][2], 2,:))
plot!(p3, pk[:k], pk[:k] .* view(pk[:auto_multipoles][2], 3,:))

plot!(p1, pk[:k], pk[:k] .* view(pk[:cross_multipoles], 1,:))
plot!(p2, pk[:k], pk[:k] .* view(pk[:cross_multipoles], 2,:))
plot!(p3, pk[:k], pk[:k] .* view(pk[:cross_multipoles], 3,:))
plot(p1, p2, p3, layout=(3,1), xscale = :log10)
savefig("test/test.png")


pk = power_spectrum((data, data), (data_w, data_w), (random, random), (random_w, random_w), "test/powspec_auto.conf", precision = :single)

plot!(p1, pk[:k], pk[:k] .* view(pk[:auto_multipoles][1], 1,:))
plot!(p2, pk[:k], pk[:k] .* view(pk[:auto_multipoles][1], 2,:))
plot!(p3, pk[:k], pk[:k] .* view(pk[:auto_multipoles][1], 3,:))

plot!(p1, pk[:k], pk[:k] .* view(pk[:auto_multipoles][2], 1,:))
plot!(p2, pk[:k], pk[:k] .* view(pk[:auto_multipoles][2], 2,:))
plot!(p3, pk[:k], pk[:k] .* view(pk[:auto_multipoles][2], 3,:))

plot!(p1, pk[:k], pk[:k] .* view(pk[:cross_multipoles], 1,:))
plot!(p2, pk[:k], pk[:k] .* view(pk[:cross_multipoles], 2,:))
plot!(p3, pk[:k], pk[:k] .* view(pk[:cross_multipoles], 3,:))
plot(p1, p2, p3, layout=(3,1), xscale = :log10, dpi=300, size = (1000,1000), legend = :topleft)
savefig("test/test.png")