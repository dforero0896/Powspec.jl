using Revise
using Plots
using CSV
using DataFrames
using Powspec


const P0 = 5f3
fkp_weights(nz, P0) = 1 / (1 + nz * P0)


data_cat_fn = "/global/cfs/projectdirs/desi/mocks/UNIT/HOD_Shadab/multiple_snapshot_lightcone/UNIT_lightcone_multibox_ELG_footprint_nz_NGC.dat"
rand_cat_fn = "/global/cfs/projectdirs/desi/mocks/UNIT/HOD_Shadab/multiple_snapshot_lightcone/UNIT_lightcone_multibox_ELG_footprint_nz_1xdata_5.ran_NGC.dat"

data_cat = DataFrame(CSV.File(data_cat_fn, delim=" ", header=[:ra, :dec, :d, :z, :nz], types = [Float64 for _ in 1:5]))
rand_cat = DataFrame(CSV.File(rand_cat_fn, delim=" ", header=[:ra, :dec, :d, :z, :nz], types = [Float64 for _ in 1:5]))

data_cat = data_cat[map(z -> ((z > 0.8) & (z < 1)), data_cat.z), :]
rand_cat = rand_cat[map(z -> ((z > 0.8) & (z < 1)), rand_cat.z), :]

data = (values(data_cat[!,:ra]), values(data_cat[!,:dec]), values(data_cat[!,:z]))
data_w = zero(values(data_cat[!,:ra])) .+ 1
data_fkp = fkp_weights.(values(data_cat[!,:nz]), Ref(P0))
data_nz = values(data_cat[!,:nz])

random = (values(rand_cat[!,:ra]), values(rand_cat[!,:dec]), values(rand_cat[!,:z]))
random_w = zero(values(rand_cat[!,:ra])) .+ 1
random_fkp = fkp_weights.(values(rand_cat[!,:nz]), Ref(P0))
random_nz = values(rand_cat[!,:nz])



pk = power_spectrum(data, data_w, data_fkp, data_nz, random, random_w, random_fkp, random_nz, "test/powspec_lc.conf"; output_file = "lc.dat", precision = :single)

p1 = plot(pk[:k], pk[:k] .* view(pk[:multipoles], 1,:))
p2 = plot(pk[:k], pk[:k] .* view(pk[:multipoles], 2,:))
p3 = plot(pk[:k], pk[:k] .* view(pk[:multipoles], 3,:))
savefig("test/lc.png")



pk = power_spectrum((data, data), (data_w, data_w), (data_fkp, data_fkp), (data_nz, data_nz), 
                    (random, random), (random_w, random_w), (random_fkp, random_fkp), (random_nz, random_nz), "test/powspec_lc_cross.conf", precision = :single)

plot!(p1, pk[:k], pk[:k] .* view(pk[:auto_multipoles][1], 1,:))
plot!(p2, pk[:k], pk[:k] .* view(pk[:auto_multipoles][1], 2,:))
plot!(p3, pk[:k], pk[:k] .* view(pk[:auto_multipoles][1], 3,:))

plot!(p1, pk[:k], pk[:k] .* view(pk[:auto_multipoles][2], 1,:)) #This gives wrong results
plot!(p2, pk[:k], pk[:k] .* view(pk[:auto_multipoles][2], 2,:))
plot!(p3, pk[:k], pk[:k] .* view(pk[:auto_multipoles][2], 3,:))

plot!(p1, pk[:k], pk[:k] .* view(pk[:cross_multipoles], 1,:))
plot!(p2, pk[:k], pk[:k] .* view(pk[:cross_multipoles], 2,:))
plot!(p3, pk[:k], pk[:k] .* view(pk[:cross_multipoles], 3,:))
plot(p1, p2, p3, layout=(3,1), dpi=300, size = (1000,1000), legend = :topleft)
savefig("test/lc.png")
