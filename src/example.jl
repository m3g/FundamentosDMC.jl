

using Plots
using Revise
using FundamentosDMC

sys = System(n=10_000)

v = init_velocities(sys,Options(),dist=:flat)
printxyz(v, sys, "test.xyz")
v0, d0 = velocity_distribution(sys,"test.xyz")
plot(v0,d0)

v = init_velocities(sys,Options(),dist=:normal)
printxyz(v, sys, "test.xyz")
v0, d0 = velocity_distribution(sys,"test.xyz")
plot!(v0,d0)

sys = System(n=100)
minimize!(sys)
md_out = md_langevin(
  sys,
  Options(lambda=0.01,nsteps=20_000,trajectory_file="md.xyz")
)
v, d = velocity_distribution(sys,"vel.xyz")
plot!(v,d)

