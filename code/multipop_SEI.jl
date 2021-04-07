#using Plots
using DataFrames
using Random
using StaticArrays
using DifferentialEquations
using ParameterizedFunctions

seasonality(t, β, ε) = (1 .+ ε .* sin(π * t)^6) .* β

# after Anderson, et al. 1981 
function sei_multi(du, u, p, t)
    r, K, β, σ, μ, ν, χ = p
    X = vec(sum(u, dims=2)) ./ K
    du[:,1] .= r .* u[:,1] .* (1 .- X) - β .* u[:,1] .* u[:,3] .+ χ' * u[:,1]
    du[:,2] .= β .* u[:,1] .* u[:,3] .- (σ .+ μ .+ r .* X) .* u[:,2] .+ χ' * u[:,2]
    du[:,3] .= σ .* u[:,2] .- (ν .+ μ .+ r .* X) .* u[:,3] .+ χ' * u[:,3]
end

function runsimulation(;tmax=1000.0, tstep=missing,
                        r=0.5, K=15.0, beta=80.0, sigma=13.0, mu=0.5, nu=73.0,
                        chi=[-0.1 0.1; 0.0 0.0], rseed=0, initvals=[4.0 0.0 1.0; 14.0 0.0 1.0])

    chi = isempty(size(chi)) ? SMatrix{1,1}(chi) : SMatrix{size(chi)...}(chi)

    params = [[length(x) > 1 ? vec(x) : x for x in [r, K, beta, sigma, mu, nu]]..., chi]

    prob = ODEProblem(sei_multi, initvals, (0.0, tmax), params)
    if ismissing(tstep)
        sol = solve(prob, alg_hints=[:stiff], callback=PositiveDomain(), abstol=1e-8, reltol=1e-8)
        #plot(sol, xlabel="Time", ylabel="Number", lw=2)
        return(sol)
    else
        sol = solve(prob, alg_hints=[:stiff], isoutofdomain=(u,p,t) -> any(x -> x < 0, u), abstol=1e-8, reltol=1e-8, maxiters=1e9, saveat=tstep)
        return(filter(row -> row[end] ∈ tstep, DataFrame(hcat([[reshape(u, (1,prod(size(u))))..., t]
                                                               for (u,t) in tuples(sol)]...)')))
    end
end
