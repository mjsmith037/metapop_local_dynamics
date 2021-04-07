using DataFrames
using Random
using StaticArrays
using DifferentialEquations

seasonality(t, β, ε) = (1 .+ ε .* sin(π * t)^6) .* β

dydt(y, z, w, λ, γ, σ, μ, ν, χ) = λ .* ((1 .- w) + (1 .- γ) .* (w - z)) - σ .* y - (μ .+ ν) .* y + χ' * y
dzdt(z, λ, μ, χ) = λ .* (1 .- z) - μ .* z + χ' * z
dwdt(w, ssλ, μ, χ) = ssλ .* (1 .- w) - μ .* w + χ' * w

function mantis!(deriv, state, params, t)
     y = @view state[:,:,1];  z = @view state[:,:,2];  w = @view state[:,:,3]
    dy = @view deriv[:,:,1]; dz = @view deriv[:,:,2]; dw = @view deriv[:,:,3]
    # adjust β to address seasonality
    βi = seasonality(t, params.β, params.ε)
    # force of infection in y,z equations
    λ = βi .* y
    # force of infection for strains that share epitopes with each strain
    # [i.e. force of infection for each virus in w equation]
    ssλ = hcat([sum(βi[:,x] .* view(y, :, params.sharedalleles[:,x]), dims=2) for x in 1:size(y, 2)]...)
    # calculate the derivatives
    dy .= dydt(y, z, w, λ, params.γ, params.σ, params.μ, params.ν, params.χ)
    dz .= dzdt(z, λ, params.μ, params.χ)
    dw .= dwdt(w, ssλ, params.μ, params.χ)
end

macro Name(arg)
   string(arg)
end

function parametersizecheck(parameter, parametername::String, npops::Int64, nstrains::Int64)
    if length(parameter) == 1
        resized_parameter = SMatrix{npops, nstrains}(fill(parameter, npops, nstrains))
    # note, if npops == nstrains, this formulation assumes differences
    # are between populations rather than between strains
    elseif length(parameter) == npops
        resized_parameter = SMatrix{npops, nstrains}(hcat(fill(reshape(parameter, npops, 1), nstrains)...))
    elseif length(parameter) == nstrains
        resized_parameter = SMatrix{npops, nstrains}(vcat(fill(reshape(parameter, 1, nstrains), npops)...))
    elseif length(parameter) == npops * nstrains
        resized_parameter = SMatrix{npops, nstrains}(reshape(parameter, npops, nstrains))
    else error("Invalid $parametername value") end
    return(resized_parameter)
end

function runMANTIS(;strainstructure=[2 2], tmax=1000.0, tstep=missing,
                   beta=4.0, gamma=0.72, epsilon=0.0, sigma=1.0, mu=0.0, nu=0.0,
                   chi=[-0.02 0.0; 0.0 -0.02], rseed=0, initvals::Any="random")
   # strainstructure=[2 2];tmax=1000;tstep=0:0.1:1000;beta=200;gamma=0.75;
   # epsilon=0;sigma=50;mu=0.02;chi=0;nu=0;initvals="random"

    chi = isempty(size(chi)) ? SMatrix{1,1}(chi) : SMatrix{size(chi)...}(chi)

    strains = reshape(collect(Iterators.product([1:x for x in strainstructure]...)), :, 1)

    npops = size(chi, 1)
    nstrains = length(strains)

    if typeof(initvals) == String && initvals == "random"
        # Random.seed!(rseed)
        initvals = rand(npops, nstrains) |> x->x./10
    end

    sharedalleles = hcat([[any(x .== y) for x in strains] for y in strains]...)

    initcond = cat(initvals, initvals,
                   min.(hcat([sum(initvals[:,sharedalleles[:,x]], dims=2) for x in 1:nstrains]...), 1),
                   dims=3)

    params = (β = parametersizecheck(beta, @Name(beta), npops, nstrains),
              σ = parametersizecheck(sigma, @Name(sigma), npops, nstrains),
              μ = parametersizecheck(mu, @Name(mu), npops, nstrains),
              ν = parametersizecheck(nu, @Name(nu), npops, nstrains),
              γ = parametersizecheck(gamma, @Name(gamma), npops, nstrains),
              ε = parametersizecheck(epsilon, @Name(epsilon), npops, nstrains),
              χ = chi, sharedalleles = sharedalleles)

    prob = ODEProblem(mantis!, initcond, (0.0, tmax), params)
    if ismissing(tstep)
        sol = solve(prob, save_everystep=false, isoutofdomain=(u,p,t) -> any(x -> (x < 0) | (x > 1), u))
        return(sol)
    else
        sol = solve(prob, saveat=tstep, abstol=1e-8, reltol=1e-8, maxiters=1e9,
                    isoutofdomain=(u,p,t) -> any(x -> (x < 0) | (x > 1), u))
        return(Dict("timeseries" => filter(row -> row[end] ∈ tstep, DataFrame(hcat([[reshape(u, (1,prod(size(u))))..., t]
                                                        for (u,t) in tuples(sol)]...)')),
                    "parameters" => params))
    end
end
