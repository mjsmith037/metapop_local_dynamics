import DelimitedFiles
import Random
import Statistics
import StatsBase
import ShiftedArrays
import CSV

include("multipop_SEI.jl")

args = tryparse.(Float64, ARGS)

function get_cycle_lengths(vector)
    maxima_locations = findall(isless.(ShiftedArrays.lag(vector), vector) .& isless.(ShiftedArrays.lead(vector), vector))
    return(maxima_locations .- ShiftedArrays.lag(maxima_locations))
end

function ave_cycle(vector)
    cycle_lengths = get_cycle_lengths(vector)
    if all(ismissing.(cycle_lengths))
         return(missing)
    end
    return(Statistics.mean(skipmissing(cycle_lengths)))
end

function var_cycle(vector)
    cycle_lengths = get_cycle_lengths(vector)
    if all(ismissing.(cycle_lengths))
         return(missing)
    end
    return(Statistics.var(skipmissing(cycle_lengths)))
end

function summarise_dynamics(vector, digits=nothing)
    # -1 = extinct, 0 = stable, 1 = unconverged, >1 = chaos/cycles
    vector = isnothing(digits) ? vector : round.(vector, digits=digits)
    # eliminate sequentially duplicate prevalences introduced due to rounding
    vector = StatsBase.rle(vector)[1]
    # if the last element is 0, the prevalence is extinct, regardless of other dynamics
    if vector[end] == 0 return(-1) end
    # pick out the local minima (use isless bc handles missing values introduced by lead/lag
    # silently: missing > all numbers)
    minima = vector[isless.(vector, ShiftedArrays.lag(vector)) .& isless.(vector, ShiftedArrays.lead(vector))]
    maxima = vector[isless.(ShiftedArrays.lag(vector), vector) .& isless.(ShiftedArrays.lead(vector), vector)]
    unique_minima = unique(minima)
    unique_maxima = unique(maxima)
    # if there are no maxima, the prevalence is constant or has a period proportional to the stepsize (unlikely)
    # if only one maxima, the vector is monotonic (implies failure to converge)
    if length(maxima) <= 1 return(length(maxima)) end
    # now check for negative monotonicity in minima => assume going extinct
    if all(((unique_minima - ShiftedArrays.lag(unique_minima)) .<= 0)[2:end]) return(-1) end
    # if positive monotonicity in minima along with negative monotonicity in maxima, then call it stable (just not yet converged)
    if all(((unique_minima - ShiftedArrays.lag(unique_minima)) .>= 0)[2:end])
        if all(((unique_maxima - ShiftedArrays.lag(unique_maxima)) .<= 0)[2:end]) return(0) else return(1) end
    end
    # if not caught by any of the earlier conditions, then the dynamics are chaotic/cyclical
    return(length(unique_maxima))
end

# NOT YET FUNCTIONAL WITH CURRENT DataFrames VERSION
# function test(vector, digits=nothing)
#     # -1 = extinct, 0 = stable, 1 = unconverged, >1 = chaos/cycles
#     vector = isnothing(digits) ? vector : round.(vector, digits=digits)
#     # eliminate sequentially duplicate prevalences introduced due to rounding
#     vector = StatsBase.rle(vector)[1]
#     # if the last element is 0, the prevalence is extinct, regardless of other dynamics
#     if vector[end] == 0 return(-1) end
#     # pick out the local minima (use isless bc handles missing values introduced by lead/lag
#     # silently: missing > all numbers)
#     minima = vector[isless.(vector, ShiftedArrays.lag(vector)) .& isless.(vector, ShiftedArrays.lead(vector))]
#     maxima_locations = findall(isless.(ShiftedArrays.lag(vector), vector) .& isless.(ShiftedArrays.lead(vector), vector))
#     maxima = vector[maxima_locations]
#     unique_minima = unique(minima)
#     unique_maxima = unique(maxima)
#     # if there are no maxima, the prevalence is constant or has a period proportional to the stepsize (unlikely)
#     # if only one maxima, the vector is monotonic (implies failure to converge)
#     if length(maxima) <= 1 return(length(maxima)) end
#     # now check for negative monotonicity in minima => assume going extinct
#     if all(((unique_minima - ShiftedArrays.lag(unique_minima)) .<= 0)[2:end]) return(-1) end
#     # if positive monotonicity in minima along with negative monotonicity in maxima, then call it stable (just not yet converged)
#     if all(((unique_minima - ShiftedArrays.lag(unique_minima)) .>= 0)[2:end])
#         if all(((unique_maxima - ShiftedArrays.lag(unique_maxima)) .<= 0)[2:end]) return(0) else return(1) end
#     end
#     # if not caught by any of the earlier conditions, then the dynamics are chaotic/cyclical

#     cycle_lengths = maxima_locations .- ShiftedArrays.lag(maxima_locations)

#     ave = all(ismissing.(cycle_lengths)) ? missing : Statistics.mean(skipmissing(cycle_lengths))
#     var = all(ismissing.(cycle_lengths)) ? missing : Statistics.var(skipmissing(cycle_lengths))

#     return([ave, var, length(unique_maxima)])
# end

networks = vcat([joinpath.(root, files) for (root, dirs, files) in
                    walkdir("../Data/synthetic_networks/")]...)[parse(Int, ARGS[1]):parse(Int, ARGS[2])]

for net in networks
    outfile = foldl(replace, ("../Data/synthetic_networks" => "../Results/LNS", "mat" => "csv"), init=net)
    if isfile(outfile) & !parse(Bool, ARGS[12])
        continue
    end
    rseed = parse(Int, replace(basename(net), ".mat" => ""))
    Random.seed!(rseed)
    
    chi = DelimitedFiles.readdlm(net)
    populations = size(chi)[1]
    initial_conditions = rand(populations, 3)
    
    timeseries = runsimulation(tstep    = Int(0.95*args[3]):args[3],
                               tmax     = args[3],
                               r        = args[4],
                               K        = rand(populations) * (args[11] - args[10]) .+ args[10],
                               beta     = args[5],
                               sigma    = args[6],
                               mu       = args[7],
                               nu       = args[8],
                               chi      = chi * args[9],
                               initvals = initial_conditions)


    rename!(timeseries, [["$(x)_$y" for x in ["S" "E" "I"] for y in 1:25]..., "time"])
    long_timeseries = stack(timeseries, Not(:time), variable_name=:variable, value_name=:prevalence)
    insertcols!(long_timeseries, [n => getindex.(split.(long_timeseries.variable, '_'), i)
                                    for (i, n) in enumerate([:equation, :population])]...)
    filter!(row -> row.equation == "I", long_timeseries)
    long_timeseries.prevalence = round.(long_timeseries.prevalence; digits=8)

    grouped = groupby(long_timeseries, :population)
    summary = combine(grouped, :prevalence => sum,
                            :prevalence => Statistics.var, 
                            :prevalence => ave_cycle, 
                            :prevalence => var_cycle, 
                            :prevalence => summarise_dynamics)

    summary[!, "net_type"] .= splitpath(net)[4]
    summary[!, "replicate"] .= rseed

    CSV.write(outfile, summary)
end

