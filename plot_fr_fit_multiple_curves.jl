
using StatsPlots
using Distributions
using Random
using MCMCChains
using Turing
using DataFrames
using LambertW
using CSV
using Tables
using Serialization

gradient = cgrad(:berlin,9,categorical = true) # create a color gradient
subcolors = [1,2,8,9]

df_data = CSV.read("copepods.csv",DataFrame)
exp_labels = ["No cue; 26C","Cue; 26C","No cue; 20C","Cue; 20C"]

# make a collection of x ranges for the different prey types
xrange = collect(range(0.0, 40, length=50))

# create an array of empty plots to fill up with the loop
allPlots = Array{Plots.Plot{Plots.GRBackend},1}(undef,4)
allPlots_post = Array{Plots.Plot{Plots.GRBackend},1}(undef,2)

include("plot_priors.jl")
figure1 = plot_priors(allPlots_post,prey_offered,prey_eaten)

########################################################################
# Now loop over the data sets to plot each fit and the posteriors
# on top of the priors
########################################################################

# how many and what are the experiments?
exps = unique(df_data.Experiment)

chain_of_h_yes = []
chain_of_h_no = []


for i in eachindex(exps)
    # grab data again
    indices = findall(df_data.Experiment .== exps[i])
    prey_eaten = df_data.eaten[indices]
    prey_offered = df_data.Density[indices]
    # open up plot canvas
    allPlots[i] = plot(ylims=[0,30],
        ylabel = "Number of prey eaten",
        xlabel = "Number of prey offered",
        size=(500,300),
        dpi=500)

    # grab the saved chain
    restored_chain = deserialize(string("chain-file_DS_",string(i),"type2.jls"))

    # pull out a sample of the specific parameter chains
    chain_of_a = sample(restored_chain["a"],100)
    chain_of_h = sample(restored_chain["h"],100)

    if i ==3 
        chain_of_h_no = chain_of_h
    elseif i == 4
        chain_of_h_yes = chain_of_h
    end

    # calculate the fits across all params in the chains
    y = xrange .- lambertw.(chain_of_a' .* chain_of_h' .* xrange .* exp.(-chain_of_a' .* (0.5 .- chain_of_h' .* xrange))) ./ (chain_of_a' .* chain_of_h')

    # plot all of those fits
    plot!(xrange, y,
        label="",
        color=gradient[subcolors[i]],
        Linewidth=1.5,
        alpha=0.4)

    fitted_params = DataFrame(summarystats(restored_chain))

    # put the mean fit on top of this
    mean_a = fitted_params[1,2]
    mean_h = fitted_params[2,2]

    y = xrange .- lambertw.(mean_a .* mean_h .* xrange .* exp.(-mean_a .* (0.5 .- mean_h .* xrange))) ./ (mean_a * mean_h)
    plot!(xrange, y,
        label=exp_labels[i],
        color=gradient[subcolors[i]],
        linewidth=2)

    # add a black line on top
    plot!(xrange, y,
        label="",
        color=:black,
        linestyle=:dash,
        linewidth=2)

    # drop the data on top
    scatter!(prey_offered,prey_eaten,label="",
        markerstrokecolor=gradient[subcolors[i]],
        markercolor=gradient[subcolors[i]],
        markersize=4)

    # plot the posteriors on the post canvas
        density!(allPlots_post[1], chain_of_a, label=exp_labels[i], color=gradient[subcolors[i]], linewidth = 2)
        density!(allPlots_post[2], chain_of_h, label=exp_labels[i], color=gradient[subcolors[i]], linewidth = 2)

end

########################################################################
# Now bring the panels together
########################################################################

figure_fits = plot(allPlots[1],allPlots[4],allPlots[2],allPlots[3], layout=(2,2), size=(500,600), dpi=1000)
savefig(figure_fits,"FR_fits_copepods.png")

figure_posteriors = plot(allPlots_post[1],allPlots_post[2])
savefig(figure_posteriors,"FR_posterios_copepods.png")

# if you want to test for differences, do the cross-comparison across all values
chain_of_h_diffs = chain_of_h_yes .- chain_of_h_no'
   
# then determine the crossing point for zero from the cumulative probability distribution
my_ecdf = ecdf(chain_of_h_diffs[:])
quantile_at_zero = my_ecdf(0.0)
