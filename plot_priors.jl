function plot_priors(allPlots_post,prey_offered,prey_eaten)

# set up the empirically determined starting points
    max_a = maximum(prey_eaten./prey_offered)
    starting_a = max_a
    starting_h = 1/maximum(df_data.eaten)

# pass these and the whole data set to the model for a common Posterior
    model_type2 = fun_res_RPE(prey_offered,prey_eaten,starting_a,max_a,starting_h)

# sample the priors from the model
    prior_chain = sample(model_type2, chain_type=MCMCChains.Chains, Prior(), 1000)

# plot the space clearance rate prior
    allPlots_post[1] = density(prior_chain[:a],label="Prior", color=:black, linestyle=:dash, linewidth = 2, xlabel="a", legend=:best, xlims=[0,20],)

# plot the handling time prior
    allPlots_post[2] = density(prior_chain[:h],label="Prior", color=:black, linestyle=:dash, linewidth = 2, xlabel="h", legend=:best, xlims=[0,0.1],)

# take a peak at the two prior plots
    figure_test = plot(allPlots_post[1],allPlots_post[2])

return figure_test

end
