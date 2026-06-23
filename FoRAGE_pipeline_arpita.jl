dir = "C:/Users/johnp/Downloads/"
cd(dir)

using Turing
using MCMCChains
using StatsPlots
using Distributions
using Random
using DataFrames, LambertW, CSV, Tables
using ParetoSmooth

gradient = cgrad(:berlin,9,categorical = true) # create a color gradient

# specify that the fitting function be accessed
#include("fun_res_HDE.jl")
#include("fun_res_RPE.jl")

df_data = CSV.read("copepods.csv",DataFrame)
exps = unique(df_data.Experiment)

# MODEL 2 - Rogers random predator
@model function fun_res_RPE(prey_offered,prey_eaten,starting_a,max_a,starting_h) 
    a ~ truncated(Normal(starting_a[1],10*starting_a[1]),0,10) # space clearance rate
    h ~ truncated(Normal(starting_h,10*starting_h),0,0.25) # handling time
	#σ ~ LogNormal(1,2)
	#σ ~ LogNormal(var(prey_eaten),2)
    σ ~ InverseGamma(1,3)
    for i in 1:length(prey_eaten)
        prey_eaten[i] ~ truncated(Normal(prey_offered[i]-lambertw(a*h*prey_offered[i]*exp(-a*(0.5-h*prey_offered[i])))/(a*h),σ),0,maximum(prey_eaten))
	end
end

# MODEL 4 - type 3 asymptotic a model, modified from RPE
@model function fun_res_RPE_AA(prey_offered,prey_eaten,starting_a,max_a,starting_h) 
    a ~ truncated(Normal(starting_a[1],10*starting_a[1]),0,1) # space clearance rate
    h ~ truncated(Normal(starting_h,10*starting_h),0,100) # handling time
    k ~ truncated(Normal(0,0.1*maximum(prey_offered)),0,maximum(prey_offered)) # refuge size
	#σ ~ LogNormal(1,2)
    σ ~ InverseGamma(2,3)
    for i in 1:length(prey_eaten)
        prey_eaten[i] ~ truncated(Normal(prey_offered[i]-lambertw((a*prey_offered[i]/(k+prey_offered[i]))*h*prey_offered[i]*exp(-(a*prey_offered[i]/(k+prey_offered[i]))*(6-h*prey_offered[i])))/((a*prey_offered[i]/(k+prey_offered[i]))*h),σ),0,maximum(prey_eaten))
	end
end

#= MODEL 6 - type 3 refuge model, modified from RPE
@model function fun_res_RPE_REF(prey_offered,prey_eaten,starting_a,max_a,starting_h) 
    a ~ truncated(Normal(starting_a[1],10*starting_a[1]),0,1) # space clearance rate
    h ~ truncated(Normal(starting_h,10*starting_h),0,100) # handling time
    k ~ truncated(Normal(0,0.1*maximum(prey_offered)),0,maximum(prey_offered)) # refuge size
	#σ ~ LogNormal(1,2)
    σ ~ InverseGamma(2,3)
    for i in 1:length(prey_eaten)
        prey_eaten[i] ~ truncated(Normal((prey_offered[i]-k)-lambertw(a*h*(prey_offered[i]-k)*exp(-a*(6-h*(prey_offered[i]-k))))/(a*h),σ),0,maximum(prey_eaten))
	end
end=#



num_studies = 4
fitted_a_type2 = zeros(Float64,1,num_studies)
fitted_h_type2 = zeros(Float64,1,num_studies)

#=fitted_a_asym = zeros(Float64,1,num_studies)
fitted_h_asym = zeros(Float64,1,num_studies)
fitted_r_asym = zeros(Float64,1,num_studies)
fitted_a_ref = zeros(Float64,1,num_studies)
fitted_h_ref = zeros(Float64,1,num_studies)
fitted_k_ref = zeros(Float64,1,num_studies)
cv_elpd_type2 = zeros(Float64,1,num_studies)
cv_elpd_type3aa = zeros(Float64,1,num_studies)
cv_elpd_type3ref = zeros(Float64,1,num_studies)
naive_lpd_type2 = zeros(Float64,1,num_studies)
naive_lpd_type3aa = zeros(Float64,1,num_studies)
naive_lpd_type3ref = zeros(Float64,1,num_studies)=#
chain_a = zeros(Float64,4,2000, 2)
chain_h = zeros(Float64,4,2000)

# make a collection of x ranges for the different prey types
xrange = collect(range(0.0, 40, length=50))

# create an array of empty plots to fill up with the loop
allPlots = Array{Plots.Plot{Plots.GRBackend},1}(undef,4)
allPlots_post = Array{Plots.Plot{Plots.GRBackend},1}(undef,2)

allPlots_post[1] = plot(
    ylabel = "Frequency",
    xlabel = "Value",
    size=(600,400),
    dpi=1000,
    title="a")

    prior_a = truncated(Normal(starting_a[1],10*starting_a[1]),0,10) # space clearance rate
    x = 0:0.01:10
    y = pdf.(prior_a,x)

    plot!(allPlots_post[1],x,y,
        label="Prior", 
        linewidth=2.5, 
        xlabel="Value (a)", 
        ylabel="Density", 
        color=:red)

allPlots_post[2] = plot(
    ylabel = "Frequency",
    xlabel = "Value",
    size=(600,400),
    dpi=1000,
    title="h")
   
    prior_h = truncated(Normal(starting_h,10*starting_h),0,0.25) # handling time
    x = 0:0.01:0.25
    y = pdf.(prior_h,x)

    plot!(allPlots_post[2],x,y,
        label="Prior", 
        linewidth=2.5, 
        xlabel="Value (h)", 
        ylabel="Density", 
        color=:red)

for i = 1:4
    which_exp = exps[i]
    println(which_exp)
    indices = findall(df_data.Experiment .== exps[i])
    prey_eaten = df_data.eaten[indices]
    prey_offered = df_data.Density[indices]

    # pull starting values from dataset
    max_a = maximum(prey_eaten./prey_offered)
    starting_a = 0.5*max_a
    starting_h = 1/maximum(prey_eaten)

    # set up model object for type 2
    model_type2 = fun_res_RPE(prey_offered,prey_eaten,starting_a,max_a,starting_h)

    # call the fitting for type 2
        fr_chain_type2 = sample(
        model_type2,
        NUTS(20000,0.65),
        MCMCSerial(),
        chain_type=MCMCChains.Chains,
        2000,
        init_params = [(starting_a,starting_h,10)],
        2)

        # pull parameters and log results
        fitted_params_type2 = DataFrame(summarystats(fr_chain_type2))
  
        fitted_a_type2[i] = median(fr_chain_type2["a"])
        fitted_h_type2[i] = median(fr_chain_type2["h"])

        # store chains of a and h to access later
        chain_a[i,:,:] = fr_chain_type2["a"]
        chain_h[i,:,:] = fr_chain_type2["h"]

        density!(allPlots_post[1],chain_a[i,:], label = "Posterior", color=:gray, linewidth = 2)
        density!(allPlots_post[2],chain_h[i,:], label = "Posterior", color=:gray, linewidth = 2)

end


plot(allPlots_post[1],allPlots_post[2])


    # set new starting values from the type 2 fit
    #starting_a = fitted_params_type2[1,2]
    #starting_h = fitted_params_type2[2,2]

    # set up model object for type 3, asymptotic a model
    #model_type3_aa = fun_res_RPE_AA(prey_offered,prey_eaten,starting_a,max_a,starting_h)

    #= call the fitting for asymptotic a type 3
        fr_chain_type3_aa = sample(
        model_type3_aa,
        NUTS(5000,0.65),
        MCMCSerial(),
        4000,
        init_params = [(starting_a,starting_h,0.0001*maximum(prey_offered),10)],
        1)
        loo_type3aa = psis_loo(model_type3_aa, fr_chain_type3_aa)
        cv_elpd_type3aa[i] = loo_type3aa.estimates[1,1]
        naive_lpd_type3aa[i] = loo_type3aa.estimates[1,2]
        
        # pull parameters and log results
        fitted_params_type3 = DataFrame(summarystats(fr_chain_type3_aa))
        # store params, adjust a back to original units
        fitted_a_asym[i] = fitted_params_type3[1,2]
        fitted_h_asym[i] = fitted_params_type3[2,2]
        fitted_r_asym[i] = fitted_params_type3[3,2]
        #fitted_σ_asym = fitted_params_type3[4,2]        

    # set up model object for type 3, refuge model
    model_type3_ref = fun_res_RPE_REF(prey_offered,prey_eaten,starting_a,max_a,starting_h)

    # call the fitting for refuge type 3
    fr_chain_type3_ref = sample(
        model_type3_ref,
        NUTS(5000,0.65),
        MCMCSerial(),
        4000,
        init_params = [(starting_a,starting_h,0.0001*maximum(prey_offered),10)],
        1)
        # put fitted parameters into a dataframe
        fitted_params_type3 = DataFrame(summarystats(fr_chain_type3_ref))
        # error check - if posterior has no variation, don't run loo step
        if fitted_params_type3[1,3] > 0
            loo_type3ref = psis_loo(model_type3_ref, fr_chain_type3_ref)
            cv_elpd_type3ref[i] = loo_type3ref.estimates[1,1]
            naive_lpd_type3ref[i] = loo_type3ref.estimates[1,1]
        end

        # pull parameters and log results
        fitted_a_ref[i] = fitted_params_type3[1,2]
        fitted_h_ref[i] = fitted_params_type3[2,2]
        fitted_k_ref[i] = fitted_params_type3[3,2] =#

    #write(string("chain-file_DS_",string(i),"type2.jls"), fr_chain_type2)
    #write(string("chain-file_DS_",string(i),"type3.jls"), fr_chain_type3)

#chn2 = read("chain-file_DS_3013.jls", Chains)
#plot(chn2)
#summarystats(fr_chain_type3_ref)


# set up a dataframe of all the fitted parameters to write out
col_names = ["Experiment"]
df = DataFrame(Tables.table(exps),col_names)
df[!,"Fitted_a_2"] .= fitted_a_type2'
df[!,"Fitted_h_2"] .= fitted_h_type2'

#=
df[!,"CV_ELPD_2"] .= cv_elpd_type2'

df[!,"Fitted_a_3_aa"] .= fitted_a_asym'
df[!,"Fitted_h_3_aa"] .= fitted_h_asym'
df[!,"Fitted_k_3_aa"] .= fitted_r_asym'
df[!,"CV_ELPD_3_aa"] .= cv_elpd_type3aa'

df[!,"Fitted_a_3_ref"] .= fitted_a_ref'
df[!,"Fitted_h_3_ref"] .= fitted_h_ref'
df[!,"Fitted_k_3_ref"] .= fitted_k_ref'
df[!,"CV_ELPD_3_ref"] .= cv_elpd_type3ref'

CSV.write("Fitted_params_and_LOO.csv",df)
=#

#writeout = [fitted_a_type2',fitted_h_type2']
#CSV.write("a_and_h.csv",Tables.table(fitted_a_type2'))
#CSV.write("r.csv",Tables.table(fitted_r_ref'))

# make a collection of x ranges for the different prey types
xrange = collect(range(0.0, 40, length=50))

gradient = cgrad(:berlin,9,categorical = true) # create a color gradient
# create an array of empty plots to fill up with the loop
allPlots = Array{Plots.Plot{Plots.GRBackend},1}(undef,4)
allPlots_post = Array{Plots.Plot{Plots.GRBackend},1}(undef,2)

allPlots_post[1] = plot(
    ylabel = "Frequency",
    xlabel = "Value",
    size=(600,400),
    dpi=1000,
    title="a")



    # sample the priors
    chain = sample(model_type2, Prior(), 100, 
        chain_type=MCMCChains.Chains,)

    density!(allPlots_post[1], chain[:a], label = "Prior", color=:blue, linewidth = 2, legend=:best)

    #, xlabel="a", xlims=(0.005, 0.18)
    # open up plot canvasylims=[0,10],
allPlots_post[2] = plot(
    ylabel = "Frequency",
    xlabel = "Value",
    size=(600,400),
    dpi=1000,
    title="h")

            #p2 = density(chain[:ρ], label = "Prior", color=:black, linewidth = 2, xlabel="ρ", xlims=(0.005, 0.18), legend=:best)

# loop over the six data sets
for i = 1:4
    # grab data again
    indices = findall(df_data.Experiment .== exps[i])
    prey_eaten = df_data.eaten[indices]
    prey_offered = df_data.Density[indices]
    # open up plot canvas
    allPlots[i] = plot(ylims=[0,30],
        ylabel = "Number of prey eaten",
        xlabel = "Number of prey offered",
        size=(600,400),
        dpi=1000,
        title=exps[i])
    # grab chains of parameters from type 2 fit
    chain_of_a = chain_a[i,:]
    chain_of_h = chain_h[i,:]
    # calculate the fits across all params in the chains
    y = xrange .- lambertw.(chain_of_a' .* chain_of_h' .* xrange .* exp.(-chain_of_a' .* (0.5 .- chain_of_h' .* xrange))) ./ (chain_of_a' .* chain_of_h')
    # plot all of those fits
    plot!(xrange, y,
        label="",
        color=gradient[2],
        alpha=0.01)
    # put the mean fit on top of this
    mean_a = df.Fitted_a_2[i]
    mean_h = df.Fitted_h_2[i]
    y = xrange .- lambertw.(mean_a .* mean_h .* xrange .* exp.(-mean_a .* (0.5 .- mean_h .* xrange))) ./ (mean_a * mean_h)
    plot!(xrange, y,
        label="",
        color=gradient[7],
        linewidth=2)
    # drop the data on top
    scatter!(prey_offered,prey_eaten,label="",
        markerstrokecolor=gradient[5],
        markercolor=gradient[5],
        markersize=5)

    # plot the posteriors on the post canvas
    density!(allPlots_post[1],chain_a[i,:], label = "Posterior", color=:gray, linewidth = 2)
    density!(allPlots_post[2],chain_h[i,:], label = "Posterior", color=:gray, linewidth = 2)

end

figure2 = plot(allPlots[1],allPlots[4],allPlots[2],
    allPlots[3], layout=(2,2), size=(500,600), dpi=1000)

display(figure2)

plot(allPlots_post[1],allPlots_post[2])

figure1 = plot_fr_fit_multiple_curves(chain_a,chain_h,exps,colors_to_use)


# sample the priors
chain = sample(model_type2, Prior(), 100)

p2 = density(chain[:ρ], label = "Prior", color=:black, linewidth = 2, xlabel="ρ", xlims=(0.005, 0.18), legend=:best)
density!(p2, chain_tpc[:ρ], label = "Posterior", color=:gray, linewidth = 2)

p3 = density(chain[:ΔT], label = "Prior", color=:black, linewidth = 2, xlabel="ΔT", xlims=(4, 8), legend=false)
density!(p3, chain_tpc[:ΔT], label = "Posterior", color=:gray, linewidth = 2)

p4 = density(chain[:λ], label = "Prior", color=:black, linewidth = 2, xlabel="λ", xlims=(-2.5, 0), legend=false)
density!(p4, chain_tpc[:λ], label = "Posterior", color=:gray, linewidth = 2)

p5 = density(chain[:tmax], label = "Prior", color=:black, linewidth = 2, xlabel="tmax", xlims=(38, 41), legend=false)
density!(p5, chain_tpc[:tmax], label = "Posterior", color=:gray, linewidth = 2)



# random predator equation plots
# type 2
y = xrange .- lambertw.(fitted_a_type2[i] .* fitted_h_type2[i] .* xrange .* exp.(-fitted_a_type2[i] .* (1 .- fitted_h_type2[i] .* xrange))) ./ (fitted_a_type2[i] * fitted_h_type2[i])
plot!(xrange, y,
    label="Type 2",
    color=:red,
    ylims=(0, 30))
# refuge model
plot!(xrange, (xrange.-fitted_r_ref[i]) .- lambertw.(fitted_a_ref[i] .* fitted_h_ref[i] .* (xrange.-fitted_r_ref[i]) .* exp.(-fitted_a_ref[i] .* (1 .- fitted_h_ref[i] .* (xrange.-fitted_r_ref[i])))) ./ (fitted_a_ref[i] * fitted_h_ref[i]),
    label="Refuge model",
    color=:green)
# asymptotic a
plot!(xrange, (fitted_a_asym[i].*xrange./(fitted_r_asym[i].+xrange)) .* (xrange) ./ (1 .+ (fitted_a_asym[i].*xrange./(fitted_r_asym[i].+xrange)) .* fitted_h_asym[i] .* (xrange)),
    label="Asymptotic a",
    color=:blue)
scatter(prey_offered,prey_eaten,label="")
hline!([1/starting_h],label="Initial guess 1/h")
hline!([1/fitted_h_type2[i]],label="Fitted 1/h",Linewidth=2)
    xlabel!("Prey offered")
    ylabel!("Prey eaten")    

# get the rsquare
SSt = sum((prey_eaten .- mean(prey_eaten)).^2)
y_type2 = prey_offered .- lambertw.(fitted_a_type2[i] .* fitted_h_type2[i] .* prey_offered .* exp.(-fitted_a_type2[i] .* (1 .- fitted_h_type2[i] .* prey_offered))) ./ (fitted_a_type2[i] * fitted_h_type2[i])
y_ref = (prey_offered.-fitted_r_ref[i]) .- lambertw.(fitted_a_ref[i] .* fitted_h_ref[i] .* (prey_offered.-fitted_r_ref[i]) .* exp.(-fitted_a_ref[i] .* (1 .- fitted_h_ref[i] .* (prey_offered.-fitted_r_ref[i])))) ./ (fitted_a_ref[i] * fitted_h_ref[i])
y_asym = (fitted_a_asym[i].*prey_offered./(fitted_r_asym[i].+prey_offered)) .* prey_offered ./ (1 .+ (fitted_a_asym[i].*prey_offered./(fitted_r_asym[i].+prey_offered)) .* fitted_h_asym[i] .* prey_offered)
SSx_2 = sum((y_type2 .- mean(prey_eaten)).^2)
SSx_ref = sum((y_ref .- mean(prey_eaten)).^2)
SSx_asym = sum((y_asym .- mean(prey_eaten)).^2)
rsquare_2 = SSx_2 / SSt
rsquare_ref = SSx_ref / SSt
rsquare_asym = SSx_asym / SSt



plot!(xrange, fitted_a_type2[i] .* (xrange) ./ (1 .+ fitted_a_type2[i] .* fitted_h_type2[i] .* (xrange)),
    label="",
    color=:blue)

plot(xrange, fitted_a_type2[i] .* xrange ./ (1 .+ fitted_a_type2[i] .* fitted_h_type2[i] .* xrange),
    label="",
    color=:blue)
    plot(xrange, xrange .- lambertw.(fitted_a_type2[i] .* fitted_h_type2[i] .* xrange .* exp.(-fitted_a_type2[i] .* (1 .- fitted_h_type2[i] .* xrange))) ./ (fitted_a_type2[i] * fitted_h_type2[i]),
    label="",
    color=:red)





###### working out the MSFR fits

# MODEL 6b -- NULL VECTOR type 3 refuge model, modified from RPE
@model function fun_res_multi_spp(null_vector,prey1_offered,prey1_eaten,prey2_offered,prey2_eaten,T,starting_a,starting_h) 
    a1 ~ truncated(Normal(starting_a[1],10*starting_a[1]),0,1) # space clearance rate
    h1 ~ truncated(Normal(starting_h,10*starting_h),0,100) # handling time
    a2 ~ truncated(Normal(starting_a[1],10*starting_a[1]),0,1) # space clearance rate
    h2 ~ truncated(Normal(starting_h,10*starting_h),0,100) # handling time
    σ ~ InverseGamma(2,3)
    for i in 1:length(prey1_eaten)

        null_vector[1,i] ~ Normal(prey1_offered[i]*(1-exp(a1*(h2*prey2_eaten[i]+h1*prey1_eaten[i]-T)))-prey1_eaten[i],σ)
        null_vector[2,i] ~ Normal(prey2_offered[i]*(1-exp(a2*(h2*prey2_eaten[i]+h1*prey1_eaten[i]-T)))-prey2_eaten[i],σ)

	end
end

# read in the data
msfr_data = CSV.read("FunctionalResponse_msfr.csv",DataFrame)
prey1_offered = msfr_data.sinella_init
prey2_offered = msfr_data.drosophila_init
prey1_eaten = msfr_data.sinella_eaten
prey2_eaten = msfr_data.drosophila_eaten

# pull starting values from dataset
max_a = maximum([prey1_eaten; prey2_eaten]./[prey1_offered; prey2_offered])
starting_a = 0.5*max_a
starting_h = 1/maximum([prey1_eaten; prey2_eaten])
null_vector = zeros(2,size(prey1_offered,1))

# set up model object for type 3, refuge model
model_MSFR = fun_res_multi_spp(null_vector,prey1_offered,prey1_eaten,prey2_offered,prey2_eaten,6,starting_a,starting_h)

# call the fitting for refuge type 3
fr_chain_MSFR = sample(
    model_MSFR,
    NUTS(5000,0.65),
    MCMCSerial(),
    4000,
    init_params = [(starting_a,starting_h,starting_a,starting_h,10)],
    2)

    plot(fr_chain_MSFR)
    
    scatter3d(prey1_offered,prey2_offered,prey1_eaten,label="Sinella")
    scatter3d!(prey1_offered,prey2_offered,prey2_eaten,label="Drosophila")
    xlabel!("Sinella offered")    
    ylabel!("Drosophila offered")
    zlabel!("Prey eaten")  

