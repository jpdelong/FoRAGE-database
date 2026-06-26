dir = "C:/Users/johnp/Downloads/"
cd(dir)

using Turing
using MCMCChains
using StatsPlots
using Distributions
using Random
using DataFrames, LambertW, CSV, Tables
using Serialization

df_data = CSV.read("copepods.csv",DataFrame)
exps = unique(df_data.Experiment)

# MODEL 2 - Rogers random predator
@model function fun_res_RPE(prey_offered,prey_eaten,starting_a,max_a,starting_h) 
    #a ~ LogNormal(log(max_a) - (0.5^2) / 2, 2*starting_a) # space clearance rate
    a ~ truncated(Normal(starting_a,10*starting_a),0,20) # space clearance rate
    h ~ truncated(Normal(starting_h,10*starting_h),0,2) # handling time
    σ ~ InverseGamma(2,3)
    for i in 1:length(prey_eaten)
        prey_eaten[i] ~ truncated(Normal(prey_offered[i]-lambertw(a*h*prey_offered[i]*exp(-a*(0.5-h*prey_offered[i])))/(a*h),σ),0,maximum(prey_eaten))
	end
end

# MODEL 6 - type 3 refuge model, modified from RPE
@model function fun_res_RPE_REF(prey_offered,prey_eaten,starting_a,max_a,starting_h) 
    a ~ truncated(Normal(starting_a,10*starting_a),0,20) # space clearance rate
    h ~ truncated(Normal(starting_h,10*starting_h),0,2) # handling time
    k ~ Normal(0,0.1*maximum(prey_offered)) # refuge size
	#σ ~ LogNormal(1,2)
    σ ~ InverseGamma(2,3)
    for i in 1:length(prey_eaten)
        prey_eaten[i] ~ truncated(Normal((prey_offered[i]-k)-lambertw(a*h*(prey_offered[i]-k)*exp(-a*(0.5-h*(prey_offered[i]-k))))/(a*h),σ),0,maximum(prey_eaten))
	end
end



#= MODEL 4 - type 3 asymptotic a model, modified from RPE
@model function fun_res_RPE_AA(prey_offered,prey_eaten,starting_a,max_a,starting_h) 
    a ~ truncated(Normal(starting_a[1],10*starting_a[1]),0,1) # space clearance rate
    h ~ truncated(Normal(starting_h,10*starting_h),0,100) # handling time
    k ~ truncated(Normal(0,0.1*maximum(prey_offered)),0,maximum(prey_offered)) # refuge size
    σ ~ InverseGamma(2,3)
    for i in 1:length(prey_eaten)
        prey_eaten[i] ~ truncated(Normal(prey_offered[i]-lambertw((a*prey_offered[i]/(k+prey_offered[i]))*h*prey_offered[i]*exp(-(a*prey_offered[i]/(k+prey_offered[i]))*(6-h*prey_offered[i])))/((a*prey_offered[i]/(k+prey_offered[i]))*h),σ),0,maximum(prey_eaten))
	end
end=#


########################################################################
# use a common prior from empirical starting points
    max_a = maximum(df_data.eaten./df_data.Density)
    starting_a = 0.5*max_a
    starting_h = 1/maximum(df_data.eaten)

# open an empty vector to log refuge model zero-crossing of k (= p value)
    k_crossing_zero = zeros(Float64,4)

########################################################################
# loop over the experiments
for i = 1:length(exps)
    which_exp = exps[i]
    println(which_exp)
    indices = findall(df_data.Experiment .== exps[i])
    prey_eaten = df_data.eaten[indices]
    prey_offered = df_data.Density[indices]

    ######################################
    # set up model object for type 2
    model_type2 = fun_res_RPE(prey_offered,prey_eaten,starting_a,max_a,starting_h)

    # call the fitting for type 2
    fr_chain_type2 = sample(
    model_type2,
    NUTS(30000,0.85),
    MCMCSerial(),
    chain_type=MCMCChains.Chains,
    3000,
    init_params = [(starting_a,starting_h,0.1)],
    4)

    serialize(string("chain-file_DS_",string(i),"type2.jls"), fr_chain_type2)

    ######################################
    # set up model object for type 3, refuge model

    # set new starting values from the type 2 fit
    fitted_params_type2 = DataFrame(summarystats(fr_chain_type2))
    starting_a = fitted_params_type2[1,2]
    starting_h = fitted_params_type2[2,2]

    model_type3_ref = fun_res_RPE_REF(prey_offered,prey_eaten,starting_a,max_a,starting_h)

    # call the fitting for refuge type 3
    fr_chain_type3_ref = sample(
        model_type3_ref,
        NUTS(30000,0.85),
        MCMCSerial(),
        chain_type=MCMCChains.Chains,
        3000,
        init_params = [(starting_a,starting_h,0.0001*maximum(prey_offered),10)],
        4)

    serialize(string("chain-file_DS_",string(i),"type3_ref.jls"), fr_chain_type3_ref)

    # find the quantile for the number of chain samples that are negative
    my_ecdf = ecdf(fr_chain_type3_ref[:k][:])
    k_crossing_zero[i] = my_ecdf(0.0)

end

restored_chain = deserialize("chain-file_DS_4type2.jls")
plot(restored_chain)
summarystats(restored_chain)

describe(fr_chain_type3_ref)
plot(fr_chain_type3_ref)



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

