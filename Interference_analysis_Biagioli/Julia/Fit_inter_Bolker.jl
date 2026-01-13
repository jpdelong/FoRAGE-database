#Set Directory
dir = "~/Interference/" 
cd(dir)

#Load Packages
using DataFrames, CSV, Statistics
using Distributions
using Turing
using StatsPlots
using Plots
using LsqFit
using Roots
using LambertW

#Confirm functions
plotlyjs()
gr()
lambertw(1)

# --- Set Output directories ---
trace_dir = "Figures/Fit Trace Plots Bolker/"
surface_dir = "Figures/FR Surfaces Bolker/"
obs_v_pred_dir = "Figures/Obs_v_Predict_Bolker/"
animated_surface = "Figures/Animated FR Surface_Bolker"
isdir(trace_dir) || mkpath(trace_dir)
isdir(surface_dir) || mkpath(surface_dir)
isdir(obs_v_pred_dir) || mkpath(obs_v_pred_dir)
isdir(animated_surface) || mkpath(animated_surface)

#Load in Data simulation function
include("simulate_data_set_INT.jl")

# --- Read in meta and curve data ---
df_meta = CSV.read("FoRAGE_V5_sources_and_meta_inter.csv", DataFrame)
dimensions = df_meta.Dim
pred_nums = df_meta.PredPerArena
data_types = df_meta.DataOrigin
Arena_size_2D = df_meta.TwoD_Arena_size_cm2
Arena_size_3D = df_meta.ThreeD_arena_size_cm3
replenished = df_meta.PreyReplaced

df_curves = CSV.read("FoRAGE_db_V5_Dec_20_2024_original_curves.csv", DataFrame)
density_units_2D = df_curves.TwoDDensityUnits
density_units_3D = df_curves.ThreeDDensityUnits

# --- Models ---
@model function fun_res_BDE_HDE_null(null_vector, prey_offered, prey_eaten, pred_numb, starting_a, starting_h, starting_w)
    a ~ truncated(Normal(starting_a[1], 10*starting_a[1]), lower = 0)
    h ~ truncated(Normal(starting_h[1], 2*starting_h[1]), lower = 0)
    w ~ truncated(Normal(starting_w[1], 10*starting_w[1]), lower = 0)
    σ ~ InverseGamma(2, 1)
    for i in 1:length(prey_eaten)
        null_vector[i] ~ Normal(a*prey_offered[i]/(1 + a*h*prey_offered[i] + w*pred_numb[i]) - prey_eaten[i], σ)
    end
end

@model function fun_res_BDE_RPE_null(null_vector, prey_offered, prey_eaten, pred_numb, T, starting_a, starting_h, starting_w)
    a ~ truncated(Normal(starting_a[1], 10*starting_a[1]), lower = 0)
    h ~ truncated(Normal(starting_h[1], 2*starting_h[1]), lower = 0)
    w ~ truncated(Normal(starting_w[1], 5*starting_w[1]), lower = 0)
    σ ~ InverseGamma(2, 1)
    for i in 1:length(prey_eaten)
        null_vector[i] ~ Normal(prey_offered[i]*(1 - exp((a*(h*prey_eaten[i] - T)) / (w*pred_numb[i] + 1))) - prey_eaten[i], σ)
    end
end


# MODEL 2 - Bolker Beddington-DeAngelis LambertW function
@model function fun_res_BBDEL(prey_offered,prey_eaten,pred_numb,starting_a,starting_h,starting_w,AS)
    a ~ truncated(Normal(starting_a[1], 10*starting_a[1]), lower = 0)
    h ~ truncated(Normal(starting_h[1], 2*starting_h[1]), lower = 0)
    w ~ truncated(Normal(starting_w[1], 10*starting_w[1]), lower = 0)
    σ ~ InverseGamma(2, 1)
    for i in 1:length(prey_eaten)
        ϕ = a/(1+w*(pred_numb[i]-1/AS)) # AS is arena size
        #ϕ = a/(1+w*(pred_numb[i])) # AS is arena size
        prey_eaten[i] ~ truncated(Normal(prey_offered[i]-LambertW.lambertw(ϕ*h*prey_offered[i]*exp(-ϕ*(1-h*prey_offered[i])))/(ϕ*h),σ),0,maximum(prey_eaten))
    end
end

# Beddington-DeAngelis analytic prediction (vectorized)
bd_predict(a, h, w, N0_vec, P_vec) = (a .* N0_vec) ./ (1 .+ a .* h .* N0_vec .+ w .* P_vec)

function rpe_predict(RE, p)
    a, h, w, RO, P = p # parameters
    RO*(1 - exp((a*(h*RE - 1)) / (w*P + 1))) - RE
end


# --- Store fitted parameters ---
fitted_a = zeros(maximum(df_meta.Inter_ID), 1)
fitted_h = zeros(maximum(df_meta.Inter_ID), 1)
fitted_w = zeros(maximum(df_meta.Inter_ID), 1)

# --- Create Vector to store prey_scale_adjust for each study
prey_scale_adjust_vec = zeros(43)

#declair number of studies, iterations, chains, and total posterior draws and build empty posterior matricies
n_studies = 43
n_draws_per_chain = 20000
n_chains = 3
n_draws = n_draws_per_chain * n_chains

W_mat = fill(NaN, n_studies, n_draws)
A_mat = fill(NaN, n_studies, n_draws)
H_mat = fill(NaN, n_studies, n_draws)

# --- Fit models across studies ---
for i = 1:43
    # --- Skip study 38 and fill outputs with NaN ---
    if i == 38
        println("Skipping study 38")

        fitted_a[i] = NaN
        fitted_h[i] = NaN
        fitted_w[i] = NaN
        prey_scale_adjust_vec[i] = NaN

        continue
    end
    
    println(i)
    #initialize empty vectors
        prey_offered_2 = Vector{Float64}()
        prey_offered_3 = Vector{Float64}()
        indices = Vector{Int}()
        preds = Vector{Float64}()

    # read in data from meta sheet
        pred_numbs = df_meta.PredPerArena[findall(df_meta.Inter_ID .== i)] # number of preds per arena
        which_ID = df_meta.FoRAGE_ID[findall(df_meta.Inter_ID .== i)] #identifyer linking Inter_ID to FoRAGE ID

    #Pull dimension and arena sizes from metadata
        dimension = dimensions[findall(df_meta.Inter_ID .== i)][1]
        arena_size_2D = Arena_size_2D[findall(df_meta.Inter_ID .== i)][1]
        arena_size_3D = Arena_size_3D[findall(df_meta.Inter_ID .== i)][1]

    # --- compute AS once per study `i` ---
    if dimension == 2
        AS = arena_size_2D
    elseif dimension == 3
        AS = arena_size_3D
    elseif dimension == 2.5
        if !ismissing(arena_size_3D)
            AS = arena_size_3D^(2.5/3)
        else
            AS = arena_size_2D^(2.5/2)
        end
    end

    #find all rows in curves dataset with matching FoRAGE ID
    for j = 1:length(which_ID) #for each set of rows with unique FoRAGE IDs 
        toappend = findall(df_curves.Column1 .== which_ID[j]) # find all columns in curves dataset with the same FoRAGE ID
        append!(indices, toappend) # concatenates all rows within the curves dataset that belong to the same Inter_ID into a single indices

        # --- Compute predator density depending on arena dimension ---
        if dimension == 2
            # Keep cm² t units
            density_list = (pred_numbs[j] ./ arena_size_2D) .* ones(length(toappend), 1) # calculates the predator density for each curves dataset row for each Inter_ID, calculating the correct predator denity for each FoRAGE ID

        elseif dimension == 3
            # Keep cm³ units
            density_list = (pred_numbs[j] ./ arena_size_3D) .* ones(length(toappend), 1) #same as above but for 3D studies

        elseif dimension == 2.5
        # approximate “2.5D” arena using geometric scaling
            #Keep cm^2.5 units
            if !ismissing(arena_size_3D)
                arena_size_2p5D = (arena_size_3D)^(2.5 / 3)
                density_list = (pred_numbs[j] ./ arena_size_2p5D) .* ones(length(toappend), 1) # same as above, but converts 2D to 2.5D
            else
                arena_size_2p5D = (arena_size_2D)^(2.5 / 2)
                density_list = (pred_numbs[j] ./ arena_size_2p5D) .* ones(length(toappend), 1) # same as above, but converts 3D to 2.5D
            end
             # Diagnostic printout for studies
            println("Study $i:")
            println("  Dim = $dimension")
            println("  Preds = ", pred_numbs[j])
            println("  Arena size raw (cm²) = ", arena_size_2D)
            println("  Arena size raw (cm3) = ", arena_size_3D)
            println("  Density example = ", density_list[1])
            println("  ---")
        end
        append!(preds, density_list) #appends values from current study i to the complete list of predator density accross all looped studies 1:i in preds
    end

    # --- Call in data from FoRAGE curve dataset ---
        prey_offered_2 = df_curves.TwoDPreyDensity[indices]
        prey_offered_3 = df_curves.ThreeDPreyDensity[indices]
        prey_eaten = df_curves.ForagingRate[indices]
        std_err = df_curves.StandardError[indices]
        samp_size = df_curves.SampleSize[indices]

    # Estimate missing std errors or sample sizes
    if ismissing(std_err[1])
        std_err = exp(-1.8886) .* prey_eaten .^ 0.95
    end
    if ismissing(samp_size[1])
        samp_size = 3 .* ones(length(prey_eaten), 1)
    end

    # Unit conversions from m to cm
    # m^2 to cm^2
    units_2D = density_units_2D[indices][1]
    if !ismissing(units_2D) && units_2D == "prey per m2"
        prey_offered_2 = prey_offered_2 ./ 1e4
    end
    # m^3 to cm^3
    units_3D = density_units_3D[indices][1]
    if !ismissing(units_3D) && units_3D == "prey per m3"
        prey_offered_3 = prey_offered_3 ./ 1e6
    end

    # Handle 2D / 2.5D / 3D data
    if dimension[] == 2
        prey_offered = prey_offered_2
    elseif dimension[] == 3
        prey_offered = prey_offered_3
    elseif dimension[] == 2.5
        if ismissing(prey_offered_3[1])
            prey_offered = prey_offered_2 .^ (2.5/2)
            preds = preds .^ (2.5/2)
        else
            prey_offered = prey_offered_3 .^ (2.5/3)
            preds = preds .^ (2.5/3)
        end
    end

    # Simulate if data are means
    data_type = data_types[findall(df_meta.Inter_ID .== i)][1]
    if data_type == "Mean"
        prey_eaten, prey_offered, preds = simulate_data_set_INT(length(prey_eaten), std_err, samp_size, prey_eaten, prey_offered, preds)
    end

    # Check and rescale eaten/offered
    #This block ensures that the number of prey offered does not exceed the number of prey available
    prey_scale_test = prey_eaten[prey_offered .> 0] ./ prey_offered[prey_offered .> 0]
    if sum(prey_scale_test .> 1) > 0
        prey_scale_adjust = 2 * maximum(prey_scale_test)
    elseif mean(prey_scale_test) < 1e-5
        prey_scale_adjust = 1e-2
    else
        prey_scale_adjust = 1
    end
        prey_offered = prey_offered .* prey_scale_adjust

    # store prey_scale_adjust for each study in empty vector
    prey_scale_adjust_vec[i] = prey_scale_adjust

    # Starting values
    max_a = maximum(prey_eaten ./ prey_offered)
    starting_a = max_a #*0.5
    starting_h = (1 / maximum(prey_eaten))

    # estimate with nonlinear inverse fit
        @. int_estimate(x, p) = 1/(p[1]*x) # model
        fit = curve_fit(int_estimate, preds, prey_eaten, [10.0])
        starting_w = fit.param[1]

    # Choose model type
    replenished = df_meta.PreyReplaced[findall(df_meta.Inter_ID .== i)][1]
    null_vector = zeros(size(prey_offered))
    if replenished == "Y"
        model_type2 = fun_res_BDE_HDE_null(null_vector, prey_offered, prey_eaten, preds, starting_a, starting_h, starting_w)
    else
        #model_type2 = fun_res_BDE_RPE_null(null_vector, prey_offered, prey_eaten, preds, 1, starting_a, starting_h, starting_w)
        model_type2 = fun_res_BBDEL(prey_offered, prey_eaten, preds, starting_a, starting_h, starting_w,AS)
    end

    # --- Fit model ---
    fr_chain = sample(
        model_type2,
        NUTS(20000, 0.8), #num of warmup iterations and target acceptance rate
        MCMCThreads(),
        20000, # number of sampling iterations
        init_params = [(
        starting_a, # inital value for a
        starting_h, # inital value for h
        starting_w, #initial value for w
        0.5) # initial value for sigma
        for _ in 1:3],
        3 # num of chains
    )

    #save posterior matricies
    df_chain = DataFrame(fr_chain)

    #pull posterior draws from df_chain for each par
    a_draws = df_chain.a
    h_draws = df_chain.h
    w_draws = df_chain.w

    #scale a draws
    scaled_a_draws = a_draws .* prey_scale_adjust

    #populate matricies with posterior draws
    A_mat[i, :] = scaled_a_draws
    H_mat[i, :] = h_draws
    W_mat[i, :] = w_draws

    # Save trace plots
    plt_trace = plot(fr_chain)
    savefig(plt_trace, joinpath(trace_dir, "Inter_ID_$(i).png"))

    # Extract parameter summaries
    fitted_params = DataFrame(summarystats(fr_chain))
    fitted_a[i] = fitted_params[1,2]
    fitted_h[i] = fitted_params[2,2]
    fitted_w[i] = fitted_params[3,2]

    a_fit, h_fit, w_fit = fitted_a[i], fitted_h[i], fitted_w[i]

    # --- 3D Surface Plot (use same model used for fitting) ---
    prey_range = range(minimum(prey_offered), stop=maximum(prey_offered), length=50)
    pred_range = range(minimum(preds), stop=maximum(preds), length=50)

    # Build Z depending on model used
    if replenished == "Y"
        # BD formula is analytic -> fast
        Z = [a_fit * p / (1 + a_fit*h_fit*p + w_fit*d) for p in prey_range, d in pred_range]
    else
        # RPE: need to solve Rogers eqn for each (p,d)
        N0_vec = collect(prey_range)
        P_vec = collect(pred_range)
        # Initialize Z
        Z = zeros(length(N0_vec), length(P_vec))
        for q = 1:length(N0_vec)
            for z = 1:length(P_vec)
                RO = N0_vec[q]
                P = P_vec[z]
                params = (a_fit, h_fit, w_fit, RO, P)
                init_guess = bd_predict(a_fit, h_fit, w_fit, RO, P)
                pforaging = find_zero(rpe_predict, init_guess, p = params)
                Z[q,z] = pforaging
            end
        end
    end

# ----- ANIMATED SURFACE PLOT -----

    # Create a surface plot 
    plt_surface = surface(
        prey_range, pred_range, Z';
        xlabel = "Prey Density",
        ylabel = "Predator Density",
        zlabel = "Prey Eaten",
        title = "Functional Response Fit (Inter_ID $(i))",
        legend = false,
        alpha = 0.7,
        zlims = (0, 1.1 * maximum(prey_eaten)),
        color = :viridis
    )

    # Overlay observed data points
    scatter3d!(plt_surface, prey_offered, preds, prey_eaten;
        markersize = 5, markercolor = :red, label = "Observed Data"
    )

    # Animate the rotation
    anim = @animate for angle in 0:2:360
        plot!(plt_surface, camera = (angle, 30))
    end

    # Save as GIF
    gif(anim, joinpath(animated_surface, "FR_surface_Inter_ID_$(i).gif"), fps = 15)

    
end

# --- Save fitted parameters ---
df_fitted = DataFrame(
    Inter_ID = 1:maximum(df_meta.Inter_ID),
    fitted_a = vec(fitted_a),
    scale_adjust_factor = prey_scale_adjust_vec,
    scaled_a = vec(fitted_a) .* prey_scale_adjust_vec,
    fitted_h = vec(fitted_h),
    fitted_w = vec(fitted_w)
)

# --- Save estimated pars in CSV ---
CSV.write("Inter_fitted_parameters_Bolker.csv", df_fitted)

CSV.write("posterior_w_matrix.csv", DataFrame(W_mat, :auto))
CSV.write("posterior_scaled_a_matrix.csv", DataFrame(A_mat, :auto))
CSV.write("posterior_h_matrix.csv", DataFrame(H_mat, :auto))
