# pick a location for your files
dir = "C:/Users/jdelong2/OneDrive - University of Nebraska-Lincoln/Projects in progress/Functional responses all species analysis/type2type3/"
cd(dir)

using Plots, Statistics, StatsBase, StatsPlots, MixedModels, Plots.PlotMeasures
using DataFrames, CSV, MAT, Distributions, KernelDensity
using GLM

df = CSV.read("FoRAGE_V4_Sept_27_2023_working_type2type3.csv",DataFrame)
Dataset = df.Dataset
df = dropmissing(df,:Dataset)
df.mediandeltaAIC = replace(df.mediandeltaAIC, NaN=>missing)
df = dropmissing(df,:mediandeltaAIC)

# restrict data we use to FRs with more than a minimum amount of data
df = df[findall(df.NumTrials .> 10),:]

# plot histogram of aic differences
f1 = histogram(df.mediandeltaAIC,bins=140,xlims=(-10, 25),ylims=(0,1500),color = "white",
  xlabel = "Type 2 AIC - Type 3 AIC",ylabel = "Frequency",label ="",
  size=(300,300),dpi=600)
# plots a vertical line at [#]
vline!([-2], linewidth=3, color="blue1", label = "-2")
vline!([0], linewidth=3, color="black", label = "0")
vline!([2], linewidth=3, linestyle=:dot, color="blue1",label = "2")
text(f1,"Figure1")
savefig(f1,"Figure1_evidence.pdf")

# designate type 3 yes or no
df.type3yesno .= 0
df.type3yesno[findall(df.mediandeltaAIC .> 2)] .= 3
df.type3yesno[findall(df.mediandeltaAIC .< -2)] .= 2

findall(df.mediandeltaAIC .< 0)
findall(df.mediandeltaAIC .> -2 .&& df.mediandeltaAIC .< 2)

# calculate the predator prey size ratio
df.PPRatio = df.PredatorMass./df.PreyMass

# figure 2 for ms
p1 = scatter(df.NumTrials,df.meandeltaAIC,markerstrokecolor="blue4",markercolor="blue4",alpha=0.2,
  xlabel = "Number of trials",ylabel = "AIC diff",legend=false,size=(200,200))
  annotate!(25, 100, ["A"])
p2 = scatter(df.FullyEaten,df.meandeltaAIC,markerstrokecolor="blue4",markercolor="blue4",alpha=0.2,
  xlabel = "Proportion fully eaten",ylabel = "",legend=false,size=(300,200))
  annotate!(0.05, 100, ["B"])
p3 = scatter(df.Temp,df.meandeltaAIC,markerstrokecolor="blue4",markercolor="blue4",alpha=0.2,
  xlabel = "Temperature",ylabel = "",legend=false,size=(200,200))
  annotate!(1, 100, ["C"])
p4 = scatter(log.(df.PPRatio),df.meandeltaAIC,markerstrokecolor="blue4",markercolor="blue4",alpha=0.2,
  xlabel = "PPRatio (ln)",ylabel = "AIC diff",legend=false,size=(200,200))
  annotate!(-3, 100, ["D"])
p5 = scatter(log.(df.PreyMass),df.meandeltaAIC,markerstrokecolor="blue4",markercolor="blue4",alpha=0.2,
  xlabel = "Prey mass (ln)",ylabel = "",legend=false,size=(200,200))
  annotate!(-20, 100, ["E"])
p6 = scatter(df.Habitat,df.meandeltaAIC,xlabel = "Habitat",markerstrokecolor="blue4",markercolor="blue4",alpha=0.2,
  ylabel = "",legend=false,size=(300,200),
  right_margin = 10px)
  annotate!(0.56, 100, ["F"])

figure2 = plot(p1,p2,p3,p4,p5,p6, layout = (2,3),size=(600,400),
  right_margin = 20px, dpi=600)
savefig(figure2,"Figure2.pdf")

# make plots of other BIOLOGICAL predictors against AIC_diff for supplement
p1 = scatter(log.(df.PredatorMass),df.meandeltaAIC,markerstrokecolor="blue4",markercolor="blue4",alpha=0.2,
  xlabel = "Predator mass",ylabel = "AIC diff",legend=false,size=(200,200))
  annotate!(-15, 100, ["A"])
p2 = scatter(df.PreyMajorGrouping,df.meandeltaAIC,markerstrokecolor="blue4",markercolor="blue4",alpha=0.2,
  xlabel = "Prey group",ylabel = "AIC diff",legend=false,size=(400,400),xrotation =90)
  annotate!(1, 100, ["B"])
p3 = scatter(df.PredMajorGrouping,df.meandeltaAIC,markerstrokecolor="blue4",markercolor="blue4",alpha=0.2,
  xlabel = "Predator group",ylabel = "AIC diff",legend=false,size=(400,400),xrotation =90)
  annotate!(1, 100, ["C"])

figureS1 = plot(p1,p2,p3, layout = (3,1),size=(400,500), dpi=600)
savefig(figureS1,"figureS1.png")

# building models
df_model = select(df,:mediandeltaAIC,:Temp,:PPRatio,:PredatorMass,:PreyMass,:PredMajorGrouping,:PreyMajorGrouping,
  :Habitat,:FullyEaten,:NumTrials)
df_model = dropmissing(df_model,:)

# use this to grab r2, only works for lm not mixed models
r2(lm_temp) 

# build single factor models for experimental stuff
lm_temp = lm(@formula(mediandeltaAIC ~ FullyEaten^2), df_model)
lm_temp = lm(@formula(mediandeltaAIC ~ NumTrials^2), df_model)
lm_temp = lm(@formula(abs(mediandeltaAIC) ~ NumTrials^2), df_model)

# build single factor linear models with statistical predictors as covariate
lm_temp = lm(@formula(mediandeltaAIC ~ Temp + FullyEaten^2 + NumTrials^2), df_model)
lm_temp = lm(@formula(mediandeltaAIC ~ log(PPRatio) + FullyEaten^2 + NumTrials^2), df_model)
lm_temp = lm(@formula(mediandeltaAIC ~ log(PredatorMass) + FullyEaten^2 + NumTrials^2), df_model)
lm_temp = lm(@formula(mediandeltaAIC ~ log(PreyMass) + FullyEaten^2 + NumTrials^2), df_model)
lm_temp = lm(@formula(mediandeltaAIC ~ Habitat + FullyEaten^2 + NumTrials^2), df_model)
lm_temp = lm(@formula(mediandeltaAIC ~ PredMajorGrouping + FullyEaten^2 + NumTrials^2), df_model)
lm_temp = lm(@formula(mediandeltaAIC ~ PreyMajorGrouping + FullyEaten^2 + NumTrials^2), df_model)

# do a check for arena size effects for 2D systems
df_arena = select(df,:mediandeltaAIC,:ArenaSize2D,:FullyEaten,:NumTrials)
df_arena = dropmissing(df_arena,:)
lm_arena = lm(@formula(mediandeltaAIC ~ log(ArenaSize2D) + FullyEaten^2 + NumTrials^2), df_arena)
# do a check for arena size effects for 3D systems
df_arena = select(df,:mediandeltaAIC,:ArenaSize3D,:FullyEaten,:NumTrials)
df_arena = dropmissing(df_arena,:)
lm_arena = lm(@formula(mediandeltaAIC ~ log(ArenaSize3D) + FullyEaten^2 + NumTrials^2), df_arena)

# build multiple regression model of all main effects
#lm_temp = lm(@formula(mediandeltaAIC ~ Temp * log(PPRatio) * Habitat + NumTrials^2 + FullyEaten^2), df_model)
#lm_temp = lm(@formula(mediandeltaAIC ~ Temp + log(PPRatio) * Habitat + NumTrials^2 + FullyEaten^2), df_model)
#lm_temp = lm(@formula(mediandeltaAIC ~ Temp * log(PPRatio) + NumTrials^2 + FullyEaten^2), df_model)

lm_temp = lm(@formula(mediandeltaAIC ~ Temp + log(PPRatio) + Habitat + NumTrials^2 + FullyEaten^2), df_model)
mainr2 = r2(lm_temp)
# sequentially drop each in turn and grab marginal r2
lm_temp = lm(@formula(mediandeltaAIC ~ log(PPRatio) + Habitat + NumTrials^2 + FullyEaten^2), df_model)
mainr2 - r2(lm_temp)
lm_temp = lm(@formula(mediandeltaAIC ~ Temp + Habitat + NumTrials^2 + FullyEaten^2), df_model)
mainr2 - r2(lm_temp)
lm_temp = lm(@formula(mediandeltaAIC ~ Temp + log(PPRatio) + NumTrials^2 + FullyEaten^2), df_model)
mainr2 - r2(lm_temp)
lm_temp = lm(@formula(mediandeltaAIC ~ Temp + log(PPRatio) + Habitat + NumTrials^2), df_model)
mainr2 - r2(lm_temp)
lm_temp = lm(@formula(mediandeltaAIC ~ Temp + log(PPRatio) + Habitat + FullyEaten^2), df_model)
mainr2 - r2(lm_temp)

# make some regression diagnostic plots
qq1 = histogram(residuals(lm_temp),ylabel = "Frequency",
  xlabel = "Residual", legend=false)
qq2 = scatter(fitted(lm_temp),residuals(lm_temp),ylabel = "Residuals",
  xlabel = "Fitted", legend=false)
qq3 = qqplot(fitted(lm_temp),residuals(lm_temp),qqline = :fit,ylabel = "Residual",
  xlabel = "Fitted", legend=false)

qq = plot(qq1,qq2,qq3)
savefig(qq,"diagnostics.png")
