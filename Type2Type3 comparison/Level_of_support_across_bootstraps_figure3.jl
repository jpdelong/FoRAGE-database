# pick a location for your files
dir = "C:/Users/jdelong2/OneDrive - University of Nebraska-Lincoln/Projects in progress/Functional responses all species analysis/type2type3/"
cd(dir)
df = CSV.read("FoRAGE_V4_Sept_27_2023_working_type2type3.csv",DataFrame)
df.type3yesno .= 0
df.medianAICdiff .= 0.0
df.list_not_three .= 0.0

# pick a location for your files
dir = "C:/Users/jdelong2/OneDrive - University of Nebraska-Lincoln/Projects in progress/Functional responses all species analysis/Fits/"
cd(dir)

# this particular code runs through each fit to gather AIC data
# from both type 2 and 3 fits. To use it, you need to have all
# the fits in separate files in the same place (>3000 files).
num_studies = 3013

for i = 1:num_studies
    file_to_open = string("DS_",string(i),"_fits_type2.mat")
    file = matopen(file_to_open)
    DS_AIC_2 = read(file, "AIC")
    close(file)

    file_to_open = string("DS_",string(i),"_fits_type3.mat")
    file = matopen(file_to_open)
    DS_AIC_3 = read(file, "AIC")
    close(file)

    DS_diffs = DS_AIC_2 .- DS_AIC_3 
    df.list_not_three[i] = sum(DS_diffs .< 2)
    df.medianAICdiff[i] = median(DS_diffs)

end

figure3 = scatter(df.medianAICdiff,df.list_not_three,markerstrokecolor="blue4",markercolor="blue4",alpha=0.2,
    xlabel = "Median ΔAIC",ylabel = "Number of ΔAIC < 2",
    legend=false,size=(400,300), dpi=600)
vline!([2], linewidth=3, color="orange3",label = "2")
hline!([5], linewidth=3, color="purple4",label = "2")
savefig(figure3,"Figure3.pdf")

