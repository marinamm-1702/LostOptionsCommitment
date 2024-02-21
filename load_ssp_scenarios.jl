#-------------------------
# Historic emission data from 1765 to 1989
#-------------------------
data_historic = readdlm("emissions/historic");
historic_yrs = data_historic[:,1]
historic_em = data_historic[:,2] + data_historic[:,3]
#-------------------------
# SSP scenarios
#-------------------------
data_ssp = readdlm("emissions/SSP/ssp");
ssp_names = data_ssp[2:end,1];
ssp_temp = [data_ssp[i,5:end]*12*1e-3/((12+32)) for i in 2:5]
ssp_time = data_ssp[1,5:end]
yearsafter = [i for i in 2101:2250]
yearsafter2 = [i for i in 2251:12000]
allyears = vcat(historic_yrs,ssp_time,yearsafter,yearsafter2)
ssp_itp = []
for i in 1:4
    em2100 = ssp_temp[i][end]
    emafter = [-em2100*(t-2100)/150 + em2100 for t in yearsafter]
    emafter2 = zeros(length(yearsafter2))
    allem= vcat(historic_em,ssp_temp[i],emafter,emafter2)
    em = LinearInterpolation(allyears,allem);
    push!(ssp_itp,em)
end
ssp_names = copy(ssp_names[[2,3,1,4]])
ssp_itp = copy(ssp_itp[[2,3,1,4]]);

# Clear temp variables
data_ssp = nothing
ssp_temp = nothing
ssp_time = nothing
historic_yrs = nothing
historic_em = nothing
yearsafter = nothing
yearsafter2 = nothing
allyears = nothing
