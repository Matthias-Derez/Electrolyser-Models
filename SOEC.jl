## Backbone to structure your codebuild_model2!
# author: Kenneth Bruninx, Sebastian Gonzato
# last update: October 26, 2020
# description: backbone to structure your code. You're not obliged to use this
# in your assignment, but you may.

## Step 0: Activate environment - ensure consistency accross computers
using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate() # If a Manifest.toml file exist in the current project, download all the packages declared in that manifest. Else, resolve a set of feasible packages from the Project.toml files and install them.
Pkg.add("LaTeXStrings")
Pkg.add("PyPlot")
##  Step 1: input data
using CSV
using DataFrames
using YAML
using Logging
using MAT

data = YAML.load_file(joinpath(@__DIR__, "data.yaml"))
# ePrice = CSV.read(joinpath(@__DIR__, "eprice.csv"), DataFrame, silencewarnings=true)
# ePrice = CSV.read(joinpath(@__DIR__, "Day-ahead Prices_April_2023.csv"), DataFrame, silencewarnings=true)
# ePrice = CSV.read(joinpath(@__DIR__, "Day-ahead Prices_2023.csv"), DataFrame, silencewarnings=true)
ePrice = CSV.read(joinpath(@__DIR__, "Day-ahead Prices_2019.csv"), DataFrame, silencewarnings=true)


coeffs2SOEC1 = CSV.read(joinpath(@__DIR__, "coeffs2SOEC1.csv"), DataFrame, silencewarnings=true)
coeffs2SOEC2 = CSV.read(joinpath(@__DIR__, "coeffs2SOEC2.csv"), DataFrame, silencewarnings=true)
coeffs2SOEC3 = CSV.read(joinpath(@__DIR__, "coeffs2SOEC3.csv"), DataFrame, silencewarnings=true)
coeffs2SOEC4 = CSV.read(joinpath(@__DIR__, "coeffs2SOEC4.csv"), DataFrame, silencewarnings=true)
coeffs2SOEC5 = CSV.read(joinpath(@__DIR__, "coeffs2SOEC5.csv"), DataFrame, silencewarnings=true)
coeffs2SOEC6 = CSV.read(joinpath(@__DIR__, "coeffs2SOEC6.csv"), DataFrame, silencewarnings=true)

## Step 2: create model & pass data to model
using JuMP
using Gurobi
m = Model(optimizer_with_attributes(Gurobi.Optimizer))

# Step 2a: create sets
function define_sets!(m::Model, data::Dict, i::Int64)
    # create dictionary to store sets
    m.ext[:sets] = Dict()

    # define the sets
    m.ext[:sets][:J] = 1:i # Timesteps 
    m.ext[:sets][:parameterset] = [i for i in keys(data)]

    # return model
    return m
end

# Step 2b: add time series
function process_time_series_data_week!(m::Model,data::Dict, i::Int64)
    # extract the relevant sets
    #J = m.ext[:sets][:J] # Time steps

    # create dictionary to store time series
    m.ext[:timeseries] = Dict()

    # example: add time series to dictionary
    m.ext[:timeseries][:π_e] = ePrice."Day-ahead Price [EUR/MWh]"[(i-1)*7*24+1:i*7*24]/1000000 #[€/Wh]


    # Both 0, very large storage tank (no problem bcs in case 1 in paper storage never full so no determening constraint)
    # p_l = 0 so just more p_u (no power to local load but everything to utility grid)
    m.ext[:timeseries][:m_LH] = 0
    m.ext[:timeseries][:p_l] = 0 
    # return model
    return m
end
# Step 2b: add time series
function process_time_series_data_day!(m::Model, data::Dict, i::Int64)
    # extract the relevant sets
    #J = m.ext[:sets][:J] # Time steps

    # create dictionary to store time series
    m.ext[:timeseries] = Dict()

    # example: add time series to dictionary
    m.ext[:timeseries][:π_e] = ePrice."Day-ahead Price [EUR/MWh]"[(i-1)*24+1:i*24]/1000000 #[€/Wh]


    # Both 0, very large storage tank (no problem bcs in case 1 in paper storage never full so no determening constraint)
    # p_l = 0 so just more p_u (no power to local load but everything to utility grid)
    m.ext[:timeseries][:m_LH] = 0
    m.ext[:timeseries][:p_l] = 0 
    # return model
    return m
end
# step 2c: process input parameters
function process_parameters!(m::Model, data::Dict)
    # extract the sets you need
    parameterset = m.ext[:sets][:parameterset]
    # generate a dictonary "parameters"
    m.ext[:parameters] = Dict()
    # create parameters
    for i in parameterset
        m.ext[:parameters][i] = data[i]
    end
    m.ext[:parameters]["A"] = coeffs2SOEC1."Area [m^2]"[1]
    m.ext[:parameters]["M_HTMax"] = m.ext[:parameters]["SOCMax"] * m.ext[:parameters]["M_tank"]
    m.ext[:parameters]["M_HTMin"] = m.ext[:parameters]["SOCMin"] * m.ext[:parameters]["M_tank"]
    m.ext[:parameters]["alfa"] = m.ext[:parameters]["R"]*m.ext[:parameters]["T_sSOEC"]/(2*(m.ext[:parameters]["gamma"]-1)*m.ext[:parameters]["η_c"])*((m.ext[:parameters]["P_out"]/m.ext[:parameters]["P_in"])^((m.ext[:parameters]["gamma"]-1)/m.ext[:parameters]["gamma"])-1)
    m.ext[:parameters]["M_ini"] = m.ext[:parameters]["SOC_ini"] * m.ext[:parameters]["M_tank"]
    # m.ext[:parameters]["p_s"] = (m.ext[:parameters]["T_sSOEC"]-m.ext[:parameters]["T_a"])/(m.ext[:parameters]["R_tSOEC"]*m.ext[:parameters]["n_c"])
    # return model
    return m
end






#########################################################
#########################################################
#########################################################
##### ALFA NOG AANPASSEN ZOALS IN THESIS ################
#########################################################
#########################################################
#########################################################


















## Step 3: construct your model


function build_model1a!(m::Model)
    # Model with power in 1 segment, heating up and evaporating water by electrical power, heating up steam 40°C with electrical power
    # With old objective

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    C_CS = m.ext[:parameters]["C_CSSOEC"]  
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2SOEC1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2SOEC1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2SOEC1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    p_u = m.ext[:variables][:p_u] =  @variable(m, [j=J], lower_bound=0, base_name="Power to utility grid") 
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t + p_u[j]*π_e[j]*delta_t for j in J) - sum(C_HS*Y_b[j-1] for j in J[2:end])
    )

    # Formulate constraints
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
        s_b[j]*p_s[j] <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N + s_b[j]*p_s[j]
    )    
    m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
         p_s[j] == (T[j] - T_a)/(R_t*n_c)
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] 
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j-1]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= p_b[j]
    )

    # Constraint (46) is about ramping rates and can be relieved

    m.ext[:constraints][:con47] = @constraint(m, [j=J],
    p_N - n_c*p[j]*(1+s_b[j]*(-1+1/η_EH)) - p_c[j] - p_u[j] -(Q_H2O[j]+Q_steam[j])/η_EH - Q_cool[j]/400 == 0
    )
    # m.ext[:constraints][:con47] = @constraint(m, [j=J],
    # p_N - n_c*p[j] -n_c*p_Q[j]*η_turb- p_l - p_c[j] - p_u[j] == 0
    # )
    m.ext[:constraints][:con48a] = @constraint(m, [j=J],
    p_uMin <= p_u[j]
    )
    m.ext[:constraints][:con48b] = @constraint(m, [j=J],
    p_u[j] <= p_uMax
    )    
    m.ext[:constraints][:con52] = @constraint(m, [j=J],
    p_c[j] == alfa*mdotH2[j]
    )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    )
    
    return m
end
function build_model2a!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
    build_model1a!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
 

    return m
end
function build_model3a!(m::Model)
    # Model with power in 9 segment, heating up and evaporating water by electrical power
    build_model1a!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:9,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Delete constraints
    # Constraint 34 was implemented in parameter earlier so doesn't have to be deleted
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:9)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*3+1:(k-1)*3+3) for k in 2:3) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))
    )

    return m
end
function build_model4a!(m::Model)
    # Model with power in 16 segment, heating up and evaporating water by electrical power
    build_model1a!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:16,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Delete constraints
    # Constraint 34 was implemented in parameter earlier so doesn't have to be deleted
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))
    )

    return m
end
function build_model5a!(m::Model)
    # Model with power in 4 segments, with heat addition from LWR for:
        # 1) heating up water from T_a to 100°C and evaporate it
        #( 2) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    build_model2a!(m)

    # Extract timeseries

    # Extract parameters
    p_N = m.ext[:parameters]["p_N"]
    η_turb = m.ext[:parameters]["η_turb"]
    η_EH = m.ext[:parameters]["η_EH"]
    n_c = m.ext[:parameters]["n_cSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    s_b = m.ext[:variables][:s_b]
    p_c = m.ext[:variables][:p_c]
    p_u = m.ext[:variables][:p_u]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    

    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con47][j])
    end

    # New constraints
    m.ext[:constraints][:con47] = @constraint(m, [j=J],
    p_N - n_c*p[j]*(1+s_b[j]*(-1+1/η_EH)) - p_c[j] - p_u[j] - Q_H2O[j]*η_turb - Q_steam[j]/η_EH - Q_cool[j]/400 == 0
    )

    return m
end
function build_model6a!(m::Model)
    # Model with power in 4 segments, with heat addition from HTGR for:
        # 1) heating up water from T_a to 100°C and evaporate it with heat 
        # 2) heating up steam last 40°C with heat
        # 3) keeping temperature constant in stand-by with high temperature heat
        # 4) adding heat directly to electrolyser
        #( 5) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    build_model2a!(m)

    # Extract timeseries

    # Extract parameters
    p_N = m.ext[:parameters]["p_N"]
    η_turb = m.ext[:parameters]["η_turb"]
    n_c = m.ext[:parameters]["n_cSOEC"]
    delta_t = m.ext[:parameters]["delta_t"]
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"]
    T_a = m.ext[:parameters]["T_a"]
    R_t = m.ext[:parameters]["R_tSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    s_b = m.ext[:variables][:s_b]
    p_c = m.ext[:variables][:p_c]
    p_u = m.ext[:variables][:p_u]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    
    
    # Create variables
    Q = m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser")

    # Delete constraints
    for j in J[2:end]
        delete(m,m.ext[:constraints][:con30][j])
    end
    for j in J
        delete(m,m.ext[:constraints][:con47][j])
    end
    
    # New constraints
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con47] = @constraint(m, [j=J],
        p_N - n_c*p[j]*(1+s_b[j]*(-1+η_turb))  - p_c[j] - p_u[j] - (Q_H2O[j]+ Q_steam[j]+Q[j])*η_turb - Q_cool[j]/400 == 0
    )
    
    # Voorlopig geen constraints op Q!
    # m.ext[:constraints][:con70] = @constraint(m, [j=J],
    # Q[j] <= n_c*p[j-1] - n_c*U_tn*I[j-1]
    # )
    
    return m
end
function build_model7a!(m::Model)
    # Model with power in 9 segments, with heat addition from HTGR for:
        # 1) heating up water from T_a to 100°C and evaporate it with heat 
        # 2) heating up steam last 40°C with heat
        # 3) keeping temperature constant in stand-by with high temperature heat
        # 4) adding heat directly to electrolyser
        #( 5) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    build_model3a!(m)

    # Extract timeseries

    # Extract parameters
    p_N = m.ext[:parameters]["p_N"]
    η_turb = m.ext[:parameters]["η_turb"]
    n_c = m.ext[:parameters]["n_cSOEC"]
    delta_t = m.ext[:parameters]["delta_t"]
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"]
    T_a = m.ext[:parameters]["T_a"]
    R_t = m.ext[:parameters]["R_tSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    s_b = m.ext[:variables][:s_b]
    p_c = m.ext[:variables][:p_c]
    p_u = m.ext[:variables][:p_u]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    
    
    # Create variables
    Q = m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser")

    # Delete constraints
    for j in J[2:end]
        delete(m,m.ext[:constraints][:con30][j])
    end
    for j in J
        delete(m,m.ext[:constraints][:con47][j])
    end
    
    # New constraints
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con47] = @constraint(m, [j=J],
    p_N - n_c*p[j]*(1+s_b[j]*(-1+η_turb)) - p_c[j] - p_u[j] - (Q_H2O[j]+ Q_steam[j]+Q[j])*η_turb - Q_cool[j]/400 == 0
    )
    
    # Voorlopig geen constraints op Q!
    # m.ext[:constraints][:con70] = @constraint(m, [j=J],
    # Q[j] <= n_c*p[j-1] - n_c*U_tn*I[j-1]
    # )
    
    return m
end

function build_model1b!(m::Model)
    # Model with power in 1 segment, heating up and evaporating water by electrical power, heating up steam 40°C with electrical power
    # New objective
    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2SOEC1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2SOEC1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2SOEC1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    p_u = m.ext[:variables][:p_u] =  @variable(m, [j=J], lower_bound=0, base_name="Power to utility grid") 
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
     


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J) - sum(C_HS*Y_b[j-1] for j in J[2:end])
    )

    # Formulate constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+1/η_EH)) + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
        s_b[j]*p_s[j] <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N + s_b[j]*p_s[j]
    )    
    m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
         p_s[j] == (T[j] - T_a)/(R_t*n_c)
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] 
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j-1]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= p_b[j]
    )

    # Constraint (46) is about ramping rates and can be relieved

    
    m.ext[:constraints][:con48a] = @constraint(m, [j=J],
    p_uMin <= p_u[j]
    )
    m.ext[:constraints][:con48b] = @constraint(m, [j=J],
    p_u[j] <= p_uMax
    )    
    m.ext[:constraints][:con52] = @constraint(m, [j=J],
    p_c[j] == alfa*mdotH2[j]
    )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    )
    
    return m
end
function build_model2b!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
    build_model1b!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
 

    return m
end
function build_model3b!(m::Model)
    # Model with power in 9 segment, heating up and evaporating water by electrical power
    build_model1b!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:9,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Delete constraints
    # Constraint 34 was implemented in parameter earlier so doesn't have to be deleted
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:9)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*3+1:(k-1)*3+3) for k in 2:3) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))
    )

    return m
end
function build_model4b!(m::Model)
    # Model with power in 16 segment, heating up and evaporating water by electrical power
    build_model1b!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:16,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Delete constraints
    # Constraint 34 was implemented in parameter earlier so doesn't have to be deleted
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))
    )

    return m
end
function build_model5b!(m::Model)
    # Model with power in 4 segments, with heat addition from LWR for:
        # 1) heating up water from T_a to 100°C and evaporate it
        #( 2) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    build_model2b!(m)

    # Extract timeseries

    # Extract parameters
    η_turb = m.ext[:parameters]["η_turb"]
    η_EH = m.ext[:parameters]["η_EH"]
    n_c = m.ext[:parameters]["n_cSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    s_b = m.ext[:variables][:s_b]
    p_c = m.ext[:variables][:p_c]
    p_E = m.ext[:variables][:p_E]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    

    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end

    # New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J], 
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+1/η_EH)) + p_c[j] +Q_H2O[j]*η_turb+Q_steam[j]/η_EH + Q_cool[j]/400 
    )

    return m
end
function build_model6b!(m::Model)
    # Model with power in 4 segments, with heat addition from HTGR for:
        # 1) heating up water from T_a to 100°C and evaporate it with heat 
        # 2) heating up steam last 40°C with heat
        # 3) keeping temperature constant in stand-by with high temperature heat
        # 4) adding heat directly to electrolyser
        #( 5) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    
    build_model2b!(m)

    # Extract timeseries

    # Extract parameters
    η_turb = m.ext[:parameters]["η_turb"]
    n_c = m.ext[:parameters]["n_cSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"]
    U_tn= m.ext[:parameters]["U_tnSOEC"]
    R_t= m.ext[:parameters]["R_tSOEC"]
    delta_t= m.ext[:parameters]["delta_t"]
    T_a= m.ext[:parameters]["T_a"]
    # Extract variables
    p = m.ext[:variables][:p]
    s_b = m.ext[:variables][:s_b]
    p_c = m.ext[:variables][:p_c]
    p_E = m.ext[:variables][:p_E]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    T = m.ext[:variables][:T]
    I = m.ext[:variables][:I]

    # Create variables
    Q = m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser")


    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end
    for j in J[2:end]
        delete(m,m.ext[:constraints][:con30][j])
    end
    # New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J], 
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+η_turb)) + p_c[j] +(Q_H2O[j]+Q_steam[j]+Q[j])*η_turb + Q_cool[j]/400 
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )


    return m
end
function build_model7b!(m::Model)
    # Model with power in 9 segments, with heat addition from HTGR for:
        # 1) heating up water from T_a to 100°C and evaporate it with heat 
        # 2) heating up steam last 40°C with heat
        # 3) keeping temperature constant in stand-by with high temperature heat
        # 4) adding heat directly to electrolyser
        #( 5) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    build_model3b!(m)

    # Extract timeseries

    # Extract parameters
    η_turb = m.ext[:parameters]["η_turb"]
    n_c = m.ext[:parameters]["n_cSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"]
    U_tn= m.ext[:parameters]["U_tnSOEC"]
    R_t= m.ext[:parameters]["R_tSOEC"]
    delta_t= m.ext[:parameters]["delta_t"]
    T_a= m.ext[:parameters]["T_a"]

    # Extract variables
    p = m.ext[:variables][:p]
    s_b = m.ext[:variables][:s_b]
    p_c = m.ext[:variables][:p_c]
    p_E = m.ext[:variables][:p_E]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    T = m.ext[:variables][:T]
    I = m.ext[:variables][:I]

    # Create variables
    Q = m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser")


    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end
    for j in J[2:end]
        delete(m,m.ext[:constraints][:con30][j])
    end
    # New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J], 
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+η_turb)) + p_c[j] +(Q_H2O[j]+Q_steam[j]+Q[j])*η_turb + Q_cool[j]/400 
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    
    return m
end

function build_model1c!(m::Model)
    # Model with power in 1 segment, heating up and evaporating water by electrical power, heating up steam 40°C with electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # EXTRA cyclic boundary constraint


    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2SOEC1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2SOEC1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2SOEC1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    p_u = m.ext[:variables][:p_u] =  @variable(m, [j=J], lower_bound=0, base_name="Power to utility grid") 
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
     


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+1/η_EH)) + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
        s_b[j]*p_s[j] <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N + s_b[j]*p_s[j]
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c) <= p_s[j]
    )
    # m.ext[:constraints][:con34b] = @constraint(m, [j=J],
    #       p_s[j] <= M*s_b[j]
    # )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] 
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j-1]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= p_b[j]
    )

    # Constraint (46) is about ramping rates and can be relieved

    
    m.ext[:constraints][:con48a] = @constraint(m, [j=J],
        p_uMin <= p_u[j]
    )
    m.ext[:constraints][:con48b] = @constraint(m, [j=J],
        p_u[j] <= p_uMax
    )    
    m.ext[:constraints][:con52] = @constraint(m, [j=J],
        p_c[j] == alfa*mdotH2[j]
    )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
        Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
        Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    return m
end
function build_model2c!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
        # EXTRA heat addition via electrical heaters

    build_model1c!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
 

    return m
end
function build_model3c!(m::Model)
    # Model with power in 9 segment, heating up and evaporating water by electrical power
        # EXTRA heat addition via electrical heaters

    build_model1c!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:9,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Delete constraints
    # Constraint 34 was implemented in parameter earlier so doesn't have to be deleted
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:9)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*3+1:(k-1)*3+3) for k in 2:3) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))
    )

    return m
end
function build_model4c!(m::Model)
    # Model with power in 16 segment, heating up and evaporating water by electrical power
        # EXTRA heat addition via electrical heaters

    build_model1c!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:16,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Delete constraints
    # Constraint 34 was implemented in parameter earlier so doesn't have to be deleted
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))
    )

    return m
end
function build_model5c!(m::Model)
    # Model with power in 4 segments, with heat addition from LWR for:
        # 1) heating up water from T_a to 100°C and evaporate it
        #( 2) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 
    
        # EXTRA heat addition via electrical heaters

    build_model2c!(m)

    # Extract timeseries

    # Extract parameters
    η_turb = m.ext[:parameters]["η_turb"]
    η_EH = m.ext[:parameters]["η_EH"]
    n_c = m.ext[:parameters]["n_cSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    s_b = m.ext[:variables][:s_b]
    p_c = m.ext[:variables][:p_c]
    p_E = m.ext[:variables][:p_E]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    Q = m.ext[:variables][:Q]

    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end

    # New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J], 
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+1/η_EH)) + p_c[j] +Q_H2O[j]*η_turb+Q_steam[j]/η_EH + Q_cool[j]/400 + Q[j]/η_EH 
    )

    return m
end
function build_model6c!(m::Model)
    # Model with power in 4 segments, with heat addition from HTGR for:
        # 1) heating up water from T_a to 100°C and evaporate it with heat 
        # 2) heating up steam last 40°C with heat
        # 3) keeping temperature constant in stand-by with high temperature heat
        # 4) adding heat directly to electrolyser
        #( 5) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    #REMARK: same as model 6b
    build_model2c!(m)

    # Extract timeseries

    # Extract parameters
    η_turb = m.ext[:parameters]["η_turb"]
    n_c = m.ext[:parameters]["n_cSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    p_c = m.ext[:variables][:p_c]
    p_E = m.ext[:variables][:p_E]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    Q = m.ext[:variables][:Q] 
    s_b = m.ext[:variables][:s_b] 


    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end
    # New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J], 
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+η_turb)) + p_c[j] +(Q_H2O[j]+Q_steam[j]+Q[j])*η_turb + Q_cool[j]/400 
    )

    return m
end
function build_model7c!(m::Model)
    # Model with power in 9 segments, with heat addition from HTGR for:
        # 1) heating up water from T_a to 100°C and evaporate it with heat 
        # 2) heating up steam last 40°C with heat
        # 3) keeping temperature constant in stand-by with high temperature heat
        # 4) adding heat directly to electrolyser
        #( 5) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 
    
        #REMARK: same as model 7b
    build_model3c!(m)

    # Extract timeseries
    # Extract parameters
    η_turb = m.ext[:parameters]["η_turb"]
    n_c = m.ext[:parameters]["n_cSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    p_c = m.ext[:variables][:p_c]
    p_E = m.ext[:variables][:p_E]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    Q = m.ext[:variables][:Q] 
    s_b = m.ext[:variables][:s_b] 



    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end
    # New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J], 
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+η_turb)) + p_c[j] +(Q_H2O[j]+Q_steam[j]+Q[j])*η_turb + Q_cool[j]/400 
    )
    
    return m
end



function build_model1e!(m::Model)
    # Model with power in 1 segment, heating up and evaporating water by electrical power, heating up steam 40°C with electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2SOEC1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2SOEC1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2SOEC1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )



    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] 
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )

    # Constraint (46) is about ramping rates and can be relieved

    m.ext[:constraints][:con52] = @constraint(m, [j=J],
        p_c[j] == alfa*mdotH2[j]
    )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
        Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
        Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    return m
end
function build_model2e!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
        # EXTRA heat addition via electrical heaters

    build_model1e!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
 

    return m
end
function build_model3e!(m::Model)
    # Model with power in 9 segment, heating up and evaporating water by electrical power
        # EXTRA heat addition via electrical heaters

    build_model1e!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:9,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Delete constraints
    # Constraint 34 was implemented in parameter earlier so doesn't have to be deleted
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:9)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*3+1:(k-1)*3+3) for k in 2:3) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))
    )

    return m
end
function build_model4e!(m::Model)
    # Model with power in 16 segment, heating up and evaporating water by electrical power
        # EXTRA heat addition via electrical heaters

    build_model1e!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:16,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC4."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC4."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC4."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC4."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC4."Current density starts"

    # Delete constraints
    # Constraint 34 was implemented in parameter earlier so doesn't have to be deleted
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))
    )

    return m
end
function build_model5e!(m::Model)
    # Model with power in 4 segments, with heat addition from LWR for:
        # 1) heating up water from T_a to 100°C and evaporate it
        #( 2) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 
    
        # EXTRA heat addition via electrical heaters

    build_model2e!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract parameters
    η_turb = m.ext[:parameters]["η_turb"]
    η_EH = m.ext[:parameters]["η_EH"]
    n_c = m.ext[:parameters]["n_cSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    p_s = m.ext[:variables][:p_s]
    p_c = m.ext[:variables][:p_c]
    p_E = m.ext[:variables][:p_E]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    Q = m.ext[:variables][:Q]

    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end

    # New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J], 
        p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +Q_H2O[j]*η_turb+Q_steam[j]/η_EH + Q_cool[j]/400 + Q[j]/η_EH 
    )
    return m
end
function build_model6e!(m::Model)
    # Model with power in 4 segments, with heat addition from HTGR for:
        # 1) heating up water from T_a to 100°C and evaporate it with heat 
        # 2) heating up steam last 40°C with heat
        # 3) keeping temperature constant in stand-by with high temperature heat
        # 4) adding heat directly to electrolyser
        #( 5) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    
    build_model2e!(m)

    # Extract sets
    J = m.ext[:sets][:J]
    # Extract timeseries

    # Extract parameters
    η_turb = m.ext[:parameters]["η_turb"]
    n_c = m.ext[:parameters]["n_cSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    p_s = m.ext[:variables][:p_s]
    p_c = m.ext[:variables][:p_c]
    p_E = m.ext[:variables][:p_E]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    Q = m.ext[:variables][:Q] 


    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end
    # New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J], 
        p_E[j] == n_c*p[j]+n_c*p_s[j]*η_turb + p_c[j] +(Q_H2O[j]+Q_steam[j]+Q[j])*η_turb + Q_cool[j]/400 
    )


    return m
end
function build_model7e!(m::Model)
    # Model with power in 9 segments, with heat addition from HTGR for:
        # 1) heating up water from T_a to 100°C and evaporate it with heat 
        # 2) heating up steam last 40°C with heat
        # 3) keeping temperature constant in stand-by with high temperature heat
        # 4) adding heat directly to electrolyser
        #( 5) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 
    
        #REMARK: same as model 7b
    build_model3e!(m)

     # Extract sets
     J = m.ext[:sets][:J]
    # Extract timeseries
    # Extract parameters
    η_turb = m.ext[:parameters]["η_turb"]
    n_c = m.ext[:parameters]["n_cSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    p_s = m.ext[:variables][:p_s]
    p_c = m.ext[:variables][:p_c]
    p_E = m.ext[:variables][:p_E]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    Q = m.ext[:variables][:Q] 


    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end
    # New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J], 
        p_E[j] == n_c*p[j]+n_c*p_s[j]*η_turb + p_c[j] +(Q_H2O[j]+Q_steam[j]+Q[j])*η_turb + Q_cool[j]/400 
    )
    
    return m
end


function build_model1e2!(m::Model)
    # Model with power in 1 segment, heating up and evaporating water by electrical power, heating up steam 40°C with electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2SOEC1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2SOEC1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2SOEC1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    # p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p = m.ext[:expressions][:p] = @expression(m, [j=J],
        a*T[j]+ b*I[j]/A + c + delta_1[j] 
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    # m.ext[:constraints][:con37] = @constraint(m, [j=J],
    #     p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] 
    # )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )

    # Constraint (46) is about ramping rates and can be relieved

   
    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    return m
end
function build_model2e2!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 

    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    # p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p = m.ext[:expressions][:p] = @expression(m, [j=J],
        sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    #New constraints
    # m.ext[:constraints][:con37] = @constraint(m, [j=J],
    #     p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    # )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )

    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    return m
end
function build_model3e2!(m::Model)
    # Model with power in 9 segment, heating up and evaporating water by electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 

    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:9,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p = m.ext[:expressions][:p] = @expression(m, [j=J],
        sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    #New constraints
    # m.ext[:constraints][:con37] = @constraint(m, [j=J],
    #     p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    # )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:9)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*3+1:(k-1)*3+3) for k in 2:3) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )

    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    return m
end
function build_model4e2!(m::Model)
    # Model with power in 16 segment, heating up and evaporating water by electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC4."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC4."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC4."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC4."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC4."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    p_u = m.ext[:variables][:p_u] =  @variable(m, [j=J], lower_bound=0, base_name="Power to utility grid") 
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:16,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
     


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )

    # Constraint (46) is about ramping rates and can be relieved

    
    m.ext[:constraints][:con48a] = @constraint(m, [j=J],
        p_uMin <= p_u[j]
    )
    m.ext[:constraints][:con48b] = @constraint(m, [j=J],
        p_u[j] <= p_uMax
    )    
    m.ext[:constraints][:con52] = @constraint(m, [j=J],
        p_c[j] == alfa*mdotH2[j]
    )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
        Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
        Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    return m
end
function build_model5e2!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 

    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]
    η_turb = m.ext[:parameters]["η_turb"]

    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    # p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p = m.ext[:expressions][:p] = @expression(m, [j=J],
        sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +Q_H2O[j]*η_turb+Q_steam[j]/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    #New constraints
    # m.ext[:constraints][:con37] = @constraint(m, [j=J],
    #     p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    # )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )

    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    return m
end
function build_model6e2!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 

    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]
    η_turb = m.ext[:parameters]["η_turb"]

    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    # p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p = m.ext[:expressions][:p] = @expression(m, [j=J],
        sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]*η_turb + p_c[j] +Q_H2O[j]*η_turb+Q_steam[j]*η_turb + Q_cool[j]/400 + Q[j]*η_turb
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    #New constraints
    # m.ext[:constraints][:con37] = @constraint(m, [j=J],
    #     p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    # )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )

    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    return m
end


### FINALE MODELLEN ###
function build_model1g!(m::Model)
    # Model with power in 1 segment, heating up and evaporating water by electrical power, heating up steam 40°C with electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2SOEC1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2SOEC1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2SOEC1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    # Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     a*T[j]+ b*I[j]/A + c + delta_1[j] 
    # )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] + Q[j])*(3600*delta_t)/ C_h
    )


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] 
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    # m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
    #     s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    # )
    # m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= s_b[j-1]
    # )
    # m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= p_b[j]
    # )

    # Constraint (46) is about ramping rates and can be relieved

   
    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    # m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    return m
end
function build_model2g!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # EXTRA same as model2e2 but T[end] not constraint to T_0!!!

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 

    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state")

    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j]
    # )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )    
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] + Q[j])*(3600*delta_t)/ C_h
    )

    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    # m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
    #     s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    # )
    # m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= s_b[j-1]
    # )
    # m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= p_b[j]
    # )

    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    # m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    return m
end
function build_model3g!(m::Model)
    # Model with power in 9 segment, heating up and evaporating water by electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 

    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    # Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:9,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j]
    # )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] + Q[j])*(3600*delta_t)/ C_h
    )

    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:9)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*3+1:(k-1)*3+3) for k in 2:3) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    # m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
    #     s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    # )
    # m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= s_b[j-1]
    # )
    # m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= p_b[j]
    # )

    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    # m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    return m
end
function build_model4g!(m::Model)
    # Model with power in 16 segment, heating up and evaporating water by electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC4."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC4."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC4."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC4."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC4."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    p_u = m.ext[:variables][:p_u] =  @variable(m, [j=J], lower_bound=0, base_name="Power to utility grid") 
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    # Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:16,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j]
    # )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] + Q[j])*(3600*delta_t)/ C_h
    )


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    # m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
    #     s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    # )
    # m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= s_b[j-1]
    # )
    # m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= p_b[j]
    # )

    # Constraint (46) is about ramping rates and can be relieved

    
  
    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
        Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
        Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    # m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    return m
end
function build_model5g!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 

    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]
    η_turb = m.ext[:parameters]["η_turb"]

    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    # Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j]
    # )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +Q_H2O[j]*η_turb+Q_steam[j]/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] + Q[j])*(3600*delta_t)/ C_h
    )

    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    # m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
    #     s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    # )
    # m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= s_b[j-1]
    # )
    # m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= p_b[j]
    # )

    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    # m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    return m
end
function build_model6g!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
    # New objective
    # No off state
    # EXTRA heat addition via electrical heaters
    # EXTRA cost hot start = 0
    # Cyclic boundary conditions

    # EXTRA andere formulatie van p en p_s die veel sneller is

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 

    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]
    η_turb = m.ext[:parameters]["η_turb"]

    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    # Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j]
    # )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]*η_turb + p_c[j] +Q_H2O[j]*η_turb+Q_steam[j]*η_turb + Q_cool[j]/400 + Q[j]*η_turb
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] + Q[j])*(3600*delta_t)/ C_h
    )

    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    # m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
    #     s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    # )
    # m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= s_b[j-1]
    # )
    # m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= p_b[j]
    # )

    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    # m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    return m
end

function build_model8g!(m::Model)
    # 1 segment maar met heat addition van HTGR
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 

    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    p_N = m.ext[:parameters]["p_N"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]
    η_turb = m.ext[:parameters]["η_turb"]

    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2SOEC1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2SOEC1."Coefficients j [W/(A/m^2)]"[1]

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    # Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    # Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")
    # p_E = m.ext[:variables][:p_E] = @variable(m,[j=J],lower_bound = 0, base_name="Electrical power")
    Q= m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser via electrical heaters")
    # t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j]
    # )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j]
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    Q_steam = m.ext[:expressions][:Q_steam] = @expression(m, [j=J], 
        mdotH2O[j]*cp_steam*delta_T
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+ n_c*p_s[j]*η_turb + p_c[j] +Q_H2O[j]*η_turb+Q_steam[j]*η_turb + Q_cool[j]/400 + Q[j]*η_turb
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] + Q[j])*(3600*delta_t)/ C_h
    )

    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J)
    )


     
    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+ n_c*p_s[j]/η_EH + p_c[j] +(Q_H2O[j]+Q_steam[j])/η_EH + Q_cool[j]/400 + Q[j]/η_EH
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] + n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
       0 <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N 
    )    
    # m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
    #      p_s[j] == (T[j] - T_a)/(R_t*n_c)*s_b[j]
    # )
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c)*s_b[j] <= p_s[j]
    )
    m.ext[:constraints][:con34b] = @constraint(m, [j=J],
          p_s[j] <= M*s_b[j]
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j]  == 1
    )
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == a*T[j] + b*I[j]/A + c + delta_1[j] 
    )
    
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    
    # Constraint (40) and (41) are replaced by (42) - (45)

    # m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
    #     s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    # )
    # m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= s_b[j-1]
    # )
    # m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
    #     Y_b[j] <= p_b[j]
    # )

    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    # )
    # #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    # m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    #     Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    # )
    m.ext[:constraints][:con54c] = @constraint(m, [j=J], 
        Q[j] <= p_b[j]*M 
    )
    # m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    return m
end

function build_model1test!(m::Model)
    # Model with power in 1 segment, heating up and evaporating water by electrical power

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price
    m_LH = m.ext[:timeseries][:m_LH] # = 0
    p_l = m.ext[:timeseries][:p_l] # = 0

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    # η_f = m.ext[:parameters]["η_f"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    C_CS = m.ext[:parameters]["C_CSSOEC"]  
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    M_ini = m.ext[:parameters]["M_ini"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    M_HTMin = m.ext[:parameters]["M_HTMin"]
    M_HTMax = m.ext[:parameters]["M_HTMax"]
    p_N = m.ext[:parameters]["p_N"]
    η_turb = m.ext[:parameters]["η_turb"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]

    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2SOEC1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2SOEC1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2SOEC1."Coefficients j [W/(A/m^2)]"[1]
    # Tborders = m.ext[:parameters]["Tborders"] = coeffs2SOEC1."Temperature borders"
    # jborders = m.ext[:parameters]["jborders"] = coeffs2SOEC1."Current density borders"
    # num_segments = m.ext[:parameters]["num_segments"] = 1

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    M_HT = m.ext[:variables][:M_HT] =  @variable(m, [j=J], lower_bound=0, base_name="Mass of stored hydrogen") #kg
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    p_u = m.ext[:variables][:p_u] =  @variable(m, [j=J], lower_bound=0, base_name="Power to utility grid") 
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # p_Q = m.ext[:variables][:p_Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to electrolyser")
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    

    # Create affine expressions (= linear combinations of variables)
    # T = m.ext[:expressions][:T] = @expression(m, [j=J],
    # if j == 1
    #     T_0
    # else
    #     sum(3600*delta_t/C_h*(n_c*p[j-i]-n_c*U_tn*I[j-i]+T_a/R_t-Q_cool[j-i])*(1-3600*delta_t/(R_t*C_h))^(i-1) for i in 1:j-1) +T_0* (1-3600*delta_t/(R_t*C_h))^(j-1)
    # end
    # )
    # with heat
    # T = m.ext[:expressions][:T] = @expression(m, [j=J[2:end]],
    #     sum(3600*delta_t/C_h*(n_c*p[j-i]-n_c*U_tn*I[j-i]+p_Q[j]+T_a/R_t-Q_cool[j-i])*(1-3600*delta_t/(R_t*C_h))^(i-1) for i in 1:j-1)+T_0*(1-3600*delta_t/(R_t*C_h))^(j-1)
    # )
    # Without heat
    T = m.ext[:expressions][:T] = @expression(m, [j=J[2:end]],
        sum(3600*delta_t/C_h*(n_c*p[j-i]-n_c*U_tn*I[j-i]+T_a/R_t-Q_cool[j-i])*(1-3600*delta_t/(R_t*C_h))^(i-1) for i in 1:j-1)+T_0*(1-3600*delta_t/(R_t*C_h))^(j-1)
    )
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    # mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
    #     M_H2*n_c*η_f*I[j]/(2*F) 
    # )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t + (p_u[j] + p_l)*π_e[j]*delta_t for j in J) - sum(C_HS*Y_b[j-1]+C_CS*Z_b[j-1] for j in J[2:end])
    )

    # Formulate constraints
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        (273+800) <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
        s_b[j]*p_s[j]+p_b[j] <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N + s_b[j]*p_s[j]
    )    
    m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
         p_s[j] == (T[j] - T_a)/(R_t*n_c)/η_EH
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    # m.ext[:constraints][:con37] = @constraint(m, [j=J],
    #     p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] + delta_2[j]
    # )

    m.ext[:constraints][:con37] = @constraint(m, [j=J[2:end]],
        p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] + delta_2[j]
    )
    m.ext[:constraints][:con37b] = @constraint(m, [j=1],
        p[j] == a*T_0+ b*I[j]/A + c  + delta_1[j] + delta_2[j]
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    m.ext[:constraints][:con38c] = @constraint(m, [j=J],
        -M*i_b[j] <= delta_2[j] 
    )
    m.ext[:constraints][:con38d] = @constraint(m, [j=J],
        delta_2[j] <= M*i_b[j]
    )
    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
    )
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j-1]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= p_b[j]
    )
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j-1]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j-1] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j-1] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved

    m.ext[:constraints][:con47] = @constraint(m, [j=J],
    p_N - n_c*p[j] - p_l - p_c[j] - p_u[j]  == 0
    )
    # m.ext[:constraints][:con47] = @constraint(m, [j=J],
    # p_N - n_c*p[j] -n_c*p_Q[j]*η_turb- p_l - p_c[j] - p_u[j] == 0
    # )
    m.ext[:constraints][:con48a] = @constraint(m, [j=J],
    p_uMin <= p_u[j]
    )
    m.ext[:constraints][:con48b] = @constraint(m, [j=J],
    p_u[j] <= p_uMax
    )
    m.ext[:constraints][:con49] = @constraint(m, [j=J[2:end]],
    M_HT[j] - M_HT[j-1] == mdotH2[j]*delta_t - m_LH*delta_t
    )
    m.ext[:constraints][:con50] = @constraint(m, [j=1],
    M_HT[j] == M_ini
    )
    m.ext[:constraints][:con51a] = @constraint(m, [j=J],
    M_HTMin <= M_HT[j]
    )
    m.ext[:constraints][:con51b] = @constraint(m, [j=J],
    M_HT[j] <= M_HTMax
    )
    m.ext[:constraints][:con52] = @constraint(m, [j=J],
    p_c[j] == alfa*mdotH2[j]
    )
    # m.ext[:constraints][:con53] = @constraint(m, [j=J],
    # p_Q[j] ==  s_b[j]*(T[j] - T_a)/(R_t*n_c)
    # )
    
    # #constraint hoeveelheid electrical power needed voor opwarmen water 20-100°C en verdampen 
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    # Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373.15-T_a) + delta_H/M_H2O*mdotH2O[j])/η_turb/η_EH
    # )
    
    return m
end
function build_model5test!(m::Model)
    # Model with power in 2 segment, with heat addition from LWR



    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price
    m_LH = m.ext[:timeseries][:m_LH] # = 0
    p_l = m.ext[:timeseries][:p_l] # = 0

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    # η_f = m.ext[:parameters]["η_f"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    C_CS = m.ext[:parameters]["C_CSSOEC"]  
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    M_ini = m.ext[:parameters]["M_ini"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    M_HTMin = m.ext[:parameters]["M_HTMin"]
    M_HTMax = m.ext[:parameters]["M_HTMax"]
    p_N = m.ext[:parameters]["p_N"]
    η_turb = m.ext[:parameters]["η_turb"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"
    
    
    # num_segments = m.ext[:parameters]["num_segments"] = 1

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    M_HT = m.ext[:variables][:M_HT] =  @variable(m, [j=J], lower_bound=0, base_name="Mass of stored hydrogen") #kg
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    p_u = m.ext[:variables][:p_u] =  @variable(m, [j=J], lower_bound=0, base_name="Power to utility grid") 
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    p_Q = m.ext[:variables][:p_Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to electrolyser")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    
    
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segment 4") #W


    # Create affine expressions (= linear combinations of variables)
    # T = m.ext[:expressions][:T] = @expression(m, [j=J[2:end]],
    #     sum(3600*delta_t/C_h*(n_c*p[j-i]-n_c*U_tn*I[j-i]+T_a/R_t-Q_cool[j-i])*(1-3600*delta_t/(R_t*C_h))^(i-1) for i in 1:j-1)+T_0*(1-3600*delta_t/(R_t*C_h))^(j-1)
    # )

    # T = m.ext[:expressions][:T] = @expression(m, [j=J],
    # if j == 1
    #     T_0
    # else
    #     sum(3600*delta_t/C_h*(n_c*p[j-i]-n_c*U_tn*I[j-i]+T_a/R_t-Q_cool[j-i])*(1-3600*delta_t/(R_t*C_h))^(i-1) for i in 1:j-1) +T_0* (1-3600*delta_t/(R_t*C_h))^(j-1)
    # end
    # )

    # With heat
    T = m.ext[:expressions][:T] = @expression(m, [j=J[2:end]],
        sum(3600*delta_t/C_h*(n_c*p[j-i]-n_c*U_tn*I[j-i]+p_Q[j]+T_a/R_t-Q_cool[j-i])*(1-3600*delta_t/(R_t*C_h))^(i-1) for i in 1:j-1)+T_0*(1-3600*delta_t/(R_t*C_h))^(j-1)
    )

    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*_I[j]/(2*F) 
    )
    # mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
    #     M_H2*n_c*1*I[j]/(2*F) 
    # )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    
    


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t + (p_u[j] + p_l)*π_e[j]*delta_t for j in J) - sum(C_HS*Y_b[j-1]+C_CS*Z_b[j-1] for j in J[2:end])
    )

    # Formulate constraints
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        (273+58)*s_b[j] <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
        s_b[j]*p_s[j]+p_b[j] <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N + s_b[j]*p_s[j]
    )    
    m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
         p_s[j] == (T[j] - T_a)/(R_t*n_c)
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j] + i_b[j] == 1
    )

    # m.ext[:constraints][:con37] = @constraint(m, [j=J],
    #     p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] + delta_2[j]
    # )
    #  m.ext[:constraints][:con37] = @constraint(m, [j=J[2:end]],
    #     p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] + delta_2[j]
    # )
    # m.ext[:constraints][:con37b] = @constraint(m, [j=1],
    #     p[j] == a*T_0+ b*I[j]/A + c + delta_1[j] + delta_2[j]
    # )

    m.ext[:constraints][:con37] = @constraint(m, [j=J[2:end]],
    p[j] == t_b[1,j]*(a[1]*T[j]+ b[1]*I[j]/A + c[1]) + t_b[2,j]*(a[2]*T[j]+ b[2]*I[j]/A + c[2])+ t_b[3,j]*(a[3]*T[j]+ b[3]*I[j]/A + c[3]) + t_b[4,j]*(a[4]*T[j]+ b[4]*I[j]/A + c[4])+delta_1[j] + delta_2[j]
    )
    m.ext[:constraints][:con37b] = @constraint(m, [j=1],
    p[j] == t_b[1,j]*(a[1]*T_0+ b[1]*I[j]/A + c[1]) + t_b[2,j]*(a[2]*T_0+ b[2]*I[j]/A + c[2])+ t_b[3,j]*(a[3]*T_0+ b[3]*I[j]/A + c[3]) + t_b[4,j]*(a[4]*T_0+ b[4]*I[j]/A + c[4])+delta_1[j] + delta_2[j]
    )

    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    m.ext[:constraints][:con38c] = @constraint(m, [j=J],
        -M*i_b[j] <= delta_2[j] 
    )
    m.ext[:constraints][:con38d] = @constraint(m, [j=J],
        delta_2[j] <= M*i_b[j]
    )
    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
    )
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j-1]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= p_b[j]
    )
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j-1]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j-1] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j-1] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved

    # m.ext[:constraints][:con47] = @constraint(m, [j=J],
    # p_N - n_c*p[j] - p_l - p_c[j] - p_u[j] == 0
    # )
    m.ext[:constraints][:con47] = @constraint(m, [j=J],
    p_N - n_c*p[j] -n_c*p_Q[j]*η_turb- p_l - p_c[j] - p_u[j] == 0
    )
    m.ext[:constraints][:con48a] = @constraint(m, [j=J],
    p_uMin <= p_u[j]
    )
    m.ext[:constraints][:con48b] = @constraint(m, [j=J],
    p_u[j] <= p_uMax
    )
    m.ext[:constraints][:con49] = @constraint(m, [j=J[2:end]],
    M_HT[j] - M_HT[j-1] == mdotH2[j]*delta_t - m_LH*delta_t
    )
    m.ext[:constraints][:con50] = @constraint(m, [j=1],
    M_HT[j] == M_ini
    )
    m.ext[:constraints][:con51a] = @constraint(m, [j=J],
    M_HTMin <= M_HT[j]
    )
    m.ext[:constraints][:con51b] = @constraint(m, [j=J],
    M_HT[j] <= M_HTMax
    )
    m.ext[:constraints][:con52] = @constraint(m, [j=J],
    p_c[j] == alfa*mdotH2[j]
    )
    # m.ext[:constraints][:con53] = @constraint(m, [j=J],
    # p_Q[j] ==  s_b[j]*(T[j] - T_a)/(R_t*n_c)
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    # Q_H2O[j] == mdotH2O[j]*cp_H2O*delta_T 
    # )
    # m.ext[:constraints][:con55] = @constraint(m, [j=J[2:end]],
    # η_f_product_I[j] == a_f*T[j] + b_f*I[j]/A + c_f
    # )
    # m.ext[:constraints][:con55b] = @constraint(m, [j=1],
    # η_f_product_I[j] == a_f*T_0 + b_f*I[j]/A + c_f
    # )

    ### POGING 1 ###

    # m.ext[:constraints][:con56] = @constraint(m, [j=J],
    # T_starts[1]*t_1b[j] <= T[j]
    # )
    # m.ext[:constraints][:con56b] = @constraint(m, [j=J],
    # T[j] <= T_starts[2]*t_1b[j]
    # )
    # m.ext[:constraints][:con57] = @constraint(m, [j=J],
    # T_starts[2]*t_2b[j] <= T[j]
    # ) #segment with highest temperature so this one constraint is enough
    # m.ext[:constraints][:con58] = @constraint(m, [j=J],
    # j_starts[1]*j_1b[j] <= I[j]/A
    # )
    # m.ext[:constraints][:con58b] = @constraint(m, [j=J],
    # I[j]/A <= j_starts[2]*j_1b[j]
    # )
    # m.ext[:constraints][:con59] = @constraint(m, [j=J],
    # j_starts[2]*j_2b[j] <= I[j]/A
    # ) #segment with highest temperature so this one constraint is enough
    # m.ext[:constraints][:con60] = @constraint(m, [j=J],
    # 1 == t_1b[j] + t_2b[j]
    # ) 
    # m.ext[:constraints][:con60b] = @constraint(m, [j=J],
    # 1 == j_1b[j] + j_2b[j]
    # )
    # m.ext[:constraints][:con61] = @constraint(m, [j=J],
    # η_f_product_I[j] == t_1b[j]*j_1b[j]*(a_f[1]*T[j] + b_f[1]*I[j]/A + c_f[1]) +t_1b[j]*j_2b[j]*(a_f[2]*T[j] + b_f[2]*I[j]/A + c_f[2]) +t_2b[j]*j_1b[j]*(a_f[3]*T[j] + b_f[3]*I[j]/A + c_f[3]) +t_2b[j]*j_2b[j]*(a_f[4]*T[j] + b_f[4]*I[j]/A + c_f[4]) 
    # )

    ### POGING 2 ### (werkt maar duurt 30min)

    # m.ext[:constraints][:con56] = @constraint(m, [j=J],
    # T_starts[1]*(t_1b[j]+t_2b[j]) <= T[j] * (1- t_3b[j] - t_4b[j])
    # )
    # m.ext[:constraints][:con56b] = @constraint(m, [j=J],
    # T[j] * (1- t_3b[j] - t_4b[j]) <= T_starts[2]*(t_1b[j]+t_2b[j])
    # )
    # m.ext[:constraints][:con57] = @constraint(m, [j=J],
    # T_starts[2]*(t_3b[j]+t_4b[j]) <= T[j] 
    # ) #segment with highest temperature so this one constraint is enough
    # m.ext[:constraints][:con58] = @constraint(m, [j=J],
    # j_starts[1]*(t_1b[j]+t_3b[j]) <= I[j]/A* (1- t_2b[j] - t_4b[j])
    # )
    # m.ext[:constraints][:con58b] = @constraint(m, [j=J],
    # I[j]/A* (1- t_2b[j] - t_4b[j]) <= j_starts[2]*(t_1b[j]+t_3b[j])
    # )
    # m.ext[:constraints][:con59] = @constraint(m, [j=J],
    # j_starts[2]*(t_2b[j]+t_4b[j]) <= I[j]/A
    # ) #segment with highest temperature so this one constraint is enough
    # m.ext[:constraints][:con60] = @constraint(m, [j=J],
    # 1 == t_1b[j] + t_2b[j] + t_3b[j] + t_4b[j]
    # ) 
    
    # m.ext[:constraints][:con61] = @constraint(m, [j=J],
    # η_f_product_I[j] == t_1b[j]*(a_f[1]*T[j] + b_f[1]*I[j]/A + c_f[1]) +t_2b[j]*(a_f[2]*T[j] + b_f[2]*I[j]/A + c_f[2]) +t_3b[j]*(a_f[3]*T[j] + b_f[3]*I[j]/A + c_f[3]) +t_4b[j]*(a_f[4]*T[j] + b_f[4]*I[j]/A + c_f[4]) 
    # )


    # ### POGING 3 ###
    # m.ext[:constraints][:con61] = @constraint(m, [j=J],
    # 1 == t_1b[j] + t_2b[j] + t_3b[j] + t_4b[j]
    # )
    # m.ext[:constraints][:con62] = @constraint(m, [j=J],
    # T_starts[2]*(t_3b[j]+t_4b[j]) <= T[j] 
    # )
    # m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    # T[j] <= T_starts[2]*(t_1b[j]+t_2b[j]) + T_s *(t_3b[j] + t_4b[j])
    # )
    # m.ext[:constraints][:con63] = @constraint(m, [j=J],
    # j_starts[2]*(t_2b[j]+t_4b[j]) <= I[j]/A 
    # )
    # m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    # I[j]/A <= j_starts[2]*(t_1b[j]+t_3b[j]) + iMax *(t_2b[j] + t_4b[j])
    # )
    #  m.ext[:constraints][:con64] = @constraint(m, [j=J],
    # η_f_product_I[j] == t_1b[j]*(a_f[1]*T[j] + b_f[1]*I[j]/A + c_f[1]) +t_2b[j]*(a_f[2]*T[j] + b_f[2]*I[j]/A + c_f[2]) +t_3b[j]*(a_f[3]*T[j] + b_f[3]*I[j]/A + c_f[3]) +t_4b[j]*(a_f[4]*T[j] + b_f[4]*I[j]/A + c_f[4]) 
    # )
    

    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )

    # m.ext[:constraints][:con62] = @constraint(m, [j=J],
    # T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    # )
    # m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    # T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    # )
    # m.ext[:constraints][:con63] = @constraint(m, [j=J],
    # j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    # )
    # m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    # I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    # )
    # m.ext[:constraints][:con64] = @constraint(m, [j=J],
    # η_f_product_I[j] == sum(t_b[i,j]*(a_f[i]*T[j] + b_f[i]*I[j]/A + c_f[i]) for i in 1:4)
    # # η_f_product_I[j] == t_1b[j]*(a_f[1]*T[j] + b_f[1]*I[j]/A + c_f[1]) +t_2b[j]*(a_f[2]*T[j] + b_f[2]*I[j]/A + c_f[2]) +t_3b[j]*(a_f[3]*T[j] + b_f[3]*I[j]/A + c_f[3]) +t_4b[j]*(a_f[4]*T[j] + b_f[4]*I[j]/A + c_f[4]) + t_5b[j]*(a_f[5]*T[j] + b_f[5]*I[j]/A + c_f[5]) +t_6b[j]*(a_f[6]*T[j] + b_f[6]*I[j]/A + c_f[6]) +t_7b[j]*(a_f[7]*T[j] + b_f[7]*I[j]/A + c_f[7]) +t_8b[j]*(a_f[8]*T[j] + b_f[8]*I[j]/A + c_f[8]) + t_9b[j]*(a_f[9]*T[j] + b_f[9]*I[j]/A + c_f[9]) +t_10b[j]*(a_f[10]*T[j] + b_f[10]*I[j]/A + c_f[10]) +t_11b[j]*(a_f[11]*T[j] + b_f[11]*I[j]/A + c_f[11]) +t_12b[j]*(a_f[12]*T[j] + b_f[12]*I[j]/A + c_f[12]) + t_13b[j]*(a_f[13]*T[j] + b_f[13]*I[j]/A + c_f[13]) +t_14b[j]*(a_f[14]*T[j] + b_f[14]*I[j]/A + c_f[14]) +t_15b[j]*(a_f[15]*T[j] + b_f[15]*I[j]/A + c_f[15]) +t_16b[j]*(a_f[16]*T[j] + b_f[16]*I[j]/A + c_f[16]) 
    # )
    m.ext[:constraints][:con62] = @constraint(m, [j=J[2:end]],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=1],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T_0 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J[2:end]],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con62d] = @constraint(m, [j=1],
    T_0 <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
       

    return m
end
function build_model1atest!(m::Model)
    # Model with power in 1 segment, heating up and evaporating water by electrical power, heating up steam 40°C with electrical power

    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data
    π_e = m.ext[:timeseries][:π_e] # Electricity price
    m_LH = m.ext[:timeseries][:m_LH] # = 0
    p_l = m.ext[:timeseries][:p_l] # = 0

    # Extract parameters
    delta_t = m.ext[:parameters]["delta_t"]
    π_H = m.ext[:parameters]["π_H"] 
    # η_f = m.ext[:parameters]["η_f"] 
    F = m.ext[:parameters]["F"]
    iMin = m.ext[:parameters]["iMinSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cSOEC"] 
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    R_t = m.ext[:parameters]["R_tSOEC"] 
    C_HS = m.ext[:parameters]["C_HSSOEC"] 
    C_CS = m.ext[:parameters]["C_CSSOEC"]  
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    M_ini = m.ext[:parameters]["M_ini"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0SOEC"]
    M_HTMin = m.ext[:parameters]["M_HTMin"]
    M_HTMax = m.ext[:parameters]["M_HTMax"]
    p_N = m.ext[:parameters]["p_N"]
    η_turb = m.ext[:parameters]["η_turb"]
    delta_T= m.ext[:parameters]["delta_T"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    delta_H= m.ext[:parameters]["delta_H"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminSOEC"]
    cp_steam = m.ext[:parameters]["cp_steam"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2SOEC1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2SOEC1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2SOEC1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    M_HT = m.ext[:variables][:M_HT] =  @variable(m, [j=J], lower_bound=0, base_name="Mass of stored hydrogen") #kg
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    p_u = m.ext[:variables][:p_u] =  @variable(m, [j=J], lower_bound=0, base_name="Power to utility grid") 
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # p_Q = m.ext[:variables][:p_Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to electrolyser")
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] = @variable(m,[j=J],lower_bound = 0, base_name="Temperature")
    Q_steam = m.ext[:variables][:Q_steam] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to steam")

    # Create affine expressions (= linear combinations of variables)
    # T = m.ext[:expressions][:T] = @expression(m, [j=J],
    # if j == 1
    #     T_0
    # else
    #     sum(3600*delta_t/C_h*(n_c*p[j-i]-n_c*U_tn*I[j-i]+T_a/R_t-Q_cool[j-i])*(1-3600*delta_t/(R_t*C_h))^(i-1) for i in 1:j-1) +T_0* (1-3600*delta_t/(R_t*C_h))^(j-1)
    # end
    # )
    # with heat
    # T = m.ext[:expressions][:T] = @expression(m, [j=J[2:end]],
    #     sum(3600*delta_t/C_h*(n_c*p[j-i]-n_c*U_tn*I[j-i]+p_Q[j]+T_a/R_t-Q_cool[j-i])*(1-3600*delta_t/(R_t*C_h))^(i-1) for i in 1:j-1)+T_0*(1-3600*delta_t/(R_t*C_h))^(j-1)
    # )
    # Without heat
    # T = m.ext[:expressions][:T] = @expression(m, [j=J[2:end]],
    #     sum(3600*delta_t/C_h*(n_c*p[j-i]-n_c*U_tn*I[j-i]+T_a/R_t-Q_cool[j-i])*(1-3600*delta_t/(R_t*C_h))^(i-1) for i in 1:j-1)+T_0*(1-3600*delta_t/(R_t*C_h))^(j-1)
    # )
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F) 
    )
    # mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
    #     M_H2*n_c*η_f*I[j]/(2*F) 
    # )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t + (p_u[j] + p_l)*π_e[j]*delta_t for j in J) - sum(C_HS*Y_b[j-1]+C_CS*Z_b[j-1] for j in J[2:end])
    )

    # Formulate constraints
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
        s_b[j]*p_s[j]+p_b[j] <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= p_b[j]*p_N + s_b[j]*p_s[j]
    )    
    m.ext[:constraints][:con34] = @constraint(m, [j=J[2:end]],
         p_s[j] == (T[j] - T_a)/(R_t*n_c)
    )
    #constraint 35a already implemented when defining variable Q_cool
    m.ext[:constraints][:con35b] = @constraint(m, [j=J],
        Q_cool[j]<=p_b[j]*M
    )
    m.ext[:constraints][:con36] = @constraint(m, [j=J],
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] + delta_2[j]
    )
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -M*s_b[j] <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= M*s_b[j]
    )
    m.ext[:constraints][:con38c] = @constraint(m, [j=J],
        -M*i_b[j] <= delta_2[j] 
    )
    m.ext[:constraints][:con38d] = @constraint(m, [j=J],
        delta_2[j] <= M*i_b[j]
    )
    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
    )
    # Constraint (40) and (41) are replaced by (42) - (45)

    m.ext[:constraints][:con42a] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + p_b[j] - 1 <= Y_b[j-1]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j-1] <= p_b[j]
    )
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j-1]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j-1] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j-1] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved

    m.ext[:constraints][:con47] = @constraint(m, [j=J],
    p_N - n_c*p[j]*(1+s_b[j]*(-1+1/η_EH)) - p_l - p_c[j] - p_u[j] -(Q_H2O[j]+Q_steam[j])/η_EH - Q_cool[j]/400 == 0
    )
    # m.ext[:constraints][:con47] = @constraint(m, [j=J],
    # p_N - n_c*p[j] -n_c*p_Q[j]*η_turb- p_l - p_c[j] - p_u[j] == 0
    # )
    m.ext[:constraints][:con48a] = @constraint(m, [j=J],
    p_uMin <= p_u[j]
    )
    m.ext[:constraints][:con48b] = @constraint(m, [j=J],
    p_u[j] <= p_uMax
    )
    m.ext[:constraints][:con49] = @constraint(m, [j=J[2:end]],
    M_HT[j] - M_HT[j-1] == mdotH2[j]*delta_t - m_LH*delta_t
    )
    m.ext[:constraints][:con50] = @constraint(m, [j=1],
    M_HT[j] == M_ini
    )
    m.ext[:constraints][:con51a] = @constraint(m, [j=J],
    M_HTMin <= M_HT[j]
    )
    m.ext[:constraints][:con51b] = @constraint(m, [j=J],
    M_HT[j] <= M_HTMax
    )
    m.ext[:constraints][:con52] = @constraint(m, [j=J],
    p_c[j] == alfa*mdotH2[j]
    )

    
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    Q_H2O[j] == (mdotH2O[j]*cp_H2O*(373-T_a) + delta_H/M_H2O*mdotH2O[j])
    )
    #constraint amount of heat needed voor heat up necessary mass rate water from T_a to 100°C and evaporate it  
    m.ext[:constraints][:con54b] = @constraint(m, [j=J], 
    Q_steam[j] == mdotH2O[j]*cp_steam*delta_T
    )
    
    return m
end
function build_model2atest!(m::Model)
    # Model with power in 4 segment, heating up and evaporating water by electrical power
    build_model1atest!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_0 = m.ext[:parameters]["T_0SOEC"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC2."Current density starts"

    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j]
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i,j] for i in 1:4)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    T_starts[2]*(t_b[3,j]+t_b[4,j]) <= T[j] 
    )
    m.ext[:constraints][:con62c] = @constraint(m, [j=J],
    T[j] <= T_starts[2]*(t_b[1,j]+t_b[2,j]) + T_s *(t_b[3,j] + t_b[4,j])
    )
    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    j_starts[2]*(t_b[2,j]+t_b[4,j]) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= j_starts[2]*(t_b[1,j]+t_b[3,j]) + iMax *(t_b[2,j] + t_b[4,j])
    )
 

    return m
end
function build_model3atest!(m::Model)
    # Model with power in 9 segment, heating up and evaporating water by electrical power
    build_model1atest!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:9,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Delete constraints
    # Constraint 34 was implemented in parameter earlier so doesn't have to be deleted
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j] + delta_2[j]
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:9)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*3+1:(k-1)*3+3) for k in 2:3) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))
    )

    return m
end
function build_model4atest!(m::Model)
    # Model with power in 16 segment, heating up and evaporating water by electrical power
    build_model1atest!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sSOEC"] 
    iMax = m.ext[:parameters]["iMaxSOEC"] 

    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:16,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 
    p = m.ext[:variables][:p] 
    T = m.ext[:variables][:T] 
 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2SOEC3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2SOEC3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2SOEC3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2SOEC3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2SOEC3."Current density starts"

    # Delete constraints
    # Constraint 34 was implemented in parameter earlier so doesn't have to be deleted
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end
    #New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] + delta_2[j]
    )
    
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))
    )

    return m
end
function build_model5atest!(m::Model)
    # Model with power in 4 segments, with heat addition from LWR for:
        # 1) heating up water from T_a to 100°C and evaporate it
        #( 2) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    build_model2atest!(m)

    # Extract timeseries
    p_l = m.ext[:timeseries][:p_l] # = 0

    # Extract parameters
    p_N = m.ext[:parameters]["p_N"]
    η_turb = m.ext[:parameters]["η_turb"]
    η_EH = m.ext[:parameters]["η_EH"]
    n_c = m.ext[:parameters]["n_cSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    s_b = m.ext[:variables][:s_b]
    p_c = m.ext[:variables][:p_c]
    p_u = m.ext[:variables][:p_u]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    

    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con47][j])
    end

    # New constraints
    m.ext[:constraints][:con47] = @constraint(m, [j=J],
    p_N - n_c*p[j]*(1+s_b[j]*(-1+1/η_EH)) - p_l - p_c[j] - p_u[j] - Q_H2O[j]*η_turb - Q_steam[j]/η_EH - Q_cool[j]/400 == 0
    )

    return m
end
function build_model6atest!(m::Model)
    # Model with power in 4 segments, with heat addition from HTGR for:
        # 1) heating up water from T_a to 100°C and evaporate it with heat 
        # 2) heating up steam last 40°C with heat
        # 3) keeping temperature constant in stand-by with high temperature heat
        # 4) adding heat directly to electrolyser
        #( 5) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    build_model2atest!(m)

    # Extract timeseries
    p_l = m.ext[:timeseries][:p_l] # = 0

    # Extract parameters
    p_N = m.ext[:parameters]["p_N"]
    η_turb = m.ext[:parameters]["η_turb"]
    n_c = m.ext[:parameters]["n_cSOEC"]
    delta_t = m.ext[:parameters]["delta_t"]
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"]
    T_a = m.ext[:parameters]["T_a"]
    R_t = m.ext[:parameters]["R_tSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    s_b = m.ext[:variables][:s_b]
    p_c = m.ext[:variables][:p_c]
    p_u = m.ext[:variables][:p_u]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    
    
    # Create variables
    Q = m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser")

    # Delete constraints
    for j in J[2:end]
        delete(m,m.ext[:constraints][:con30][j])
    end
    for j in J
        delete(m,m.ext[:constraints][:con47][j])
    end
    
    # New constraints
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j]
    )
    m.ext[:constraints][:con47] = @constraint(m, [j=J],
    p_N - n_c*p[j]*(1+s_b[j]*(-1+η_turb)) - p_l - p_c[j] - p_u[j] - (Q_H2O[j]+ Q_steam[j]+Q[j])*η_turb - Q_cool[j]/400 == 0
    )
    
    # Voorlopig geen constraints op Q!
    # m.ext[:constraints][:con70] = @constraint(m, [j=J],
    # Q[j] <= n_c*p[j-1] - n_c*U_tn*I[j-1]
    # )
    
    return m
end
function build_model7atest!(m::Model)
    # Model with power in 9 segments, with heat addition from HTGR for:
        # 1) heating up water from T_a to 100°C and evaporate it with heat 
        # 2) heating up steam last 40°C with heat
        # 3) keeping temperature constant in stand-by with high temperature heat
        # 4) adding heat directly to electrolyser
        #( 5) Partly heat up electrolyser during cold start, just changes cold start up cost, NOT INTEGRATED IN MODEL!!!) 

    build_model3atest!(m)

    # Extract timeseries
    p_l = m.ext[:timeseries][:p_l] # = 0

    # Extract parameters
    p_N = m.ext[:parameters]["p_N"]
    η_turb = m.ext[:parameters]["η_turb"]
    n_c = m.ext[:parameters]["n_cSOEC"]
    delta_t = m.ext[:parameters]["delta_t"]
    U_tn = m.ext[:parameters]["U_tnSOEC"]
    C_h = m.ext[:parameters]["C_hSOEC"]
    T_a = m.ext[:parameters]["T_a"]
    R_t = m.ext[:parameters]["R_tSOEC"]

    # Extract variables
    p = m.ext[:variables][:p]
    s_b = m.ext[:variables][:s_b]
    p_c = m.ext[:variables][:p_c]
    p_u = m.ext[:variables][:p_u]
    Q_H2O = m.ext[:variables][:Q_H2O]
    Q_cool = m.ext[:variables][:Q_cool]
    Q_steam = m.ext[:variables][:Q_steam]
    
    
    # Create variables
    Q = m.ext[:variables][:Q] = @variable(m,[j=J],lower_bound = 0, base_name="Heat directly to electrolyser")

    # Delete constraints
    for j in J[2:end]
        delete(m,m.ext[:constraints][:con30][j])
    end
    for j in J
        delete(m,m.ext[:constraints][:con47][j])
    end
    
    # New constraints
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t)== n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] + Q[j]
    )
    m.ext[:constraints][:con47] = @constraint(m, [j=J],
    p_N - n_c*p[j]*(1+s_b[j]*(-1+η_turb)) - p_l - p_c[j] - p_u[j] - (Q_H2O[j]+ Q_steam[j]+Q[j])*η_turb - Q_cool[j]/400 == 0
    )
    
    # Voorlopig geen constraints op Q!
    # m.ext[:constraints][:con70] = @constraint(m, [j=J],
    # Q[j] <= n_c*p[j-1] - n_c*U_tn*I[j-1]
    # )
    
    return m
end


### DAGEN ###

objective_values = Matrix{Float64}(undef, 367*24, 12)
define_sets!(m, data,24)
process_parameters!(m, data)
for j in 1:1
    operating_hours = 0
    average_eff = 0
    m.ext[:parameters]["π_H"] =2.5 + (j - 1) 
    T_end = 1073
    for i in 55:55
        m.ext[:parameters]["T_0SOEC"] = T_end
        m.ext[:parameters]["C_hSOEC"] = m.ext[:parameters]["C_hSOEC"]
        m.ext[:parameters]["R_tSOEC"] = m.ext[:parameters]["R_tSOEC"]
        process_time_series_data_day!(m, data, i)
        build_model2g!(m)
        # set_optimizer_attribute(m, "mip_gap", 0.001)
        optimize!(m)
        objective_values[i, j] = value.(m.ext[:objective])
        operating_hours = operating_hours + value.(m.ext[:expressions][:p_bSUM]) 
        objective_values[(i-1)*24+1:i*24, 4+j] = [value.(m.ext[:expressions][:mdotH2])[k]*value.(m.ext[:variables][:p_b])[k] for k in 1:24]
        objective_values[(i-1)*24+1:i*24, 8+j] = [value.(m.ext[:expressions][:p_E])[k] for k in 1:24]
        T_end = value.(first(m.ext[:expressions][:T25])) #geen constraint op T[24]!!
        print(
            """
            hprice = $j
            dag = $i
            T_end = $T_end
            """
        )
    end
    objective_values[366, j] = operating_hours
    objective_values[367, j] = average_eff
end
matfile = matopen("objective_values_SOEC.mat", "w")
header = "objective_values"
write(matfile, header, objective_values)
close(matfile)

# ### WEKEN ### 
# objective_values = Matrix{Float64}(undef, 367*24, 12)
# define_sets!(m, data,7*24)
# process_parameters!(m, data)
# for j in 1:4
#     operating_hours = 0
#     average_eff = 0
#     m.ext[:parameters]["π_H"] = 2.5 + (j - 1) 
#     for i in 1:52
#         process_time_series_data_week!(m, data, i)
#         build_model2e2!(m)
#         set_optimizer_attribute(m, "mip_gap", 0.001)
#         optimize!(m)
#         objective_values[i, j] = value.(m.ext[:objective])
#         operating_hours = operating_hours + value.(m.ext[:expressions][:p_bSUM]) 
#         objective_values[(i-1)*7*24+1:i*7*24, 4+j] = [value.(m.ext[:expressions][:mdotH2])[k]*value.(m.ext[:variables][:p_b])[k] for k in 1:7*24]
#         objective_values[(i-1)*7*24+1:i*7*24, 8+j] = [value.(m.ext[:expressions][:p_E])[k] for k in 1:7*24]
#     end
#     objective_values[53, j] = operating_hours
#     objective_values[53, j] = average_eff
# end
# matfile = matopen("objective_values_SOEC.mat", "w")
# header = "objective_values"
# write(matfile, header, objective_values)
# close(matfile)



### TESTEN ###

# function process_time_series_data!(m::Model, data::Dict)
#     # extract the relevant sets
#     #J = m.ext[:sets][:J] # Time steps

#     # create dictionary to store time series
#     m.ext[:timeseries] = Dict()

#     # example: add time series to dictionary
#     m.ext[:timeseries][:π_e] = ePrice."Day-ahead Price [EUR/MWh]"[20*7*24+1:3*7*24]/1000000 #[€/Wh]


#     # Both 0, very large storage tank (no problem bcs in case 1 in paper storage never full so no determening constraint)
#     # p_l = 0 so just more p_u (no power to local load but everything to utility grid)
#     m.ext[:timeseries][:m_LH] = 0
#     m.ext[:timeseries][:p_l] = 0 
#     # return model
#     return m
# end
# process_time_series_data!(m, data)
# define_sets!(m, data,7*24)
# process_parameters!(m, data)
# build_model2e3!(m)
# set_optimizer_attribute(m, "mip_gap", 0.001)

# ## Step 4: solve
# optimize!(m)

# # check termination status
# print(
#     """

#     Termination status: $(termination_status(m))

#     """
# )


# # print some output
# @show value(m.ext[:objective])


# Step 5: interpretation
# access data

# sets
J = m.ext[:sets][:J]
# # parameters
eprice = value.(m.ext[:timeseries][:π_e])
n_c = value.(m.ext[:parameters]["n_cSOEC"])
T_0 = value.(m.ext[:parameters]["T_0SOEC"])-273 #°C
p_N = value.(m.ext[:parameters]["p_N"])
# variables/expressions
p = value.(m.ext[:variables][:p])
s_b = value.(m.ext[:variables][:s_b])
p_b = value.(m.ext[:variables][:p_b])
# Y_b = value.(m.ext[:variables][:Y_b])
Q_cool = value.(m.ext[:variables][:Q_cool])
Q = value.(m.ext[:variables][:Q])
# Q_H2O = value.(m.ext[:variables][:Q_H2O])
# p_c = value.(m.ext[:variables][:p_c])
I = value.(m.ext[:variables][:I])
delta_1 = value.(m.ext[:variables][:delta_1])

# T = value.(m.ext[:expressions][:T]) #K
T = value.(m.ext[:variables][:T]) #K

mdotH2 = value.(m.ext[:expressions][:mdotH2])
# mdotH2 = value.(m.ext[:variables][:mdotH2])
p_s = value.(m.ext[:variables][:p_s])
p_svec = [p_s[j] for j in J] #[MW]


# create arrays for plotting
pvec = [p[j]*n_c/1000000 for j in J] #*n_c for whole Electrolyser /1000000 so in [MW]
epricevec = [eprice[j]*1000000 for j in J] #[€/MWh]
p_Nvec =  [p_N/1000000 for j in J] #[MW]
Tvec = [T[j]-273 for j in J] #°C
s_bvec = [round(s_b[j]) for j in J]
p_bvec = [round(p_b[j]) for j in J]
# Y_bvec = [Y_b[j] for j in J]
Q_coolvec = [Q_cool[j] for j in J]
Qvec = [Q[j] for j in J]
# p_uvec = [p_u[j]/1000000 for j in J] #[MW]
# p_cvec = [p_c[j]/1000000 for j in J] 
Ivec = [I[j] for j in J] 
mdotH2vec = [mdotH2[j] for j in J] 

using Plots
using LaTeXStrings
using StatsPlots
using Plots.PlotMeasures


pyplot()
plot(pvec, label = "Electrolyser Power", xlabel="Time [h]",ylabel = "Electrolyser power [MW]", 
        legend = :topleft, left_margin = 5Plots.mm, right_margin = 15Plots.mm)
plot!(twinx(), epricevec, label = "Electricity price", ylabel="Electricity Price [€/MWh]", legend = :topright,color = :red)
savefig("power_price.png")

bar(J.-1/2, pvec, label = "Electrolyser Power", xlabel = "Time [h]", ylabel = "Electrolyser power [MW]",
    legend = :topleft, left_margin = 5Plots.mm, right_margin = 15Plots.mm, bar_width = 1.0, fillalpha = 0.5, linealpha = 0.1)

# Create a line plot for epricevec on the secondary y-axis
plot!(twinx(), J.-1/2, epricevec, label = "Electricity price", ylabel = "Electricity Price [€/MWh]",
      line = :path, legend = :topright, color = :red)

# Save the plot to a file
savefig("power_price_bars.png")

plot(J,Qvec, xlabel = "Time [h]", ylabel = "Heat directly to electrolyser [MW]",legend = false, left_margin = 5Plots.mm, right_margin = 15Plots.mm)
savefig("heat directly to electrolyser.png")

plot(epricevec, xlabel = "Time [h]", ylabel="Electricity Price [€/MWh]", legend=false)
savefig("eprice.png")

plot(J.-1/2, [Tvec], xlabel = "Time [h]", ylabel = "Electrolyser Temperature [°C]",legend = false, left_margin = 5Plots.mm, right_margin = 15Plots.mm)
savefig("Temperature.png")

bar(J.-1/2, [Tvec], xlabel = "Time [h]", ylabel = "Electrolyser Temperature [°C]", legend = false, left_margin = 5Plots.mm, right_margin = 15Plots.mm, bar_width = 1.0, fillalpha = 0.5, linealpha = 0.1)
savefig("Temperaturebars.png")


plot(J, [p_bvec s_bvec], label = ["Production state"  "Standby state" "Off state"], xlabel = "Time [h]", legend = :right,left_margin = 5Plots.mm, right_margin = 15Plots.mm)
savefig("Binary variables.png")
bar(J.+ 1/2, [p_bvec s_bvec], label = ["Production state"  "Standby state" "Off state"], xlabel = "Time [h]", legend = :right,left_margin = 5Plots.mm, right_margin = 15Plots.mm, bar_width = 1.0, fillalpha = 0.5, linealpha = 0.1)
savefig("Binary variablesbars.png")

# plot(J, [Y_bvec], label = ["Hot start"], xlabel = "Time [h]", left_margin = 5Plots.mm, right_margin = 15Plots.mm )
# savefig("Binary variables hot start.png")

plot(J, Q_coolvec, xlabel = "Time [h]", ylabel = "Cooling Power [W]",legend = false,left_margin = 5Plots.mm, right_margin = 15Plots.mm )
savefig("Cooling power.png")

# plot(J, p_uvec, xlabel = "Time [h]", ylabel = "Power to utility grid [W]",left_margin = 5Plots.mm, right_margin = 15Plots.mm )
# savefig("Power to utility grid.png")

# plot([p_uvec pvec p_Nvec p_cvec], label = ["Power to utility grid" "Electrolyser power" "Wind power" "Compressor power"], xlabel="Time [h]",ylabel = "Power [MW]", legend = :topleft, left_margin = 5Plots.mm, right_margin = 15Plots.mm )
# savefig("powers.png")

plot(J, Ivec, xlabel = "Time [h]", ylabel = "Electrolyser current [A]",left_margin = 5Plots.mm, right_margin = 15Plots.mm )
savefig("Electrolyser current.png")

plot(J, mdotH2vec, xlabel = "Time [h]", ylabel = "Hydrogen production rate [kg/s]",left_margin = 5Plots.mm, right_margin = 15Plots.mm )
savefig("Hydrogen production rate.png")

colours = ifelse.(p_bvec.==1, :green, :orange)
states = ["Production state"  "Standby state"]
# Get the unique colors in the order they appear
unique_colours = unique(colours)

# Create a blank plot
plot()


# Iterate over each unique color
for color in unique_colours
    # Get the indices of bars with the current color
    indices = findall(colours .== color)
    index = findall(unique_colours .== color)[1]
    # Get the corresponding bar heights
    bar_heights = [Tvec[indices]]
    
    # Get the corresponding x-axis positions
    bar_positions = indices .- 0.5
    
    # Create a bar plot for the current color
    bar!(bar_positions, bar_heights, fillcolor = color, bar_width = 1.0,
         fillalpha = 0.5, linealpha = 0.1, label = states[index], legend = :bottomright)
end



# Customize the plot
xlabel!("Time [h]")
ylabel!("Electrolyser Temperature [°C]")


# Save the plot to a file
savefig("Temperaturebarscolours.png")