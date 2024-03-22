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
ePrice = CSV.read(joinpath(@__DIR__, "Day-ahead Prices_2019.csv"), DataFrame, silencewarnings=true)

coeffs2PEM1 = CSV.read(joinpath(@__DIR__, "coeffs2PEM1.csv"), DataFrame, silencewarnings=true)
coeffs2PEM2 = CSV.read(joinpath(@__DIR__, "coeffs2PEM2.csv"), DataFrame, silencewarnings=true)
coeffs2PEM3 = CSV.read(joinpath(@__DIR__, "coeffs2PEM3.csv"), DataFrame, silencewarnings=true)
coeffs2PEM4 = CSV.read(joinpath(@__DIR__, "coeffs2PEM4.csv"), DataFrame, silencewarnings=true)
coeffs2PEM5 = CSV.read(joinpath(@__DIR__, "coeffs2PEM5.csv"), DataFrame, silencewarnings=true)
coeffs2PEM6 = CSV.read(joinpath(@__DIR__, "coeffs2PEM6.csv"), DataFrame, silencewarnings=true)

## Step 2: create model & pass data to model
using JuMP
using Gurobi
m = Model(optimizer_with_attributes(Gurobi.Optimizer))

# Step 2a: create sets
function define_sets!(m::Model, data::Dict,i::Int64)
    # create dictionary to store sets
    m.ext[:sets] = Dict()

    # define the sets
    m.ext[:sets][:J] = 1:i# Timesteps 
    m.ext[:sets][:parameterset] = [i for i in keys(data)]

    # return model
    return m
end

# Step 2b: add time series
function process_time_series_data!(m::Model, data::Dict, i::Int64)
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
    m.ext[:parameters]["A"] = coeffs2PEM1."Area [m^2]"[1]
    m.ext[:parameters]["M_HTMax"] = m.ext[:parameters]["SOCMax"] * m.ext[:parameters]["M_tank"]
    m.ext[:parameters]["M_HTMin"] = m.ext[:parameters]["SOCMin"] * m.ext[:parameters]["M_tank"]
    m.ext[:parameters]["alfa"] = m.ext[:parameters]["R"]*m.ext[:parameters]["T_sPEM"]/(2*(m.ext[:parameters]["gamma"]-1)*m.ext[:parameters]["η_c"])*((m.ext[:parameters]["P_out"]/m.ext[:parameters]["P_in"])^((m.ext[:parameters]["gamma"]-1)/m.ext[:parameters]["gamma"])-1)
    m.ext[:parameters]["M_ini"] = m.ext[:parameters]["SOC_ini"] * m.ext[:parameters]["M_tank"]
    # m.ext[:parameters]["p_s"] = (m.ext[:parameters]["T_sPEM"]-m.ext[:parameters]["T_a"])/(m.ext[:parameters]["R_tPEM"]*m.ext[:parameters]["n_c"])
    # return model
    return m
end

# # call functions
# define_sets!(m, data)
# process_parameters!(m, data)

## Step 3: construct your model

function build_model1a!(m::Model)
    # Model with power in 1 segment
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
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
    iMin = m.ext[:parameters]["iMinPEM"] 
    iMax = m.ext[:parameters]["iMaxPEM"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cPEM"] 
    U_tn = m.ext[:parameters]["U_tnPEM"]
    C_h = m.ext[:parameters]["C_hPEM"] 
    T_s = m.ext[:parameters]["T_sPEM"] 
    R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    C_CS = m.ext[:parameters]["C_CSPEM"]  
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0PEM"]
    p_N = m.ext[:parameters]["p_N"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminPEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2PEM1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2PEM1."Coefficients j [W/(A/m^2)]"[1]


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
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")


    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    
    


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t + p_u[j]*π_e[j]*delta_t for j in J) - sum(C_HS*Y_b[j-1]+C_CS*Z_b[j-1] for j in J[2:end])
    )

    # Formulate constraints
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin*s_b[j] <= T[j] 
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
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
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
        p_N - n_c*p[j]*(1+s_b[j]*(-1+1/η_EH))  - p_c[j] - p_u[j] - Q_cool[j]/400 - Q_H2O[j]/η_EH == 0
    )
    m.ext[:constraints][:con48a] = @constraint(m, [j=J],
        p_uMin <= p_u[j]
    )
    m.ext[:constraints][:con48b] = @constraint(m, [j=J],
        p_u[j] <= p_uMax
    )
      m.ext[:constraints][:con52] = @constraint(m, [j=J],
        p_c[j] == alfa*mdotH2[j]
    )
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )
    return m
end
function build_model2a!(m::Model)
    # Model with power in 4 segments
    build_model1a!(m)
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    T = m.ext[:variables][:T]
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2]
    p = m.ext[:variables][:p]
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM2."Current density starts"

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
function build_model3a!(m::Model)
    # Model with power in 9 segments
    build_model1a!(m)
    
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] = @variable(m,[i=1:9,j=J], binary=true, base_name="Binary variable for segment")

    # Extract variables
    I = m.ext[:variables][:I] 
    T = m.ext[:variables][:T] 
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 
    p = m.ext[:variables][:p]
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM3."Current density starts"
    

    # Delete constraints
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
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )

    
    return m
end
function build_model4a!(m::Model)
    # Model with power in 16 segments, T as a variable
    build_model1a!(m)
        
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] = @variable(m,[i=1:16,j=J], binary=true, base_name="Binary variable for segment")

    # Extract variables
    I = m.ext[:variables][:I] 
    p = m.ext[:variables][:p]
    T = m.ext[:variables][:T] 

    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM4."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM4."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM4."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM4."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM4."Current density starts"
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 


    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end

    # New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] + delta_2[j]
    )
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )

    return m
end
function build_model5a!(m::Model)
    # Model with eff_farad*current in 4 segments, power in 4 segment
    # Nuclear heat used for:
        # 1) heating up water from 46.8°C to 100°C
        # 2) keeping temperature constant in standby
        # 3) No effect on cold start up cost because minimal operating temperature is 20°C!!!

    build_model2a!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 

    # Extract parameters
    p_N = m.ext[:parameters]["p_N"] 
    n_c = m.ext[:parameters]["n_cPEM"] 
    η_turb = m.ext[:parameters]["η_turb"] 
    # Extract variables
    s_b = m.ext[:variables][:s_b] 
    p_c = m.ext[:variables][:p_c] 
    p_u = m.ext[:variables][:p_u] 
    Q_cool = m.ext[:variables][:Q_cool] 
    Q_H2O = m.ext[:variables][:Q_H2O] 
    p = m.ext[:variables][:p]
   
    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con47][j])
    end
   
    #New constraints
    m.ext[:constraints][:con47] = @constraint(m, [j=J],
    p_N - n_c*p[j]*(1+s_b[j]*(-1+η_turb)) - p_c[j] - p_u[j] - Q_cool[j]/400 - Q_H2O[j]*η_turb == 0
    )
    
    return m
end


function build_model1b!(m::Model)
    # Model with power in 1 segment
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # EXTRA With new objective

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
    iMin = m.ext[:parameters]["iMinPEM"] 
    iMax = m.ext[:parameters]["iMaxPEM"] 
    A = m.ext[:parameters]["A"]   
    n_c = m.ext[:parameters]["n_cPEM"] 
    U_tn = m.ext[:parameters]["U_tnPEM"]
    C_h = m.ext[:parameters]["C_hPEM"] 
    T_s = m.ext[:parameters]["T_sPEM"] 
    R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    C_CS = m.ext[:parameters]["C_CSPEM"]  
    p_uMin = m.ext[:parameters]["p_uMin"]
    p_uMax = m.ext[:parameters]["p_uMax"] 
    T_a = m.ext[:parameters]["T_a"]
    M_H2 = m.ext[:parameters]["M_H2"]
    M = m.ext[:parameters]["M"]
    alfa = m.ext[:parameters]["alfa"]
    T_0 = m.ext[:parameters]["T_0PEM"]
    p_N = m.ext[:parameters]["p_N"]
    cp_H2O= m.ext[:parameters]["cp_H2O"]
    M_H2O= m.ext[:parameters]["M_H2O"]
    M_H2= m.ext[:parameters]["M_H2"]
    η_EH = m.ext[:parameters]["η_EH"]
    Tmin = m.ext[:parameters]["TminPEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2PEM1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2PEM1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name ="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    p_u = m.ext[:variables][:p_u] =  @variable(m, [j=J], lower_bound=0, base_name="Power to utility grid") 
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")


    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    
    


    # Formulate objective
     m.ext[:objective] = @objective(m, Max,
        sum(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t for j in J) - sum(C_HS*Y_b[j-1]+C_CS*Z_b[j-1] for j in J[2:end])
    )
    # Formulate constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+1/η_EH)) + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin*s_b[j] <= T[j] 
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
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
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

    m.ext[:constraints][:con48a] = @constraint(m, [j=J],
        p_uMin <= p_u[j]
    )
    m.ext[:constraints][:con48b] = @constraint(m, [j=J],
        p_u[j] <= p_uMax
    )
    m.ext[:constraints][:con52] = @constraint(m, [j=J],
        p_c[j] == alfa*mdotH2[j]
    )
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
        Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )
    return m
end
function build_model2b!(m::Model)
    # Model with power in 4 segments
    build_model1b!(m)
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    T = m.ext[:variables][:T]
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2]
    p = m.ext[:variables][:p]
    
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM2."Current density starts"

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
function build_model3b!(m::Model)
    # Model with power in 9 segments
    build_model1b!(m)
    
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] = @variable(m,[i=1:9,j=J], binary=true, base_name="Binary variable for segment")

    # Extract variables
    I = m.ext[:variables][:I] 
    T = m.ext[:variables][:T] 
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 
    p = m.ext[:variables][:p]
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM3."Current density starts"
    

    # Delete constraints
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
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )

    
    return m
end
function build_model4b!(m::Model)
    # Model with power in 16 segments, T as a variable
    build_model1b!(m)
        
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] = @variable(m,[i=1:16,j=J], binary=true, base_name="Binary variable for segment")

    # Extract variables
    I = m.ext[:variables][:I] 
    p = m.ext[:variables][:p]
    T = m.ext[:variables][:T] 

    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM4."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM4."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM4."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM4."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM4."Current density starts"
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 


    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end

    # New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] + delta_2[j]
    )
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )

    return m
end
function build_model5b!(m::Model)
    # Model with power in 4 segment
    # Nuclear heat used for:
        # 1) heating up water from 46.8°C to 100°C
        # 2) keeping temperature constant in standby
        # 3) No effect on cold start up cost because minimal operating temperature is 20°C!!!

    build_model2b!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 

    # Extract parameters
    n_c = m.ext[:parameters]["n_cPEM"] 
    η_turb = m.ext[:parameters]["η_turb"] 
    # Extract variables
    s_b = m.ext[:variables][:s_b] 
    p_c = m.ext[:variables][:p_c] 
    p_E = m.ext[:variables][:p_E] 
    Q_cool = m.ext[:variables][:Q_cool] 
    Q_H2O = m.ext[:variables][:Q_H2O] 
    p = m.ext[:variables][:p]
   
    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end
   
    #New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+η_turb)) + p_c[j] + Q_H2O[j]*η_turb + Q_cool[j]/400 
    )
    
    
    return m
end

function build_model1c!(m::Model)
    # Model with power in 1 segment
   # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
   # EXTRA With new objective
   # EXTRA With cyclic boundary conditions
   
   build_model1b!(m)

   # Extract sets
   J = m.ext[:sets][:J]

   # Extract parameters
   T_0 = m.ext[:parameters]["T_0PEM"]
   # Extract variables
   T = m.ext[:variables][:T] 

   #Extra constraints
   m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
       T[j] == T_0
   )
   
   return m
end
function build_model2c!(m::Model)
   # Model with power in 4 segments    
   # Electrical energy used to heat up water 47°C and keep temperature constant in standby
   
   build_model2b!(m)

   # Extract sets
   J = m.ext[:sets][:J]

   # Extract parameters
   T_0 = m.ext[:parameters]["T_0PEM"]
   # Extract variables
   T = m.ext[:variables][:T] 

   #Extra constraints
   m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
       T[j] == T_0
   )
   
   return m
end
function build_model3c!(m::Model)
   # Model with power in 9 segments
   # Electrical energy used to heat up water 47°C and keep temperature constant in standby

   
   build_model3b!(m)

   # Extract sets
   J = m.ext[:sets][:J]

   # Extract parameters
   T_0 = m.ext[:parameters]["T_0PEM"]
   # Extract variables
   T = m.ext[:variables][:T] 

   #Extra constraints
   m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
       T[j] == T_0
   )
   
   return m
end
function build_model4c!(m::Model)
   # Model with power in 16 segments, T as a variable
   # Electrical energy used to heat up water 47°C and keep temperature constant in standby

   
   build_model4b!(m)

   # Extract sets
   J = m.ext[:sets][:J]

   # Extract parameters
   T_0 = m.ext[:parameters]["T_0PEM"]
   # Extract variables
   T = m.ext[:variables][:T] 

   #Extra constraints
   m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
       T[j] == T_0
   )
   
   return m
end
function build_model5c!(m::Model)
   # Model with power in 4 segment
   # Nuclear heat used for:
       # 1) heating up water from 46.8°C to 100°C
       # 2) keeping temperature constant in standby
       # 3) No effect on cold start up cost because minimal operating temperature is 20°C!!!

   
   build_model5b!(m)

   # Extract sets
   J = m.ext[:sets][:J]

   # Extract parameters
   T_0 = m.ext[:parameters]["T_0PEM"]
   # Extract variables
   T = m.ext[:variables][:T] 

   #Extra constraints
   m.ext[:constraints][:con101] = @constraint(m, [j=J[end]],
       T[j] == T_0
   )
   
   return m
end



# PROBLEEM: segment 1 gaat tot 333K, op bepaald pu  
function build_model1e!(m::Model)
    # Model with eff_farad*current in 1 segment, power in 1 segment
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2PEM1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2PEM1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    delta_3 = m.ext[:variables][:delta_3] = @variable(m, [j=J], lower_bound=-M, base_name="delta 3") #slack variable 
    delta_4 = m.ext[:variables][:delta_4] = @variable(m, [j=J], lower_bound=-M, base_name="delta 4") #slack variable 
    delta_5 = m.ext[:variables][:delta_5] = @variable(m, [j=J], lower_bound=-M, base_name="delta 5")
    delta_6 = m.ext[:variables][:delta_6] = @variable(m, [j=J], lower_bound=-M, base_name="delta 6")
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
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    C_CS = m.ext[:variables][:C_HS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
  
    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    )

    
    # Formulate constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+1/η_EH)) + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == (n_c*p[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1])*(1-SUTime/3600*Z_b[j-1])
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin*s_b[j] <= T[j] 
    )
    m.ext[:constraints][:con31c] = @constraint(m, [j=J[2:end]],
        T_a <= T[j]
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
    )
    m.ext[:constraints][:con32b] = @constraint(m, [j=J],
        I[j] <= A*p_b[j]*iMax
    )
    m.ext[:constraints][:con33a] = @constraint(m, [j=J],
        s_b[j]*p_s[j] <= p[j]
    )
    m.ext[:constraints][:con33b] = @constraint(m, [j=J],
        p[j] <= s_b[j]*p_s[j]+ p_b[j]*p_N
    )    
    m.ext[:constraints][:con34] = @constraint(m, [j=J],
         (T[j] - T_a)/(R_t*n_c) <= p_s[j]
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
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
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
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
        Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con206] = @constraint(m, [j=J],
       C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    )
    m.ext[:constraints][:con207] = @constraint(m, [j=J[2:end]],
        i_b[j] + s_b[j-1] <= 1
    )
    return m
end
function build_model2e!(m::Model)
    # Model with power in 4 segments
    build_model1e!(m)
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    T = m.ext[:variables][:T]
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2]
    p = m.ext[:variables][:p]
    
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM2."Current density starts"

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
function build_model3e!(m::Model)
    # Model with power in 9 segments
    build_model1e!(m)
    
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] = @variable(m,[i=1:9,j=J], binary=true, base_name="Binary variable for segment")

    # Extract variables
    I = m.ext[:variables][:I] 
    T = m.ext[:variables][:T] 
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 
    p = m.ext[:variables][:p]
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM3."Current density starts"
    

    # Delete constraints
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
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )

    
    return m
end
function build_model4e!(m::Model)
    # Model with power in 16 segments, T as a variable
    build_model1e!(m)
        
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] = @variable(m,[i=1:16,j=J], binary=true, base_name="Binary variable for segment")

    # Extract variables
    I = m.ext[:variables][:I] 
    p = m.ext[:variables][:p]
    T = m.ext[:variables][:T] 

    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM4."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM4."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM4."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM4."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM4."Current density starts"
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 


    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end

    # New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] + delta_2[j]
    )
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )

    return m
end
function build_model5e!(m::Model)
    # Model with power in 4 segment
    # Nuclear heat used for:
        # 1) heating up water from 46.8°C to 100°C
        # 2) keeping temperature constant in standby
        # 3) No effect on cold start up cost because minimal operating temperature is 20°C!!!

    build_model2e!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 

    # Extract parameters
    n_c = m.ext[:parameters]["n_cPEM"] 
    η_turb = m.ext[:parameters]["η_turb"] 
    # Extract variables
    s_b = m.ext[:variables][:s_b] 
    p_c = m.ext[:variables][:p_c] 
    p_E = m.ext[:variables][:p_E] 
    Q_cool = m.ext[:variables][:Q_cool] 
    Q_H2O = m.ext[:variables][:Q_H2O] 
    p = m.ext[:variables][:p]
   
    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end
   
    #New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]*(1+s_b[j]*(-1+η_turb)) + p_c[j] + Q_H2O[j]*η_turb + Q_cool[j]/400 
    )
    
    
    return m
end


function build_model1f!(m::Model)
    # Model with eff_farad*current in 1 segment, power in 1 segment
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2PEM1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2PEM1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    delta_3 = m.ext[:variables][:delta_3] = @variable(m, [j=J], lower_bound=-M, base_name="delta 3") #slack variable 
    delta_4 = m.ext[:variables][:delta_4] = @variable(m, [j=J], lower_bound=-M, base_name="delta 4") #slack variable 
    delta_5 = m.ext[:variables][:delta_5] = @variable(m, [j=J], lower_bound=-M, base_name="delta 5")
    delta_6 = m.ext[:variables][:delta_6] = @variable(m, [j=J], lower_bound=-M, base_name="delta 6")
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
    Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
  
    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    )

    
    # Formulate constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == (n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1])*(1-SUTime/3600*Z_b[j-1])
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin*s_b[j] <= T[j] 
    )
    m.ext[:constraints][:con31c] = @constraint(m, [j=J[2:end]],
        T_a <= T[j]
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved

    m.ext[:constraints][:con52] = @constraint(m, [j=J],
        p_c[j] == alfa*mdotH2[j]
    )
    m.ext[:constraints][:con54] = @constraint(m, [j=J], 
        Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con206] = @constraint(m, [j=J],
       C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    )
    m.ext[:constraints][:con207] = @constraint(m, [j=J[2:end]],
        i_b[j] + s_b[j-1] <= 1
    )
    return m
end
function build_model2f!(m::Model)
    # Model with power in 4 segments
    build_model1f!(m)
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Extract variables
    I = m.ext[:variables][:I] 
    T = m.ext[:variables][:T]
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2]
    p = m.ext[:variables][:p]
    
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM2."Current density starts"

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
function build_model3f!(m::Model)
    # Model with power in 9 segments
    build_model1f!(m)
    
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] = @variable(m,[i=1:9,j=J], binary=true, base_name="Binary variable for segment")

    # Extract variables
    I = m.ext[:variables][:I] 
    T = m.ext[:variables][:T] 
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 
    p = m.ext[:variables][:p]
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM3."Current density starts"
    

    # Delete constraints
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
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )

    
    return m
end
function build_model4f!(m::Model)
    # Model with power in 16 segments, T as a variable
    build_model1f!(m)
        
    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 


    # Extract parameters
    A = m.ext[:parameters]["A"] 
    T_s = m.ext[:parameters]["T_sPEM"]
    iMax = m.ext[:parameters]["iMaxPEM"]
    # Create variables
    t_b = m.ext[:variables][:t_b] = @variable(m,[i=1:16,j=J], binary=true, base_name="Binary variable for segment")

    # Extract variables
    I = m.ext[:variables][:I] 
    p = m.ext[:variables][:p]
    T = m.ext[:variables][:T] 
    delta_1 = m.ext[:variables][:delta_1] 
    delta_2 = m.ext[:variables][:delta_2] 
    # Create parameters
    c = m.ext[:parameters]["c"] = coeffs2PEM4."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM4."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM4."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM4."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM4."Current density starts"



    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con37][j])
    end

    # New constraints
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
    p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:16) + delta_1[j] + delta_2[j]
    )
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:16)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*4+1:(k-1)*4+4) for k in 2:4) <= T[j]
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*4+1:(k-2)*4+4) for k in 2:4)+ T_s*sum(t_b[i, j] for i in 13:16)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+4*i,j] for i in 0:3) for k in 2:4) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+4*i,j] for i in 0:3) for k in 2:4) + iMax*(sum(t_b[4*i,j] for i in 1:4))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )

    return m
end
function build_model5f!(m::Model)
    # Model with power in 4 segment
    # Nuclear heat used for:
        # 1) heating up water from 46.8°C to 100°C
        # 2) keeping temperature constant in standby
        # 3) No effect on cold start up cost because minimal operating temperature is 20°C!!!

    build_model2f!(m)

    # Extract sets
    J = m.ext[:sets][:J]

    # Extract time series data 

    # Extract parameters
    n_c = m.ext[:parameters]["n_cPEM"] 
    η_turb = m.ext[:parameters]["η_turb"] 
    # Extract variables
    p_c = m.ext[:variables][:p_c] 
    p_E = m.ext[:variables][:p_E] 
    Q_cool = m.ext[:variables][:Q_cool] 
    Q_H2O = m.ext[:variables][:Q_H2O] 
    p = m.ext[:variables][:p]
    p_s = m.ext[:variables][:p_s]
   
    # Delete constraints
    for j in J
        delete(m,m.ext[:constraints][:con29][j])
    end
   
    #New constraints
    m.ext[:constraints][:con29] = @constraint(m, [j=J],
        p_E[j] == n_c*p[j]+n_c*p_s[j]*η_turb + p_c[j] + Q_H2O[j]*η_turb + Q_cool[j]/400 
    )

    
    return m
end

function build_model1f2!(m::Model)
    # power in 1 segment
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2PEM1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2PEM1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    # p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    p = m.ext[:expressions][:p] = @expression(m, [j=J],
         a*T[j]+ b*I[j]/A + c + delta_1[j] + delta_2[j] 
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
       SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    )
    

    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == (n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1])*(1-SUTime/3600*Z_b[j-1])
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    # m.ext[:constraints][:con31c] = @constraint(m, [j=J[2:end]],
    #     T_a <= T[j]
    # )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    # m.ext[:constraints][:con37] = @constraint(m, [j=J],
    #     p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] + delta_2[j] 
    # )
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
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    m.ext[:constraints][:con207] = @constraint(m, [j=J[2:end]],
        i_b[j] + s_b[j-1] <= 1
    )
    return m
end
function build_model2f2!(m::Model)
    # power in 4 segments
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM2."Current density starts"


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j] 
    # )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
       SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    )
    

    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == (n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1])*(1-SUTime/3600*Z_b[j-1])
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin*s_b[j] <= T[j] 
    )
    m.ext[:constraints][:con31c] = @constraint(m, [j=J[2:end]],
        T_a <= T[j]
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j]
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
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    m.ext[:constraints][:con207] = @constraint(m, [j=J[2:end]],
        i_b[j] + s_b[j-1] <= 1
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
function build_model3f2!(m::Model)
    # power in 9 segments
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM3."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:9,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j] 
    # )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
       SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    )
    

    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == (n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1])*(1-SUTime/3600*Z_b[j-1])
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin*s_b[j] <= T[j] 
    )
    m.ext[:constraints][:con31c] = @constraint(m, [j=J[2:end]],
        T_a <= T[j]
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j] + delta_2[j]
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
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    m.ext[:constraints][:con207] = @constraint(m, [j=J[2:end]],
        i_b[j] + s_b[j-1] <= 1
    )
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:9)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*3+1:(k-1)*3+3) for k in 2:3) <= T[j]
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )
    return m
end
function build_model5f2!(m::Model)
    # Model with power in 4 segment
    # Nuclear heat used for:
        # 1) heating up water from 46.8°C to 100°C
        # 2) keeping temperature constant in standby
        # 3) No effect on cold start up cost because minimal operating temperature is 20°C!!!

    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM2."Current density starts"


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    delta_3 = m.ext[:variables][:delta_3] = @variable(m, [j=J], lower_bound=-M, base_name="delta 3") #slack variable 
    delta_4 = m.ext[:variables][:delta_4] = @variable(m, [j=J], lower_bound=-M, base_name="delta 4") #slack variable 
    delta_5 = m.ext[:variables][:delta_5] = @variable(m, [j=J], lower_bound=-M, base_name="delta 5")
    delta_6 = m.ext[:variables][:delta_6] = @variable(m, [j=J], lower_bound=-M, base_name="delta 6")
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j] 
    # )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]*η_turb + p_c[j] + Q_H2O[j]*η_turb + Q_cool[j]/400 
    )
    C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
       SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    )
    

    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == (n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1])*(1-SUTime/3600*Z_b[j-1])
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin*s_b[j] <= T[j] 
    )
    m.ext[:constraints][:con31c] = @constraint(m, [j=J[2:end]],
        T_a <= T[j]
    )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j]
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
        s_b[j-1] + p_b[j] - 1 <= Y_b[j]
    )
    m.ext[:constraints][:con42b] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= s_b[j-1]
    )
    m.ext[:constraints][:con43] = @constraint(m, [j=J[2:end]],
        Y_b[j] <= p_b[j]
    )
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    m.ext[:constraints][:con207] = @constraint(m, [j=J[2:end]],
        i_b[j] + s_b[j-1] <= 1
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


function build_model1f5!(m::Model)
    # power in 1 segment
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2PEM1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2PEM1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    # p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    p = m.ext[:expressions][:p] = @expression(m, [j=J],
         a*T[j]+ b*I[j]/A + c + delta_1[j]
    )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    # C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
    #    SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    # Formulate objective
    # m.ext[:objective] = @objective(m, Max,
    #     (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    # )
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)*(1-SUTime/3600*Z_b[j]) for j in J))
    )

    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == (n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1])*(1-SUTime/3600*Z_b[j-1])
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    # m.ext[:constraints][:con31c] = @constraint(m, [j=J[2:end]],
    #     T_a <= T[j]
    # )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    # m.ext[:constraints][:con37] = @constraint(m, [j=J],
    #     p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j] + delta_2[j] 
    # )
    #power kan bij j = 0 liggen tussen 483 en -923
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -1000*(s_b[j] + i_b[j]) <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= 2000*(s_b[j] + i_b[j])
    )

    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
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
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    m.ext[:constraints][:con207] = @constraint(m, [j=J[2:end]],
        i_b[j] + s_b[j-1] <= 1
    )
    return m
end
function build_model2f5!(m::Model)
    # power in 4 segments
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM2."Current density starts"


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j] 
    # )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    # C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
    #    SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    # Formulate objective
    # m.ext[:objective] = @objective(m, Max,
    #     (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    # )
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)*(1-SUTime/3600*Z_b[j]) for j in J))
    )
    

    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == (n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1])*(1-SUTime/3600*Z_b[j-1])
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    # for j = 0 power can vary between 1 and -294
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
       -100*(s_b[j] + i_b[j]) <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= 500*(s_b[j] + i_b[j])
    )
    
    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
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
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    m.ext[:constraints][:con207] = @constraint(m, [j=J[2:end]],
        i_b[j] + s_b[j-1] <= 1
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
function build_model3f5!(m::Model)
    # power in 9 segments
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM3."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:9,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j] 
    # )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    # C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
    #    SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    # Formulate objective
    # m.ext[:objective] = @objective(m, Max,
    #     (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    # )
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)*(1-SUTime/3600*Z_b[j]) for j in J))
    )
    

    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == (n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1])*(1-SUTime/3600*Z_b[j-1])
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j] 
    )
    # power tussen -9 en -173
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        0 <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= 300*(s_b[j]+i_b[j])
    )

    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
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
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    m.ext[:constraints][:con207] = @constraint(m, [j=J[2:end]],
        i_b[j] + s_b[j-1] <= 1
    )
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:9)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*3+1:(k-1)*3+3) for k in 2:3) <= T[j]
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )
    return m
end
# extra constraint toegevoegd voor standby boven 30gradne
# extra constraint i_b = 0
function build_model5f5!(m::Model)
    # Model with power in 4 segment
    # Nuclear heat used for:
        # 1) heating up water from 46.8°C to 100°C
        # 2) keeping temperature constant in standby
        # 3) No effect on cold start up cost because minimal operating temperature is 20°C!!!

    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]


    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM2."Current density starts"


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j] 
    # )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]*η_turb + p_c[j] + Q_H2O[j]*η_turb + Q_cool[j]/400 
    )
    # C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
    #    SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    # Formulate objective
    # m.ext[:objective] = @objective(m, Max,
    #     (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    # )
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)*(1-SUTime/3600*Z_b[j]) for j in J))
    )
    

    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == (n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1])*(1-SUTime/3600*Z_b[j-1])
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
        T[j] == T_0
    )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )




    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin<= T[j] 
    )
    # m.ext[:constraints][:con31b] = @constraint(m, [j=J],
    #    s_b[j] == 0
    # )





    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    # for j = 0 power can vary between 35 and -294
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
       -100*(s_b[j] + i_b[j]) <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= 500*(s_b[j] + i_b[j])
    )
    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
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
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    m.ext[:constraints][:con207] = @constraint(m, [j=J[2:end]],
        i_b[j] + s_b[j-1] <= 1
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



function build_model1g!(m::Model)
    # power in 1 segment
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]
   p_b0 = m.ext[:parameters]["p_b0"]
   i_b0 = m.ext[:parameters]["i_b0"]
   s_b0 = m.ext[:parameters]["s_b0"]



    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM1."Coefficients constant [W]"[1]
    a = m.ext[:parameters]["a"] = coeffs2PEM1."Coefficients T [W/K]"[1]
    b = m.ext[:parameters]["b"] = coeffs2PEM1."Coefficients j [W/(A/m^2)]"[1]


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    # Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #      a*T[j]+ b*I[j]/A + c + delta_1[j]
    # )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    # C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
    #    SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
       SUTime/3600*(3600*0.0291*π_H*delta_t - 5903400*π_e[j]*delta_t)
    )
    Q_CS = m.ext[:expressions][:C_HS] = @expression(m, [j=J],
        SUTime/3600*(n_c*668.717-n_c*U_tn*1500*A) # zie matlabbestand hotstartPEM voor getallen
    )# assumptie opstart op 90°C want hot and minimal power (j=1000 A/m²)

    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] - Q_CS[j]*Z_b[j])*(3600*delta_t)/ C_h
    )


    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    )

    # m.ext[:objective] = @objective(m, Max,
    #     (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)*(1-SUTime/3600*Z_b[j]) for j in J))
    # )

    
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con1] = @constraint(m, [j=1],
       p_b[j] == p_b0
    )
    m.ext[:constraints][:con2] = @constraint(m, [j=1],
       i_b[j] == i_b0
    )
    m.ext[:constraints][:con3] = @constraint(m, [j=1],
       s_b[j] == s_b0
    )
    m.ext[:constraints][:con5] = @constraint(m, [j=J[end]],
       T25[j] <= T_s
    )
    m.ext[:constraints][:con6] = @constraint(m, [j=J[end]],
       T_a <= T25[j]
    )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] - Q_CS[j-1]*Z_b[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    # m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    # m.ext[:constraints][:con31c] = @constraint(m, [j=J[2:end]],
    #     T_a <= T[j]
    # )
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == a*T[j]+ b*I[j]/A + c + delta_1[j]
    )
    #power kan bij j = 0 liggen tussen 483 en -923
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -2000*(s_b[j] + i_b[j]) <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= 4000*(s_b[j] + i_b[j])
    )

    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
    )
    m.ext[:constraints][:con40] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + i_b[j] <= 1
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
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
   
    return m
end
function build_model2g!(m::Model)
    # power in 4 segments
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]
   p_b0 = m.ext[:parameters]["p_b0"]
   i_b0 = m.ext[:parameters]["i_b0"]
   s_b0 = m.ext[:parameters]["s_b0"]



    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM2."Current density starts"


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j] 
    # )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
       SUTime/3600*(3600*0.0291*π_H*delta_t - 5903400*π_e[j]*delta_t)
    )
    Q_CS = m.ext[:expressions][:C_HS] = @expression(m, [j=J],
        SUTime/3600*(n_c*668.717-n_c*U_tn*1500*A) # zie matlabbestand hotstartPEM voor getallen
    )# assumptie opstart op 90°C want hot and minimal power (j=1000 A/m²)

    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] - Q_CS[j]*Z_b[j])*(3600*delta_t)/ C_h
    )



    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    )

    # m.ext[:objective] = @objective(m, Max,
    #     (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)*(1-SUTime/3600*Z_b[j]) for j in J))
    # )
    

    m.ext[:constraints][:con1] = @constraint(m, [j=1],
       p_b[j] == p_b0
    )
    m.ext[:constraints][:con2] = @constraint(m, [j=1],
       i_b[j] == i_b0
    )
    m.ext[:constraints][:con3] = @constraint(m, [j=1],
       s_b[j] == s_b0
    )
    m.ext[:constraints][:con5] = @constraint(m, [j=J[end]],
        T25[j] <= T_s
    )
    m.ext[:constraints][:con6] = @constraint(m, [j=J[end]],
       T_a <= T25[j]
    )
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] - Q_CS[j-1]*Z_b[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    # m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    # for j = 0 power can vary between 38 and -299
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
       -20000*(s_b[j] + i_b[j]) <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= 20000*(s_b[j] + i_b[j])
    )
    
    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
    )
    m.ext[:constraints][:con40] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + i_b[j] <= 1
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
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
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
    return m
end
function build_model3g!(m::Model)
    # power in 9 segments
    # Electrical energy used to heat up water 53.2°C and keep temperature constant in standby
    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
    C_HS = m.ext[:parameters]["C_HSPEM"] 
    #    C_CS = m.ext[:parameters]["C_CSPEM"]  
   p_uMin = m.ext[:parameters]["p_uMin"]
   p_uMax = m.ext[:parameters]["p_uMax"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]
   p_b0 = m.ext[:parameters]["p_b0"]
   i_b0 = m.ext[:parameters]["i_b0"]
   s_b0 = m.ext[:parameters]["s_b0"]



    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM3."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM3."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM3."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM3."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM3."Current density starts"

    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    delta_2 = m.ext[:variables][:delta_2] = @variable(m, [j=J], lower_bound=-M, base_name="delta 2")
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Y_b = m.ext[:variables][:Y_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for hot start")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:9,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j] 
    # )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    )
    # C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
    #    SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
       SUTime/3600*(3600*0.0291*π_H*delta_t - 5903400*π_e[j]*delta_t)
    )
    Q_CS = m.ext[:expressions][:C_HS] = @expression(m, [j=J],
        SUTime/3600*(n_c*668.717-n_c*U_tn*1500*A) # zie matlabbestand hotstartPEM voor getallen
    )# assumptie opstart op 90°C want hot and minimal power (j=1000 A/m²)

    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] - Q_CS[j]*Z_b[j])*(3600*delta_t)/ C_h
    )



    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    )

    # m.ext[:objective] = @objective(m, Max,
    #     (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)*(1-SUTime/3600*Z_b[j]) for j in J))
    # )
    

    m.ext[:constraints][:con1] = @constraint(m, [j=1],
       p_b[j] == p_b0
    )
    m.ext[:constraints][:con2] = @constraint(m, [j=1],
       i_b[j] == i_b0
    )
    m.ext[:constraints][:con3] = @constraint(m, [j=1],
       s_b[j] == s_b0
    )
    m.ext[:constraints][:con5] = @constraint(m, [j=J[end]],
       T25[j] <= T_s
    )
    m.ext[:constraints][:con6] = @constraint(m, [j=J[end]],
       T_a <= T25[j]
    )
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] - Q_CS[j-1]*Z_b[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    # m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )
    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin <= T[j] 
    )
    
    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:9) + delta_1[j] 
    )
    # power tussen -9 en -173
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
        -100*(s_b[j]+i_b[j])<= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= 300*(s_b[j]+i_b[j])
    )

    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
    )
    m.ext[:constraints][:con40] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + i_b[j] <= 1
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
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
    # )
    m.ext[:constraints][:con61] = @constraint(m, [j=J],
    1 == sum(t_b[i, j] for i in 1:9)
    )
    m.ext[:constraints][:con62] = @constraint(m, [j=J],
    sum(T_starts[k]*sum(t_b[i, j] for i in (k-1)*3+1:(k-1)*3+3) for k in 2:3) <= T[j]
    # T_starts[2]*sum(t_b[i, j] for i in 5:8) + T_starts[3]*sum(t_b[i, j] for i in 9:12) + T_starts[4]*sum(t_b[i, j] for i in 13:16)   <= T[j] 
    )
    m.ext[:constraints][:con62b] = @constraint(m, [j=J],
    T[j] <= sum(T_starts[k]*sum(t_b[i, j] for i in (k-2)*3+1:(k-2)*3+3) for k in 2:3)+ T_s*sum(t_b[i, j] for i in 7:9)
    # T[j] <= T_starts[2]*sum(t_b[i, j] for i in 1:4) + T_starts[3]*sum(t_b[i, j] for i in 5:8) + T_starts[4]*sum(t_b[i, j] for i in 9:12) + T_s*sum(t_b[i, j] for i in 13:16)
    )

    m.ext[:constraints][:con63] = @constraint(m, [j=J],
    sum(j_starts[k]*sum(t_b[k+3*i,j] for i in 0:2) for k in 2:3) <= I[j]/A 
    # j_starts[2]*(t_b[2,j] + t_b[6,j] + t_b[10,j] + t_b[14,j]) + j_starts[3]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + j_starts[4]*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])   <= I[j]/A 
    )
    m.ext[:constraints][:con63b] = @constraint(m, [j=J],
    I[j]/A <= sum(j_starts[k]*sum(t_b[k-1+3*i,j] for i in 0:2) for k in 2:3) + iMax*(sum(t_b[3*i,j] for i in 1:3))

    # I[j]/A <= j_starts[2]*(t_1b[j]+t_5b[j] + t_9b[j] + t_13b[j]) + j_starts[3]*(t_2b[j] + t_6b[j] + t_10b[j] + t_14b[j]) + j_starts[4]*(t_3b[j] + t_7b[j] + t_11b[j] + t_15b[j]) + iMax*(t_4b[j] + t_8b[j] + t_12b[j] + t_16b[j])
    )
    return m
end
function build_model5g!(m::Model)
    # Model with power in 4 segment
    # Nuclear heat used for:
        # 1) heating up water from 46.8°C to 100°C
        # 2) keeping temperature constant in standby
        # 3) No effect on cold start up cost because minimal operating temperature is 20°C!!!

    # With new objective
    # with new C_CS expression
    # With p_s variable!!
    # With cyclic boundary conditions

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
   iMin = m.ext[:parameters]["iMinPEM"] 
   iMax = m.ext[:parameters]["iMaxPEM"] 
   A = m.ext[:parameters]["A"]   
   n_c = m.ext[:parameters]["n_cPEM"] 
   M_H2 = m.ext[:parameters]["M_H2"]
   F = m.ext[:parameters]["F"]
   U_tn = m.ext[:parameters]["U_tnPEM"]
   C_h = m.ext[:parameters]["C_hPEM"] 
   T_s = m.ext[:parameters]["T_sPEM"] 
   R_t = m.ext[:parameters]["R_tPEM"] 
   T_a = m.ext[:parameters]["T_a"]   
   M = m.ext[:parameters]["M"]
   alfa = m.ext[:parameters]["alfa"]
   T_0 = m.ext[:parameters]["T_0PEM"]
   p_N = m.ext[:parameters]["p_N"]
   cp_H2O= m.ext[:parameters]["cp_H2O"]
   M_H2O= m.ext[:parameters]["M_H2O"]
   M_H2= m.ext[:parameters]["M_H2"]
   η_EH = m.ext[:parameters]["η_EH"]
   η_turb = m.ext[:parameters]["η_turb"]
   Tmin = m.ext[:parameters]["TminPEM"]
   SUTime = m.ext[:parameters]["SUTimePEM"]
   p_b0 = m.ext[:parameters]["p_b0"]
   i_b0 = m.ext[:parameters]["i_b0"]
   s_b0 = m.ext[:parameters]["s_b0"]



    # Create parameters
    # Created here and not in seperate function because these are specific to each model
    c = m.ext[:parameters]["c"] = coeffs2PEM2."Coefficients constant [W]"
    a = m.ext[:parameters]["a"] = coeffs2PEM2."Coefficients T [W/K]"
    b = m.ext[:parameters]["b"] = coeffs2PEM2."Coefficients j [W/(A/m^2)]"
    T_starts = m.ext[:parameters]["T_starts"] = coeffs2PEM2."Temperature starts"
    j_starts = m.ext[:parameters]["j_starts"] = coeffs2PEM2."Current density starts"


    # Create variables
    delta_1 = m.ext[:variables][:delta_1] = @variable(m, [j=J], lower_bound=-M, base_name="delta 1") #slack variable 
    Q_cool = m.ext[:variables][:Q_cool] =  @variable(m, [j=J], lower_bound=0, base_name="Cooling power") #W
    I = m.ext[:variables][:I] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser cell current") #A
    i_b = m.ext[:variables][:i_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for off state") 
    p = m.ext[:variables][:p] =  @variable(m, [j=J], lower_bound=0, base_name="Electrolyser power") #W
    p_b = m.ext[:variables][:p_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for production state") #W
    # p_c = m.ext[:variables][:p_c] =  @variable(m, [j=J], lower_bound=0, base_name="Compressor power")
    s_b = m.ext[:variables][:s_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for standby state")
    Z_b = m.ext[:variables][:Z_b] =  @variable(m, [j=J], binary=true, base_name="Binary variable for cold start")
    p_s = m.ext[:variables][:p_s] =  @variable(m, [j=J], lower_bound = 0, base_name="Power consumption in stand by state")
    # Q_H2O = m.ext[:variables][:Q_H2O] = @variable(m,[j=J],lower_bound = 0, base_name="Heat to water")
    # η_f_product_I = m.ext[:variables][:η_f_product_I] = @variable(m, [j=J], lower_bound = 0, base_name="Product Faraday efficiency and current")
    T = m.ext[:variables][:T] =  @variable(m, [j=J], lower_bound = 0, base_name="Temperature")
    # p_E = m.ext[:variables][:p_E] =  @variable(m, [j=J], lower_bound = 0, base_name="Electrical power")
    # Q_CS = m.ext[:variables][:Q_CS] =  @variable(m, [j=J],lower_bound = 0, base_name="Heat cold start")
    # C_CS = m.ext[:variables][:C_CS] =  @variable(m, [j=J],lower_bound = -M, base_name="Cost hot start")
    t_b = m.ext[:variables][:t_b] =  @variable(m, [i = 1:4,j=J], binary=true, base_name="Binary variable for segments") #W

    # Create affine expressions (= linear combinations of variables)
    mdotH2 = m.ext[:expressions][:mdotH2] = @expression(m, [j=J], 
        M_H2*n_c*1*I[j]/(2*F)
    )
    mdotH2O = m.ext[:expressions][:mdotH2O] = @expression(m, [j=J], 
        1/M_H2*1*M_H2O*mdotH2[j]
    )
    p_c = m.ext[:expressions][:p_c] = @expression(m, [j=J],
        alfa*mdotH2[j]
    )
    Q_H2O = m.ext[:expressions][:Q_H2O] = @expression(m, [j=J], 
        mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    )    
    # p = m.ext[:expressions][:p] = @expression(m, [j=J],
    #     sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] + delta_2[j] 
    # )
    p_E = m.ext[:expressions][:p_E] = @expression(m, [j=J],
        n_c*p[j]+n_c*p_s[j]*η_turb + p_c[j] + Q_H2O[j]*η_turb + Q_cool[j]/400 
    )
    C_CS = m.ext[:expressions][:C_CS] = @expression(m, [j=J],
       SUTime/3600*(3600*0.0291*π_H*delta_t - 5903400*π_e[j]*delta_t)
    )
    Q_CS = m.ext[:expressions][:C_HS] = @expression(m, [j=J],
        SUTime/3600*(n_c*668.717-n_c*U_tn*1500*A) # zie matlabbestand hotstartPEM voor getallen
    )# assumptie opstart op 90°C want hot and minimal power (j=1000 A/m²)

    p_bSUM = m.ext[:expressions][:p_bSUM] = @expression(m, 
        sum(p_b[j] for j in J)
    )
    T25 = m.ext[:expressions][:T25] = @expression(m, [j=J[end]],
        T[j] + (n_c*p[j] + n_c*p_s[j] - n_c*U_tn*I[j] - (T[j]-T_a)/R_t - Q_cool[j] - Q_CS[j]*Z_b[j])*(3600*delta_t)/ C_h
    )



    # Formulate objective
    m.ext[:objective] = @objective(m, Max,
        (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t) for j in J)) - sum(C_CS[j]*Z_b[j] for j in J[2:end])
    )
    # m.ext[:objective] = @objective(m, Max,
    #     (sum((3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)*(1-SUTime/3600*Z_b[j]) for j in J))
    # )
    

    m.ext[:constraints][:con1] = @constraint(m, [j=1],
       p_b[j] == p_b0
    )
    m.ext[:constraints][:con2] = @constraint(m, [j=1],
       i_b[j] == i_b0
    )
    m.ext[:constraints][:con3] = @constraint(m, [j=1],
       s_b[j] == s_b0
    )
    m.ext[:constraints][:con5] = @constraint(m, [j=J[end]],
       T25[j] <= T_s
    )
    m.ext[:constraints][:con6] = @constraint(m, [j=J[end]],
       T_a <= T25[j]
    )
    # Formulate constraints
    # m.ext[:constraints][:con29] = @constraint(m, [j=J],
    #     p_E[j] == n_c*p[j]+n_c*p_s[j]/η_EH + p_c[j] + Q_H2O[j]/η_EH + Q_cool[j]/400 
    # )
    m.ext[:constraints][:con30] = @constraint(m, [j=J[2:end]],
        C_h*(T[j] -T[j-1])/(3600*delta_t) == n_c*p[j-1] +n_c*p_s[j-1] - n_c*U_tn*I[j-1] - (T[j-1]-T_a)/R_t - Q_cool[j-1] - Q_CS[j-1]*Z_b[j-1]
    )
    m.ext[:constraints][:con30b] = @constraint(m, [j=1],
        T[j] == T_0
    )
    # m.ext[:constraints][:con30c] = @constraint(m, [j=J[end]],
    #     T[j] == T_0
    # )
    m.ext[:constraints][:con31] = @constraint(m, [j=J[2:end]],
        T[j] <= T_s
    )



    m.ext[:constraints][:con31b] = @constraint(m, [j=J[2:end]],
        Tmin<= T[j] 
    )
    # m.ext[:constraints][:con31b] = @constraint(m, [j=J],
    #    s_b[j] == 0
    # )

    m.ext[:constraints][:con32a] = @constraint(m, [j=J],
        A*p_b[j]*iMin  <= I[j]
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
        p_b[j] + s_b[j] + i_b[j] == 1
    )
    m.ext[:constraints][:con37] = @constraint(m, [j=J],
        p[j] == sum(t_b[i,j]*(a[i]*T[j] + b[i]*I[j]/A + c[i]) for i in 1:4) + delta_1[j] 
    )
    # for j = 0 power can vary between 35 and -294
    m.ext[:constraints][:con38a] = @constraint(m, [j=J],
       -20000*(s_b[j] + i_b[j]) <= delta_1[j]
    )
    m.ext[:constraints][:con38b] = @constraint(m, [j=J],
        delta_1[j] <= 20000*(s_b[j] + i_b[j])
    )
    m.ext[:constraints][:con39] = @constraint(m, [j=J[2:end]],
        s_b[j] + i_b[j-1] <= 1
    )
    m.ext[:constraints][:con40] = @constraint(m, [j=J[2:end]],
        s_b[j-1] + i_b[j] <= 1
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
    m.ext[:constraints][:con44a] = @constraint(m, [j=J[2:end]],
        i_b[j-1] + p_b[j] - 1 <= Z_b[j]
    )
    m.ext[:constraints][:con44b] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= i_b[j-1]
    )
    m.ext[:constraints][:con45] = @constraint(m, [j=J[2:end]],
        Z_b[j] <= p_b[j]
    )
    # Constraint (46) is about ramping rates and can be relieved


    # m.ext[:constraints][:con52] = @constraint(m, [j=J],
    #     p_c[j] == alfa*mdotH2[j]
    # )
    # m.ext[:constraints][:con54] = @constraint(m, [j=J], 
    #     Q_H2O[j] == mdotH2O[j]*cp_H2O*53.2 # benadering om het lineair te houden!!! zie notities!! 
    # )    
    #The following constraint is added to prevent cycling: on off on off, bcs at high temperatures
    #  during a long period where you want to be in off state, 
    #  it does not matter if you perform a cold start up or not because you can choose not to add heat
    #  and the hour without production also doesn't matter if you want to stay 
    #  in the standby the next couple of hours
    m.ext[:constraints][:con204] = @constraint(m, [j=J[3:end]], 
       i_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    m.ext[:constraints][:con205] = @constraint(m, [j=J[3:end]], 
       s_b[j] + p_b[j-1] + i_b[j-2] <= 2 
    )
    # m.ext[:constraints][:con206] = @constraint(m, [j=J],
    #    C_CS[j] == SUTime/3600*(3600*mdotH2[j]*π_H*delta_t - p_E[j]*π_e[j]*delta_t)
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
    return m
end

objective_values = Matrix{Float64}(undef, 367*24, 12)
define_sets!(m, data,24)
process_parameters!(m, data)
for j in 1:1
    operating_hours = 0
    average_eff = 0
    m.ext[:parameters]["π_H"] = 3.4+ (j - 1) 
    T_end = 293
    p_b0 = 1
    i_b0 = 0
    s_b0 = 0
    for i in 5:5
        if T_end < 293
            m.ext[:parameters]["T_0PEM"] = T_end + 10^-4
        elseif T_end > 373
            m.ext[:parameters]["T_0PEM"] = T_end - 10^-4
        else 
            m.ext[:parameters]["T_0PEM"] = T_end
        end
        m.ext[:parameters]["C_hPEM"] = m.ext[:parameters]["C_hPEM"]
        m.ext[:parameters]["R_tPEM"] = m.ext[:parameters]["R_tPEM"]
        m.ext[:parameters]["p_b0"] = p_b0
        m.ext[:parameters]["i_b0"] = i_b0
        m.ext[:parameters]["s_b0"] = s_b0
        process_time_series_data!(m, data, i)
        build_model1g!(m)
        # set_optimizer_attribute(m, "mip_gap", 0.001)
        optimize!(m)
        objective_values[i, j] = value.(m.ext[:objective])
        operating_hours = operating_hours + value.(m.ext[:expressions][:p_bSUM]) 
        objective_values[(i-1)*24+1:i*24, 4+j] = [value.(m.ext[:expressions][:mdotH2])[k]*value.(m.ext[:variables][:p_b])[k] for k in 1:24]
        objective_values[(i-1)*24+1:i*24, 8+j] = [value.(m.ext[:expressions][:p_E])[k] for k in 1:24]
        T_end = value.(first(m.ext[:expressions][:T25])) #geen constraint op T[24]!!
        p_b0 = round(value.(m.ext[:variables][:p_b][24]))
        i_b0 = round(value.(m.ext[:variables][:i_b][24]))
        s_b0 = round(value.(m.ext[:variables][:s_b][24]))
        # # check termination status
        print(
            """
            hprice = $j
            dag = $i
            T_end = $T_end
            p_b0 = $p_b0
            i_b0 = $i_b0
            s_b0 = $s_b0
            """
        )
    end
    objective_values[366, j] = operating_hours
    objective_values[367, j] = average_eff
end
matfile = matopen("objective_values_PEM.mat", "w")
header = "objective_values"
write(matfile, header, objective_values)
close(matfile)


# objective_values = Matrix{Float64}(undef, 367*24, 12)
# define_sets!(m, data,24)
# process_parameters!(m, data)
# for j in 1:1
#     operating_hours = 0
#     average_eff = 0
#     m.ext[:parameters]["π_H"] =2.5+ (j - 1) 
#     T_end = 293
#     p_b0 = 1
#     i_b0 = 0
#     s_b0 = 0
#     for i in 1:365
#         if T_end < 293
#             m.ext[:parameters]["T_0PEM"] = T_end + 10^-4
#         elseif T_end > 373
#             m.ext[:parameters]["T_0PEM"] = T_end - 10^-4
#         else 
#             m.ext[:parameters]["T_0PEM"] = T_end
#         end
#         m.ext[:parameters]["T_0PEM"] = T_end
#         m.ext[:parameters]["C_hPEM"] = m.ext[:parameters]["C_hPEM"]/2
#         m.ext[:parameters]["R_tPEM"] = m.ext[:parameters]["R_tPEM"]
#         m.ext[:parameters]["p_b0"] = p_b0
#         m.ext[:parameters]["i_b0"] = i_b0
#         m.ext[:parameters]["s_b0"] = s_b0
#         process_time_series_data!(m, data, i)
#         build_model2g!(m)
#         # set_optimizer_attribute(m, "mip_gap", 0.001)
#         optimize!(m)
#         objective_values[i, j] = value.(m.ext[:objective])
#         operating_hours = operating_hours + value.(m.ext[:expressions][:p_bSUM]) 
#         objective_values[(i-1)*24+1:i*24, 4+j] = [value.(m.ext[:expressions][:mdotH2])[k]*value.(m.ext[:variables][:p_b])[k] for k in 1:24]
#         objective_values[(i-1)*24+1:i*24, 8+j] = [value.(m.ext[:expressions][:p_E])[k] for k in 1:24]
#         T_end = value.(first(m.ext[:expressions][:T25])) #geen constraint op T[24]!!
#         p_b0 = round(value.(m.ext[:variables][:p_b][24]))
#         i_b0 = round(value.(m.ext[:variables][:i_b][24]))
#         s_b0 = round(value.(m.ext[:variables][:s_b][24]))
#         # # check termination status
#         print(
#             """
#             hprice = $j
#             dag = $i
#             T_end = $T_end
#             p_b0 = $p_b0
#             i_b0 = $i_b0
#             s_b0 = $s_b0
#             """
#         )
#     end
#     objective_values[366, j] = operating_hours
#     objective_values[367, j] = average_eff
# end
# matfile = matopen("objective_values_PEM6.mat", "w")
# header = "objective_values"
# write(matfile, header, objective_values)
# close(matfile)


# objective_values = Matrix{Float64}(undef, 367*24, 12)
# define_sets!(m, data,24)
# process_parameters!(m, data)
# for j in 1:1
#     operating_hours = 0
#     average_eff = 0
#     m.ext[:parameters]["π_H"] =2.5+ (j - 1) 
#     T_end = 293
#     p_b0 = 1
#     i_b0 = 0
#     s_b0 = 0
#     for i in 1:365
#         if T_end < 293
#             m.ext[:parameters]["T_0PEM"] = T_end + 10^-4
#         elseif T_end > 373
#             m.ext[:parameters]["T_0PEM"] = T_end - 10^-4
#         else 
#             m.ext[:parameters]["T_0PEM"] = T_end
#         end
#         m.ext[:parameters]["T_0PEM"] = T_end
#         m.ext[:parameters]["C_hPEM"] = m.ext[:parameters]["C_hPEM"]*2
#         m.ext[:parameters]["R_tPEM"] = m.ext[:parameters]["R_tPEM"]
#         m.ext[:parameters]["p_b0"] = p_b0
#         m.ext[:parameters]["i_b0"] = i_b0
#         m.ext[:parameters]["s_b0"] = s_b0
#         process_time_series_data!(m, data, i)
#         build_model2g!(m)
#         # set_optimizer_attribute(m, "mip_gap", 0.001)
#         optimize!(m)
#         objective_values[i, j] = value.(m.ext[:objective])
#         operating_hours = operating_hours + value.(m.ext[:expressions][:p_bSUM]) 
#         objective_values[(i-1)*24+1:i*24, 4+j] = [value.(m.ext[:expressions][:mdotH2])[k]*value.(m.ext[:variables][:p_b])[k] for k in 1:24]
#         objective_values[(i-1)*24+1:i*24, 8+j] = [value.(m.ext[:expressions][:p_E])[k] for k in 1:24]
#         T_end = value.(first(m.ext[:expressions][:T25])) #geen constraint op T[24]!!
#         p_b0 = round(value.(m.ext[:variables][:p_b][24]))
#         i_b0 = round(value.(m.ext[:variables][:i_b][24]))
#         s_b0 = round(value.(m.ext[:variables][:s_b][24]))
#         # # check termination status
#         print(
#             """
#             hprice = $j
#             dag = $i
#             T_end = $T_end
#             p_b0 = $p_b0
#             i_b0 = $i_b0
#             s_b0 = $s_b0
#             """
#         )
#     end
#     objective_values[366, j] = operating_hours
#     objective_values[367, j] = average_eff
# end
# matfile = matopen("objective_values_PEM7.mat", "w")
# header = "objective_values"
# write(matfile, header, objective_values)
# close(matfile)


# objective_values = Matrix{Float64}(undef, 367*24, 12)
# define_sets!(m, data,24)
# process_parameters!(m, data)
# for j in 1:1
#     operating_hours = 0
#     average_eff = 0
#     m.ext[:parameters]["π_H"] =2.5+ (j - 1) 
#     T_end = 293
#     p_b0 = 1
#     i_b0 = 0
#     s_b0 = 0
#     for i in 1:365
#         if T_end < 293
#             m.ext[:parameters]["T_0PEM"] = T_end + 10^-4
#         elseif T_end > 373
#             m.ext[:parameters]["T_0PEM"] = T_end - 10^-4
#         else 
#             m.ext[:parameters]["T_0PEM"] = T_end
#         end
#         m.ext[:parameters]["T_0PEM"] = T_end
#         m.ext[:parameters]["C_hPEM"] = m.ext[:parameters]["C_hPEM"]
#         m.ext[:parameters]["R_tPEM"] = m.ext[:parameters]["R_tPEM"]/2
#         m.ext[:parameters]["p_b0"] = p_b0
#         m.ext[:parameters]["i_b0"] = i_b0
#         m.ext[:parameters]["s_b0"] = s_b0
#         process_time_series_data!(m, data, i)
#         build_model2g!(m)
#         # set_optimizer_attribute(m, "mip_gap", 0.001)
#         optimize!(m)
#         objective_values[i, j] = value.(m.ext[:objective])
#         operating_hours = operating_hours + value.(m.ext[:expressions][:p_bSUM]) 
#         objective_values[(i-1)*24+1:i*24, 4+j] = [value.(m.ext[:expressions][:mdotH2])[k]*value.(m.ext[:variables][:p_b])[k] for k in 1:24]
#         objective_values[(i-1)*24+1:i*24, 8+j] = [value.(m.ext[:expressions][:p_E])[k] for k in 1:24]
#         T_end = value.(first(m.ext[:expressions][:T25])) #geen constraint op T[24]!!
#         p_b0 = round(value.(m.ext[:variables][:p_b][24]))
#         i_b0 = round(value.(m.ext[:variables][:i_b][24]))
#         s_b0 = round(value.(m.ext[:variables][:s_b][24]))
#         # # check termination status
#         print(
#             """
#             hprice = $j
#             dag = $i
#             T_end = $T_end
#             p_b0 = $p_b0
#             i_b0 = $i_b0
#             s_b0 = $s_b0
#             """
#         )
#     end
#     objective_values[366, j] = operating_hours
#     objective_values[367, j] = average_eff
# end
# matfile = matopen("objective_values_PEM8.mat", "w")
# header = "objective_values"
# write(matfile, header, objective_values)
# close(matfile)



# objective_values = Matrix{Float64}(undef, 367*24, 12)
# define_sets!(m, data,24)
# process_parameters!(m, data)
# for j in 1:1
#     operating_hours = 0
#     average_eff = 0
#     m.ext[:parameters]["π_H"] =2.5+ (j - 1) 
#     T_end = 293
#     p_b0 = 1
#     i_b0 = 0
#     s_b0 = 0
#     for i in 1:365
#         if T_end < 293
#             m.ext[:parameters]["T_0PEM"] = T_end + 10^-4
#         elseif T_end > 373
#             m.ext[:parameters]["T_0PEM"] = T_end - 10^-4
#         else 
#             m.ext[:parameters]["T_0PEM"] = T_end
#         end
#         m.ext[:parameters]["T_0PEM"] = T_end
#         m.ext[:parameters]["C_hPEM"] = m.ext[:parameters]["C_hPEM"]
#         m.ext[:parameters]["R_tPEM"] = m.ext[:parameters]["R_tPEM"]*2
#         m.ext[:parameters]["p_b0"] = p_b0
#         m.ext[:parameters]["i_b0"] = i_b0
#         m.ext[:parameters]["s_b0"] = s_b0
#         process_time_series_data!(m, data, i)
#         build_model2g!(m)
#         # set_optimizer_attribute(m, "mip_gap", 0.001)
#         optimize!(m)
#         objective_values[i, j] = value.(m.ext[:objective])
#         operating_hours = operating_hours + value.(m.ext[:expressions][:p_bSUM]) 
#         objective_values[(i-1)*24+1:i*24, 4+j] = [value.(m.ext[:expressions][:mdotH2])[k]*value.(m.ext[:variables][:p_b])[k] for k in 1:24]
#         objective_values[(i-1)*24+1:i*24, 8+j] = [value.(m.ext[:expressions][:p_E])[k] for k in 1:24]
#         T_end = value.(first(m.ext[:expressions][:T25])) #geen constraint op T[24]!!
#         p_b0 = round(value.(m.ext[:variables][:p_b][24]))
#         i_b0 = round(value.(m.ext[:variables][:i_b][24]))
#         s_b0 = round(value.(m.ext[:variables][:s_b][24]))
#         # # check termination status
#         print( 
#             """
#             hprice = $j
#             dag = $i
#             T_end = $T_end
#             p_b0 = $p_b0
#             i_b0 = $i_b0
#             s_b0 = $s_b0
#             """
#         )
#     end
#     objective_values[366, j] = operating_hours
#     objective_values[367, j] = average_eff
# end
# matfile = matopen("objective_values_PEM9.mat", "w")
# header = "objective_values"
# write(matfile, header, objective_values)
# close(matfile)

# objective_values = Matrix{Float64}(undef, 367*24, 12)
# define_sets!(m, data,24)
# process_parameters!(m, data)
# for j in 1:4
#     operating_hours = 0
#     average_eff = 0
#     m.ext[:parameters]["π_H"] = 2.5 + (j - 1)
#     for i in 1:365
#         process_time_series_data!(m, data, i)
#         build_model5f2!(m)
#         optimize!(m)
#         objective_values[i, j] = value.(m.ext[:objective])
#         operating_hours = operating_hours + value.(m.ext[:expressions][:p_bSUM]) 
#         objective_values[(i-1)*24+1:i*24, 4+j] = [value.(m.ext[:expressions][:mdotH2])[k]*value.(m.ext[:variables][:p_b])[k] for k in 1:24]
#         objective_values[(i-1)*24+1:i*24, 8+j] = [value.(m.ext[:expressions][:p_E])[k] for k in 1:24]
#     end
#     objective_values[366, j] = operating_hours
#     objective_values[367, j] = average_eff
# end
# matfile = matopen("objective_values_PEM.mat", "w")
# header = "objective_values"
# write(matfile, header, objective_values)
# close(matfile)


# function process_time_series_data!(m::Model, data::Dict)
#     # extract the relevant sets
#     #J = m.ext[:sets][:J] # Time steps

#     # create dictionary to store time series
#     m.ext[:timeseries] = Dict()

#     # example: add time series to dictionary
#     m.ext[:timeseries][:π_e] = ePrice."Day-ahead Price [EUR/MWh]"[101*24+1:102*24]/1000000 #[€/Wh]


#     # Both 0, very large storage tank (no problem bcs in case 1 in paper storage never full so no determening constraint)
#     # p_l = 0 so just more p_u (no power to local load but everything to utility grid)
#     m.ext[:timeseries][:m_LH] = 0
#     m.ext[:timeseries][:p_l] = 0 
#     # return model
#     return m
# end
# process_time_series_data!(m, data)
# define_sets!(m, data,24)
# process_parameters!(m, data)
# m.ext[:parameters]["T_0PEM"] = 293
# m.ext[:parameters]["p_b0"] = 1
# m.ext[:parameters]["i_b0"] = 0
# m.ext[:parameters]["s_b0"] = 0
# build_model1g!(m)

# ## Step 4: solve
# optimize!(m)

# # check termination status
# print(
#     """

#     Termination status: $(termination_status(m))

#     """
# )


# print some output
@show value(m.ext[:objective])


## Step 5: interpretation
# access data

# sets
J = m.ext[:sets][:J]
# # parameters
eprice = value.(m.ext[:timeseries][:π_e])
n_c = value.(m.ext[:parameters]["n_cPEM"])
T_0 = value.(m.ext[:parameters]["T_0PEM"])-273 #°C
p_N = value.(m.ext[:parameters]["p_N"])
# variables/expressions
p = value.(m.ext[:variables][:p])
i_b = value.(m.ext[:variables][:i_b])
s_b = value.(m.ext[:variables][:s_b])
p_b = value.(m.ext[:variables][:p_b])
Z_b = value.(m.ext[:variables][:Z_b])
Q_cool = value.(m.ext[:variables][:Q_cool])
p_c = value.(m.ext[:expressions][:p_c])
I = value.(m.ext[:variables][:I])
delta_1 = value.(m.ext[:variables][:delta_1])
# delta_2 = value.(m.ext[:variables][:delta_2])

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
i_bvec = [round(i_b[j]) for j in J]
s_bvec = [round(s_b[j]) for j in J]
p_bvec = [round(p_b[j]) for j in J]
Z_bvec = [Z_b[j] for j in J]
Q_coolvec = [Q_cool[j] for j in J]
p_cvec = [p_c[j]/1000000 for j in J] 
Ivec = [I[j] for j in J] 
mdotH2vec = [mdotH2[j] for j in J] 

using Plots
using LaTeXStrings
using StatsPlots
using Plots.PlotMeasures


pyplot()
plot(pvec, label = "Electrolyser Power", xlabel="Time [h]",ylabel = "Electrolyser power [MW]", color = :red, 
        legend = :topleft, left_margin = 5Plots.mm, right_margin = 15Plots.mm)
plot!(twinx(), epricevec, label = "Electricity price", ylabel="Electricity Price [€/MWh]", legend = :topright)
savefig("power_price.png")

plot(J,p_Nvec, xlabel = "Time [h]", ylabel = "Wind Power [MW]",legend = false, left_margin = 5Plots.mm, right_margin = 15Plots.mm)
savefig("Windpower.png")

plot(epricevec, xlabel = "Time [h]", ylabel="Electricity Price [€/MWh]", legend=false)
savefig("eprice.png")

plot(J.+1/2, [Tvec], xlabel = "Time [h]", ylabel = "Electrolyser Temperature [°C]",legend = false, left_margin = 5Plots.mm, right_margin = 15Plots.mm)
savefig("Temperature.png")

bar(J.+1/2, [Tvec], xlabel = "Time [h]", ylabel = "Electrolyser Temperature [°C]", legend = false, left_margin = 5Plots.mm, right_margin = 15Plots.mm, bar_width = 1.0, fillalpha = 0.5, linealpha = 0.1)
savefig("Temperaturebars.png")


plot(J, [p_bvec s_bvec i_bvec], label = ["Production state"  "Standby state" "Off state"], xlabel = "Time [h]", legend = :right,left_margin = 5Plots.mm, right_margin = 15Plots.mm)
savefig("Binary variables.png")
bar(J.+ 1/2, [p_bvec s_bvec i_bvec], label = ["Production state"  "Standby state" "Off state"], xlabel = "Time [h]", legend = :right,left_margin = 5Plots.mm, right_margin = 15Plots.mm, bar_width = 1.0, fillalpha = 0.5, linealpha = 0.1)
savefig("Binary variablesbars.png")

plot(J, [Z_bvec ], label = ["Cold start" ], xlabel = "Time [h]", left_margin = 5Plots.mm, right_margin = 15Plots.mm )
savefig("Binary variables cold hot start.png")

plot(J, Q_coolvec, xlabel = "Time [h]", ylabel = "Cooling Power [W]",legend = false,left_margin = 5Plots.mm, right_margin = 15Plots.mm )
savefig("Cooling power.png")


# plot([p_uvec pvec p_Nvec p_cvec], label = ["Power to utility grid" "Electrolyser power" "Wind power" "Compressor power"], xlabel="Time [h]",ylabel = "Power [MW]", legend = :topleft, left_margin = 5Plots.mm, right_margin = 15Plots.mm )
# savefig("powers.png")

plot(J, Ivec, xlabel = "Time [h]", ylabel = "Electrolyser current [A]",left_margin = 5Plots.mm, right_margin = 15Plots.mm )
savefig("Electrolyser current.png")

plot(J, mdotH2vec, xlabel = "Time [h]", ylabel = "Hydrogen production rate [kg/s]",left_margin = 5Plots.mm, right_margin = 15Plots.mm )
savefig("Hydrogen production rate.png")

bar(J.-1/2, pvec, label = "Electrolyser Power", xlabel = "Time [h]", ylabel = "Electrolyser power [MW]",
    legend = :topleft, left_margin = 5Plots.mm, right_margin = 15Plots.mm, bar_width = 1.0, fillalpha = 0.5, linealpha = 0.1)

# Create a line plot for epricevec on the secondary y-axis
plot!(twinx(), J.-1/2, epricevec, label = "Electricity price", ylabel = "Electricity Price [€/MWh]",
      line = :path, legend = :topright, color = :red)

# Save the plot to a file
savefig("power_price_bars.png")

colours = ifelse.(p_bvec.==1, :green, (ifelse.(s_bvec .== 1, :orange, :red)))
states = ["Production state"  "Off state" "Standby state"]
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