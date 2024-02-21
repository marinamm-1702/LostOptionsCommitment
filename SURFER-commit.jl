#-----------------------------------------------------
# fixed CarbonCycle+Aerosol parameters
#-----------------------------------------------------
k_AU = 1/4              # 1/years
h_U = 150             # m. Thickness of upper ocean layer
h_D =  h_U*20         # m. Thickness of lower ocean layer
δ = h_D/h_U           # h_L/h_U dimensionless 
δDIC = 1.15           # dimensionless DIC_L/DIC_U at pre industrial
k_UD = (δ*δDIC)/1000    # 1/years
K0 = 0.0314806         # mol/(kg atm) Henry CO2 solubility
K1 = 1.32648E-6        # mol/kg
K2 = 9.19803E-10     # mol/kg
mA = 1.727E20         # moles in atmosphere
mO = 7.8E22           # moles in ocean
GtCtoppm(M_A) = (M_A * 1E6 * 1E15)/(12*mA) # ppm in atmosphere/ PgC in atmosphere
ppmtoGtC(conc) = conc*12*mA/(1E21)
W_U = mO*18E-3*h_U/(h_U+h_D) # whole upper ocean mass kg
W_D = mO*18E-3*h_D/(h_U+h_D) # whole lower ocean mass kg
pH_PI = 8.17
H_PI = 10^(-pH_PI) # mol/kg
CO2conc_a_PI = 280 # ppm
pCO2_a(CO2conc_a) = CO2conc_a*1E-6 # atm
Q = (K1/H_PI + 2*K1*K2/H_PI^2)*K0*pCO2_a(CO2conc_a_PI) # mol/kg
Qm = Q*W_U*12E-3*1E-12 # PgC
k_AL = 1/40            #1/yr
β_L = 1.7             # dimensionless [0.5 - 2.3]
#-----------------------------------------------------
# pre industrial initial conditions for the carbon reservoirs
#-----------------------------------------------------
M_A_PI = ppmtoGtC(CO2conc_a_PI) # PgC
M_U_PI = M_A_PI*(1 + K1/H_PI + K1*K2/H_PI^2)*W_U*K0/mA     # PgC
M_D_PI = M_U_PI*δ*δDIC    # PgC
M_L_PI = 2200 # PgC
#-----------------------------------------------------
# parameters entering temperature computation (greenhouse and heat transport parameters)
#-----------------------------------------------------
ECS = 3.5#2.0#5.0           # Equilibrium Climate Sensitivity. Celcius / doubling of CO2
TCR = 2.0#1.2#2.4           # Transient Climate Response. Celcius / doubling of CO2
F2X = 3.9            # Watts / m^2. Forcing due to a doubling of CO2 concentration
β = F2X/ECS          # Watts / (m^2 Celcius). Inverse equilibrium cliM_Ae sensitivity.
γ = F2X/TCR - β      # Watts / (m^2 Celcius). Thermal conductivity between layers
cvol = 0.13          # (Watts year) / (m^3 Celcius). Volumetric heat capacity of seawater.
# TCR = F2X/(β+γ)      # Celsius per doubling of CO2 concentration
#-----------------------------------------------------
# parameters entering temperature computation (aerosol parameters taken from Helwegen2019, changed η)
#-----------------------------------------------------
#η = 0.742         # dimensionless
η = 1.0
αSO2 = 65        # Watts / m^2
βSO2 = 2246      # Mt of S / year
γSO2 = 0.23;     # dimensionless

#----------------------------------------------------
# a couple of auxiliary functions 
#----------------------------------------------------

function B(M_U) #dimensionles
    (sqrt( K1*( Qm*(K1-4K2)*(Qm-2M_U) + K1*M_U^2 ) ) - ( Qm*(K1-4K2) + M_U*(-K1+8K2) ) )/(2M_U*(K1-4K2))
end

function F_CO2(M_A)
    return F2X*log2(M_A/M_A_PI)
end

function F_SO2(I)
    return - η*αSO2*exp(-(βSO2/I)^γSO2)
end

function F(M_A,I) # antropogenic forcing
    if I > 0.0
        return  F_CO2(M_A) + F_SO2(I)
    else
        return F_CO2(M_A)
    end
end

function SO2_needed_withD(δTeq,M_A,δT_D)
    return βSO2*(-log(max(0, ( γ*(δT_D - δTeq) - β*δTeq + F_CO2(M_A) )) / (αSO2) ))^(-1/γSO2)
end

function SO2_needed_withD_SDP(δT_U,M_A,δT_D; δT_goal = 1.5)
    if δT_U > δT_goal
        return βSO2*(-log(max(0, ( β*(δT_U - δT_goal)+ γ*(δT_D - δT_U) - β*δT_U + F_CO2(M_A) )) / (αSO2) ))^(-1/γSO2)
    else
        return 0
    end
end

function SO2_needed_old(δTeq,M_A)
    return βSO2*(-log(max(0, ( - β*δTeq + F_CO2(M_A) )) / (αSO2) ))^(-1/γSO2)
end

function CDR_needed(M_A,M_U,M_L)
     CDR_rate = + k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) + k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))
    return CDR_rate
end

function H(M_U) # Proton (hydrogen ion) concentration. mol/kg
    H = ( -Qm*K1 + K1*M_U + sqrt( K1*(Qm^2*(K1 - 4*K2) - 2*Qm*(K1 - 4*K2)*M_U + K1*M_U^2) ) )/(2*Qm)
    return H
end

function pH(M_U) # pH global scale
    return -log10(H(M_U))
end

function DIC(M_U) # μmol/kg
    return 1E6(M_U/12E-15)/(W_U)
end

function DIC_D(M_D) # μmol/kg
    return 1E6(M_D/12E-15)/(W_D)
end

function H2CO3(M_U) # μmol/kg
    H2CO3 = DIC(M_U)/(1 + K1/H(M_U) + K1*K2/H(M_U)^2)
    return H2CO3
end

function HCO3(M_U) # μmol/kg
    HCO3 = K1*H2CO3(M_U)/H(M_U)
    return HCO3
end

function CO3(M_U) # μmol/kg
    CO3 = K2*HCO3(M_U)/H(M_U)
    return CO3
end

CO3sat = CO3(M_U_PI)/3.44 

function Ω(M_U)
    return CO3(M_U)/CO3sat
end

Ω_PB = 0.8*Ω(M_U_PI)
#-----------------------------------------------------
# ICE parameters and functions
#-----------------------------------------------------

function Vmcons(model_parameters)
    
    Tp, Tm, Vp, Vm, τmelt, τfreeze = model_parameters
    
    x = ((-Tm + Tp)/(Tm + Tp + 2*sqrt(Tm*Tp)))^(1/3)
    Vm = ( -2 + Vp*(1 + x + 1/x ) )/( -1 + x + 1/x )
    
    return Vm
end

Greenland_params = [1.52, 0.3, 0.77, 0.3526554620064224, 470.0, 5000.0]#[Tp, Tm, Vp, Vm, tau_melt, tau_freeze] 
Greenland_params[4] = Vmcons(Greenland_params)
Antarctica_params = [6.8, 4, 0.44, 0.07857839308355193, 3000.0, 5500.0]#[Tp, Tm, Vp, Vm, tau_melt, tau_freeze] 
Antarctica_params[4] = Vmcons(Antarctica_params)
# sea level rise potential
SLRpotentialG = 7.4
SLRpotentialA = 55

function coeffs(model_parameters)
    Tp, Tm, Vp, Vm, τmelt, τfreeze = model_parameters
    
    a = 3*(Vm + Vp)/2
    b = -3*Vm*Vp
    c = (Vp - Vm)^3/(2*(Tm - Tp))
    d = ( Tp*Vm^2*(Vm-3Vp) - Tm*Vp^2*(Vp-3Vm) )/(2*(Tm - Tp))

    return a,b,c,d
end

Greenland_coeffs = coeffs(Greenland_params)
Antarctica_coeffs = coeffs(Antarctica_params)

function dV_dt_v2(V, Tf, model_parameters, coeffs)

    Tp, Tm, Vp, Vm, τmelt, τfreeze = model_parameters
    a,b,c,d = coeffs
    
    aux = (- V^3 + a*V^2 + b*V + c*Tf + d) 
    
    if aux > 0
        μ = 1/τfreeze
    else
        if V < 1.0e-4
            μ = 0
        else
            μ = 1/τmelt
        end
    end
       
    return μ*aux
end

function dV_dt(V, Tf, model_parameters)

    Tp, Tm, Vp, Vm, τmelt, τfreeze = model_parameters
    
    a = 3*(Vm + Vp)/2
    b = -3*Vm*Vp
    c = (Vp - Vm)^3/(2*(Tm - Tp))
    d = ( Tp*Vm^2*(Vm-3Vp) - Tm*Vp^2*(Vp-3Vm) )/(2*(Tm - Tp))
    
    function μ(V,Tf)
        
        if (- V^3 + a*V^2 + b*V + c*Tf + d) > 0
            return 1/τfreeze
        else
            if V < 1.0e-4
                return 0.0
            else
                return 1/τmelt
            end
        end
    end
       
    return μ(V,Tf)*(- V^3 + a*V^2 + b*V + c*Tf + d)
end


dV_dtA(V,Tf) = dV_dt_v2(V,Tf,Antarctica_params,Antarctica_coeffs)
dV_dtG(V,Tf) = dV_dt_v2(V,Tf,Greenland_params,Greenland_coeffs) 


#-----------------------------------------------------------------------------------
# Sea level rise
#-----------------------------------------------------------------------------------
# SLR ice (m)
SLR_G(VG) = SLRpotentialG*(1-VG)
SLR_A(VA) = SLRpotentialA*(1-VA)
SLRice(VG,VA) =  SLR_G(VG) + SLR_A(VA)
SLRrate_ice(VG,VA,δT_U) = - SLRpotentialG*(dV_dtG(VG,δT_U)) - SLRpotentialA*(dV_dtA(VA,δT_U)) # m/yr

# SLR thermal (m)
αU = 2.3e-4 #K^-1
αD = 1.3e-4 #K^-1
SLRthermal(δT_U,δT_D;αU = αU, αD = αD) = αU*δT_U*h_U + αD*δT_D*h_D
dδT_U_dt(M_A,I,δT_U,δT_D) = ( F(M_A,I) - β*δT_U - γ*(δT_U - δT_D) )/(cvol*h_U)
dδT_D_dt(δT_U,δT_D) = γ*(δT_U - δT_D) /(cvol*h_D)
SLRrate_thermal(M_A,I,δT_U,δT_D) = αU*dδT_U_dt(M_A,I,δT_U,δT_D)*h_U + αD*dδT_D_dt(δT_U,δT_D)*h_D


# SLR glacier (m)
τ_gl = 200 # years
SLRpotential_gl = 0.5 # meters
ζ = 2 # Related to equlibrium SLRglacier sensitivity to temperature change (Celsius)
SLRgl_eq(δT_U) = SLRpotential_gl*tanh(δT_U/ζ)
SLRglacier(S) = S
SLRrate_gl(δT_U,S) = (1/τ_gl)*(SLRgl_eq(δT_U) - S)


# SLRtotal(δT_U,δT_D,VG,VA,S) (m)
SLRtotal(δT_U,δT_D,VG,VA,S) = SLRthermal(δT_U,δT_D) + SLRice(VG,VA) + SLRglacier(S)

# SLRrate(M_A,I,δT_U,δT_D,VG,VA,S) (m/yr)
SLRrate(M_A,I,δT_U,δT_D,VG,VA,S) = SLRrate_thermal(M_A,I,δT_U,δT_D) + SLRrate_gl(δT_U,S) + SLRrate_ice(VG,VA,δT_U) 


#-----------------------------------------------------------------------------------
# SURFER model
#-----------------------------------------------------------------------------------

function surfer!(du,u,p,t)
    M_A, M_U, M_D, M_L, δT_U, δT_D, VG, VA, S = u
    
    Emissions, Injections = p
    
    du[1] = dM_A = Emissions(t) - k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) - k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))
    du[2] = dM_U = k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) - k_UD*(M_U - M_D/(δ*δDIC))
    du[3] = dM_D = k_UD*(M_U - M_D/(δ*δDIC))
    du[4] = dM_L = k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))
    du[5] = dδT_U = ( F(M_A,Injections(t)) - β*δT_U - γ*(δT_U - δT_D) )/(cvol*h_U)
    du[6] = dδT_D = γ*(δT_U - δT_D)/(cvol*h_D)
    du[7] = dVG = dV_dtG(VG,δT_U)
    du[8] = dVA = dV_dtA(VA,δT_U)
    du[9] = dS = (1/τ_gl)*(SLRgl_eq(δT_U) - S)

end

#-----------------------------------------------------------------------------------
# SURFER model with forcings depending also in Mat for smarter "options" in commitment experiments
#-----------------------------------------------------------------------------------
function surfer_commit!(du,u,p,t) 
    # surfer model with emissions and injections depending also in carbon reservoir masses and thermal reservoirs respectively to accomodate the carbon dioxide removal and SRM strategies
    M_A, M_U, M_D, M_L, δT_U, δT_D, VG, VA, S = u
    
    Emissions, Injections = p
    
    du[1] = dM_A = Emissions(t, M_A,M_U,M_L) - k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) - k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))
    du[2] = dM_U = k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) - k_UD*(M_U - M_D/(δ*δDIC))
    du[3] = dM_D = k_UD*(M_U - M_D/(δ*δDIC))
    du[4] = dM_L = k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))
    du[5] = dδT_U = ( F(M_A,Injections(t,M_A,δT_U,δT_D)) - β*δT_U - γ*(δT_U - δT_D) )/(cvol*h_U)
    du[6] = dδT_D = γ*(δT_U - δT_D)/(cvol*h_D)
    du[7] = dVG = dV_dtG(VG,δT_U)
    du[8] = dVA = dV_dtA(VA,δT_U)
    du[9] = dS = (1/τ_gl)*(SLRgl_eq(δT_U) - S)
        
end


