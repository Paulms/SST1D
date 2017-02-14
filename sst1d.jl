# Constantes:
const vo = 3.47 #m/h
const rv = 0.37 #m^3/kg
const α = 4.00 #Pa
const β = 4.00 #Kg/m³
const ρs = 1050 #Kg/m³
const Δρ = 52 #Kg/m³
const grav = 9.81 #m/s²
const Cc = 6.00 #Kg/m³
const Cmax = 20 #Kg/m³
# Dimensiones del Tanque
const H = 1 #m
const B = 3 #m
const A = 400 #m²
#Discretización
const N = 90
const Δt = 0.001
const Tmax = 800 #h

# Flujos volumétricos
Qf(t) = 250 #m³/h
Qu(t) = 80 #m³/h
Qe(t) = Qf(t) - Qu(t)

#Se van a dejar variar
const α1 = 0 #m-1
const α2 = 0 #m-1

# Leemos funciones auxiliares
include("sst1dNum.jl")

# Iteración temporal
function simular(tiempo, Ct, Cf)
  #Inicializamos sistema
  C = Ct[:,end]
  Δz = (B+H)/N
  jf = ceil(H/Δz)
  ΔC, Dtilde = ComputeD()
  t = 2
  #Iteramos sobre el tiempo
  while tiempo <= Tmax
    tiempo = tiempo + Δt
    if tiempo%50 < Δt
      println("Tiempo: $tiempo")
    end
    #Update Layers
    for j = -1:(N+2)
      Cl = j < 0 ? 0.0 : C[j+1]
      Cr = j>(N+1) ? 0.0 : C[j+3]
      C[j+2] = update_layer(j, Cl, C[j+2], Cr, Δz, jf, ΔC, Dtilde, Cf, tiempo)
    end
    #Save C each hour
    if tiempo%1 < Δt
      Ct[:,t] = C
      t = t + 1
    end
  end
end

#Llegar primero a estado estacionario
Cf(t) = 4.0
C = ones(N+4, Tmax+1)*4.0
@time simular(0.0, C, Cf)
#Ahora simulamos ejemplo
Cf2(t) = begin
  if (t < 50)
    4.0
  elseif t < 250
    3.7
  else
    4.1
  end
end
@time simular(0.0, C, Cf2)

#Profile
#Profile.init(delay=0.001)
#Profile.clear()
#C = ones(N+4, Tmax+1)*4.0
#@profile simular(0.0, C, Cf)
#using ProfileView
#ProfileView.view()

#Graficar
#using Plots
#Z = (-1:(N+2))*(B+H)/N-H
#plot(Z,C[:,Tmax+1])

#using PyPlot
#fig = figure(figsize=(8,6))
#ax = fig[:gca](projection="3d")
#grid_a = [i for i in Z, j in 1:(Tmax+1)]
#grid_b = [j for i in Z, j in 1:(Tmax+1)]
#ax[:plot_surface](grid_a, grid_b,C, rstride=2, cstride=2, cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
#fig[:show]()
