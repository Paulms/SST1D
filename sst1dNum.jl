#Kynch batch density function
fbk(C) = C*vo*exp(-rv*C)
#Encontramos el máximo
const Chat = 1/rv

#Gudonov flux for fbk
function Gj(Cl, Cr)
  if Cl < Cr
    min(fbk(Cr), fbk(Cl))
  else
    if (Chat - Cl)*(Chat - Cr) < 0
      fbk(Chat)
    else
      max(fbk(Cl), fbk(Cr))
    end
  end
end

# Numerical Flux
function Fj(Cl, Cr, j, jf, t)
  if j == -2 || j == -1
    -(Qe(t)/A)*Cr
  elseif j >= 0 && j <= (jf-1)
    -(Qe(t)/A)*Cr + Gj(Cl, Cr)
  elseif j >= jf && j <= N
    (Qu(t)/A)*Cl + Gj(Cl, Cr)
  elseif j == N+1 || j == N+2
    (Qu(t)/A)*Cl
  else
    error("Error celda no definida: $j")
  end
end

#Dispersion function
function ddisp(z,t=0)
  if abs(z) < α2*Qf(t)
    α1*Qf(t)*exp((-z^2/(α2*Qf(t))^2)/(1-abs(z)/(α2*Qf(t))))
  else
    0.0
  end
end

# Compresion
function dcomp(C)
  if C >= Cc
    ρs*α*vo*exp(-rv*C)/(grav*Δρ*(β+C-Cc))
  else
    0.0
  end
end

#Precompute of Di (Alg 3.2)
function ComputeD()
  M = N^2
  ΔC = (Cmax-Cc)/M
  D = zeros(M+1)
  di = zeros(M+1)
  D[1] = 0
  di[1] = dcomp(Cc)
  for i = 2:(M+1)
    di[i] = dcomp(Cc + i*ΔC)
    D[i] = D[i-1] + ΔC/2*(di[i-1]+di[i])
  end
  return ΔC, D
end

#Dj numérico (Nota: Dtilde se numera desde 1 no desde 0
#como en el paper)
function Djnum(Cj, ΔC, Dtilde)
  if Cj <= Cc
    0.0
  else
    i = Int(floor((Cj-Cc)/ΔC))
    ((Cj-Cc)/ΔC-i)*Dtilde[i+2]+(i+1-(Cj-Cc)/ΔC)*Dtilde[i+1]
  end
end

#Method of lines with explicit euler
function update_layer(j, Cl, Ca, Cr, Δz, jf, ΔC, Dtilde, Cf, t)
  zj = j*Δz-H
  zj1 = (j-1)*Δz-H
  if j == -1
    Ca+Δt/Δz*(Qe(t)/A*(Cr-Ca))
  elseif j == 0
    Ca+Δt/Δz*(Qe(t)/A*(Cr-Ca)-Gj(Ca, Cr)+(Djnum(Ca,ΔC,Dtilde)+Djnum(Cr,ΔC,Dtilde))/Δz)
  elseif j == 1
    Ca+Δt/Δz*(Qe(t)/A*(Cr-Ca)-(Gj(Ca, Cr)-Gj(Cl, Ca))+(ddisp(zj)*(Cr-Ca)+
    Djnum(Cr,ΔC,Dtilde)-2*Djnum(Ca,ΔC,Dtilde)+Djnum(Cl,ΔC,Dtilde))/Δz)
  elseif j >= 2 && j <= (jf-1)
    Ca+Δt/Δz*(Qe(t)/A*(Cr-Ca)-(Gj(Ca, Cr)-Gj(Cl, Ca))+(ddisp(zj)*(Cr-Ca)-ddisp(zj1)*(Ca-Cl)+
    Djnum(Cr,ΔC,Dtilde)-2*Djnum(Ca,ΔC,Dtilde)+Djnum(Cl,ΔC,Dtilde))/Δz)
  elseif j == jf
    Ca-Δt/Δz*((Qe(t)+Qu(t))/A*Ca+(Gj(Ca, Cr)-Gj(Cl, Ca))-(ddisp(zj)*(Cr-Ca)-ddisp(zj1)*(Ca-Cl)+
    Djnum(Cr,ΔC,Dtilde)-2*Djnum(Ca,ΔC,Dtilde)+Djnum(Cl,ΔC,Dtilde))/Δz-Qf(t)*Cf(t)/A)
  elseif j >= (jf+1) && j <= (N-1)
    Ca-Δt/Δz*(Qu(t)/A*(Ca-Cl)+(Gj(Ca, Cr)-Gj(Cl, Ca))-(ddisp(zj)*(Cr-Ca)-ddisp(zj1)*(Ca-Cl)+
    Djnum(Cr,ΔC,Dtilde)-2*Djnum(Ca,ΔC,Dtilde)+Djnum(Cl,ΔC,Dtilde))/Δz)
  elseif j == N
    Ca-Δt/Δz*(Qu(t)/A*(Ca-Cl)+(Gj(Ca, Cr)-Gj(Cl, Ca))-(-ddisp(zj1)*(Ca-Cl)+
    Djnum(Cr,ΔC,Dtilde)-2*Djnum(Ca,ΔC,Dtilde)+Djnum(Cl,ΔC,Dtilde))/Δz)
  elseif j == N+1
    Ca-Δt/Δz*(Qu(t)/A*(Ca-Cl)-Gj(Cl, Ca)+(Djnum(Ca,ΔC,Dtilde)-Djnum(Cl,ΔC,Dtilde))/Δz)
  elseif j == N+2
    Ca-Δt/Δz*(Qu(t)/A*(Ca-Cl))
  else
    error("Undefined layer: $j")
  end
end
