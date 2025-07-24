struct Components 
    R1 :: Float64
    R2 :: Float64
    R3 :: Float64
    R4 :: Float64
    C1 :: Float64
    C2 :: Float64
    C3 :: Float64
end 

mutable struct Constants 
    beta11 :: Float64
    beta12 :: Float64
    beta13 :: Float64
    beta14 :: Float64
    beta21 :: Float64
    beta22 :: Float64
    beta23 :: Float64
    beta24 :: Float64
    beta25 :: Float64
    beta26 :: Float64
    beta31 :: Float64
    beta32 :: Float64
    beta33 :: Float64
    beta34 :: Float64
    beta35 :: Float64
    beta36 :: Float64
    alpha11 :: Float64
    alpha12 :: Float64
    alpha13 :: Float64
    alpha21 :: Float64
    alpha22 :: Float64
    alpha23 :: Float64
    alpha24 :: Float64
    alpha25 :: Float64
    alpha31 :: Float64
    alpha32 :: Float64
    alpha33 :: Float64
    alpha34 :: Float64
    alpha35 :: Float64
    
    Constants() = new()

end 

Engl = Components(250e3, 1e6, 20e3, 47e3, 0.47e-9, 47e-9, 22e-9)
JCM = Components(220e3, 1e6, 22e3, 33e3, 0.47e-9, 22e-9, 22e-9)
Soldano = Components(250e3, 1e6, 25e3, 47e3, 0.47e-9, 20e-9, 20e-9)


comps = Soldano

ctes = Constants()

ctes.beta11 = comps.C1*comps.R1;
ctes.beta12 = comps.C3*comps.R3;
ctes.beta13 = comps.C1*comps.R2 + comps.C2*comps.R2;
ctes.beta14 = comps.C1*comps.R3 + comps.C2*comps.R3;

ctes.beta21 = comps.C1*comps.C2*comps.R1*comps.R4 + comps.C1*comps.C3*comps.R1*comps.R4;
ctes.beta22 = comps.C1*comps.C3*comps.R3*comps.R3 + comps.C2*comps.C3*comps.R3*comps.R3;
ctes.beta23 = comps.C1*comps.C3*comps.R1*comps.R3 + comps.C1*comps.C3*comps.R3*comps.R3 + comps.C2*comps.C3*comps.R3*comps.R3;
ctes.beta24 = comps.C1*comps.C2*comps.R1*comps.R2 + comps.C1*comps.C2*comps.R2*comps.R4 + comps.C1*comps.C3*comps.R2*comps.R4;
ctes.beta25 = comps.C1*comps.C3*comps.R2*comps.R3 + comps.C2*comps.C3*comps.R2*comps.R3;
ctes.beta26 = comps.C1*comps.C2*comps.R1*comps.R3 + comps.C1*comps.C2*comps.R3*comps.R4 + comps.C1*comps.C3*comps.R3*comps.R4;

ctes.beta31 = comps.C1*comps.C2*comps.C3*comps.R1*comps.R2*comps.R3 + comps.C1*comps.C2*comps.C3*comps.R2*comps.R3*comps.R4;
ctes.beta32 = comps.C1*comps.C2*comps.C3*comps.R1*comps.R3*comps.R3 + comps.C1*comps.C2*comps.C3*comps.R3*comps.R3*comps.R4;
ctes.beta33 = comps.C1*comps.C2*comps.C3*comps.R1*comps.R3*comps.R3 + comps.C1*comps.C2*comps.C3*comps.R3*comps.R3*comps.R4;
ctes.beta34 = comps.C1*comps.C2*comps.C3*comps.R1*comps.R3*comps.R4;
ctes.beta35 = comps.C1*comps.C2*comps.C3*comps.R1*comps.R3*comps.R4;
ctes.beta36 = comps.C1*comps.C2*comps.C3*comps.R1*comps.R2*comps.R4;

ctes.alpha11 = comps.C1*comps.R1 + comps.C1*comps.R3 + comps.C2*comps.R3 + comps.C2*comps.R4 + comps.C3*comps.R4;
ctes.alpha12 = comps.C3*comps.R3;
ctes.alpha13 = comps.C1*comps.R2 + comps.C2*comps.R2;

ctes.alpha21 = comps.C1*comps.C3*comps.R1*comps.R3 
             - comps.C2*comps.C3*comps.R3*comps.R4 
             + comps.C1*comps.C3*comps.R3*comps.R3 
             + comps.C2*comps.C3*comps.R3*comps.R3;

ctes.alpha22 = comps.C1*comps.C3*comps.R2*comps.R3 + comps.C2*comps.C3*comps.R2*comps.R3;
ctes.alpha23 = comps.C1*comps.C3*comps.R3*comps.R3 + comps.C2*comps.C3*comps.R3*comps.R3;
ctes.alpha24 = comps.C1*comps.C2*comps.R2*comps.R4 
             + comps.C1*comps.C2*comps.R1*comps.R2 
             + comps.C1*comps.C3*comps.R2*comps.R4 
             + comps.C2*comps.C3*comps.R2*comps.R4;
             
ctes.alpha25 = comps.C1*comps.C2*comps.R1*comps.R4 
             + comps.C1*comps.C3*comps.R1*comps.R4 
             + comps.C1*comps.C2*comps.R3*comps.R4 
             + comps.C1*comps.C2*comps.R1*comps.R3 
             + comps.C1*comps.C3*comps.R3*comps.R4 
             + comps.C2*comps.C3*comps.R3*comps.R4;

ctes.alpha31 = comps.C1*comps.C2*comps.C3*comps.R1*comps.R2*comps.R3 + comps.C1*comps.C2*comps.C3*comps.R2*comps.R3*comps.R4;
ctes.alpha32 = comps.C1*comps.C2*comps.C3*comps.R1*comps.R3*comps.R3 + comps.C1*comps.C2*comps.C3*comps.R3*comps.R3*comps.R4;
ctes.alpha33 = comps.C1*comps.C2*comps.C3*comps.R3*comps.R3*comps.R4 
             + comps.C1*comps.C2*comps.C3*comps.R1*comps.R3*comps.R3 
             - comps.C1*comps.C2*comps.C3*comps.R1*comps.R3*comps.R4;
ctes.alpha34 = comps.C1*comps.C2*comps.C3*comps.R1*comps.R2*comps.R4;
ctes.alpha35 = comps.C1*comps.C2*comps.C3*comps.R1*comps.R3*comps.R4;


@show ctes
