%% Transfer Matrix Method Summary
%   plots reflactance and transmittance over wavelength range
%   5 layer dielectric. Layer 3 is Al2O3, Layer 5 is TiO2
%% Can change the following parameters

% wavelength range [micrometers] (should be within range of imported data)
wavelength_min = 0.3;
wavelength_max = 2;

% angle in incidence [degrees]
incidence = 60;

% polarization ["s" for TE / "p" for TM]
polarization = "p";

%refractive index of layers
    % n3 is fixed to Al2O3
    % n5 is fixed to TiO2
nclad = 1.3;
n1 = 2.0;
n2 = 2.5;
n4 = 3.5;
nsub = 1.3;

% thickness of each layers [nanometers]
d1 = 150;
d2 = 200;
d3 = 200;
d4 = 150;
d5 = 150;


%% The algorithm

a = importdata('Al2O3.csv');
b = importdata('Siefke.csv'); % make sure wavelength range is within data
na = a.data;
nb = b.data;

d1 = d1*10^-9;
d2 = d2*10^-9;
d3 = d3*10^-9;
d4 = d4*10^-9;
d5 = d5*10^-9;

E0 = 8.854*10^-12;
U0 = 4*pi*10^-7;

Ec = (nclad.^2) * E0;
E1 = (n1.^2) * E0;
E2 = (n2.^2) * E0;
E4 = (n4.^2) * E0;
Es = (nsub.^2) * E0;
U  = (1) * U0; % Ur = 1 for all layers

index = 0;
for i=wavelength_min:0.001:wavelength_max
    index = index+1;
    wavelength = i;

    n3 = interp1(na(:,1), na(:,2), wavelength);
    n5 = interp1(nb(:,1), nb(:,2), wavelength);
    
    
    
    E3 = (n3.^2) * E0;
    E5 = (n5.^2) * E0;
    
    k0 = 2*pi/(wavelength*10^-6);

    theta1 = asind((nclad/n1) * sind(incidence));
    theta2 = asind((n1/n2) * sind(theta1));
    theta3 = asind((n2/n3) * sind(theta2));
    theta4 = asind((n3/n4) * sind(theta3));
    theta5 = asind((n4/n5) * sind(theta4));
    refraction = asind((n5/nsub) * sind(theta5));

    phi1 = k0*n1*d1*cosd(theta1);
    phi2 = k0*n2*d2*cosd(theta2);
    phi3 = k0*n3*d3*cosd(theta3);
    phi4 = k0*n4*d4*cosd(theta4);
    phi5 = k0*n5*d5*cosd(theta5);

    P1 = [ exp(1i*phi1) 0    ;   0 exp(-1i*phi1) ];
    P2 = [ exp(1i*phi2) 0    ;   0 exp(-1i*phi2) ];
    P3 = [ exp(1i*phi3) 0    ;   0 exp(-1i*phi3) ];
    P4 = [ exp(1i*phi4) 0    ;   0 exp(-1i*phi4) ];
    P5 = [ exp(1i*phi5) 0    ;   0 exp(-1i*phi5) ];

    if polarization == "s"
        Dc = [ 1 1 ; sqrt(Ec/U)*cosd(incidence)   -sqrt(Ec/U)*cosd(incidence) ];
        D1 = [ 1 1 ; sqrt(E1/U)*cosd(theta1)      -sqrt(E1/U)*cosd(theta1)    ];
        D2 = [ 1 1 ; sqrt(E2/U)*cosd(theta2)      -sqrt(E2/U)*cosd(theta2)    ];
        D3 = [ 1 1 ; sqrt(E3/U)*cosd(theta3)      -sqrt(E3/U)*cosd(theta3)    ];
        D4 = [ 1 1 ; sqrt(E4/U)*cosd(theta4)      -sqrt(E4/U)*cosd(theta4)    ];
        D5 = [ 1 1 ; sqrt(E5/U)*cosd(theta5)      -sqrt(E5/U)*cosd(theta5)    ];
        Ds = [ 1 1 ; sqrt(Es/U)*cosd(refraction)  -sqrt(Es/U)*cosd(refraction)];

    elseif polarization == "p"
        Dc = [ cosd(incidence)  cosd(incidence)   ; sqrt(Ec/U)    -sqrt(Ec/U) ];
        D1 = [ cosd(theta1)     cosd(theta1)      ; sqrt(E1/U)    -sqrt(E1/U) ];
        D2 = [ cosd(theta2)     cosd(theta2)      ; sqrt(E2/U)    -sqrt(E2/U) ];
        D3 = [ cosd(theta3)     cosd(theta3)      ; sqrt(E3/U)    -sqrt(E3/U) ];
        D4 = [ cosd(theta4)     cosd(theta4)      ; sqrt(E4/U)    -sqrt(E4/U) ];
        D5 = [ cosd(theta5)     cosd(theta5)      ; sqrt(E5/U)    -sqrt(E5/U) ];
        Ds = [ cosd(refraction) cosd(refraction)  ; sqrt(Es/U)    -sqrt(Es/U) ];

    else
        disp('ERROR: use either "s" or "p" for polarization')
    end

    M = (Dc\D1)*P1*(D1\D2)*P2*(D2\D3)*P3*(D3\D4)*P4*(D4\D5)*P5*(D5\Ds);

    r = M(2,1)/M(1,1);
    t = 1/M(1,1);

    R(index) = abs(r).^2;
    T(index) = abs(t).^2;
    wl(index) = i;

end

    subplot(2,1,1)
    plot(wl,R)
    ylabel('Reflectance')
    xlabel('Wavelength [micrometers]')
    
    subplot(2,1,2)
    plot(wl,T)
    ylabel('Transmittance')
    xlabel('Wavelength [micrometers]')
    
    
    