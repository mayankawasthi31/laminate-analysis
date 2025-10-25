clear all; clc;

% Assumptions:
%  No voids are considered
%  Plane stress condition

% --------------------------------------------------------------PHASE-1------------------------------------------------------------------


% Taking fiber and matrix properties, and the angle of the lamina as input :

Ef = 257e9;  % Young's modulus of fiber (in Pa)  
Em = 3e9;    % Young's modulus of matrix (in Pa)  
vf = 0.25;   % Poisson’s ratio of fiber  
vm = 0.35;    % Poisson’s ratio of matrix  
Gf = 103e9;   % Shear modulus of fiber (in Pa)  
Gm = 1.1e9;    % Shear modulus of matrix (in Pa)  
Vf = 0.70;   % Volume fraction of fiber  
Vm = 1 - Vf; % Volume fraction of matrix  
theta = 60;  % Fiber orientation angle (in degrees, measured anticlockwise from the global x-axis)

% Given global strains :
epsilonxy = [-7.31e-3; 3.18e-3; 7.82e-4];

% Ultimate strengths
F1t = 1500e6;  % (in Pa)
F1c = 1500e6;  % (in Pa)
F2t = 40e6;    % (in Pa)
F2c = 246e6;   % (in Pa)
F12 = 88e6;    % (in Pa)

%  Laminate definition
nLayersi = 4;
stack_seqi = [45, 60, -45, -60];  % (in degrees)            
ti = [0.003, 0.002, 0.003, 0.002]; % (in m)

% loads (in N/m and N·m/m)
Nx  = 1000e3;   Ny  = 200e3;  Nxy = 0;
Mx  = 0;    My  = 0;      Mxy = 0;

% Computing the engineering constants for the lamina :
E1 = Vf * Ef + Vm * Em; % (in Pa)
E2 = 1 / ( (Vf/Ef) + (Vm/Em) ); % (in Pa)
G12 = 1 / ( (Vf/Gf) + (Vm/Gm) ); % (in Pa)
v12 = Vf * vf + Vm * vm;
v21 = (E2 * v12) / E1;

% Engineering constants:
fprintf("\nComputed engineering constants :\n");
fprintf("E1 = %.2e Pa\n", E1);
fprintf("E2 = %.2e Pa\n", E2);
fprintf("G12 = %.2e Pa\n", G12);
fprintf("v12 = %.5f\n", v12);
fprintf("v21 = %.5f\n", v21);

% Computing the reduced compliance and stiffness matrices : 
S12 = [1/E1, -v12/E1, 0;
       -v21/E2, 1/E2, 0;
       0, 0, 1/G12]; % (in Pa^-1)

Q12 = inv(S12); % (in Pa)

fprintf("\nComputed compliance and stiffness matrices:\n");
fprintf("\nReduced compliance matrix S12 (in Pa^-1):\n");
disp(S12);

fprintf("Reduced stiffness matrix Q12 (in Pa):\n");
disp(Q12);

% Obtaining transformed reduced compliance and stiffness matrices : 
c = cosd(theta);
s = sind(theta);

T = [c^2, s^2, 2*c*s;
     s^2, c^2, -2*c*s;
     -c*s, c*s, c^2 - s^2];

R = [1, 0, 0;
     0, 1, 0;
     0, 0, 2];

Sxy = R * inv(T) * inv(R) * S12 * T; % (in Pa^-1)
Qxy = inv(Sxy); % (in Pa)

fprintf("\nObtained transformed compliance and stiffness matrices:\n");
fprintf("\nTransformed reduced compliance matrix Sxy (in Pa^-1):\n");
disp(Sxy);

fprintf("Transformed reduced stiffness matrix Qxy (in Pa):\n");
disp(Qxy);

% Computing the local strains, global and local stresses in the lamina :
epsilon12 = R * T * inv(R) * epsilonxy;
sigmaxy = Qxy * epsilonxy; % (in Pa)
sigma12 = T * sigmaxy; % (in Pa)

fprintf("\nComputed strains and stresses :\n");
fprintf("\nLocal Strain:\n");
disp(epsilon12);

fprintf("Global stress (in Pa):\n");
disp(sigmaxy);

fprintf("Local stress (in Pa):\n");
disp(sigma12);


% --------------------------------------------------------------PHASE-2------------------------------------------------------------------


% Implementation of Failure Theories: (1) Maximum Stress    (2) Maximum Strain    (3) Tsai–Hill    (4) Tsai–Wu  to  identify the failure mode.
% If the lamina is safe, compute Strength Ratio.

% 1) Maximum Stress Failure Theory

sigma1 = sigma12(1); % Local stresses (in Pa)
sigma2 = sigma12(2);                % (in Pa)
tau12 = sigma12(3);                 % (in Pa) 

failModes_stress = {};
ratios_stress = [];

% Longitudinal stress check
if sigma1 >= 0
    if sigma1 >= F1t
        failModes_stress{end+1} = 'Longitudinal tensile stress failure';
        ratios_stress(end+1) = sigma1 / F1t;
    else
        SR_stress1 = F1t / sigma1;
    end
elseif sigma1 < 0
    if abs(sigma1) >= F1c
        failModes_stress{end+1} = 'Longitudinal compressive stress failure';
        ratios_stress(end+1) = abs(sigma1) / F1c;
    else
        SR_stress1 = F1c / abs(sigma1);
    end
end

% Transverse stress check
if sigma2 >= 0
    if sigma2 >= F2t
        failModes_stress{end+1} = 'Transverse tensile stress failure';
        ratios_stress(end+1) = sigma2 / F2t;
    else
        SR_stress2 = F2t / sigma2;
    end
elseif sigma2 < 0
    if abs(sigma2) >= F2c
        failModes_stress{end+1} = 'Transverse compressive stress failure';
        ratios_stress(end+1) = abs(sigma2) / F2c;
    else
        SR_stress2 = F2c / abs(sigma2);
    end
end

% In-plane stress check
if tau12 >= 0
    if tau12 >= F12
        failModes_stress{end+1} = 'In-plane shear stress failure';
        ratios_stress(end+1) = tau12 / F12;
    else
        SR_stress12 = F12 / tau12;
    end
elseif tau12 < 0
    if abs(tau12) >= F12
        failModes_stress{end+1} = 'In-plane shear stress failure';
        ratios_stress(end+1) = abs(tau12) / F12;
    else
        SR_stress12 = F12 / abs(tau12);
    end
end

if ~isempty(failModes_stress)
    fprintf('\n** Maximum Stress: FAILED **\n');
    for i = 1:length(failModes_stress)
        fprintf('%s, ratio = %.3f\n', failModes_stress{i}, ratios_stress(i));
    end
    [maxRatio_stress, index_stress] = max(ratios_stress);
    fprintf('Dominant mode: %s (ratio = %.3f)\n', failModes_stress{index_stress}, maxRatio_stress);
    fprintf('No SR since lamina fails here.\n\n');
else
    fprintf('\nMaximum Stress: Lamina is SAFE.\n');
    SR_maxStress = min([SR_stress1, SR_stress2, SR_stress12]);
    fprintf('Strength Ratio (SR) from Maximum Stress: %.3f\n\n', SR_maxStress);
end



% 2) Maximum Strain Failure Theory

epsilon1 = epsilon12(1); % Local strains
epsilon2 = epsilon12(2);
gamma12 = epsilon12(3);


eps1_tens = F1t / E1; % Ultimate strains
eps1_comp = -F1c / E1;
eps2_tens = F2t / E2;
eps2_comp = -F2c / E2;
gamma12_ult = F12 / G12;

failModes_strain = {};
ratios_strain = [];

% Longitudinal strain check
if epsilon1 >= 0
    if epsilon1 >= eps1_tens
        failModes_strain{end+1} = 'Longitudinal tensile strain failure';
        ratios_strain(end+1) = epsilon1 / eps1_tens;
    else
        SR_strain1 = eps1_tens / epsilon1;
    end
else
    if epsilon1 <= eps1_comp
        failModes_strain{end+1} = 'Longitudinal compressive strain failure';
        ratios_strain(end+1) = abs(epsilon1) / abs(eps1_comp);
    else
        SR_strain1 = abs(eps1_comp) / abs(epsilon1);
    end
end

% Transverse strain check
if epsilon2 >= 0
    if epsilon2 >= eps2_tens
        failModes_strain{end+1} = 'Transverse tensile strain failure';
        ratios_strain(end+1) = epsilon2 / eps2_tens;
    else
        SR_strain2 = eps2_tens / epsilon2;
    end
else
    if epsilon2 <= eps2_comp
        failModes_strain{end+1} = 'Transverse compressive strain failure';
        ratios_strain(end+1) = abs(epsilon2) / abs(eps2_comp);
    else
        SR_strain2 = abs(eps2_comp) / abs(epsilon2);
    end
end

% In-plane strain check
if abs(gamma12) >= gamma12_ult
    failModes_strain{end+1} = 'In-plane shear strain failure';
    ratios_strain(end+1) = abs(gamma12) / gamma12_ult;
else
    SR_strain12 = gamma12_ult / abs(gamma12);
end

if ~isempty(failModes_strain)
    fprintf('** Maximum Strain: FAILED **\n');
    for i = 1:length(failModes_strain)
        fprintf('%s, ratio = %.3f\n', failModes_strain{i}, ratios_strain(i));
    end
    [maxRatio_strain, index_strain] = max(ratios_strain);
    fprintf('Dominant mode: %s (ratio = %.3f)\n', failModes_strain{index_strain}, maxRatio_strain);
    fprintf('No SR since lamina fails here.\n\n');
else
    SR_strain = min([SR_strain1, SR_strain2, SR_strain12]);
    fprintf('Maximum Strain: Lamina is SAFE.\n');
    fprintf('Strength Ratio (SR) from Maximum Strain: %.3f\n\n', SR_strain);
end


% 3)  Tsai–Hill Failure Theory

if sigma1 >= 0
    X1 = F1t;   % (in Pa)
else
    X1 = F1c;  % (in Pa)
end

if sigma2 >= 0
    X2 = F2t;  % (in Pa)
else
    X2 = F2c;  % (in Pa)
end

S = F12;       % (in Pa)

tsaiHillVal = (sigma1/X1)^2 + (sigma2/X2)^2 - (sigma1*sigma2)/(X1^2) + (tau12/S)^2;
fprintf('Tsai–Hill index = %.3f\n', tsaiHillVal);

if tsaiHillVal >= 1
    fprintf('** Tsai–Hill: FAILED **\n');
    if sigma1 >= 0
        ratio_TH1 = sigma1 / F1t;  mode1 = 'Longitudinal tensile failure';
    else
        ratio_TH1 = abs(sigma1) / F1c;  mode1 = 'Longitudinal compressive failure';
    end
    if sigma2 >= 0
        ratio_TH2 = sigma2 / F2t;  mode2 = 'Transverse tensile failure';
    else
        ratio_TH2 = abs(sigma2) / F2c;  mode2 = 'Transverse compressive failure';
    end
    ratio_TH3 = abs(tau12) / F12;  mode3 = 'In-plane shear failure';
    
    [maxRatio_TH, index_TH] = max([ratio_TH1, ratio_TH2, ratio_TH3]);
    if index_TH == 1
        dominantMode_TH = mode1;
    elseif index_TH == 2
        dominantMode_TH = mode2;
    else
        dominantMode_TH = mode3;
    end
    fprintf('Ratios: Fiber = %.3f, Matrix = %.3f, Shear = %.3f\n', ratio_TH1, ratio_TH2, ratio_TH3);
    fprintf('Dominant mode: %s (ratio = %.3f)\n', dominantMode_TH, maxRatio_TH);
    fprintf('No SR since lamina fails here.\n\n');
else
    fprintf('Tsai–Hill: Lamina is SAFE.\n');
    denominator = (sigma1^2/X1^2) + (sigma2^2/X2^2) + (tau12^2/S^2) - (sigma1*sigma2)/(X1^2);
    if denominator <= 0
        fprintf('Warning: denominator <= 0, cannot compute SR.\n\n');
    else
        SR_TsaiHill = sqrt(1/denominator);
        fprintf('Strength Ratio (SR) from Tsai–Hill: %.3f\n\n', SR_TsaiHill);
    end
end


% 4)  Tsai–Wu Failure Theory

H1  = 1/F1t - 1/F1c;                 % (in Pa^-1)
H2  = 1/F2t - 1/F2c;                 % (in Pa^-1)
H6  = 0;                             % (in Pa^-1)
H11 = 1/(F1t*F1c);                   % (in Pa^-2)
H22 = 1/(F2t*F2c);                   % (in Pa^-2)
H66 = 1/(F12^2);                     % (in Pa^-2)
H12 = -0.5 * sqrt(H11 * H22);        % (in Pa^-2)

TW = H1*sigma1 + H2*sigma2 + H6*tau12 + H11*sigma1^2 + H22*sigma2^2 + H66*tau12^2 + 2*H12*sigma1*sigma2;
fprintf('Tsai–Wu index = %.3f\n', TW);

if TW >= 1
    fprintf('** Tsai–Wu: FAILED **\n');
    if sigma1 >= 0
        ratio_TW1 = sigma1 / F1t;  mode1 = 'Longitudinal tensile failure';
    else
        ratio_TW1 = abs(sigma1) / F1c;  mode1 = 'Longitudinal compressive failure';
    end
    if sigma2 >= 0
        ratio_TW2 = sigma2 / F2t;  mode2 = 'Transverse tensile failure';
    else
        ratio_TW2 = abs(sigma2) / F2c;  mode2 = 'Transverse compressive failure';
    end
    ratio_TW3 = abs(tau12) / F12;  mode3 = 'In-plane shear failure';
    [maxRatio_TW, index_TW] = max([ratio_TW1, ratio_TW2, ratio_TW3]);
    if index_TW == 1
        dominantMode_TW = mode1;
    elseif index_TW == 2
        dominantMode_TW = mode2;
    else
        dominantMode_TW = mode3;
    end
    fprintf('Ratios: Fiber = %.3f, Matrix = %.3f, Shear = %.3f\n', ratio_TH1, ratio_TH2, ratio_TH3);
    fprintf('Dominant mode: %s (ratio = %.3f)\n', dominantMode_TW, maxRatio_TW);
    fprintf('No SR since lamina fails here.\n\n');
else
    fprintf('Tsai–Wu: Lamina is SAFE.\n');
    A = H11*sigma1^2 + H22*sigma2^2 + H66*tau12^2 + 2*H12*sigma1*sigma2;
    B = H1*sigma1 + H2*sigma2 + H6*tau12;
    if A == 0
        if B ~= 0
            SR_TW = 1/B;
        else
            fprintf('Cannot compute SR, A and B both zero.\n');
            SR_TW = NaN;
        end
    else
        discriminant = sqrt(B^2 + 4*A);
        SR_stress_TW1 = (-B + discriminant)/(2*A);
        SR_stress_TW2 = (-B - discriminant)/(2*A);
        if SR_stress_TW1 > 0 && SR_stress_TW2 > 0
            SR_TW = min(SR_stress_TW1, SR_stress_TW2);
        elseif SR_stress_TW1 > 0
            SR_TW = SR_stress_TW1;
        elseif SR_stress_TW2 > 0
            SR_TW = SR_stress_TW2;
        else
            fprintf('No positive SR found.\n');
            SR_TW = NaN;
        end
    end
    fprintf('Strength Ratio (SR) from Tsai–Wu: %.3f\n\n', SR_TW);
end

% --------------------------------------------------------------PHASE-3-----------------------------------------------------------------

% Lamina properties and stiffness matrix S12 are defined in Phase-1 and used here directly.

% Validate inputs
if numel(stack_seqi)~=nLayersi || numel(ti)~=nLayersi
    error('nLayersi (%d) must equal length(stack_seqi) and length(ti).', nLayersi);
end

% Symmetry check
isSym = isequal(stack_seqi, fliplr(stack_seqi)) && isequal(ti, fliplr(ti));

% Mid-plane split 
H = sum(ti);
halfH = H/2;
cum_t = cumsum(ti);
k_mid = find(cum_t >= halfH, 1);

if isSym || cum_t(k_mid)==halfH
    stack_seq = stack_seqi;  
    t_k = ti;
else
    t_top = cum_t(k_mid) - halfH;
    t_bot = ti(k_mid) - t_top;
    stack_seq = [ stack_seqi(1:k_mid-1), stack_seqi(k_mid), stack_seqi(k_mid), stack_seqi(k_mid+1:end) ];
    t_k = [ ti(1:k_mid-1), t_top, t_bot, ti(k_mid+1:end) ];   % (in m)
end
nLayers = numel(t_k);

% Computing interface heights about mid-plane
h = zeros(nLayers+1,1);  % (in m)
h(1) = -halfH;
for j = 1:nLayers
    h(j+1) = h(j) + t_k(j);
end

% Assemble A, B, D and store per-ply stiffness Qbar
A = zeros(3); % (in Pa.m)
B = zeros(3); % (in Pa.m^2)
D = zeros(3); % (in Pa.m^3)
Qbar_all = zeros(3,3,nLayers);

for k = 1:nLayers
    thetak = stack_seq(k);
    c_k = cosd(thetak);  s_k = sind(thetak);
    T_k = [ c_k^2, s_k^2, 2*c_k*s_k;
          s_k^2, c_k^2, -2*c_k*s_k;
         -c_k*s_k, c_k*s_k, c_k^2-s_k^2 ];
    R_k = diag([1 1 2]);
    S_xy = ((R_k / T_k) / R_k) * S12 * T_k;
    Qbar = inv(S_xy);
    Qbar_all(:,:,k) = Qbar;
    A = A + Qbar * t_k(k);
    B = B + Qbar * ((h(k+1)^2 - h(k)^2)/2);
    D = D + Qbar * ((h(k+1)^3 - h(k)^3)/3);
end

% enforce exact zero coupling if symmetric
if isSym
    B = zeros(3);
end

% display A, B, D
fprintf('A matrix (Pa·m):\n'); disp(A);
fprintf('B matrix (Pa·m^2):\n'); disp(B);
fprintf('D matrix (Pa·m^3):\n'); disp(D);

% Solve for mid-plane strains ε0 and curvatures κ
ABD   = [A, B; B, D];
loads = [Nx; Ny; Nxy; Mx; My; Mxy];
sol   = ABD \ loads;
epsilon0 = sol(1:3);
kappa    = sol(4:6);
fprintf('Mid-plane strains: [%g, %g, %g], Curvatures: [%g, %g, %g]\n', sol(1:3), sol(4:6));

% Through-thickness sampling of strains & stresses
nPts = 200;
zz   = linspace(h(1), h(end), nPts);
ex   = zeros(1,nPts);  ey=zeros(1,nPts);  gxy=zeros(1,nPts);
sx   = zeros(1,nPts);  sy=zeros(1,nPts);  txy=zeros(1,nPts);

for m = 1:nPts
    zi = zz(m);
    ply = find(zi >= h(1:end-1) & zi <= h(2:end), 1);
    epsv = epsilon0 + zi * kappa;
    ex(m)  = epsv(1);  ey(m) = epsv(2);  gxy(m) = epsv(3);
    sig   = Qbar_all(:,:,ply) * epsv;
    sx(m)  = sig(1);  sy(m) = sig(2);  txy(m) = sig(3);
end

% Plot strains (strain on x-axis, z on y-axis)
figure; plot(ex, zz, 'LineWidth',1.2);
xlabel('\epsilon_x'); ylabel('z (m)');
title('\epsilon_x through thickness');

figure; plot(ey, zz, 'LineWidth',1.2);
xlabel('\epsilon_y'); ylabel('z (m)');
title('\epsilon_y through thickness');

figure; plot(gxy, zz, 'LineWidth',1.2);
xlabel('\gamma_{xy}'); ylabel('z (m)');
title('\gamma_{xy} through thickness');

% Plot stresses (stress on x-axis, z on y-axis)
figure; plot(sx, zz, 'LineWidth',1.2);
xlabel('\sigma_x (Pa)'); ylabel('z (m)');
title('\sigma_x through thickness');

figure; plot(sy, zz, 'LineWidth',1.2);
xlabel('\sigma_y (Pa)'); ylabel('z (m)');
title('\sigma_y through thickness');

figure; plot(txy, zz, 'LineWidth',1.2);
xlabel('\tau_{xy} (Pa)'); ylabel('z (m)');
title('\tau_{xy} through thickness');
