clear; clc; close all;

%% IM7 / HexPly 8552 — RTD lamina properties

E1   = 162.0;    % GPa (fiber direction)
E2   = 8.96;     % GPa (transverse)
G12  = 4.69;     % GPa (in-plane shear)
nu12 = 0.316;

Xt = 2.558;      % GPa (longitudinal tension)
Xc = 1.689;      % GPa (longitudinal compression)
Yt = 0.0640;     % GPa (transverse tension)
Yc = 0.2856;     % GPa (transverse compression)
S  = 0.0911;     % GPa (in-plane shear)

%% --- Geometry ---
plyThickness = 0.18;          
numPlies = 12;
totalThickness = numPlies * plyThickness;  

specimenThickness = 2.2;
supportSpan = 60 * specimenThickness;
specimenLength = 1.2 * supportSpan;
specimenWidth  = 13;

targetFlexuralModulus = 106;

%% --- Layup ---
%layup = [0 90 0 -45 90 45 45 90 -45 0 90 0];
layup = [27.5 0 0 0 0 0 0 0 0 0 0 27.5];

fprintf('=== CLT ANALYSIS FOR LAYUP ===\n');
fprintf('Layup: [ %s]\n\n', sprintf('%g ', layup));

%% --- Reduced Stiffness Matrix Q ---
nu21  = nu12 * E2 / E1;
denom = 1 - nu12 * nu21;

Q = [E1/denom, nu12*E2/denom, 0;
     nu12*E2/denom, E2/denom, 0;
     0, 0, G12];

Q11 = Q(1,1); Q12 = Q(1,2); Q22 = Q(2,2); Q66 = Q(3,3);

%% --- Compute z-coordinates ---
t = plyThickness;
theta_full = layup;
nply = numel(theta_full);

Ttot = nply * t;
z = zeros(1, nply+1);
z(1) = -Ttot/2;
for i = 1:nply, z(i+1) = z(i) + t; end

%% --- ABD Matrices ---
A = zeros(3,3); B = zeros(3,3); D = zeros(3,3);

for k = 1:nply
    th = theta_full(k);
    c = cosd(th); s = sind(th);
    c2=c^2; s2=s^2; c3=c^3; s3=s^3; c4=c^4; s4=s^4;

    % Qbar (transformed stiffness)
    Q11bar = Q11*c4 + Q22*s4 + 2*(Q12+2*Q66)*s2*c2;
    Q12bar = (Q11+Q22-4*Q66)*s2*c2 + Q12*(c4+s4);
    Q22bar = Q11*s4 + Q22*c4 + 2*(Q12+2*Q66)*s2*c2;
    Q16bar = (Q11-Q12-2*Q66)*c3*s - (Q22-Q12-2*Q66)*c*s3;
    Q26bar = (Q11-Q12-2*Q66)*c*s3 - (Q22-Q12-2*Q66)*c3*s;
    Q66bar = (Q11+Q22-2*Q12-2*Q66)*s2*c2 + Q66*(s4+c4);

    Qbar = [Q11bar Q12bar Q16bar;
            Q12bar Q22bar Q26bar;
            Q16bar Q26bar Q66bar];

    zk1 = z(k);
    zk  = z(k+1);

    % A, B, D accumulation
    A = A + Qbar * (zk - zk1);
    B = B + 0.5 * Qbar * (zk^2 - zk1^2);
    D = D + (1/3) * Qbar * (zk^3 - zk1^3);
end
%%
Astar = inv(A);          % A*
Bstar = -(A\B);          % -A^{-1}B
Cstar = B/A;             % BA^{-1}
Dstar = D - B*(A\B);     % D - BA^{-1}B

Di = inv(Dstar);                 % D' = (D*)^{-1}
Ai = Astar - Bstar*Di*Cstar;     % A'
Bi = Bstar*Di;                   % B'
Ci = -Di*Cstar;                  % C'
%% --- Flexural Modulus ---
Ef = 12 / (Di(1,1)*totalThickness^3);
Ef_GPa = Ef;
errorPercent = ((Ef_GPa - targetFlexuralModulus)/targetFlexuralModulus) * 100;

fprintf('Flexural Modulus (Ef): %.2f GPa  (error vs target: %.1f%%)\n', Ef_GPa, errorPercent);

%% --- Tensile Modulus ---
Ex = 1/(Ai(1,1)*totalThickness);
Ex_GPa = Ex;
fprintf('Tensile Modulus (Ex): %.2f GPa\n\n', Ex_GPa);

%% --- First Ply Failure (Tsai-Hill) ---
strainRange = 0.0001:0.0001:0.02;
ultimateStrain = 0;

for strain = strainRange
    epsilon = [strain; 0; 0];
    sigma_global = A * epsilon / totalThickness;

    failed = false;

    for k = 1:numPlies
        theta = deg2rad(layup(k));
        c = cos(theta); s = sin(theta);

        sig1  = sigma_global(1)*c^2 + sigma_global(2)*s^2 + 2*sigma_global(3)*s*c;
        sig2  = sigma_global(1)*s^2 + sigma_global(2)*c^2 - 2*sigma_global(3)*s*c;
        tau12 = -sigma_global(1)*s*c + sigma_global(2)*s*c + sigma_global(3)*(c^2-s^2);

        if sig1 >= 0, X = Xt; else, X = Xc; end
        if sig2 >= 0, Y = Yt; else, Y = Yc; end

        F = (sig1/X)^2 - (sig1*sig2)/(X^2) + (sig2/Y)^2 + (tau12/S)^2;

        if F >= 1, failed = true; break; end
    end

    if failed, ultimateStrain = strain; break; end
end

if ultimateStrain == 0
    warning('No ply failure detected within strain range. Increase strainRange upper limit.');
end

ultimateStrength_MPa = Ex_GPa * ultimateStrain * 1000;

fprintf('Ultimate Tensile Strength: %.1f MPa\n', ultimateStrength_MPa);

%% --- Progressive Ply Failure ---
fprintf('\nPROGRESSIVE FAILURE ANALYSIS (1%% strain):\n');

activePlies = layup;
iteration = 1;
appliedStrain = 0.01;

while ~isempty(activePlies) && iteration <= 20
    [A_current,~,~] = calculateABD(activePlies, Q, plyThickness);
    h_current = length(activePlies) * plyThickness;

    epsilon = [appliedStrain; 0; 0];
    sigma_global = A_current * epsilon / h_current;

    maxF = 0; failedIndex = 0;

    for k = 1:length(activePlies)
        theta = deg2rad(activePlies(k));
        c=cos(theta); s=sin(theta);

        sig1  = sigma_global(1)*c^2 + sigma_global(2)*s^2 + 2*sigma_global(3)*s*c;
        sig2  = sigma_global(1)*s^2 + sigma_global(2)*c^2 - 2*sigma_global(3)*s*c;
        tau12 = -sigma_global(1)*s*c + sigma_global(2)*s*c + sigma_global(3)*(c^2-s^2);

        if sig1 < 0, X = Xc; else, X = Xt; end
        if sig2 < 0, Y = Yc; else, Y = Yt; end

        F = (sig1/X)^2 - (sig1*sig2)/(X^2) + (sig2/Y)^2 + (tau12/S)^2;

        if F > maxF
            maxF = F;
            if F >= 1, failedIndex = k; end
        end
    end

    if failedIndex > 0
        fprintf('Iteration %d: %g° ply failed\n', iteration, activePlies(failedIndex));
        activePlies(failedIndex) = [];
        iteration = iteration + 1;
    else
        break;
    end
end

fprintf('Remaining plies: %d\n\n', length(activePlies));

%% --- Prepreg Area (angle-dependent bounding box) ---
% Panel holds 4 specimens side-by-side along the width direction,
% with a 0.5 in (12.7 mm) margin on every edge.
margin_mm      = 12.7;
panelWidth_mm  = 4 * specimenWidth + 2*margin_mm;   % 4 strips across
panelLength_mm = specimenLength    + 2*margin_mm;

% For each ply at angle theta, cutting from a 0-deg prepreg roll requires
% a bounding rectangle larger than the panel itself. Sum over all plies.
prepregArea_mm2 = 0;
for k = 1:numPlies
    c = abs(cosd(layup(k)));
    s = abs(sind(layup(k)));
    bb_w = panelWidth_mm * c + panelLength_mm * s;
    bb_l = panelWidth_mm * s + panelLength_mm * c;
    prepregArea_mm2 = prepregArea_mm2 + bb_w * bb_l;
end
prepregArea_cm2 = prepregArea_mm2 / 100;

fprintf('Prepreg Area Required: %.1f cm²\n', prepregArea_cm2);

%% --- Summary ---
fprintf('\n=== SUMMARY ===\n');
fprintf('Flexural Modulus: %.2f GPa  (target: %g GPa, error: %.1f%%)\n', Ef_GPa, targetFlexuralModulus, errorPercent);
fprintf('Tensile Modulus:  %.2f GPa\n', Ex_GPa);
fprintf('Ultimate Strength: %.1f MPa\n', ultimateStrength_MPa);
fprintf('Prepreg Area:     %.1f cm²\n', prepregArea_cm2);


%% --- ABD helper ---
function [A,B,D] = calculateABD(layup,Q,t)

    n = numel(layup);
    z = linspace(-n*t/2, n*t/2, n+1);

    Q11=Q(1,1); Q12=Q(1,2); Q22=Q(2,2); Q66=Q(3,3);

    A=zeros(3,3); B=zeros(3,3); D=zeros(3,3);

    for k=1:n
        th = layup(k);
        c=cosd(th); s=sind(th);
        c2=c^2; s2=s^2; c3=c^3; s3=s^3; c4=c^4; s4=s^4;

        Q11b = Q11*c4 + Q22*s4 + 2*(Q12+2*Q66)*s2*c2;
        Q12b = (Q11+Q22-4*Q66)*s2*c2 + Q12*(c4+s4);
        Q22b = Q11*s4 + Q22*c4 + 2*(Q12+2*Q66)*s2*c2;
        Q16b = (Q11-Q12-2*Q66)*c3*s - (Q22-Q12-2*Q66)*c*s3;
        Q26b = (Q11-Q12-2*Q66)*c*s3 - (Q22-Q12-2*Q66)*c3*s;
        Q66b = (Q11+Q22-2*Q12-2*Q66)*s2*c2 + Q66*(s4+c4);

        Qbar=[Q11b Q12b Q16b; Q12b Q22b Q26b; Q16b Q26b Q66b];

        zk=z(k+1); zk1=z(k);

        A=A + Qbar*(zk-zk1);
        B=B + 0.5*Qbar*(zk^2 - zk1^2);
        D=D + (1/3)*Qbar*(zk^3 - zk1^3);
    end
end