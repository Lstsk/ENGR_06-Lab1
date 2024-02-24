% E6 Lab 1, 2/14/2024; William Grenville, Carson Lin, Victoria Mann,
% Charlie Schuetz, Katherine Wang
clear all
close all
clc
%% approximated scenario
% Pulley location
Pulley2 = [-95.5, 64, 116.61] % inches from origin
Pulley3 = [-95.75, -68, 116.96]
Pulley4 = [94.25, -7.25, 116.73]
Centerweight = [2, -10.5, 50.5; -36.5, -12.75, 27; -24, -12, 35.5]

% Distance vectors
DistVectors_Pulley2 = [Pulley2 - Centerweight(1,:); Pulley2 - Centerweight(2,:); Pulley2 - Centerweight(3,:)] % this will give us the distance vector for Pulley 2 in each trial as its rows
DistVectors_Pulley3 = [Pulley3 - Centerweight(1,:); Pulley3 - Centerweight(2,:); Pulley3 - Centerweight(3,:)]
DistVectors_Pulley4 = [Pulley4 - Centerweight(1,:); Pulley4 - Centerweight(2,:); Pulley4 - Centerweight(3,:)]
% transpose because sum is weird
DistVectors_Pulley2 = DistVectors_Pulley2';
DistVectors_Pulley3 = DistVectors_Pulley3';
DistVectors_Pulley4 = DistVectors_Pulley4';

% Magnitude of the distance vectors
Mag_DistVectors_Pulley2 = [(sum(DistVectors_Pulley2.^2).^(1/2))];
Mag_DistVectors_Pulley3 = [(sum(DistVectors_Pulley3.^2).^(1/2))];
Mag_DistVectors_Pulley4 = [(sum(DistVectors_Pulley4.^2).^(1/2))];
% transpose back
DistVectors_Pulley2 = DistVectors_Pulley2';
DistVectors_Pulley3 = DistVectors_Pulley3';
DistVectors_Pulley4 = DistVectors_Pulley4';

% Lambda of the distance vectors
lambdasPulley2 = [DistVectors_Pulley2(1,:)/Mag_DistVectors_Pulley2(1);DistVectors_Pulley2(2,:)/Mag_DistVectors_Pulley2(2);DistVectors_Pulley2(3,:)/Mag_DistVectors_Pulley2(3)]
lambdasPulley3 = [DistVectors_Pulley3(1,:)/Mag_DistVectors_Pulley3(1);DistVectors_Pulley3(2,:)/Mag_DistVectors_Pulley3(2);DistVectors_Pulley3(3,:)/Mag_DistVectors_Pulley3(3)]
lambdasPulley4 = [DistVectors_Pulley4(1,:)/Mag_DistVectors_Pulley4(1);DistVectors_Pulley4(2,:)/Mag_DistVectors_Pulley4(2);DistVectors_Pulley4(3,:)/Mag_DistVectors_Pulley4(3)]

% Tension
syms T2 T3 T4
rectx = lambdasPulley2 * T2
recty = lambdasPulley3 * T3
rectz = lambdasPulley4 * T4
for i = 1:3
  Rx =  rectx(i,1) + recty(i,1) + rectz(i,1) == 0
  Ry =  rectx(i,2) + recty(i,2) + rectz(i,2) == 0
  Rz =  rectx(i,3) + recty(i,3) + rectz(i,3) == 3.02
  soln = solve(Rx, Ry, Rz, T2, T3, T4)
  trial = [double(soln.T2), double(soln.T3),double(soln.T4)]
  Tension(i,:) = trial % pounds; [trial 1; trial 2; trial 3]
end

%% non-approximated scenario
% pulley coordinates
Pulleys = [Pulley2; Pulley3; Pulley4]
r_pulley = 0.75 + 0.1875 % in
Pulley2_xy = [-95.5, 64]; % inches from origin
Pulley3_xy = [-95.75, -68];
Pulley4_xy = [94.25, -7.25];
Centerweight = [2, -10.5, 50.5; -36.5, -12.75, 27; -24, -12, 35.5] % [trial 1; trial 2; trial 3]

% pulley 2
AC2 =  [Pulley2_xy; Pulley2_xy; Pulley2_xy] - [Centerweight(:,1), Centerweight(:,2)] % [AC pulley2 trial 1; AC pulley2 trial 2; AC pulley2 trial 3]
AC2_mag = [sqrt(sum(AC2(1,:).^2)); sqrt(sum(AC2(2,:).^2)); sqrt(sum(AC2(3,:).^2))]
AE2_mag = AC2_mag - r_pulley
for i = 1:3
   alpha2 = atand(Pulley2(:,3)./AE2_mag(i))
   alpha2_array(i) = alpha2 % [pulley2 trial 1, pulley2 trial 2, pulley2 trial 3]
end
AF2_mag = sqrt(AE2_mag.^2 + (Pulley2(:,3)).^2) % [AF pulley2 trial 1; AF pulley2 trial 2; AF pulley2 trial 3]
for i = 1:3
   beta2 = asind(r_pulley/AF2_mag(i))
   beta2_array(i) = beta2 % [pulley2 trial 1, pulley2 trial 2, pulley2 trial 3]
end

for i = 1:3
   Pulley2_new_z = AC2_mag(i)*tand(alpha2_array(i) + beta2_array(i))
   Pulley2_new_z_array(i) = Pulley2_new_z  % [pulley2 trial 1, pulley2 trial 2, pulley2 trial 3]
end

Pulley2_new = [Pulley2(1,1), Pulley2(1,2), Pulley2_new_z_array(1); Pulley2(1,1), Pulley2(1,2), Pulley2_new_z_array(2); Pulley2(1,1), Pulley2(1,2), Pulley2_new_z_array(3)]
DistVectors_Pulley2_new = Pulley2_new - Centerweight
DistVectors_Pulley2_new = DistVectors_Pulley2_new'
Mag_DistVectors_Pulley2_new = [(sum(DistVectors_Pulley2_new.^2).^(1/2))] % [trial 1, trial 2, trial 3]
DistVectors_Pulley2_new = DistVectors_Pulley2_new'
lambdasPulley2_new = [DistVectors_Pulley2_new(1,:)/Mag_DistVectors_Pulley2_new(1);DistVectors_Pulley2_new(2,:)/Mag_DistVectors_Pulley2_new(2);DistVectors_Pulley2_new(3,:)/Mag_DistVectors_Pulley2_new(3)]

% pulley 3
AC3 =  [Pulley3_xy; Pulley3_xy; Pulley3_xy] - [Centerweight(:,1), Centerweight(:,2)] % [AC pulley3 trial 1; AC pulley3 trial 2; AC pulley3 trial 3]
AC3_mag = [sqrt(sum(AC3(1,:).^2)); sqrt(sum(AC3(2,:).^2)); sqrt(sum(AC3(3,:).^2))]
AE3_mag = AC3_mag - r_pulley
for i = 1:3
   alpha3 = atand(Pulley3(:,3)./AE3_mag(i))
   alpha3_array(i) = alpha3 % [pulley3 trial 1, pulley3 trial 2, pulley3 trial 3]
end

AF3_mag = sqrt(AE3_mag.^2 + (Pulley3(:,3)).^2) % [AF pulley3 trial 1; pulley3 trial 2; AF pulley3 trial 3]
for i = 1:3
   beta3 = asind(r_pulley/AF3_mag(i))
   beta3_array(i) = beta3 % [pulley3 trial 1, pulley3 trial 2, pulley3 trial 3]
end

for i = 1:3
   Pulley3_new_z = AC3_mag(i)*tand(alpha3_array(i) + beta3_array(i))
   Pulley3_new_z_array(i) = Pulley3_new_z  % [pulley3 trial 1, pulley3 trial 2, pulley3 trial 3]
end

Pulley3_new = [Pulley3(1,1), Pulley3(1,2), Pulley3_new_z_array(1); Pulley3(1,1), Pulley3(1,2), Pulley3_new_z_array(2); Pulley3(1,1), Pulley3(1,2), Pulley3_new_z_array(3)]
DistVectors_Pulley3_new = Pulley3_new - Centerweight
DistVectors_Pulley3_new = DistVectors_Pulley3_new'
Mag_DistVectors_Pulley3_new = [(sum(DistVectors_Pulley3_new.^2).^(1/2))] % [trial 1, trial 2, trial 3]
DistVectors_Pulley3_new = DistVectors_Pulley3_new'
lambdasPulley3_new = [DistVectors_Pulley3_new(1,:)/Mag_DistVectors_Pulley3_new(1);DistVectors_Pulley3_new(2,:)/Mag_DistVectors_Pulley3_new(2);DistVectors_Pulley3_new(3,:)/Mag_DistVectors_Pulley3_new(3)]

% pulley 4
AC4 =  [Pulley4_xy; Pulley4_xy; Pulley4_xy] - [Centerweight(:,1), Centerweight(:,2)] % [AC pulley4 trial 1; AC pulley4 trial 2; AC pulley4 trial 3]
AC4_mag = [sqrt(sum(AC4(1,:).^2)); sqrt(sum(AC4(2,:).^2)); sqrt(sum(AC4(3,:).^2))]
AE4_mag = AC4_mag - r_pulley
for i = 1:3
   alpha4 = atand(Pulley4(:,3)./AE4_mag(i))
   alpha4_array(i) = alpha4 % [pulley4 trial 1, pulley4 trial 2, pulley4 trial 3]
end

AF4_mag = sqrt(AE4_mag.^2 + (Pulley4(:,3)).^2) % [AF pulley4 trial 1; pulley4 trial 2; AF pulley4 trial 3]
for i = 1:3
   beta4 = asind(r_pulley/AF4_mag(i))
   beta4_array(i) = beta4 % [pulley4 trial 1, pulley4 trial 2, pulley4 trial 3]
end

for i = 1:3
   Pulley4_new_z = AC4_mag(i)*tand(alpha4_array(i) + beta4_array(i))
   Pulley4_new_z_array(i) = Pulley4_new_z  % [pulley4 trial 1, pulley4 trial 2, pulley4 trial 3]
end

Pulley4_new = [Pulley4(1,1), Pulley4(1,2), Pulley4_new_z_array(1); Pulley4(1,1), Pulley4(1,2), Pulley4_new_z_array(2); Pulley4(1,1), Pulley4(1,2), Pulley4_new_z_array(3)]
DistVectors_Pulley4_new = Pulley4_new - Centerweight
DistVectors_Pulley4_new = DistVectors_Pulley4_new'
Mag_DistVectors_Pulley4_new = [(sum(DistVectors_Pulley4_new.^2).^(1/2))] % [trial 1, trial 2, trial 3]
DistVectors_Pulley4_new = DistVectors_Pulley4_new'
lambdasPulley4_new = [DistVectors_Pulley4_new(1,:)/Mag_DistVectors_Pulley4_new(1);DistVectors_Pulley4_new(2,:)/Mag_DistVectors_Pulley4_new(2);DistVectors_Pulley4_new(3,:)/Mag_DistVectors_Pulley4_new(3)]

% Tension
syms T2_new T3_new T4_new
rectx_new = lambdasPulley2_new * T2_new
recty_new = lambdasPulley3_new * T3_new
rectz_new = lambdasPulley4_new * T4_new
for i = 1:3
  Rx_new =  rectx_new(i,1) + recty_new(i,1) + rectz_new(i,1) == 0
  Ry_new =  rectx_new(i,2) + recty_new(i,2) + rectz_new(i,2) == 0
  Rz_new =  rectx_new(i,3) + recty_new(i,3) + rectz_new(i,3) == 3.02
  soln_new = solve(Rx_new, Ry_new, Rz_new)
  trial_new = [double(soln_new.T2_new), double(soln_new.T3_new), double(soln_new.T4_new)]
  Tension_new(i,:) = trial_new % pounds; [trial 1; trial 2; trial 3]
end
