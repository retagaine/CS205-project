function [pot] = get_pot(x,locs,chg) 
% define elementary charge (e), dielectric constant (er), and permittivity
% of free space (e0)
e = 1.60217662e-19;
er = sqrt(10.03*9.66);
e0 = 8.854187817e-12;

% compute distances between particles and then the potential
inv = sqrt(sum([(x(1) - locs(1,1,:)).^2,(x(2) - locs(1,2,:)).^2,(x(3) - locs(1,3,:)).^2]))+eps;
pot = e/4/pi/e0/er*sum(chg./inv);