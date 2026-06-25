% Fluids - Flow over Diamond Airfoil
% homework 4 diamond part

LE  = [0, 0];
Top = [0.5, 0.25];
TE  = [3, 0];
Bot = [0.5, -0.25];
 
N_vortices = 30;
Nstep      = 100;
 
x_grid = linspace(-1, 4, Nstep);
y_grid = linspace(-1, 1, Nstep);
[X, Y] = meshgrid(x_grid, y_grid);
 
x_poly = [LE(1), Top(1), TE(1), Bot(1), LE(1)];
y_poly = [LE(2), Top(2), TE(2), Bot(2), LE(2)];
 
U_inf_mag = 155;         % m/s 
alpha     = deg2rad(2);  % rad 
 
rho      = 1.006;   % kg/m3  air density at 2000m from standard atmosphere
rho_2000 = rho;     
 
U_inf_x = U_inf_mag * cos(alpha);
U_inf_y = U_inf_mag * sin(alpha);
 
fprintf('Free stream velocity : %.2f m/s\n', U_inf_mag);
fprintf('Angle of attack      : %.2f deg\n', rad2deg(alpha));
fprintf('Air density at 2000m : %.4f kg/m3\n', rho_2000);
 
xlocations = linspace(LE(1)+0.01, TE(1)-0.01, N_vortices);
y_vort     = zeros(1, N_vortices);
thickness  = zeros(1, N_vortices);
 
for i = 1:length(xlocations)
    xi = xlocations(i);
    if xi <= 0.5
        top_y = 0.5 * xi;
    else
        top_y = 0.25 - 0.1 * (xi - 0.5);
    end
    thickness(i) = 2 * top_y;
end
 
st_base = 15;
st      = st_base * (thickness ./ max(thickness));
st = st + 10 * alpha;
 
% Compute velocities
U_vort = zeros(Nstep, Nstep);
V_vort = zeros(Nstep, Nstep);
 
for i = 1:N_vortices
    X_shift = X - xlocations(i);
    Y_shift = Y - y_vort(i);
    R2      = X_shift.^2 + Y_shift.^2 + 1e-6;
    U_vort  = U_vort + (st(i)/(2*pi)) * ( Y_shift ./ R2);
    V_vort  = V_vort + (st(i)/(2*pi)) * (-X_shift ./ R2);
end
 
U_tot = U_vort + U_inf_x;
V_tot = V_vort + U_inf_y;
 
% layers are here
inside        = inpolygon(X, Y, x_poly, y_poly);
U_tot(inside) = NaN;
V_tot(inside) = NaN;
 
V_mag      = sqrt(U_tot.^2 + V_tot.^2);
lift_force = 0;
dx         = 3.0 / N_vortices;
 
for i = 1:N_vortices
    xi = xlocations(i);
    if xi <= 0.5
        y_surf = 0.5 * xi;
    else
        y_surf = 0.25 - 0.1 * (xi - 0.5);
    end
 
    % Find nearest grid point above upper surface and below lower surface
    V_top = interp2(X, Y, V_mag, xi,  y_surf + 0.05);
    V_bot = interp2(X, Y, V_mag, xi, -y_surf - 0.05);
 
    if ~isnan(V_top) && ~isnan(V_bot)
        % P + 0.5*rho*V^2 = const  ->  dP = 0.5*rho*(V_top^2 - V_bot^2)
        dP         = 0.5 * rho * (V_top^2 - V_bot^2);
        lift_force = lift_force + dP * dx;
    end
end
 
fprintf('Calculated Lift Force: %.2f N/m\n', lift_force);
 
% Plot Vector Field
figure('Name','Airfoil Flow','Color','black');
quiver(X, Y, U_tot, V_tot, 1.5, 'b');
hold on;
plot(x_poly, y_poly, 'k-', 'LineWidth', 2);
axis equal;
xlim([-0.5, 3.5]); ylim([-0.5, 0.5]);
grid on;
xlabel('x (m)'); 
ylabel('y (m)');
title(sprintf('Flow Field at 2000m  (V = %.1f m/s,  \\alpha = %.1f deg)', ...
    U_inf_mag, rad2deg(alpha)));
 
Cl_6DoF            = 0.2 + 1.2 * alpha;                              % Cl0 + Cl_alph*alpha from params_air
Lift_6DoF_per_span = 0.5 * rho_2000 * U_inf_mag^2 * (2.0) * Cl_6DoF; % c = 2.0 m
 
fprintf('Reach Idea: 6DoF Lift = %.2f N/m | Fluids Lift = %.2f N/m\n', ...
    Lift_6DoF_per_span, lift_force);
