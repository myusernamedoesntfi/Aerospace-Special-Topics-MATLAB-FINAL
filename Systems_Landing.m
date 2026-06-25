%% Systems - Aircraft Landing
% heavily modified homework 5

clear all; 
close all; 
clc;
 
fprintf('[INIT] calling Aircraft Landing ...\n');
 
fprintf('[ATMOS] Building standard atmosphere layers...\n');
 
hBase_km = [0 11 25 47 53 79 90 105];
a_lapse  = [-0.0065 0 0.003 0 -0.0045 0 0.004];
T0 = 288.16; P0 = 101325; R = 287; g_atm = 9.81;
 
[Ti, Pi, rhoi] = buildAtmosphere(hBase_km * 1000, a_lapse, T0, P0, R, g_atm);
 
fprintf('[ATMOS] Done. %d base heights, %d lapse rates.\n', length(hBase_km), length(a_lapse));
 
% there is a problem and I cant find it
for alt_test = [0, 5, 11, 20]
    [Tt, ~, rhot] = stdAtm(alt_test, struct('a',a_lapse,'hBase',hBase_km*1000,...
        'Ti',Ti,'Pi',Pi,'rhoi',rhoi,'R',R,'g',g_atm));
    fprintf('[ATMOS]   %3d km: T=%.1fK  rho=%.4f kg/m3\n', alt_test, Tt, rhot);
end
 
fprintf(' params are here...\n');
 
params_air.g       = 9.81;
params_air.mass    = 1000;
params_air.inertia = diag([300, 350, 280]);
params_air.c       = 2;
params_air.S       = 20;
 
% Aerodynamic derivatives
params_air.Cl0     = 0.2;    
params_air.Cl_alph = 1.2;   
params_air.Cl_dele = 0.8;
params_air.Cm0     = 0;      
params_air.Cm_alph = -0.5;  
params_air.Cm_dele = -1;
params_air.Cm_q    = -15;
params_air.Cd0     = 0.02;   
params_air.Cd_alph = 0.4;   
params_air.Cd_dele = 0.2;
 
% Crosswind / lateral derivatives
params_air.Cy_beta = -0.8;   
params_air.Cy_delr = 0.4;
params_air.Cn_beta =  0.15;  
params_air.Cn_delr = -0.1;  
params_air.Cn_r = -0.2;
 
% Engine & termination
params_air.T      = 0;            % coasting
params_air.Talt   = 100;          % stop at 100 m
params_air.V_wind = [0; -10; 0];  % 10 m/s lateral crosswind (inertial frame)
 
% Atmosphere sub-struct
params_air.atm.a     = a_lapse;
params_air.atm.hBase = hBase_km * 1000;
params_air.atm.Ti    = Ti;
params_air.atm.Pi    = Pi;
params_air.atm.rhoi  = rhoi;
params_air.atm.R     = R;
params_air.atm.g     = g_atm;
  
%  STEP 3: Initial conditions + EoM sanity check before ODE

% FIX: initialise theta at the desired descent angle so the controller
%      starts with zero error -- prevents the large initial torque spike
%      that forced ode45 into microsecond timesteps.
theta_init = deg2rad(-5);
X0_air = [0; 0; 10000; 150; 0; 0; 0; theta_init; 0; 0; 0; 0; 0; 0];
 
fprintf('[IC] Alt=%.0fm  u=%.1fm/s  theta=%.1fdeg\n', ...
    X0_air(3), X0_air(4), rad2deg(X0_air(8)));
 
fprintf('[DEBUG] Evaluating EoM at t=0 to check for NaN/Inf...\n');
Xdot0 = Flight_EoM_6DoF(0, X0_air, params_air);
if any(~isfinite(Xdot0))
    bad = find(~isfinite(Xdot0));
    fprintf('[ERROR] Xdot NAN: %s\n', mat2str(bad));
    error('EoM returned NaN/Inf at t=0. Fix before running ODE.');
else
    fprintf('[DEBUG] Xdot0 = [');
    fprintf(' %7.3f', Xdot0);
    fprintf(' ]\n');
    fprintf('[DEBUG] dz/dt = %.2f m/s  (negative = descending  OK)\n', Xdot0(3));
end

%  STEP 4: ODE solve
opts_air = odeset( ...
    'Events',    @(t,X) Flight_Termination(t, X, params_air), ...
    'RelTol',    1e-5, ...
    'AbsTol',    1e-5, ...
    'OutputFcn', @(t,X,flag) odeProgress(t, X, flag));
 
fprintf('[ODE] Starting ode45 (tspan 0..1000 s)...\n');
[t_air, X_air] = ode45( ...
    @(t,X) Flight_EoM_6DoF(t, X, params_air), ...
    [0 1000], X0_air, opts_air);
 
fprintf('\n[ODE] Done.  Steps=%d  t_final=%.1fs  alt_final=%.1fm\n', ...
    length(t_air), t_air(end), X_air(end,3));
 
%  STEP 5: Extract 2000 m conditions

idx_2000 = find(X_air(:,3) <= 2000, 1, 'first');
if isempty(idx_2000)
    fprintf('[POST] Aircraft did not reach 2000 m.\n');
else
    fprintf('\n--- Conditions at 2000 m ---\n');
    fprintf('  Time        : %.1f s\n',   t_air(idx_2000));
    fprintf('  Velocity u  : %.2f m/s\n', X_air(idx_2000, 4));
    fprintf('  Pitch angle : %.2f deg\n', rad2deg(X_air(idx_2000, 8)));
    fprintf('  Yaw angle   : %.2f deg\n', rad2deg(X_air(idx_2000, 9)));
    fprintf('  Lateral Y   : %.2f m\n',   X_air(idx_2000, 2));
end
 
%  STEP 6: Plots
fprintf('[PLOT] Generating figures...\n');
 
figure('Name','Aircraft Landing Analysis','Color','black');
subplot(3,1,1);
    plot3(X_air(:,1), X_air(:,2), X_air(:,3), 'b', 'LineWidth', 1.5);
    grid on;
    title('Flight Trajectory'); 
    xlabel('X (m)'); 
    ylabel('Y (m)'); 
    zlabel('Alt (m)');
subplot(3,1,2);
    plot(t_air, rad2deg(X_air(:,13)), 'LineWidth', 1.5);
    grid on; 
    yline(0,'k--'); 
    yline(15,'r:'); 
    yline(-15,'r:');
    title('Elevator Deflection'); 
    xlabel('Time (s)'); 
    ylabel('deg');
subplot(3,1,3);
    plot(t_air, rad2deg(X_air(:,14)), 'LineWidth', 1.5);
    grid on; yline(0,'k--'); 
    yline(20,'r:'); 
    yline(-20,'r:');
    title('Rudder Deflection'); % this one is goofy 
    xlabel('Time (s)'); 
    ylabel('deg');
 
figure('Name','Stability Angles','Color','black');
subplot(2,1,1);
    plot(t_air, rad2deg(X_air(:,8)), 'b', 'LineWidth', 1.5);
    hold on; 
    yline(-5, 'r--', 'Target -5 deg');
    grid on; 
    title('Pitch (\theta) vs Time'); 
    xlabel('Time (s)'); 
    ylabel('deg');

subplot(2,1,2);
    plot(t_air, rad2deg(X_air(:,9)), 'r', 'LineWidth', 1.5);
    grid on; 
    title('Yaw (\psi) vs Time'); 
    xlabel('Time (s)'); 
    ylabel('deg');
 
figure('Name','Lateral Tracking','Color','black');
subplot(2,1,1);
    plot(t_air, X_air(:,2), 'LineWidth', 1.5);
    yline(0,'k--','Runway CL');
    grid on; 
    title('Lateral Position Y'); 
    xlabel('Time (s)'); 
    ylabel('m');

subplot(2,1,2);
    plot(t_air, X_air(:,5), 'LineWidth', 1.5);
    grid on; 
    title('Lateral Velocity v'); 
    xlabel('Time (s)'); 
    ylabel('m/s');
 
fprintf('[DONE] All figures generated.\n');
 
%% Meat and potatoes (cause claude cooked here)
function Xdot = Flight_EoM_6DoF(~, X, params)
    % Unpack state vector
    Pos      = X(1:3);    % inertial position [x;y;z]  (m)
    Vel      = X(4:6);    % body-frame velocity [u;v;w] (m/s)
    Rot      = X(7:9);    % Euler angles [phi;theta;psi] (rad)
    angVel   = X(10:12);  % body angular rates [p;q;r]  (rad/s)
    elev_ang = X(13);     % actuator elevator state (rad)
    rud_ang  = X(14);     % the rudder
 
    phi   = Rot(1);
    theta = Rot(2);
    psi   = Rot(3);
 
    alt_m  = max(Pos(3), 0);
    alt_km = min(alt_m / 1000, params.atm.hBase(end) / 1000);
    [~, ~, rho] = stdAtm(alt_km, params.atm);
    rho = max(rho, 1e-6);
 
    Cbn          = Cbn_matrix(phi, theta, psi); 
    V_wind_body  = Cbn * params.V_wind;
    V_rel        = Vel - V_wind_body;
    Vmag         = max(norm(V_rel), 1e-3);
 
    alpha = atan2(V_rel(3), V_rel(1));
    beta  = asin(max(min(V_rel(2)/Vmag, 1.0), -1.0));
 
    theta_des = deg2rad(-5);
    dele_cmd  = -0.5 * (theta - theta_des) - 0.3 * angVel(2);
    delr_cmd  =  0.03 * Pos(2) + 0.08 * Vel(2) - 0.3 * angVel(3);
 
    dele_cmd = max(min(dele_cmd, deg2rad(15)), deg2rad(-15));
    delr_cmd = max(min(delr_cmd, deg2rad(20)), deg2rad(-20));
 

    Q  = 0.5 * rho * Vmag^2;
 
    Cl = params.Cl0 + params.Cl_alph * alpha      + params.Cl_dele * elev_ang;
    Cd = params.Cd0 + params.Cd_alph * abs(alpha) + params.Cd_dele * abs(elev_ang);
    Cy = params.Cy_beta * beta  + params.Cy_delr * rud_ang;
    Cm = params.Cm0 + params.Cm_alph * alpha + params.Cm_dele * elev_ang ...
         + params.Cm_q * (params.c / (2*Vmag)) * angVel(2);
    Cn = params.Cn_beta * beta + params.Cn_delr * rud_ang ...
         + params.Cn_r  * (params.c / (2*Vmag)) * angVel(3);
 
    L_aero = Q * params.S * Cl;
    D_aero = Q * params.S * Cd;
    Y_aero = Q * params.S * Cy;
 
    T = params.T;
    Fx_aero =  (T - D_aero) * cos(alpha) + L_aero * sin(alpha);
    Fy_aero =  Y_aero;
    Fz_aero =  (T - D_aero) * sin(alpha) - L_aero * cos(alpha);

    % Gravity vector in body frame  (standard: g acts in inertial -Z)
    g_body = params.g * [-sin(theta);
                          cos(theta)*sin(phi);
                          cos(theta)*cos(phi)];
 
    % Total acceleration in body frame (Coriolis term from rotating frame)
    Vdot = [Fx_aero; Fy_aero; Fz_aero] / params.mass ...
           + g_body ...
           - cross(angVel, Vel);
 
    I  = params.inertia;
 
    My = Q * params.S * params.c * Cm;
    Mz = Q * params.S * params.c * Cn;
    Mx = -(I(2,2) - I(3,3)) * angVel(2) * angVel(3);
    Mz = Mz - (I(1,1) - I(2,2)) * angVel(1) * angVel(2);
 
    angAccel = I \ [Mx; My; Mz];

    cos_theta = cos(theta);
    if abs(cos_theta) < 1e-4
        cos_theta = 1e-4 * sign(cos_theta + 1e-12);
    end
    tan_theta = sin(theta) / cos_theta;
 
    T_euler = [1,  sin(phi)*tan_theta,  cos(phi)*tan_theta;
               0,  cos(phi),           -sin(phi);
               0,  sin(phi)/cos_theta,  cos(phi)/cos_theta];
    Euler_dot = T_euler * angVel;
 
    Pos_dot = Cbn' * Vel;   % Cbn' = Cnb (body->inertial)
 
    elev_rate = (dele_cmd - elev_ang) / 0.2;
    rud_rate  = (delr_cmd - rud_ang)  / 0.2;
 
    Xdot = [Pos_dot; Vdot; Euler_dot; angAccel; elev_rate; rud_rate];
end
 
% goof
function R = Cbn_matrix(phi, theta, psi)
    cp = cos(phi);   sp = sin(phi);
    ct = cos(theta); st = sin(theta);
    cs = cos(psi);   ss = sin(psi);
    R = [ ct*cs,            ct*ss,           -st;
          sp*st*cs-cp*ss,   sp*st*ss+cp*cs,   sp*ct;
          cp*st*cs+sp*ss,   cp*st*ss-sp*cs,   cp*ct];
end
 
% kill attempt 1
function [v, t, d] = Flight_Termination(~, X, params)
    v = X(3) - params.Talt;
    t = 1;
    d = -1;
end
 
% rename it so i can know whats going on
function [Ti, Pi, rhoi] = buildAtmosphere(hBase, a, T0, P0, R, g)
    n    = length(a);
    npts = length(hBase);
    if npts ~= n + 1
        error('buildAtmosphere: length(hBase) must be length(a)+1. Got %d vs %d.', npts, n);
    end
    Ti   = zeros(1, npts);
    Pi   = zeros(1, npts);
    rhoi = zeros(1, npts);
    Ti(1)   = T0;
    Pi(1)   = P0;
    rhoi(1) = P0 / (R * T0);
    fprintf('[ATMOS]   Layer 1 (sea level): T=%.2fK  P=%.0fPa  rho=%.4f\n', Ti(1), Pi(1), rhoi(1));
    for i = 1:n
        dh = hBase(i+1) - hBase(i);
        if a(i) ~= 0
            Ti(i+1)   = max(Ti(i) + a(i)*dh, 1);
            Pi(i+1)   = Pi(i) * (Ti(i+1)/Ti(i))^(-g/(a(i)*R));
            rhoi(i+1) = Pi(i+1) / (R * Ti(i+1));
        else
            Ti(i+1)   = Ti(i);
            Pi(i+1)   = Pi(i) * exp(-g/(R*Ti(i)) * dh);
            rhoi(i+1) = Pi(i+1) / (R * Ti(i+1));
        end
        fprintf('[ATMOS]   Layer %d  h=%gkm: T=%.2fK  P=%.0fPa  rho=%.4f\n', ...
            i+1, hBase(i+1)/1000, Ti(i+1), Pi(i+1), rhoi(i+1));
    end
end
 

function [T, P, rho] = stdAtm(alt_km, atm)
    h     = max(alt_km*1000, atm.hBase(1));
    h     = min(h, atm.hBase(end));
    idx   = find(h >= atm.hBase, 1, 'last');
    if isempty(idx), idx = 1; end
    a_idx = min(idx, length(atm.a));
    dh    = h - atm.hBase(idx);
    if atm.a(a_idx) ~= 0
        T   = max(atm.Ti(idx) + atm.a(a_idx)*dh, 1);
        P   = atm.Pi(idx)   * (T/atm.Ti(idx))^(-atm.g/(atm.a(a_idx)*atm.R));
        rho = atm.rhoi(idx) * (T/atm.Ti(idx))^(-atm.g/(atm.a(a_idx)*atm.R)-1);
    else
        T   = atm.Ti(idx);
        P   = atm.Pi(idx)   * exp(-atm.g/(atm.R*atm.Ti(idx))*dh);
        rho = atm.rhoi(idx) * exp(-atm.g/(atm.R*atm.Ti(idx))*dh);
    end
end
 
% kill attempt 2
function status = odeProgress(t, X, flag)
    % kill ODE with fire
    persistent last_alt
    status = 0;
    switch flag
        case 'init'
            last_alt = X(3);        % FIX: seed with actual starting altitude
        case 'done'
            % nothing needed
        otherwise
            if ~isempty(t) && ~isempty(X)
                alt = X(3, end);
                % FIX: trigger when altitude decreases by 1000 m from last report
                if (last_alt - alt) >= 1000
                    fprintf('[ODE]   t=%6.1fs  alt=%6.0fm  u=%5.1fm/s  theta=%5.1fdeg  Y=%6.1fm\n', ...
                        t(end), alt, X(4,end), rad2deg(X(8,end)), X(2,end));
                    last_alt = alt;
                end
            end
    end
end
 
