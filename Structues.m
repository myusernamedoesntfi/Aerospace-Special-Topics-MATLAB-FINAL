%% Structure

clear all
close all
clc

LHS = 1e3 .* [-5.333;-8;5.333]; % forces and moments defined on left hand side
K = 1e6 .* [42, -7.875, 10.5;
-7.875, 3.9375, -7.875;
10.5, -7.875, 21]; % stiffness matrix

syms PhB vC PhC

disp_vars = [PhB; vC; PhC]; % disp is not Disp

RHS = mtimes(K,disp_vars);

eq1 = RHS(1,:) == LHS(1);
eq2 = RHS(2,:) == LHS(2);
eq3 = RHS(3,:) == LHS(3);

% Solve system
[Y,Z] = equationsToMatrix([eq1 eq2 eq3],[PhB vC PhC]);
X = linsolve(Y, Z);

Disp_Results = vpa(X); % [PhiB (rad); vC (m); PhiC (rad)]

E = 210e9; % Modulus (Pa)
A = 1e-4;  % cross-sectional area (m^2)

% I am guessing on these coordinates Im not really sure how else
% youre supposed to come up with these
nodes = [
    0.0,  0.000; % Node 1: Leading edge tip
    0.25, 0.125; % Node 2: Top inner 
    0.25,-0.125; % Node 3: Bottom inner
    0.5,  0.250; % Node 4: Top max thickness
    0.5,  0.000; % Node 5: Mid max thickness
    0.5, -0.250; % Node 6: Bottom max thickness
    1.0,  0.200; % Node 7
    1.0, -0.200; % Node 8
    1.5,  0.150; % Node 9
    1.5, -0.150; % Node 10
    2.0,  0.100; % Node 11
    2.0, -0.100; % Node 12
    2.5,  0.050; % Node 13
    2.5, -0.050; % Node 14
    3.0,  0.000  % Node 15: Trailing edge tip
];

% Elements (connecty lines)
elements = [
    1,2; 1,3; 2,3; % Front triangle
    2,4; 3,6; 2,5; 3,5; 4,5; 5,6; % second triangle transition
    4,7; 6,8; 4,8; 6,7; 7,8; % up and lower diagonals and verticals
    7,9; 8,10; 7,10; 8,9; 9,10; % next midle diamond
    9,11; 10,12; 9,12; 10,11; 11,12; % next one
    11,13; 12,14; 11,14; 12,13; 13,14; % same deal but the skinnyiest one
    13,15; 14,15 % the tip on the end
];

numNodes = size(nodes,1);
DoF = 2 * numNodes;

K = zeros(DoF);

for i = 1:size(elements, 1)
    n1 = elements(i,1);
    n2 = elements(i,2);
    
    x1 = nodes(n1,1); 
    y1 = nodes(n1,2);
    x2 = nodes(n2,1); 
    y2 = nodes(n2,2);
    
    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    c = (x2-x1)/L;
    s = (y2-y1)/L;
    
    % Local stiffness matrix
    k_loc = (E*A/L) * [ c^2,  c*s, -c^2, -c*s;
                        c*s,  s^2, -c*s, -s^2;
                       -c^2, -c*s,  c^2,  c*s;
                       -c*s, -s^2,  c*s,  s^2];
                   
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    K(dofs, dofs) = K(dofs, dofs) + k_loc;
end

%Force Array
F = zeros(DoF, 1);
F(2*13) = -1000; % I dont thing this is at node 4 anymore

% This is probably not the only fixed part
fixDoF = [2,3,4]; 

freDoF = setdiff(1:DoF, fixDoF);

K_ff = K(freDoF, freDoF);
F_f = F(freDoF);

u = zeros(DoF, 1);
u(freDoF) = K_ff\F_f;

el_F = zeros(size(elements, 1), 1);

% plot it
figure('Color', 'black'); 
hold on; 
grid on;
title('Airfoil Truss Structural Analysis');
xlabel('X (m)'); 
ylabel('Y (m)');

for i = 1:size(elements, 1)
    n1 = elements(i,1); 
    n2 = elements(i,2);
    x = nodes([n1, n2], 1); 
    y = nodes([n1, n2], 2);
    
    % Calculate Element Force
    L = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    c = (x(2)-x(1))/L;
    s = (y(2)-y(1))/L;
   
    dofs = [2*n1-1,2*n1,2*n2-1,2*n2];
    u_e = u([2*n1-1, 2*n1, 2*n2-1, 2*n2]);
    
    T = [-c -s c s];
    el_F(i) = (E*A/L)*T*u_e;

    % Plot element
    plot(x, y, 'k-', 'LineWidth', 1.5);

end

plot(nodes(:,1), nodes(:,2), 'ro', 'MarkerFaceColor', 'r');
quiver(nodes(4,1), nodes(4,2), 0, -0.2, 'b', 'LineWidth', 2, 'MaxHeadSize', 1);
axis equal;

Elemental_Forces = el_F;