% Structural Analysis
% last lecture

E  = 70e9;  % Modulus for Aluminum (Pa)
A1 = 0.100; % Orange internal truss elements (m^2)
A2 = 0.075; % Blue outer shell elements (m^2)

% Node positions [x, y] — row index = node number
nodes = [
    0.0,  0.000; % Node 1:  (pinned)
    0.5,  0.250; % Node 2:  
    0.5, -0.250; % Node 3:  
    0.5,  0.000; % Node 4:  
    1.0,  0.200; % Node 5
    1.0, -0.200; % Node 6
    1.5,  0.150; % Node 7
    1.5, -0.150; % Node 8
    2.0,  0.100; % Node 9
    2.0, -0.100; % Node 10
    2.5,  0.050; % Node 11
    2.5, -0.050; % Node 12
    3.0,  0.000  % Node 13: (Pinned)
];

elements = [
% Outer Shell (A2 = 0.075 m^2, Blue)
    1,  2, 2;   
    1,  3, 2;
    2,  5, 2;  
    3,  6, 2;
    5,  7, 2;   
    6,  8, 2;
    7,  9, 2;   
    8, 10, 2;
    9, 11, 2;  
    10, 12, 2;
    11, 13, 2;  
    12, 13, 2;
% Internal Truss (A1 = 0.1 m^2, Orange)
    2,  4, 1;   
    3,  4, 1;
    2,  6, 1;   
    3,  5, 1;
    4,  5, 1;   
    4,  6, 1; 
    5,  6, 1;   
    5,  8, 1;
    6,  7, 1;   
    7,  8, 1;
    7, 10, 1;   
    8,  9, 1;
    9, 10, 1;   
    9, 12, 1;
   10, 11, 1;  
   11, 12, 1;
];

% Stiffness Matrix

numNodes = size(nodes, 1);
DoF = 2 * numNodes;
K = zeros(DoF);

for i = 1:size(elements, 1)
    n1   = elements(i, 1);
    n2   = elements(i, 2);
    type = elements(i, 3);

    % Select cross-sectional area based on element type
    if type == 1
        A_val = A1;
    else
        A_val = A2;
    end

    x1 = nodes(n1, 1);  
    y1 = nodes(n1, 2);
    x2 = nodes(n2, 1);  
    y2 = nodes(n2, 2);

    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    c = (x2-x1) / L; % cos(theta)
    s = (y2-y1) / L; % sin(theta)

    % Local stiffness matrix for truss element
    k_loc = (E*A_val/L) * [ c^2,  c*s, -c^2, -c*s;
                             c*s,  s^2, -c*s, -s^2;
                            -c^2, -c*s,  c^2,  c*s;
                            -c*s, -s^2,  c*s,  s^2];

    % Global DoF indices for this element
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];

    % Assemble into global K
    K(dofs, dofs) = K(dofs, dofs) + k_loc;
end

% Boundary Conditions
fixDoF = [1, 2, 2*13-1, 2*13];
freDoF = setdiff(1:DoF, fixDoF);

total_lift = 5000; % (N)

F = zeros(DoF, 1);
top_nodes = [2, 5, 7, 9, 11]; 

for i = 1:length(top_nodes)
    F(2*top_nodes(i)) = total_lift / length(top_nodes); 
end

K_ff = K(freDoF, freDoF);
F_f  = F(freDoF);

u = zeros(DoF, 1);
u(freDoF) = K_ff \ F_f; 

% Diagnostic: check which DoFs have zero diagonal stiffness
zero_stiff = find(diag(K) == 0);
fprintf('Zero stiffness at global DoF indices: ');
disp(zero_stiff')

el_F = zeros(size(elements, 1), 1);

for i = 1:size(elements, 1)
    n1   = elements(i, 1);
    n2   = elements(i, 2);
    type = elements(i, 3);

    if type == 1
        A_val = A1;
    else
        A_val = A2;
    end

    x1 = nodes(n1, 1);  y1 = nodes(n1, 2);
    x2 = nodes(n2, 1);  y2 = nodes(n2, 2);

    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    c = (x2-x1) / L;
    s = (y2-y1) / L;

    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    u_e  = u(dofs); 

    T = [-c -s c s];
    el_F(i) = (E*A_val/L) * T * u_e;
end

fprintf('\nElement Load Table:\n');
fprintf('Elem | Nodes | Force (N)  | State\n');
fprintf('--------------------------------------\n');
for i = 1:length(el_F)
    if el_F(i) >= 0
        state = 'Tension';
    else
        state = 'Compression';
    end
    fprintf('%3d  | %2d-%2d  | %10.2f | %s\n', ...
        i, elements(i,1), elements(i,2), abs(el_F(i)), state);
end

Elemental_Forces = el_F

% plot
figure;
hold on;
title('Airfoil Truss Free Body Diagram');
xlabel('X (m)');
ylabel('Y (m)');

for i = 1:size(elements, 1)
    n1 = elements(i, 1);
    n2 = elements(i, 2);
    if elements(i, 3) == 1
        plot([nodes(n1,1), nodes(n2,1)], [nodes(n1,2), nodes(n2,2)], ...
            '-', 'Color', [0.85, 0.33, 0.1], 'LineWidth', 2); % orange internal
    else
        plot([nodes(n1,1), nodes(n2,1)], [nodes(n1,2), nodes(n2,2)], ...
            '-', 'Color', [0.1, 0.3, 0.5], 'LineWidth', 3);   % blue outer shell
    end
end

% Plot pinned supports at LE and TE
plot(nodes([1,13], 1), nodes([1,13], 2), 'k^', ...
    'MarkerFaceColor', 'g', 'MarkerSize', 10);

% Plot lift force arrows at top nodes
for i = 1:length(top_nodes)
    quiver(nodes(top_nodes(i), 1), nodes(top_nodes(i), 2), 0, 0.1, ...
        'r', 'LineWidth', 2, 'MaxHeadSize', 2);
end

% Plot node markers
plot(nodes(:,1), nodes(:,2), 'ro', 'MarkerSize', 8);

axis equal;

%% The hard part
% find which has the highest stress
% same as above with a lot of loops

rho_al  = 2700;   % kg/m^3 aluminum (google)
sigma_y = 276e6;  % Pa, 6061-T6 yield stress

triangles = [2, 3, 4, 5, 6, 7, 8];

results_nelems = zeros(size(triangles));
results_stress = zeros(size(triangles));
results_weight = zeros(size(triangles));

for b = 1:length(triangles)
    nb = triangles(b);
    x_bays = linspace(0, 3, nb+2);

    % Build node list explicitly into a growing matrix
    nd = zeros(1 + nb*3 + 1, 2);  % preallocate: LE + (top,bot,mid)*nb + TE
    nd(1,:) = [0, 0];              % Node 1: LE
    row = 2;
    for k = 2:nb+1
        xi = x_bays(k);
        if xi <= 0.5
            hy = 0.5 * xi;
        else
            hy = 0.25 - 0.1*(xi - 0.5);
        end
        nd(row,  :) = [xi,  hy];   % top
        nd(row+1,:) = [xi, -hy];   % bottom
        nd(row+2,:) = [xi,  0 ];   % mid
        row = row + 3;
    end
    nd(end,:) = [3, 0];            % TE

    nN = size(nd, 1);

    top_idx = @(k) 2 + (k-1)*3;
    bot_idx = @(k) 3 + (k-1)*3;
    mid_idx = @(k) 4 + (k-1)*3;

    % Build element list
    el = [1, top_idx(1), 2;
          1, bot_idx(1), 2];

    for k = 1:nb-1
        t1=top_idx(k); b1=bot_idx(k); m1=mid_idx(k);
        t2=top_idx(k+1); b2=bot_idx(k+1);
        el = [el;
              t1, t2, 2;    % outer shell top
              b1, b2, 2;    % outer shell bot
              t1, m1, 1;    % internal
              b1, m1, 1;
              m1, t2, 1;
              m1, b2, 1;
              t1, b2, 1;
              b1, t2, 1;
              t1, b1, 1];
    end

    % Last bay to TE
    k = nb;
    el = [el;
          top_idx(k), nN,        2;
          bot_idx(k), nN,        2;
          top_idx(k), mid_idx(k), 1;
          bot_idx(k), mid_idx(k), 1;
          mid_idx(k), nN,        1];

    nE     = size(el, 1);
    DoF_b  = 2 * nN;
    K_b    = zeros(DoF_b);

    % Assemble stiffness
    for i = 1:nE
        n1 = el(i,1); n2 = el(i,2);
        A_val = A1;
        
        if el(i,3)==2, A_val = A2; end
        x1=nd(n1,1); 
        y1=nd(n1,2);
        x2=nd(n2,1); 
        y2=nd(n2,2);
      
        L = sqrt((x2-x1)^2+(y2-y1)^2);
        c=(x2-x1)/L; s=(y2-y1)/L;
       
        k_loc=(E*A_val/L)*[c^2,c*s,-c^2,-c*s;
                            c*s,s^2,-c*s,-s^2;
                           -c^2,-c*s,c^2,c*s;
                           -c*s,-s^2,c*s,s^2];
        
        dofs=[2*n1-1,2*n1,2*n2-1,2*n2];
        K_b(dofs,dofs)=K_b(dofs,dofs)+k_loc;
    end

    fixDoF_b = [1, 2, 2*nN-1, 2*nN];
    freDoF_b = setdiff(1:DoF_b, fixDoF_b);

    top_nodes_b = top_idx(1:nb);   
    F_b = zeros(DoF_b, 1);
   
    for i = 1:length(top_nodes_b)
        F_b(2*top_nodes_b(i)) = total_lift / length(top_nodes_b);
    end

    u_b = zeros(DoF_b,1);
    u_b(freDoF_b) = K_b(freDoF_b,freDoF_b) \ F_b(freDoF_b);

    % Recover forces and stresses
    sigma_max = 0;
    total_vol = 0;
    for i = 1:nE
        n1=el(i,1); 
        n2=el(i,2);
        A_val=A1; 
        
        if el(i,3)==2, A_val=A2; 
        end
        x1=nd(n1,1); 
        y1=nd(n1,2);
        x2=nd(n2,1); 
        y2=nd(n2,2);
        
        L=sqrt((x2-x1)^2+(y2-y1)^2);
        c=(x2-x1)/L; s=(y2-y1)/L;
        
        dofs=[2*n1-1,2*n1,2*n2-1,2*n2];
        F_elem=(E*A_val/L)*[-c -s c s]*u_b(dofs);
        
        sigma=abs(F_elem)/A_val;
        
        if sigma > sigma_max
            sigma_max = sigma;
        end
       
        total_vol = total_vol + A_val*L;
    end

    results_nelems(b) = nE;
    results_stress(b) = sigma_max / 1e6;
    results_weight(b) = total_vol * rho_al;

    fprintf('bays=%d  elems=%d  sigma_max=%.4f MPa  weight=%.1f kg\n', ...
        nb, nE, results_stress(b), results_weight(b));
end

% Plot
figure;
yyaxis left
plot(results_nelems, results_stress, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
yline(sigma_y/1e6, 'r--', 'Yield (276 MPa)', 'LineWidth', 1.5);
ylabel('Max Element Stress (MPa)');

yyaxis right
plot(results_nelems, results_weight, 'k^--', 'LineWidth', 1.5, 'MarkerSize', 7);
ylabel('Structural Weight (kg)');

xlabel('Number of Elements');
title('Truss Design: Stress and Weight vs. Element Count');
grid on;

[~, cur_idx] = min(abs(triangles - 6));
hold on;
yyaxis left
plot(results_nelems(cur_idx), results_stress(cur_idx), 'r*', 'MarkerSize', 14);
legend('Max Stress','Yield Limit','Structural Weight','Current Design','Location','best');