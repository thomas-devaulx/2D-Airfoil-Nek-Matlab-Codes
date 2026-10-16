%% NACA 4-Digit Airfoil & Nek5000 .geo Generator
clear; clc;

% User Inputs
NACA  = '0012';     
aoa   = 10;         
c_f   = 0.25;       
N_p   = 200;        
t_max = str2double(NACA(3:4))/100;

R_domain = 5.0;    
X_wake   = 40.0;    

% Resolution Targets 
lc_far   = 0.5;     % outer boundary
lc_wake  = 0.03;    % wake
lc_surf  = 0.015;   % near-surface resolution
lc_LE    = 0.005;   % boundary layer
lc_TE    = 0.005;   % boundary layer

% Airfoil Math & Rotation
m = str2double(NACA(1))/100;
p = str2double(NACA(2))/10;

k = linspace(0, 1, N_p)';
x = k.^2;   % squared spacing, more points to define high curvature at leading edge

yt = 5 * t_max * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1036*x.^4);

yc = zeros(size(x)); dyc_dx = zeros(size(x));
if m > 0
    idx1 = x <= p; idx2 = x > p;
    yc(idx1) = (m/p^2) * (2*p*x(idx1) - x(idx1).^2);
    dyc_dx(idx1) = (2*m/p^2) * (p - x(idx1));
    yc(idx2) = (m/(1-p)^2) * ((1-2*p) + 2*p*x(idx2) - x(idx2).^2);
    dyc_dx(idx2) = (2*m/(1-p)^2) * (p - x(idx2));
end
theta = atan(dyc_dx);

xu = x - yt.*sin(theta); yu = yc + yt.*cos(theta);
xl = x + yt.*sin(theta); yl = yc - yt.*cos(theta);

A = -aoa * pi/180;
R = [cos(A), -sin(A); sin(A), cos(A)];

coords_u = R * [(xu - c_f)'; yu'];
coords_l = R * [(xl - c_f)'; yl'];

% Split Top and Bottom Surface
up_x = coords_u(1, :);
up_y = coords_u(2, :);
low_x = fliplr(coords_l(1, 2:end-1)); 
low_y = fliplr(coords_l(2, 2:end-1));

x_airfoil = [up_x, low_x];
y_airfoil = [up_y, low_y];

LE_ID = 1;
TE_ID = length(up_x); 
% plot
close all; clc;
gemColors = orderedcolors("gem"); 
figure('Position',[100, 100, 500, 150]);
plot(up_x+0.25,up_y,'Color',gemColors(3,:),'LineWidth',1.5); hold on;
plot(low_x+0.25,low_y,'Color',gemColors(3,:),'LineWidth',1.5); 
plot(up_x+0.25,yc,'--k','LineWidth',1.5);
xlim([0 1]); ylim([-0.2 0.2]); axis equal; grid on;
title('NACA 0012'); xlabel('x/c'); ylabel('y/c');

%% 3. WRITE GMSH .GEO FILE
filename = sprintf('NACA%s_A%d_BL.geo', NACA, aoa);
fid = fopen(filename, 'w');

fprintf(fid, '// Nek5000 Physics-Aware C-Grid\n');
fprintf(fid, 'SetFactory("OpenCASCADE");\n\n');

% generates the points for the airfoil, with more or less resolution near certain regions (TE & LE)
fprintf(fid, '// --- Airfoil Geometry (Split Topology) ---\n');
for i = 1:length(x_airfoil)
    if x_airfoil(i) < -0.20
        lc_pt = lc_LE/3;
    else
        lc_pt = lc_TE*3;
    end
    fprintf(fid, 'Point(%d) = {%f, %f, 0, %f};\n', i, x_airfoil(i), y_airfoil(i), lc_pt);
end

% Links the points together to create the airfoil surface
fprintf(fid, 'Spline(1) = {1:%d};\n', TE_ID);                            
fprintf(fid, 'Spline(2) = {%d, %d:%d, 1};\n\n', TE_ID, TE_ID+1, length(x_airfoil)); 

% Creates outer domain
fprintf(fid, '// --- Outer Domain ---\n');
pt_out = length(x_airfoil) + 1;
fprintf(fid, 'Point(%d) = {%f, %f, 0, %f};\n', pt_out,   X_wake,  R_domain, lc_far); 
fprintf(fid, 'Point(%d) = {%f, %f, 0, %f};\n', pt_out+1, 0,       R_domain, lc_far); 
fprintf(fid, 'Point(%d) = {%f, %f, 0, %f};\n', pt_out+2, -R_domain, 0,      lc_far); 
fprintf(fid, 'Point(%d) = {%f, %f, 0, %f};\n', pt_out+3, 0,      -R_domain, lc_far); 
fprintf(fid, 'Point(%d) = {%f, %f, 0, %f};\n', pt_out+4, X_wake, -R_domain, lc_far); 
fprintf(fid, 'Point(%d) = {0, 0, 0, %f};\n', pt_out+5, lc_far); 

fprintf(fid, 'Line(3) = {%d, %d};\n', pt_out, pt_out+1);           
fprintf(fid, 'Circle(4) = {%d, %d, %d};\n', pt_out+1, pt_out+5, pt_out+2); 
fprintf(fid, 'Circle(5) = {%d, %d, %d};\n', pt_out+2, pt_out+5, pt_out+3); 
fprintf(fid, 'Line(6) = {%d, %d};\n', pt_out+3, pt_out+4);         
fprintf(fid, 'Line(7) = {%d, %d};\n', pt_out+4, pt_out);           

% Links the lines together to make the 2D surface
fprintf(fid, 'Curve Loop(1) = {3, 4, 5, 6, 7};\n'); 
fprintf(fid, 'Curve Loop(2) = {1, 2};\n');          
fprintf(fid, 'Plane Surface(1) = {1, 2};\n\n');

fprintf(fid, '// --- Refinement Fields ---\n');
fprintf(fid, '// Long Box\n'); % controls refinement of the far wake
fprintf(fid, 'Field[1] = Box;\n');
fprintf(fid, 'Field[1].VIn = %f; Field[1].VOut = %f;\n', lc_wake, lc_far);
fprintf(fid, 'Field[1].XMin = -0.55; Field[1].XMax = %f;\n', X_wake);
fprintf(fid, 'Field[1].YMin = -1.0; Field[1].YMax = 1.0;\n');
fprintf(fid, 'Field[1].ZMin = -1.0; Field[1].ZMax = 1.0;\n');
fprintf(fid, 'Field[1].Thickness = 1.0;\n\n'); 

fprintf(fid, '// Short Box\n'); % controls refinement of the near-wake
fprintf(fid, 'Field[2] = Box;\n');
fprintf(fid, 'Field[2].VIn = %f; Field[2].VOut = %f;\n', lc_surf, lc_far);
fprintf(fid, 'Field[2].XMin = -0.55; Field[2].XMax = 3.75;\n');
fprintf(fid, 'Field[2].YMin = -0.25; Field[2].YMax = 0.25;\n');
fprintf(fid, 'Field[2].ZMin = -1.0; Field[2].ZMax = 1.0;\n');
fprintf(fid, 'Field[2].Thickness = 0.5;\n\n'); 

fprintf(fid, '// Boundary Layer Field\n'); % creates a boundary layer surrounding the airfoil, fanned at TE
fprintf(fid, 'Field[3] = BoundaryLayer;\n');
fprintf(fid, 'Field[3].CurvesList = {1, 2};\n');
fprintf(fid, 'Field[3].Ratio = 1.1;\n');
fprintf(fid, 'Field[3].Thickness = 0.2;\n');
fprintf(fid, 'Field[3].Size = %f;\n',lc_TE);
fprintf(fid, 'Field[3].SizeFar = %f;\n',lc_surf);

fprintf(fid, 'Mesh.BoundaryLayerFanElements = 8;\n');
fprintf(fid, 'Field[3].FanPointsList = {200};\n');

fprintf(fid, '// --- Combine and Apply the Fields ---\n');
fprintf(fid, 'Field[6] = Min;\n');
fprintf(fid, 'Field[6].FieldsList = {1, 2};\n');
fprintf(fid, 'Background Field = 6;\n\n');
fprintf(fid, 'BoundaryLayer Field = 3;\n');

fprintf(fid, 'Mesh.MeshSizeExtendFromBoundary = 0;\n'); 
fprintf(fid, 'Mesh.MeshSizeFromPoints = 1;\n');         % Forces gmsh to follow the refinement fields
fprintf(fid, 'Mesh.MeshSizeFromCurvature = 0;\n\n');

fprintf(fid, 'Recombine Surface {1};\n');
fprintf(fid, 'Mesh.Algorithm = 8;\n');              % Frontal-Delaunay for Quads
fprintf(fid, 'Mesh.RecombinationAlgorithm = 2;\n'); % Simple Full Quad
fprintf(fid, 'Mesh.RecombineAll = 1;\n\n');         

% Creates physical meaning for Nek5000
fprintf(fid, '// --- Physical Tags ---\n');
fprintf(fid, 'Physical Curve("W") = {1, 2};\n');           
fprintf(fid, 'Physical Curve("v") = {3, 4, 5, 6};\n');     
fprintf(fid, 'Physical Curve("O") = {7};\n');           
fprintf(fid, 'Physical Surface("Fluid") = {1};\n');

fclose(fid);
fprintf('Successfully generated %s\n', filename);