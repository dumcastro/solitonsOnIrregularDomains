clear all, clc, close all

addpath(['/Users/eduardocastro/Library/' ...
    'CloudStorage/GoogleDrive-eduardomdecastro@gmail.com/Meu Drive/ROOT/WORK/Files/03 Resources/sc-toolbox-3.1.3'])

% REQUIRES TOBY DRISCOLL SC TOOLBOX

%Experiment 2: crinkled walls
%{
We consider the linearized shallow-water system:
h_t + u_x + v_y = 0
u_t + h_x = 0
v_t + h_y = 0
%}

anim = 1; %1 for physical, 2 for canonical animation
travel_distance = 8;% in lambda_e units
%-----------

Ly = 1; %width of channel (idea is fix to 1 and let lambda_e be variable)
lambda_e = [4, 11, 18, 25]; 
index = 1; % 1, 2, 3 or 4, for the possibilities above

kappa = 1./lambda_e; %(width/wavelength regime parameter)

dxi = 0.006*lambda_e(index);
dzeta = 0.02; 

Lx = lambda_e(index)*(travel_distance+3);% Length of the domain (gets longer if wavelength is longer...)

%---Setting domain angle---

%---Extending a bit the domain: set a gap to avoid Jacobian singularity
ep = 0.005; 

Gamma = 1/Lx; % aspect ratio parameter
Lx = 1/Gamma; % Perimeter/lengh of full domain
b = 1; %width of channel (must be 1)
crink_depth = 0.1*b;

number_of_pts = 15;

ver = [0, Lx, 1i*b + Lx, 1i*b + Lx - Lx/4 ];

for i = 1:number_of_pts+1
    last = ver(end);
    ver = [ver, last - (Lx/2)/(number_of_pts+1)];
end

ver = [ver, 1i*b];

noise = 1i*random_array(number_of_pts, -crink_depth, crink_depth);

ver(5:end-2) = ver(5:end-2) + noise;

sang = [0.5, 0.5, 0.5, 1, ones(1,number_of_pts), 1, 0.5];

%sang = [0.5000, 1, 0.5000, 0.5000, 1, 0.5000]; %angles of typical rectangle
%--------
%Pep = polygon(verep);
P = polygon(ver);  % Gets poly

f_tilde = crrectmap(P, sang); % Gets the map
%Cep_tilde = evalinv(f_tilde,Pep);
C_tilde = evalinv(f_tilde,P); % Compute the inverse mapping to obtain the canonical domain

%---Getting the dilation paratemeter---
%Checking the boundaries of rectangle C
xi_lowlim=min(real(vertex(C_tilde)));
xi_highlim=max(real(vertex(C_tilde)));

Lxi_tilde = (xi_highlim)-(xi_lowlim); % Length of the canonical domain
%Lzeta= (zeta_highlim)-(zeta_lowlim);

alpha = Lxi_tilde/Lx; % How much the domain gets strecthed/shrinked

C = (1/alpha)*C_tilde;
%Cep = (1/alpha)*Cep_tilde;

%if ~(vertex(P) == vertex(C))
%  error('P e C não estão batendo') 
%end

%
figure;
subplot(1, 3, 1);
plot(P, 'b', 'LineWidth', 2,'k');
title('Physical Region');

subplot(1, 3, 2);
plot(C_tilde, 'r', 'LineWidth', 2);
title('Numerical canonical domain');

subplot(1, 3, 3);
plot(C, 'y', 'LineWidth', 2);
title('Canonical domain width = 1');
%}

%% Getting the canonical grid
%Checking the boundaries of rectangle C
xi_lowlim=min(real(vertex(C)));
xi_highlim=max(real(vertex(C)));

zeta_lowlim=min(imag(vertex(C)))+ep;
zeta_highlim=max(imag(vertex(C)));

Lxi = (xi_highlim) - (xi_lowlim); % Length of the canonical domain
Lzeta = (zeta_highlim) - (zeta_lowlim);

Nxi = length(dxi); % Number of grid points
Nzeta = length(dzeta);

xi = xi_lowlim:dxi:xi_highlim;
zeta = zeta_lowlim:dzeta:zeta_highlim;

[Xi, Zeta] = meshgrid(xi,zeta);
w = Xi + 1i * Zeta; % Grid in canonical variables
%w = w - ep*1i;

z = eval(f_tilde, alpha*w); % Forward map gets "physical grid"
%z = z + ep*1i;

%{
[rows, cols] = size(z);

% Generate random row and column indices
row_indices = randi(rows, 1000, 1);
col_indices = randi(cols, 1000, 1);

% Extract the 10 random points
random_w = w(sub2ind(size(w), row_indices, col_indices));
random_z = z(sub2ind(size(w), row_indices, col_indices));

% Scatter plot of the selected points
figure;
scatter(real(random_w), imag(random_w), 100, 'b', 'Filled'), hold on
scatter(real(random_z), imag(random_z), 50, 'y', 'Filled')
legend('w points', 'z points')
title('Scatter Plot of Randomly Selected Points');
grid on;
%}

%% Getting the Jacobian
dz = evaldiff(f_tilde, alpha*w); % Evaluate the derivative of the SC map across the entire grid
J = (alpha^2)*abs(dz).^2;

%J_const = mode(mode(round(J,4))); %This gets the 'scaling factor' of f (equal to 1/alpha^2)
%J = J*(1/J_const)*b; % This corrects the scaling factor to be equal to the width b

%---Plot the Jacobian determinant as a 3D surface plot---

jump_xi = 12; jump_zeta = 1;

figure;
surf(xi(1:jump_xi:end-20), zeta(1:jump_zeta:end), J(1:jump_zeta:end,1:jump_xi:end-20));
xlabel('Real Axis');
ylabel('Imaginary Axis');
%title(['alpha= ', num2str(alpha)]);

%---Checking aspect ratio---
verC = vertex(C);
heightC = abs(verC(6) - verC(1));
lengthC = abs(verC(3) - verC(1));
AR_can = heightC/lengthC; %aspect ratio of canonical domain
AR_physical = b/Lx; %this is just Gamma

% AR_physical and AR_can should always be close

%plot(P), hold on, plot(C) %Uncomment to visualize scaling
%{
%---Save data---
if want_save == 1
    th_degree=rad2deg(theta);
    save(['L_angle= ', num2str(th_degree), ' width= ', num2str(b), ' AR= ', num2str(Gamma*100)])
end

%----------------------------------------------------------------
%----------------------------------------------------------------
%}

%dxi = 0.005*lambda_e(index)*alpha;
%dy = .5;
%x = 0:dxi:Lx;
%y = 0:dy:Ly;

% Parameters
dt = dxi/3; %Time step size

CC = sqrt(-2*log(1/100)); % parameter for getting sigma as a function of lambda_e

% Define pulse parameters
%sigma = lambda_e(index)/(2*CC);    % Width of the pulse

x0 = 3*lambda_e(index)/2;     % Center of the pulse
a = 0.1;
sigma = lambda_e(index)/(2*asech(sqrt(0.01*a)));
%h0 = a*exp(-(Xi-x0).^2/ (2 * sigma^2));
h0 = a*sech((Xi-x0)/sigma).^2;
h = h0;
u = 0;        % Initial velocity
v = h;

X = real(z);
Y = imag(z);

% Main loop
x0f = x0; 
dist = 0; %distance between current center of pulse and initial center of pulse
iter = 0;

while dist < travel_distance*lambda_e(index)
% RK4 time-stepping
    k1_u = -dt * (circshift(h, [ -1 0]) - circshift(h, [ 1 0])) / (2 * dxi);
    k1_v = -dt * (circshift(h, [ 0 -1]) - circshift(h, [0 1])) / (2 * dzeta);
    k1_h = -dt * ((circshift(u, [ -1 0]) - circshift(u, [ 1 0])) / (2 * dxi) + (circshift(v, [ 0 -1]) - circshift(v, [ 0 1])) / (2 * dzeta))./J;
    
    k2_u = -dt * (circshift(h + 0.5 * k1_h, [ -1 0]) - circshift(h + 0.5 * k1_h, [ 1 0])) / (2 * dxi);
    k2_v = -dt * (circshift(h + 0.5 * k1_h, [ 0 -1]) - circshift(h + 0.5 * k1_h, [ 0 1])) / (2 * dzeta);
    k2_h = -dt * ((circshift(u + 0.5 * k1_u, [ -1 0]) - circshift(u + 0.5 * k1_u, [ 1 0])) / (2 * dxi) + (circshift(v + 0.5 * k1_v, [ 0 -1]) - circshift(v + 0.5 * k1_v, [ 0 1])) / (2 * dzeta))./J;
    
    k3_u = -dt * (circshift(h + 0.5 * k2_h, [ -1 0]) - circshift(h + 0.5 * k2_h, [ 1 0])) / (2 * dxi);
    k3_v = -dt * (circshift(h + 0.5 * k2_h, [ 0 -1]) - circshift(h + 0.5 * k2_h, [ 0 1])) / (2 * dzeta);
    k3_h = -dt * ((circshift(u + 0.5 * k2_u, [ -1 0]) - circshift(u + 0.5 * k2_u, [ 1 0])) / (2 * dxi) + (circshift(v + 0.5 * k2_v, [ 0 -1]) - circshift(v + 0.5 * k2_v, [ 0 1])) / (2 * dzeta))./J;
    
    k4_u = -dt * (circshift(h + k3_h, [ -1 0]) - circshift(h + k3_h, [ 1 0])) / (2 * dxi);
    k4_v = -dt * (circshift(h + k3_h, [0 -1]) - circshift(h + k3_h, [ 0 1])) / (2 * dzeta);
    k4_h = -dt * ((circshift(u + k3_u, [ -1 0]) - circshift(u + k3_u, [ 1 0])) / (2 * dxi) + (circshift(v + k3_v, [ 0 -1]) - circshift(v + k3_v, [ 0 1])) / (2 * dzeta))./J;
    
    u = u + (1/6) * (k1_u + 2 * k2_u + 2 * k3_u + k4_u);
    v = v + (1/6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
    h = h + (1/6) * (k1_h + 2 * k2_h + 2 * k3_h + k4_h);
    
    % Apply reflecting boundary conditions
    h(1,:) = h(2,:); % left boundary
    h(end,:) = h(end-1,:); % right boundary
    h(:,1) = h(:,2); % bottom boundary
    h(:,end) = h(:,end-1); % top boundary
    
    v(:,1) = 0;
    v(:,end) = 0;
    
    if anim == 1 && mod(iter, 1/dxi) == 0 %PHYSICAL
        jump_xi = 5; jump_zeta = 1;

        XX = X(1:jump_zeta:end,1:jump_xi:end);
        YY = Y(1:jump_zeta:end,1:jump_xi:end);
        hh = h(1:jump_zeta:end,1:jump_xi:end);
        
        mesh(XX, YY, hh,'edgecolor', 'k'); 
        zlim([a*(-0.2),a*1.5])
        xlabel('X');
        ylabel('Y');
        zlabel('h');
        title(['iterations = ',num2str(iter), ' kappa= ', num2str(kappa(index))]);
        drawnow, pause(0.1)
    elseif anim == 2 && mod(iter, 100) == 0 %CANONICAL
        mesh(Xi, Zeta, h,'edgecolor', 'k'); 
        zlim([a*(-0.2),a*1.5])
        xlabel('Xi');
        ylabel('Zeta');
        zlabel('h');
        title(['Number of iterations = ',num2str(iter)]);
        drawnow
    elseif anim ==3 && mod(iter, 10/dxi) == 0 %BOTH
       %...
    end
     
hh = h(floor(end/2),:);   
[~,loc] = findpeaks(hh,'MinPeakHeight', 0.7*a); % finds index for which final wave prof peaks (center of final gaussian)
x0f = xi(loc);
dist = abs(x0 - x0f); %gets distance between current and initial peak positions
iter = iter + 1;
end


%% Comparing final against initial wave profile
figure;
plot(xi,h0(floor(end/2),:),'r'), hold on % middle cross-section of init prof
%hh = h(floor(end/2),:);

%[~,loc] = findpeaks(hh,'MinPeakHeight', 0.9*a); % finds index for which final wave prof peaks (center of final gaussian)
loc0 = find(round(xi,4)==x0); % gets index for inital center of gaussian pulse 
plot(xi,circshift(hh,-(loc-loc0)),'--k') %plot initial prof against final prof translated for matching peaks
legend('Initial','Final')


function rand_array = random_array(n, min_val, max_val)
    rand_array = min_val + (max_val - min_val) * rand(1, n);
end




