clear all, clc, close all
%% Defining space discretization

lambda_e = 13;
%travel_distance = 10;% in lambda_e units

%Lx = lambda_e*(travel_distance+3);% Length of the domain
Lx = 50;
Ly = 5;

dx = 0.2;
dy = 0.2;

x = 0:dx:Lx;
y = 0:dy:Ly;

Nx = length(x);
Ny = length(y);

%% Load Jacobian data (WIP)
%J = 1; %constant J for now
J = ones(Ny,Nx);

%% Defining time discretization
T = 10;
dt = 0.5*dx;
t = 0:dt:T;
Nt = length(t);
speed = 2;

%% Check CFL (WIP)
CFL = dt / sqrt(dx^2 + dy^2);
if CFL > 1
    error('CFL condition not satisfied. Stopping execution.');
end

%% Initial data
% Parameters
alpha = 0.3;
beta = sqrt(3 * alpha / (4 * (1 + 0.68 * alpha)));
c = sqrt(6 * (1 + alpha)^2 / (alpha^2 * (3 + 2 * alpha)) * ((1 + alpha) * log(1 + alpha) - alpha));

% Define eta_0 as a function of x, with default x_0 = 0 and t = 0
eta_0 = @(x, x_0, t) alpha ./ cosh(beta * (x - x_0 - c * t)).^2 ./ (1 + alpha * tanh(beta * (x - x_0 - c * t)).^2);

% Parameter d for potential phi
d = c * alpha / (beta * (1 + alpha));

% Define phi_0 as a function of x, with default x_0 = 0 and t = 0
phi_0 = @(x, x_0, t) c * alpha / (beta * (1 + alpha)) * tanh(beta * (x - x_0 - c * t));

x_0 = 3*lambda_e/2;     % Center of the pulse

eta = zeros(Nt,Ny,Nx);
phi = zeros(Nt,Ny,Nx);

eta(1,:,:) = repmat(eta_0(x, x_0, 0),Ny,1);
phi(1,:,:) = repmat(phi_0(x, x_0, 0),Ny,1);

%% Fin dif funcs
m = @(f,n) f(n,2:end-1, 2:end-1);

Dx = @(f,n) (f(n,2:end-1, 3:end) - f(n,2:end-1, 1:end-2))/(2*dx);
 
Dy = @(f,n) (f(n,3:end, 2:end-1) - f(n,1:end-2,2:end-1))/(2*dy);

DDx = @(f,n) (f(n,2:end-1, 3:end) - 2*f(n,2:end-1, 2:end-1) + f(n,2:end-1, 1:end-2))/(dx^2);

DDy = @(f,n) (f(n,3:end ,2:end-1) - 2*f(n,2:end-1,2:end-1) + f(n,1:end-2,2:end-1))/(dy^2);

%{
eta_x = Dx(eta,dx); eta_y = Dy(eta, dy);
phi_x = Dx(phi,dx); phi_y = Dy(phi, dy);

phi_x2 = DDx(phi,dx); phi_y2 = DDy(phi, dy);
%
%% Aux functions
%E = @(n) -((1+m(eta,n)).*(DDx(phi,n)+DDy(phi,n)) + Dx(eta,n).*Dx(phi,n) + Dy(eta,n).*Dy(phi,n))./J;
%F = @(n) m(phi,n) -(DDx(phi,n)+DDy(phi,n))./(3*J);
%G = @(n) - m(eta,n) - (Dx(phi,n).^2 + Dy(phi,n).^2)./(2*J);

%E = -((1+m(eta)).*(phi_x2 + phi_y2) + eta_x.*phi_x + eta_y.*phi_y)./J;
%F = m(phi) - (phi_x2 + phi_y2)./(3*J);
%G = - m(eta) - (phi_x.^2 + phi_y.^2)./(2*J);
%}
%% Constructing the matrix A
% Mapping from 2D index (i, j) to 1D index k
index = @(i, j) (i - 1) * Nx + j;

% Discretization coefficients
r = (1/3)/(dx^2);
s = (1/3)/(dy^2);

% Initialize sparse matrix A and right-hand side vector f
N = Nx * Ny;

%
A = sparse(N, N);  % Sparse matrix initialization in MATLAB
% Construct the matrix A
for i = 2:Ny-1
    for j = 2:Nx-1
        k = index(i, j);  % 1D index for grid point (i, j)

        % Elliptical equation in 2D: u - (u_xx + u_yy) = 3F
        A(k, k) = 1 + 2*(r + s)/J(i,j);  % Center point: u_{i,j}

        % x-direction neighbors: u_{i+1,j} and u_{i-1,j}
        A(k, index(i+1, j)) = -r/J(i,j);   % u_{i+1, j}
        A(k, index(i-1, j)) = -r/J(i,j);   % u_{i-1, j}

        % y-direction neighbors: u_{i,j+1} and u_{i,j-1}
        A(k, index(i, j+1)) = -s/J(i,j);   % u_{i, j+1}
        A(k, index(i, j-1)) = -s/J(i,j);   % u_{i, j-1}
    end
end

% Left boundary (x = 0)
for j = 1:Ny
    k = index(j, 1);
    A(k, :) = 0;
    A(k, index(j, 1)) = 1;
    A(k, index(j, 2)) = -1;  % u_{0,j} = u_{1,j}
end

% Right boundary (x = Lx)
for j = 1:Ny
    k = index(j, Nx);
    A(k, :) = 0;
    A(k, index(j, Nx)) = 1;
    A(k, index(j, Nx-1)) = -1;  % u_{Nx-1,j} = u_{Nx-2,j}
end

% Bottom boundary (y = 0)
for i = 1:Nx
    k = index(1, i);
    A(k, :) = 0;
    A(k, index(1, i)) = 1;
    A(k, index(2, i)) = -1;  % u_{i,0} = u_{i,1}
end

% Top boundary (y = Ly)
for i = 1:Nx
    k = index(Ny, i);
    A(k, :) = 0;
    A(k, index(Ny, i)) = 1;
    A(k, index(Ny-1, i)) = -1;  % u_{i,Ny-1} = u_{i,Ny-2}
end

clear i j
%save('myA.mat','A')
%}

%% Loading A_python for comparison
%load('myA.mat')
%load('A_python')

%A_dense = csvread('A_python.csv');  % For older MATLAB versions
%A_dense = readmatrix('A_python.csv'); % For newer MATLAB versions (R2019b or later)

%A_python = sparse(A_dense);

[L, U] = lu(A);

EE = zeros(Nt, Ny, Nx);
FF = zeros(Nt, Ny, Nx);
GG = zeros(Nt, Ny, Nx);

%% Main loop
innerIterCap = 3;
J = reshape(J, [1,size(J)]);
jump_x = 1; %grid spacing for better visualization
dist = 0;

outerIter = 1;
n = outerIter;
outerIterCap = 60;
%while dist < travel_distance*lambda_e
while outerIter < min(outerIterCap, Nt)
    
	%% Predictor initial guess
    EE(n,2:end-1,2:end-1) = ...
        -((1+m(eta,n)).*(DDx(phi,n)+DDy(phi,n)) +...
        Dx(eta,n).*Dx(phi,n) + ...
        Dy(eta,n).*Dy(phi,n))./J(1,2:end-1,2:end-1);
    FF(n,2:end-1,2:end-1) = ...
        m(phi,n) -(DDx(phi,n)+...
        DDy(phi,n))./(3*J(1,2:end-1,2:end-1));
    GG(n,2:end-1,2:end-1) = ...
        -m(eta,n) - (Dx(phi,n).^2 + Dy(phi,n).^2)./(2*J(1,2:end-1,2:end-1));
    eta(n+1,2:end-1,2:end-1) = m(eta,n) + dt*m(EE,n); 
    FF(n+1,2:end-1,2:end-1) = m(FF,n) + dt*m(GG,n);
    FF(n+1,:,:) = neumann_correction(squeeze(FF(n+1,:,:)));
    eta(n+1,:,:) = neumann_correction(squeeze(eta(n+1,:,:)));
    aux = solve(squeeze(FF(n+1,:,:)),L,U);
    %aux2 = solve(squeeze(FF(n+1,:,:)),A_python);
    phi(n+1,:,:) = aux; 
    %{
    %% Debug (prova real)
    phir = squeeze(phi(1,:,:));
    phir = set_f(phir);
    %phir = reshape(phir,N,1);
    sbF = A*phir;
    sbF = reshape_new(sbF,Ny,Nx);
    
    figure(1)
    plot(x(2:end), squeeze(FF(1,5,2:end))), hold on
    plot(x(2:end), squeeze(sbF(5,2:end)))
    %plot(x(2:end), squeeze(phi(1,5,2:end)))
    
    legend('This if F', 'This should be F', 'phi')
    %{
    p = 34;
    
    sbphi = solve(squeeze(FF(1,:,:)),L,U); %should be phi(1)
    figure(2)
    plot(x, squeeze(phi(1,5,:))), hold on
    plot(x, squeeze(sbphi(5,:)))
    
    pp = 23;
    %}
    
    %plot(x(2:end-1), squeeze(FF(1,5,2:end-1))), hold on
    %plot(x(2:end-1), squeeze(FF(2,5,2:end-1)))
    
    %figure(2)
    %plot(x(2:end-1), squeeze(phi(1,5,2:end-1))), hold on
    %plot(x(2:end-1), squeeze(phi(2,5,2:end-1)))
    %}
      
    %% Corrector
    innerIter = 0; 
    while innerIter < innerIterCap
        %EE(n+1,2:end-1,2:end-1) = E(n+1); %EE(n+1,:,:) = neumann_correction(squeeze(EE(n+1,:,:)));
        %FF(n+1,2:end-1,2:end-1) = F(n+1);
        %GG(n+1,2:end-1,2:end-1) = G(n+1);
        
        EE(n+1,2:end-1,2:end-1) = ...
            -((1+m(eta,n+1)).*(DDx(phi,n+1)+DDy(phi,n+1)) + Dx(eta,n+1).*Dx(phi,n+1) + Dy(eta,n+1).*Dy(phi,n+1))./J(1,2:end-1,2:end-1);
        FF(n+1,2:end-1,2:end-1) = ...
            m(phi,n+1) -(DDx(phi,n+1)+DDy(phi,n+1))./(3*J(1,2:end-1,2:end-1));
        GG(n+1,2:end-1,2:end-1) = ...
            - m(eta,n+1) - (Dx(phi,n+1).^2 + Dy(phi,n+1).^2)./(2*J(1,2:end-1,2:end-1));
           
        eta(n+1,2:end-1,2:end-1) = m(eta,n) + (dt/2)*(m(EE,n+1)+m(EE,n));
        eta(n+1,:,:) = neumann_correction(squeeze(eta(n+1,:,:)));
        
        FF(n+1,2:end-1,2:end-1) = m(FF,n) + (dt/2)*(m(GG,n+1)+m(GG,n));
        phi(n+1,:,:) = solve(squeeze(FF(n+1,:,:)),L,U);
        
        innerIter = innerIter + 1;
    end
    
    %
    if mod(n,speed)==0
        
        %{        
        figure(1)
        plot(x(1:jump_x:end), squeeze(eta(n,5,1:jump_x:end)),'k');
        title('eta')
        drawnow, pause(0.001)
        %}
        
        figure(1)
        %plot(x(1:jump_x:end), squeeze(eta(n,5,1:jump_x:end)),'k');
        eeta = squeeze(eta(n,:,:));
        surf(x,y,eeta)
        title('eta')
        drawnow, pause(0.001)
        %}
        
        %{        
        figure(2)
        plot(x(1:jump_x:end), squeeze(phi(n,5,1:jump_x:end)), 'k');
        title('phi')
        
        figure(3)
        plot(x(1:jump_x:end), squeeze(EE(n,5,1:jump_x:end)), 'k');
        title('E')
        
        figure(4)
        plot(x(1:jump_x:end), squeeze(FF(n,5,1:jump_x:end)), 'k');
        title('F')
        
        figure(5)
        plot(x(1:jump_x:end), squeeze(GG(n,5,1:jump_x:end)), 'k');
        title('G')
        
        
        p = 32;
        
        figure(1)
        plot(x(1:jump_x:end), squeeze(eta(n+1,5,1:jump_x:end)),'k');
        title('eta')
        
        figure(2)
        plot(x(1:jump_x:end), squeeze(phi(n+1,5,1:jump_x:end)), 'k');
        title('phi')
        
        figure(3)
        plot(x(1:jump_x:end), squeeze(EE(n+1,5,1:jump_x:end)), 'k');
        title('E')
        
        figure(4)
        plot(x(1:jump_x:end), squeeze(FF(n+1,5,1:jump_x:end)), 'k');
        title('F')
        
        figure(5)
        plot(x(1:jump_x:end), squeeze(GG(n+1,5,1:jump_x:end)), 'k');
        title('G')
         %}
        %p = 33;
      
    end
    %}
    
    %hh = squeeze(eta(n,5,:));   
    %[~,loc] = findpeaks(hh,'MinPeakHeight', 0.7*alpha); % finds index for which final wave prof peaks (center of final gaussian)
    %x0f = x(loc);
    %dist = abs(x_0 - x0f)
    
    %n = n + 1;
    outerIter = outerIter + 1
    Nt
    n = outerIter;
end


%% Finite dif and extra functions
%{
function result = Dx(f,dx)
    result = (f(3:end, 2:end-1) - f(1:end-2,2:end-1)) / (2 * dx);
end
 
function result = Dy(f,dy)
    result = (f(2:end-1, 3:end) - f(2:end-1, 1:end-2)) / (2 * dy);
end

function result = DDx(f,dx)
    result = (f(3:end ,2:end-1) - 2*f(2:end-1,2:end-1) + f(1:end-2,2:end-1))/(dx^2);
end

function result = DDy(f,dy)
    result = (f(2:end-1, 3:end) - 2*f(2:end-1, 2:end-1) + f(2:end-1, 1:end-2))/(dy^2);
end

function result = m(f)
    result = f(2:end-1,2:end-1);
end
%}

function f = set_f(FF)
    % Initialize right-hand side vector f
    aux = size(FF);
    Ny = aux(1); Nx = aux(2);
    N = prod(aux);
    f = zeros(N, 1);  % Column vector in MATLAB
    
    % Mapping from 2D index (i, j) to 1D index k
    index = @(i, j) (i - 1)*Nx + j;
    
    % Construct the source term in the interior of the domain
    for i = 2:Ny-1
        for j = 2:Nx-1
            k = index(i, j);  % 1D index for grid point (i, j)

            % Source term f(x, y) at the point (i, j)
            f(k) = FF(i, j);
        end
    end

    % Apply Dirichlet boundary conditions (u = 0) on the left and right boundaries
    for j = 1:Ny
        % Left boundary (x = 0)
        f(index(j, 1)) = 0;

        % Right boundary (x = Lx)
        f(index(j, Nx)) = 0;
    end

    % Apply Neumann boundary conditions (no flux) on the bottom and top boundaries
    for i = 1:Nx
        % Bottom boundary (y = 0)
        f(index(1, i)) = 0;  % No flux means no contribution to f

        % Top boundary (y = Ly)
        f(index(Ny, i)) = 0;  % No flux means no contribution to f
    end
end

%--------

function h = neumann_correction(h)
    h(:, 1) = h(:, 2); % Left boundary
    h(:, end) = h(:, end-1); % Right boundary
    h(1, :) = h(2, :); % Bottom boundary
    h(end, :) = h(end-1, :); % Top boundary
end

function result = solve(FF, L, U)
    aux = size(FF);
    Ny = aux(1); Nx = aux(2);
    N = prod(aux);

    f = set_f(FF);
    
    %f = reshape(FF, N, 1);
    
    u = U\(L\f);
    
    result = reshape_new(u, Ny, Nx);
end





