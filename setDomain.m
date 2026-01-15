function [w,J] = setDomain(b, geometry,options)
% Creates the SC computational domain of a complex channel

%% Step 1: set physical pob
%ver = [0, Lx, 1i*b + Lx, 1i*b + rugosityInterval(2)];

%ver = [0, Lx/3:crinkleLength:2*Lx/3, Lx, b*1i + Lx, b*1i + tmp, b*1i];

%{
for i = 1:options.numberOfPoints+1

    last = ver(end);
    ver = [ver, last - diff(rugosityInterval)/(options.numberOfPoints+1)];
end
%}


switch geometry
    case 'crinkled'
        Lx = options.hallLength+options.corridorLength+options.outdoorsLength;
        tmp = options.hallLength+options.corridorLength:-options.crinkleLength:options.hallLength;
        switch options.type
            case 'periodic'

                noise = zeros(1,length(tmp)-2);

            for i = 1:length(tmp)-2
               noise(i) = options.rugosityDepth*1i*(-1)^i;
            end

            case 'random'
            noise = 1i*randomArray(length(tmp)-2, -options.rugosityDepth,options.rugosityDepth);
    
        end

        switch options.symmetry
            case 'oneSideOnly'
                ver = [0, Lx, 1i*b + Lx, 1i*b + tmp, 1i*b];
                ver(5:end-2) = ver(5:end-2) + noise;

            sang = [0.5, 0.5, 0.5, 1, ones(1,length(tmp)-2), 1, 0.5];
            case 'mirrored'
                ver = [0, Lx/3:options.crinkleLength:2*Lx/3, Lx, 1i*b + Lx, 1i*b + tmp, 1i*b];
        
                ver(3:length(tmp)) = ver(3:length(tmp)) - noise;
        
                ver((end-length(tmp)+1):(end-2)) = ver((end-length(tmp)+1):(end-2)) + noise;
        
        
                sang = [0.5, ones(1,length(tmp)), 0.5, 0.5, ones(1,length(tmp)), 0.5];

        end

    case 'funil'

        gamma = 1i*(b-options.outdoorsWidth)/2;

        ver = [0, options.hallLength,...
            options.hallLength + options.corridorLength + gamma];
        ver = [ver, ver(end) + options.outdoorsLength];
        aux = flip(conj(ver));
        ver = [ver, aux + b*1i];

        %ver = [ver, ver(end) +  options.sharpness - gamma];
        %ver = [ver, ver(end) +  options.outdoorsLength];
        %aux = flip(conj(ver));
        %ver = [ver, aux + b*1i];

        sang = [0.5,1,1,0.5,0.5,1,1,0.5];

    case 'asymFunil'

        gamma = 1i*(b-options.gapWidth)/2;
        
        ver = [0, options.hallLength,...
            options.hallLength + options.sharpness + gamma];
        ver = [ver, ver(end) + options.corridorLength];
        ver = [ver, ver(end) +  options.sharpness - gamma];
        ver = [ver, ver(end) +  options.outdoorsLength];
        aux = flip(conj(ver));
        ver = [ver, aux + b*1i];
        
        sang = [0.5,1,1,1,1,0.5,0.5,1,1,1,1,0.5];





    case 'hourglass'
        gamma = 1i*(b-options.gapWidth)/2;

        ver = [0, options.hallLength,...
            options.hallLength + options.sharpness + gamma];
        ver = [ver, ver(end) + options.corridorLength];
        ver = [ver, ver(end) +  options.sharpness - gamma];
        ver = [ver, ver(end) +  options.outdoorsLength];
        aux = flip(conj(ver));
        ver = [ver, aux + b*1i];

        sang = [0.5,1,1,1,1,0.5,0.5,1,1,1,1,0.5];

    case 'asymHourglass'
        gamma = 1i*(b-options.gapWidth);

        ver = [0, options.hallLength,...
            options.hallLength + options.sharpness + gamma];
        ver = [ver, ver(end) + options.corridorLength];
        ver = [ver, ver(end) +  options.sharpness - gamma];
        ver = [ver, ver(end) +  options.outdoorsLength];
        ver = [ver, ver(end) + 1i*b, 1i*b];

        sang = [0.5,1,1,1,1,0.5,0.5,0.5];

    case 'swift'

        ver = [0, options.hallLength, options.hallLength + options.swiftness*1i];
        ver = [ver, ver(end)+options.outdoorsLength];
        ver = [ver, ver(end) + 1i*b];
        ver = [ver, ver(end) - options.outdoorsLength];
        ver = [ver, ver(end) - 1i*options.swiftness];
        ver = [ver, 1i*b];

        sang = [0.5,1,1,0.5,0.5,1,1,0.5];

end

P = polygon(ver);  % Gets poly
%P = polyedit(P) % uncomment to visualize poly before proceding

%% Step 2: calculate SC map to rectangle

disp('computing SC rect map...')
f_tilde = crrectmap(P, sang); % Gets the map
%Cep_tilde = evalinv(f_tilde,Pep);
C_tilde = evalinv(f_tilde,P); % Compute the inverse mapping to obtain the canonical domain

%% Step 3: calculate scaling factor alpha
zeta_lims = [min(imag(vertex(C_tilde))), max(imag(vertex(C_tilde)))];
Lzeta_tilde = diff(zeta_lims);
alpha = Lzeta_tilde/b; % This is the scaling factor

%% Step 4: get C by dilation
C = (1/alpha)*C_tilde;

%% Step 5: create computational grid
xi_lims = [min(real(vertex(C))), max(real(vertex(C)))];
zeta_lims = [min(imag(vertex(C))), max(imag(vertex(C)))];

dzeta = (zeta_lims(2)-zeta_lims(1))/(options.Nzeta-1); %dxi, dzeta are obtained from Nzeta choice
dxi = dzeta;

xi = xi_lims(1):dxi:xi_lims(2);
zeta = zeta_lims(1):dzeta:zeta_lims(2);

[Xi, Zeta] = meshgrid(xi, zeta);
w = Xi + 1i * Zeta;

w = w(2:end-1, :);

%% Final step: Get physical domain and Jacobian
disp('getting physical mesh...')
z = eval(f_tilde,alpha*w);
disp('calculating Jacobian of SC map...')
dz = evaldiff(f_tilde, alpha*w); % Evaluate the derivative of the SC map across the entire grid
J = (alpha^2)*abs(dz).^2;


%% Setting the mean Jacobian Jm
tmp = size(J);

Jm = zeros(1,tmp(2));

for i = 1:tmp(2)

Jslice = cumtrapz(J(:,i));   
Jm(i) = Jslice(end);

end

Jm = Jm/b;

[~, domainName] = standardNaming(geometry, 1);

save(domainName, 'w', 'J', 'z','Jm','P')
end