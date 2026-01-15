function [] = processWaveData(kappa, b, geometry, options)

[waveName, domainName] = standardNaming(geometry, kappa);

load(waveName,'H', 'h', 'h0', 'x0')
load(domainName, 'w','z','J')

%surf(real(w),imag(w),J)

X = real(z);
Y = imag(z);

%% Comparing final against initial wave profile
if options.simpleComp
    xi = real(w);
    xi = xi(1,:);
    
    h0 = squeeze(H(1,floor(end/2),:)); %init profile
    h = squeeze(H(end,floor(end/2),:)); %final profile
    
    figure;
    plot(xi,h,'r'), hold on % middle cross-section of init prof
    
    [~,loc] = findpeaks(h,'MinPeakHeight', 0.07); % finds index for which final wave prof peaks (center of final gaussian)
    [~,loc0] = findpeaks(h0,'MinPeakHeight', 0.07);
    %loc0 = find(round(xi,4)==x0, 0.1, 'last'); % gets index for inital center of gaussian pulse 
    
    shift = loc-loc0;
    
    h0_shift = circshift(h0,shift);
    
    plot(xi,h0_shift,'--k') %plot initial prof against final prof translated for matching peaks
    %plot(xi,h0,'--k')
    %legend('Initial','Final')
end


%% 3D Animation
if options.playMoviePhys
    figure
    tmp = size(H);
    for i = 1:tmp(1)
        
        h = squeeze(H(i,:,:));        
        %h = reshape(H(:,i),size(z));

        mesh(X,Y,h)

        %view(0,90);
        %zlim([-0.02,.12])
        xlabel('X'); ylabel('Y'); zlabel('h','Rotation', 0);
        title(['Time evolution of wave profile = ', num2str(i)]);
        %title(mytitle)

        set(gca, 'FontSize', 16)   % makes axis numbers larger
  
        pause(0.03)

        drawnow;

    end
end

if options.playMovieCan
    figure
    tmp = size(H);
    for i = 1:tmp(2)

        h = reshape(H(:,i),size(z));

        mesh(real(w),imag(w),h)
        %view(0,90);
        %zlim([-0.02,.12])
        xlabel('\xi'); ylabel('\zeta'); zlabel('h','Rotation', 0);
        %title(['Time evolution of wave profile = ',num2str(t)]);
        %title(mytitle)

        set(gca, 'FontSize', 16)   % makes axis numbers larger
  
        pause(0.01)

        drawnow;

    end
end

end