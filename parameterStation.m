% Secondary parameters

addpath('/home/eduardocastro/Desktop/repo lab nazareth/sc-toolbox-3.1.3')

%% Set domain options

domainOptions = struct();
domainOptions.Nzeta = 35;
domainOptions.hallLength = 2*lambda;
domainOptions.corridorLength = 5*lambda;
domainOptions.outdoorsLength = 2*lambda;

switch geometry

    case 'crinkled'
        
        domainOptions.type = 'periodic'; %Choose: 'periodic' or 'random'
        domainOptions.symmetry = 'oneSideOnly'; %Choose: 'oneSideOnly' or 'mirrored' ('asymmetric' possibly in the future)
        %domainOptions.numberOfPoints = 30;
        domainOptions.rugosityDepth = 0.0*b;
        domainOptions.rugosityInterval = [domainOptions.hallLength,...
            domainOptions.hallLength+domainOptions.corridorLength];
        domainOptions.crinkleLength = 0.05*lambda;

        %{
        switch domainOptions.type

            case 'periodic'
        
                %noise = zeros(1,domainOptions.numberOfPoints-2);
        
                %for i = 1:domainOptions.numberOfPoints-2
                   %noise(i) = domainOptions.rugosityDepth*1i*(-1)^i;
                %end
        
            case 'random'
                %noise = 1i*randomArray(domainOptions.numberOfPoints-2,...
                    %-domainOptions.rugosityDepth, domainOptions.rugosityDepth);
        end
        %}

    case 'funil'

        domainOptions.hallWidth = 3*b;
        domainOptions.outdoorsWidth = b;


    case 'hourglass'
        
        domainOptions.sharpness = 0.05*domainOptions.hallLength;
        domainOptions.gapWidth = 1.7*b;
        domainOptions.gapLength = 0.1*domainOptions.hallLength;

    case 'asymHourglass'
        
        domainOptions.sharpness = 0.05*domainOptions.hallLength;
        domainOptions.gapWidth = 1.7*b;
        domainOptions.gapLength = 0.1*domainOptions.hallLength;

    case 'swift'
    
    domainOptions.swiftness = 0.2*b;

end





%% Set solver options
waveOptions = struct();
waveOptions.frameRate = 150;
waveOptions.finalTimeCap = 30000;
waveOptions.plotFlag = false;

%% Set wave view options
waveViewOptions = struct();

waveViewOptions.playMoviePhys = true;
waveViewOptions.playMovieCan = false;
waveViewOptions.simpleComp = false;



