function [waveName, domainName] = standardNaming(geometry, kappa)

domainName = ['DomainData/geometry= ', geometry, '.mat'];

waveName = ['WaveData/kappa=',mat2str(kappa),'geometry= ', geometry,  '.mat'];

end