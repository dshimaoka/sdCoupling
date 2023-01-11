
function W = getConnectivity(allTgtPositions, sourcePosition, sigma, nConnections, tgtPosition)
%  W = getConnectivity(alltgtPositions, sourcePosition, sigma, nConnections)
%  creates connectivity around the source position with sigma
%
%  W = getConnectivity(alltgtPositions, sourcePosition, sigma, nConnections, tgtPosition)
%  creates connectivity around the tgt position with sigma
%
% OUTPUT: (target position) x (source position)

if nargin<5
    f = @(x)exp(-(x - sourcePosition).^2/2/sigma^2);% .* (1-dirac(x-sourcePosition));
else
    f = @(x)exp(-(x - tgtPosition).^2/2/sigma^2);
end

Wori = slicesample(1,2*nConnections,'pdf',f,'thin',5,'burnin',1000);
Wori(Wori<min(allTgtPositions)) = Wori(Wori<min(allTgtPositions)) + max(allTgtPositions); 
Wori(Wori>max(allTgtPositions)) = Wori(Wori>max(allTgtPositions)) - max(allTgtPositions); 
[~,idx_all] = arrayfun(@(a)(min(abs(a-allTgtPositions))), Wori); %make it discrete
idx_noAuto = find(allTgtPositions(idx_all)~=sourcePosition);
idx = idx_all(idx_noAuto);%omit autapses
idx = idx(1:nConnections);

dx = mean(diff(allTgtPositions));
edges = [allTgtPositions allTgtPositions(end)+dx]-0.5*dx;
W = histcounts(allTgtPositions(idx), edges)';
end
%plot(tgtPositions,W);
%vline(sourcePosition);
