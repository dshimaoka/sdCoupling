function p = getLRconnectivity(p, pLR)
%reconnect pLR*Ne neurons so they projct to a randomly assigned position

%p.Ne = 1024;
%pLR=0.1;%fraction of neurons with long-range projection against all excitatory neurons

%% 1. find p*Ne neurons as long-range projection neurons
allIdx = randperm(p.Ne);
nLRneurons = round(pLR*p.Ne);
LRIdx = allIdx(1:nLRneurons);
sourcePositions = p.pos_e(LRIdx);

%% 2. assign distance from a uniform distribution
allIdx2 = randperm(p.Ne);
distIdx = allIdx2(1:nLRneurons);
tgtPositions = p.pos_e(distIdx);

%% 3. find p.Ceed connections drawn from a gaussian distribution

for ii = 1:nLRneurons
    p.WEEd(:,LRIdx(ii)) = getConnectivity(p.pos_e, sourcePositions(ii), ...
        p.sigma_e, p.Ceed, tgtPositions(ii));
end

%%TODO
%add another input (tgtposition) to getConnectivity
%omit autopsy
%check Weed dimension tgt x src
