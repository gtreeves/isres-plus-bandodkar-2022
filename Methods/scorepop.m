function score = scorepop(I,vstat,lambda)

sc = fliplr(1:lambda);
sc = sc.^10;
% sc = sc;
%plot(sc)

[~,sreco]  = intersect(I, find((vstat == 1)));
[~,smuta]  = intersect(I, find((vstat == 2)));
[~,slin]   = intersect(I, find((vstat == 3)));
[~,snewt]  = intersect(I, find((vstat == 4)));

if isempty(sreco), score(1) = 0; else, score(1) = max(sc(sreco)); end
if isempty(smuta), score(2) = 0; else, score(2) = max(sc(smuta)); end
if isempty(slin),  score(3) = 0; else, score(3) = max(sc(slin));  end
if isempty(snewt), score(4) = 0; else, score(4) = max(sc(snewt)); end

score   = score/sum(score);

end