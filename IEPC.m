clear all;

load Flame2d2c;
Z = flame(:,1:2);
TL = flame(:,3)';
col = size(Z,2);
cnum = 2;
rows = size(Z,1);
dps = [];
dp_len = rows;
nmi = [];
ari = [];
for iter = 5:50
    dps = [];
for k = iter
    [ind,dy,rho,dypro] = IEPks(Z,k);
    dps = [dps;ind(1:dp_len)'];
    rho = rho(1:dp_len);

maxcossim = zeros(dp_len,1);
cosdis = zeros(dp_len,dp_len);

dynew = zeros(dp_len,k-1);
for i = 1 :dp_len
    dynew(i,:) = dy(dps(i),:);
end

A = dynew;
B = dynew';
normA = sqrt(sum(A .^ 2, 2));
normB = sqrt(sum(B .^ 2, 1));
cosdis = bsxfun(@rdivide, bsxfun(@rdivide, A * B, normA), normB);

for i = 1 : dp_len
    if i ~= 1
        maxcossim(i) = max(cosdis(i,1:i-1));
    end     
end


minnordis = zeros(dp_len,1);

Ztemp = zeros(dp_len,col);
for i = 1 : dp_len
    Ztemp(i,:) = Z(dps(i),:);
end
nordis = L2_dist(Ztemp,Ztemp);
for i = 1 : dp_len   
    if i ~= 1
        minnordis(i) = min(nordis(i,1:i-1));
    end     
end
minnordis(1) = max(nordis(1,:));
cosdis = 1 - cosdis;
mincosdis = 1 - maxcossim;

rho1 = normalize(rho,'zscore');
mincosdis1 = normalize(mincosdis,'zscore');
minnordis1 = normalize(minnordis,'zscore');

rhocosnor = rho1.*mincosdis1.*minnordis1;
[gamma, rd3] = sort(rhocosnor,'descend');
ngam = normalize(gamma);
scatter3(rho1,mincosdis1,minnordis1);
clustering = NNClu(Ztemp,cnum,(1:rows)',rd3(1:cnum));


label = zeros(1,rows);
for tt = 1:cnum
    tep = clustering{tt,1};
    label(ind(tep)) = tt;
end
disp('NMI:');
nmi = [nmi;NMI(TL,label)]
disp('ARI:')
ari = [ari;RandIndex(TL,label)]
end
max(nmi)
max(ari)