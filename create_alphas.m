indus=[];
for i=1:566
    indus=[indus;allstocks(i).industrylist(1).industry];    
end
indusList = unique(cellstr(indus));
rho = zeros(566,40);

for i=1:566
    [~,idx]=ismember(cellstr(allstocks(i).industrylist(1).industry),indusList);
    rho(i,idx) =1;
end

ctry=[];
for i=1:566
    ind = find(cellfun('isempty',{allstocks(i).indexlist.index})==0);
    ctry=[ctry;allstocks(i).indexlist(ind(1)).index];    
end
ctryList = unique(cellstr(ctry));
phi = zeros(566,11);

for i=1:566
    idx = find(isactivenow(i, :)==1);
    ind = find(cellfun('isempty',{allstocks(i).indexlist.index})==0); 
    [~,idx]=ismember(cellstr(allstocks(i).indexlist(ind(1)).index),ctryList);
    phi(i,idx) =1;
end


% Alpha 1: find arithmetic return
tri_T = size(tri,1);
retMat = tri(2:tri_T,:)./tri(1:(tri_T-1),:) - 1;

% short-term contrarian
T = size(retMat, 1);
n = size(retMat, 2);
alpharev_raw = zeros(tri_T, n);

triDecay = zeros(1,21);
for i = 0:20 
    triDecay(1,i+1) = (1/11) - (1/231) * i;
end

triDecay = flip(triDecay);

% for each column
for i = 1:n
    % for each row starting from 21
    window = 20;
    for j = 21:T
        alpharev_raw((j+1),i) = -1 * sum(triDecay * retMat((j-20):j,i) , 'omitnan');
    end
end

alpharev_raw = alpharev_raw * (eye(n) - rho * inv(transpose(rho) * rho) * transpose(rho));


% demean 
alpharevcxmean = mean(alpharev_raw, 2,'omitnan');
X_rev = alpharev_raw - alpharevcxmean;

% standardize 
alpharev = X_rev ./ std(X_rev, 0, 2, 'omitnan');
alpharev(alpharev > 3) = 3;
alpharev(alpharev < -3) = -3;


% Alpha2: SHORT-TERM PROCYCLICAL
alpharec_raw = zeros(tri_T, n);
triDecayAlpharec = zeros(1,45);

for i = 0:44 
    triDecayAlpharec(1,i+1) = (1/23) - (1/1035) * i;
end

triDecayAlpharec = flip(triDecayAlpharec);

% for each column
for i = 1:n
    % for each row starting from 21
    window = 44;
    for j = 45:T
        alpharec_raw((j+1),i) = -1 * sum(triDecayAlpharec * retMat((j- 44):j,i) , 'omitnan');
    end
end

% demean 
alpharecCXmean = mean(alpharec_raw, 2,'omitnan');
X_rec = alpharec_raw - alpharecCXmean;

% standardize
alpharec = X_rec ./ std(X_rec, 0, 2, 'omitnan');
alpharec(alpharec > 3) = 3;
alpharec(alpharec < -3) = -3;



% Alpha 3: Long-term Contrarian
% 1/mtbv: book to market ratio
alphaval = ones(1504, 566)./mtbv;
alphaval(isinf(alphaval)) = nan;
% winsorization
alphaval = alphaval - mean(alphaval, 2, 'omitnan');
alphaval = alphaval./std(alphaval, 0, 2, 'omitnan');
alphaval(alphaval > 3) = 3;
alphaval(alphaval < -3) = -3;

% Alpha 4: 
% Assuming 252 trading day per year, 21 trading day per month
% The first available date is day 253
% T-12 month ~ T-1 month return would be return from day 1 to day 231 
alphamom = tri(232:1483, :) ./ tri(1:1252, :)-1;
alphamom = alphamom - mean(alphamom, 2, 'omitnan');
alphamom = alphamom./std(alphamom, 0, 2, 'omitnan');
alphamom(alphamom > 3) = 3;
alphamom(alphamom < -3) = -3;
% To construct T*n matrix, backfill first year's alpha_4 with 0
alphamom = [zeros(252,566); alphamom];



% alpha_sum = 0.50 .* a1_windsorized + 0.25 .* a2_windsorized + 0.15 .*
% a3_windsorized + 0.10 .* a4_windsorized;

alpah_sum = 0.50 .* alpharev + 0.25 * alpharec + 0.15 * alphaval + 0.1 * alphamom;

% demean 
alpah_sum_mean = mean(alpah_sum, 2,'omitnan');
X_blend = alpah_sum - alpah_sum_mean;

% standardize
alphablend = X_blend ./ std(X_blend, 0, 2, 'omitnan');
alphablend(alphablend > 3) = 3;
alphablend(alphablend < -3) = -3;


alphablend(isnan(alphablend)) = 0;










