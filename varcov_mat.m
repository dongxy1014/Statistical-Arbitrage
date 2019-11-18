
% part 2 risk model 
y = datetime(myday,'InputFormat','dd-MMM-yyyy');
mon = month(y);
lagmon = lagmatrix(mon,1);
firstday = mon - lagmon;
index = find(firstday == 1 | firstday == -11 | isnan(firstday));

startIdx = index(13);
yrInterval = 12;
shrink = zeros(60,1);
shrinkIdx = 1;

for i = 13:length(index)
    active_now = find(isactivenow(index(i),:) == 1);
    % find out the active stock based on the first day of the current month
    past_1y = tri(index(i - yrInterval):index(i), active_now);
    
    % calculate arithmetic return
    T = size(past_1y,1);
    ret_1y = past_1y(2:T,:) ./ past_1y(1:(T-1),:) - 1;
    ret_1y(isnan(ret_1y)) = 0;
    
    dim = size(ret_1y);
    t = dim(1);
    n = dim(2);
    
    total_sum = zeros(n,n);
    
    for j = 1:t
        Xt = transpose(ret_1y(j,:));
        total_sum = total_sum + (Xt * transpose(Xt));
    end 
    % sampel variance matrix
    S = (1/t) * total_sum;
    
    sigmabar = sum(diag(S))/n;
    shrinkagetarget = sigmabar * eye(n);
    
    total_mse = 0;
    
    for k = 1:t
        Xt = transpose(ret_1y(k,:));
        total_mse = total_mse + sum((Xt * transpose(Xt) - S).^2);
    end
    betahat = 1 - ((1/(t^2 - t)) * sum(total_mse))/ sum(sum((S - shrinkagetarget).^2));
    shrink(shrinkIdx) = betahat;
    shrinkIdx = shrinkIdx + 1;
end


