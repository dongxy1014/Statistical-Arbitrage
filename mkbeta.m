%equally-weighted market portfolio
ew_mkt=mean(retMat,2,'omitnan');

dt = datetime(myday,'InputFormat','dd-MMM-yyyy');
mon = month(dt);
firDay = lagmatrix(mon,1);
firstday = mon - firDay;
firInd = find(firstday ~= 0);

yr = 12;
mktbeta = [];
for i = 13:length(firInd)
    % find out the active stock based on the first day of the current month
    pastyr = retMat(firInd(i - yr):firInd(i),:);
    pastyr(isnan(pastyr)) = 0;
    
    row = size(pastyr, 1);
    col = size(pastyr, 2);
    
    x = [ones(row,1) ew_mkt(firInd(i - yr):firInd(i),:)];
    beta = inv(transpose(x) * x) * transpose(x) * pastyr;
    mktbeta = [mktbeta;beta(2,:)];
end









