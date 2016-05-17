filroot1 = ['AMF_USPFa_'] ;
filroot2 = ['_L2_GF_V008.nc'] ;
yr1 = 1995 ;
yr2 = 2014 ;
yrtot = yr1:yr2 ;
nyr = length(yrtot) ;

TA = nan(365*nyr, 1) ;
VPD = TA ;
CumPRCP = TA ;

varn = {'TA', 'VPD', 'PREC'} ;

for varind = 1:3 ;

for i =  1:length(yrtot);
    filin = [filroot1 num2str(yrtot(i)) filroot2] ;
    tem = ncread(filin, varn{varind}) ;
    if (mod(yrtot(i), 4) == 0) ;
        tem2 = nan(366, 1) ;
        for j = 1:366 ;
            ind = 24*(j-1)+[1:24] ;
            tem2(j) = mean(tem(ind)) ;
        end
        tem2 = tem2([1:59 61:366]') ;
    else
        tem2 = nan(365, 1) ;
        for j = 1:365 ;
            ind = 24*(j-1)+[1:24] ;
            tem2(j) = mean(tem(ind)) ;
        end
    end
    eval([varn{varind} '(365*(i-1)+[1:365], 1) = tem2 ;']) ;
end

end

% Gap Fill
ind = find(TA < -50) ;
TA(ind) = NaN ;
for i = 1:length(ind) ;
    ind2 = ind(i) ;
    TA(ind2) = nanmean([TA(ind2-2) TA(ind2-1) TA(ind2+1) TA(ind2+2)]) ;
end

%
ind = find(VPD < 0) ;
VPD(ind) = 0 ;
    
%
ind = find(PREC < 0) ;
PREC(ind) = 0 ;

%% Put the data together
data = [TA VPD PREC] ;

dpy = [31 28 31 30 31 30 31 31 30 31 30 31] ;
cdpy = cumsum(dpy) ;

datJJA = nan(nyr, 122, 3) ;
datJJA_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = (cdpy(5)+1):cdpy(9) ;
    datJJA(i,:,:) = data(ind+365*(i-1),:); 
    cumJJA(i,:,:) = sum(data(ind+365*(i-1),:));
    datJJA_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end
%%
%  Get data for dry years this is for problem #2
[tem, ind] = sort(datJJA_mean(:,3)) 
driest5years = ind(1:5) ;
wettest5years = ind(16:20) ;

driest10years = ind(1:10)
wettest10years = ind(11:20)

driestdailyprec = datJJA(driest5years, :, 3) ;
driestdailyprec = driestdailyprec(:) ;
wettestdailyprec = datJJA(wettest5years, :, 3) ;
wettestdailyprec = wettestdailyprec(:) ;
%Temperature
[tem, ind] = sort(datJJA_mean(:,1))
hottest5years = ind(1:5)
coolest5years = ind(16:20)
hottest10years = ind(1:10)
coolest10years = ind(11:20)

hottestdailytemp = datJJA(hottest5years, :, 1)
hottestdailytemp = hottestdailytemp(:) ;
coolestdailytemp = datJJA(coolest5years, :, 1)
coolestdailytemp = coolestdailytemp(:) ;

%%
%basisp = find(datJJA_mean(:,3) > mean(datJJA_mean(:,3)))
%basisn = find(datJJA_mean(:,3) < mean(datJJA_mean(:,3)))
%Np = length(basisp);
%Nn = length(basisn);

%  Plot results
%figure(3);
%clf;
%orient tall;


%  Set up a subplot in the figure window:
%subplot(1,1,1);

%  Plot a line, and +'s and o's for positive and negative ENSO events.
%pp = plot(yrtot, ind, 'k', ...    
scatter(ind, datJJA_mean(:,3) > mean(datJJA_mean(:,3)), 'filled')
%axis([1 20 -5 5])
%set(gca, 'XTick', 1:1:20, 'fontsize', 10);
%grid on
%%

%  Calculate correlation matrix for summer means (annually resolved)
C = corrcoef(datJJA_mean)

%  Calculate daily correlation coefficients for summer
tem = permute(datJJA, [3, 1, 2]) ;
tem = reshape(tem, 3, 92*20) ;
tem = tem' ;
corrcoef(tem)

%  For EOF analysis - use the 'tem' variable above for daily. Or, use the
%  annually resolved datJJA_mean for annual data this is for number 3

ntim = 20*92 ;
datEOF = (tem - ones(ntim,1)*mean(tem)) ./ (ones(ntim, 1)*std(tem)) ;

[u,s,v] = svd(datEOF) ;
lam = diag(s).^2 ;
per = lam./sum(lam) ;
eofs = v(:,1:3) ;
pcs = u(:,1:3) ;

%%
for i=1:3;
    y(i)=(lam(i,1)./sum(lam))*100; % calculate percent variance
    x(i)=i; % create an array 1:10 for plotting later
    e(i)= (sqrt(2/length(lam))*lam(i,1))./2; %compute error bars by "North's Rule"
end

subplot(1,1,1)
eigplot = errorbar(x,y,(e / (1E5)),'ok');
set(eigplot,'Color','b')
axis([0.5 3.5 0 60]);
title('Precipitation Eigenspectrum of First 3 Modes','FontSize',20);
xlabel('Eigenvalue','FontSize',16);
ylabel('% Variance Explained','FontSize',16);


%%

% Pull out precip
p = datJJA(:, :, 3) ;
p = p(:) ;

% Get histogram for a particular year this is for number 1
ind = [1:20] ; %9(2003) 10(2004) 15(2009) 13(2007) 12(2006)  Drought Years
%ind = 1:20 ;
p = datJJA(ind, :, 3) ;
%p = p(p>0)
%histfit(p,8,'gamma') ;

% Pull out VPD
v = datJJA(:, :, 2) ;
v = v(:) ;

% Get histogram for a particular year
ind = [1:20] ; %
%ind = 1:20 ;
v = datJJA(ind, :, 2) ;
%histfit(v,8,'gamma') ;

% Pull out Temp
t = datJJA(:, :, 1) ;
t = t(:) ;

% Get histogram for a particular year
ind = [1:20] ; 
%ind = 1:20 ;
t = datJJA(ind, :, 1) ;
%histfit(t,8,'gamma') ;

%scatter(v,p,'filled')
