filroot1 = ['EFDC_L2_Flx_CHCha_'] ;
filroot2 = ['_L2_GF_V008.nc'] ;
yr1 = 1995 ;
yr2 = 2014 ;
yrtot = yr1:yr2 ;
nyr = length(yrtot) ;
mtot = 12 ; 
nym = 1:12 ;

TA = nan(365*nyr, 1) ;
VPD = TA ;
PREC = TA ;

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

%% January
datJan = nan(nyr, 31, 3) ;
datJan_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 1:31 ;
    datJan(i,:,:) = data(ind+365*(i-1),:); 
    cumJan(i,:,:) = sum(data(ind+365*(i-1),:));
    datJan_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    JanT(i) = mean(datJan(i,:,1))
end
JanTmean = mean(JanT)

for i = 1:20
    JanV(i) = mean(datJan(i,:,2))
end
JanVmean = mean(JanV)

for i = 1:20
    JanP(i) = mean(datJan(i,:,3))
end
JanPmean = mean(JanP)

%% February Data
datFeb = nan(nyr, 28, 3) ;
datFeb_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 32:59 ;
    datFeb(i,:,:) = data(ind+365*(i-1),:); 
    cumFeb(i,:,:) = sum(data(ind+365*(i-1),:));
    datFeb_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    FebT(i) = mean(datFeb(i,:,1))
end
FebTmean = mean(FebT)

for i = 1:20
    FebV(i) = mean(datFeb(i,:,2))
end
FebVmean = mean(FebV)

for i = 1:20
    FebP(i) = mean(datFeb(i,:,3))
end
FebPmean = mean(FebP)

%% March Data
datMar = nan(nyr, 31, 3) ;
datMar_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 60:90 ;
    datMar(i,:,:) = data(ind+365*(i-1),:); 
    cumMar(i,:,:) = sum(data(ind+365*(i-1),:));
    datMar_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    MarT(i) = mean(datMar(i,:,1))
end
MarTmean = mean(MarT)

for i = 1:20
    MarV(i) = mean(datMar(i,:,2))
end
MarVmean = mean(MarV)

for i = 1:20
    MarP(i) = mean(datMar(i,:,3))
end
MarPmean = mean(MarP)

%% April Data

datApr = nan(nyr, 30, 3) ;
datApr_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 91:120 ;
    datApr(i,:,:) = data(ind+365*(i-1),:); 
    cumApr(i,:,:) = sum(data(ind+365*(i-1),:));
    datApr_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    AprT(i) = mean(datApr(i,:,1))
end
AprTmean = mean(AprT)

for i = 1:20
    AprV(i) = mean(datApr(i,:,2))
end
AprVmean = mean(AprV)

for i = 1:20
    AprP(i) = mean(datApr(i,:,3))
end
AprPmean = mean(AprP)

%% May Data

datMay = nan(nyr, 31, 3) ;
datMay_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 121:151 ;
    datMay(i,:,:) = data(ind+365*(i-1),:); 
    cumMay(i,:,:) = sum(data(ind+365*(i-1),:));
    datMay_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    MayT(i) = mean(datMay(i,:,1))
end
MayTmean = mean(MayT)

for i = 1:20
    MayV(i) = mean(datMay(i,:,2))
end
MayVmean = mean(MayV)

for i = 1:20
    MayP(i) = mean(datMay(i,:,3))
end
MayPmean = mean(MayP)

%% June Data

datJun = nan(nyr, 30, 3) ;
datJun_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 152:181 ;
    datJun(i,:,:) = data(ind+365*(i-1),:); 
    cumJun(i,:,:) = sum(data(ind+365*(i-1),:));
    datJun_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    JunT(i) = mean(datJun(i,:,1))
end
JunTmean = mean(JunT)

for i = 1:20
    JunV(i) = mean(datJun(i,:,2))
end
JunVmean = mean(JunV)

for i = 1:20
    JunP(i) = mean(datJun(i,:,3))
end
JunPmean = mean(JunP)

%% July Data

datJul = nan(nyr, 31, 3) ;
datJul_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 182:212 ;
    datJul(i,:,:) = data(ind+365*(i-1),:); 
    cumJul(i,:,:) = sum(data(ind+365*(i-1),:));
    datJul_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    JulT(i) = mean(datJul(i,:,1))
end
JulTmean = mean(JulT);

for i = 1:20
    JulV(i) = mean(datJul(i,:,2))
end
JulVmean = mean(JulV)

for i = 1:20
    JulP(i) = mean(datJul(i,:,3))
end
JulPmean = mean(JulP)

%% August Data

datAug = nan(nyr, 31, 3) ;
datAug_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 213:243 ;
    datAug(i,:,:) = data(ind+365*(i-1),:); 
    cumAug(i,:,:) = sum(data(ind+365*(i-1),:));
    datAug_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    AugT(i) = mean(datAug(i,:,1))
end
AugTmean = mean(AugT)

for i = 1:20
    AugV(i) = mean(datAug(i,:,2))
end
AugVmean = mean(AugV)

for i = 1:20
    AugP(i) = mean(datAug(i,:,3))
end
AugPmean = mean(AugP)

%% September Data

datSep = nan(nyr, 30, 3) ;
datSep_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 244:273 ;
    datSep(i,:,:) = data(ind+365*(i-1),:); 
    cumSep(i,:,:) = sum(data(ind+365*(i-1),:));
    datSep_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    SepT(i) = mean(datSep(i,:,1))
end
SepTmean = mean(SepT)

for i = 1:20
    SepV(i) = mean(datSep(i,:,2))
end
SepVmean = mean(SepV)

for i = 1:20
    SepP(i) = mean(datSep(i,:,3))
end
SepPmean = mean(SepP)

%% October Data

datOct = nan(nyr, 31, 3) ;
datOct_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 274:304 ;
    datOct(i,:,:) = data(ind+365*(i-1),:); 
    cumOct(i,:,:) = sum(data(ind+365*(i-1),:));
    datOct_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    OctT(i) = mean(datOct(i,:,1))
end
OctTmean = mean(OctT)

for i = 1:20
    OctV(i) = mean(datOct(i,:,2))
end
OctVmean = mean(OctV)

for i = 1:20
    OctP(i) = mean(datOct(i,:,3))
end
OctPmean = mean(OctP)

%% November Data

datNov = nan(nyr, 30, 3) ;
datNov_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 305:334 ;
    datNov(i,:,:) = data(ind+365*(i-1),:); 
    cumNov(i,:,:) = sum(data(ind+365*(i-1),:));
    datNov_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    NovT(i) = mean(datNov(i,:,1))
end
NovTmean = mean(NovT)

for i = 1:20
    NovV(i) = mean(datNov(i,:,2))
end
NovVmean = mean(NovV)

for i = 1:20
    NovP(i) = mean(datNov(i,:,3))
end
NovPmean = mean(NovP)

%% December Data

datDec = nan(nyr, 31, 3) ;
datDec_mean = nan(nyr, 3) ;
for i = 1:nyr ;
    ind = 335:365 ;
    datDec(i,:,:) = data(ind+365*(i-1),:); 
    cumDec(i,:,:) = sum(data(ind+365*(i-1),:));
    datDec_mean(i,:) = mean(data(ind+365*(i-1),:)) ;
end

for i = 1:20
    DecT(i) = mean(datDec(i,:,1))
end
DecTmean = mean(DecT)

for i = 1:20
    DecV(i) = mean(datDec(i,:,2))
end
DecVmean = mean(DecV)

for i = 1:20
    DecP(i) = mean(datDec(i,:,3))
end
DecPmean = mean(DecP)
%% Lowest 5 and Highest 5 years for each variable. lv = low vpd, hv = high vpd
%Jan Prec
[tem, ind] = sort(JanP(:)) ;
Jandry = ind(1:5) ;
Janwet = ind(16:20) ;

%Jan Temperature
[tem, ind] = sort(JanT(:)) ;
Jancold = ind(1:5)
Janhot = ind(16:20)

%Jan VPD
[tem, ind] = sort(JanV(:))
Janlv = ind(1:5)
Janhv = ind(16:20)

%Feb Prec
[tem, ind] = sort(FebP(:)) 
Febdry = ind(1:5) ;
Febwet = ind(16:20) ;

%Feb Temperature
[tem, ind] = sort(FebT(:))
Febcold = ind(1:5)
Febhot = ind(16:20)

%Feb VPD
[tem, ind] = sort(FebV(:))
Feblv = ind(1:5)
Febhv = ind(16:20)

%Mar Prec
[tem, ind] = sort(MarP(:)) 
Mardry = ind(1:5) ;
Marwet = ind(16:20) ;

%Mar Temperature
[tem, ind] = sort(MarT(:))
Marcold = ind(1:5)
Marhot = ind(16:20)

%Mar VPD
[tem, ind] = sort(MarV(:))
Marlv = ind(1:5)
Marhv = ind(16:20)

%Apr Prec
[tem, ind] = sort(AprP(:)) 
Aprdry = ind(1:5) ;
Aprwet = ind(16:20) ;

%Apr Temperature
[tem, ind] = sort(AprT(:))
Aprcold = ind(1:5)
Aprhot = ind(16:20)

%Apr VPD
[tem, ind] = sort(AprV(:))
Aprlv = ind(1:5)
Aprhv = ind(16:20)

%May Prec
[tem, ind] = sort(MayP(:)) 
Maydry = ind(1:5) ;
Maywet = ind(16:20) ;

%May Temperature
[tem, ind] = sort(MayT(:))
Maycold = ind(1:5)
Mayhot = ind(16:20)

%May VPD
[tem, ind] = sort(MayV(:))
Maylv = ind(1:5)
Mayhv = ind(16:20)

%Jun Prec
[tem, ind] = sort(JunP(:)) 
Jundry = ind(1:5) ;
Junwet = ind(16:20) ;

%Jun Temperature
[tem, ind] = sort(JunT(:))
Juncold = ind(1:5)
Junhot = ind(16:20)

%Jun VPD
[tem, ind] = sort(JunV(:))
Junlv = ind(1:5)
Junhv = ind(16:20)

%Jul Prec
[tem, ind] = sort(JulP(:)) 
Juldry = ind(1:5) ;
Julwet = ind(16:20) ;

%Jul Temperature
[tem, ind] = sort(JulT(:))
Julcold = ind(1:5)
Julhot = ind(16:20)

%Jul VPD
[tem, ind] = sort(JulV(:))
Jullv = ind(1:5)
Julhv = ind(16:20)

%Aug Prec
[tem, ind] = sort(AugP(:)) 
Augdry = ind(1:5) ;
Augwet = ind(16:20) ;

%Aug Temperature
[tem, ind] = sort(AugT(:))
Augcold = ind(1:5)
Aughot = ind(16:20)

%Aug VPD
[tem, ind] = sort(AugV(:))
Auglv = ind(1:5)
Aughv = ind(16:20)

%Sep Prec
[tem, ind] = sort(SepP(:)) 
Sepdry = ind(1:5) ;
Sepwet = ind(16:20) ;

%Sep Temperature
[tem, ind] = sort(SepT(:))
Sepcold = ind(1:5)
Sephot = ind(16:20)

%Sep VPD
[tem, ind] = sort(SepV(:))
Seplv = ind(1:5)
Sephv = ind(16:20)

%Oct Prec
[tem, ind] = sort(OctP(:)) 
Octdry = ind(1:5) ;
Octwet = ind(16:20) ;

%Oct Temperature
[tem, ind] = sort(OctT(:))
Octcold = ind(1:5)
Octhot = ind(16:20)

%Oct VPD
[tem, ind] = sort(OctV(:))
Octlv = ind(1:5)
Octhv = ind(16:20)

%Nov Prec
[tem, ind] = sort(NovP(:)) 
Novdry = ind(1:5) ;
Novwet = ind(16:20) ;

%Nov Temperature
[tem, ind] = sort(NovT(:))
Novcold = ind(1:5)
Novhot = ind(16:20)

%Nov VPD
[tem, ind] = sort(NovV(:))
Novlv = ind(1:5)
Novhv = ind(16:20)

%Dec Prec
[tem, ind] = sort(DecP(:)) 
Decdry = ind(1:5) ;
Decwet = ind(16:20) ;

%Dec Temperature
[tem, ind] = sort(DecT(:))
Deccold = ind(1:5)
Dechot = ind(16:20)

%Dec VPD
[tem, ind] = sort(DecV(:))
Declv = ind(1:5)
Dechv = ind(16:20)

%% Plotting

%20 year monthly averaged variable distribution
boxplot( [JanT(:), FebT(:), MarT(:), AprT(:), MayT(:), JunT(:), JulT(:), AugT(:), SepT(:), OctT(:), NovT(:), DecT(:)],'notch', 'on', 'color', 'rb')
boxplot( [JanP(:), FebP(:), MarP(:), AprP(:), MayP(:), JunP(:), JulP(:), AugP(:), SepP(:), OctP(:), NovP(:), DecP(:)],'notch', 'on', 'color', 'rb')
boxplot( [JanV(:), FebV(:), MarV(:), AprV(:), MayV(:), JunV(:), JulV(:), AugV(:), SepV(:), OctV(:), NovV(:), DecV(:)],'notch', 'on', 'color', 'rb')

%June Spread
subplot(3,1,1)
plot(JunT(1:20), 'r')
title('June Temperatures')
xlabel('Years since 1994')
ylabel('Degrees C')
subplot(3,1,2)
plot(JunP(1:20), 'g')
title('June Precipitation')
xlabel('Years since 1994')
ylabel('mm')
subplot(3,1,3)
plot(JunV(1:20), 'b')
title('June VPD')
xlabel('Years since 1994')
ylabel('kPa')

%June Temperature with Mean reference line
plot(JunT(:), 'b')
hline = refline([0 JunTmean])
hline.Color = 'r'

