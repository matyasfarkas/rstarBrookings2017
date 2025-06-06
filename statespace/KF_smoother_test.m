%% kalman_smooth.m
clear; clc;

% 1. Load model
TTT = readmatrix('TTT.csv');      % 91×91
RRR = readmatrix('RRR.csv');      % 91×29
CCC = readmatrix('CCC.csv')';      % 91×1
ZZ  = readmatrix('ZZ.csv');       % 20×91
DD  = readmatrix('DD.csv')';       % 20×1
QQ   = readmatrix('QQ.csv');% 29×29

Q= eye(size(QQ))*10^-10; 
Q([2:5 10 22:23], [2:5 10 22:23])=QQ([2:5 10 22:23], [2:5 10 22:23]);

% 2. Load data
TBL = readtable('dataset.csv');      
dates = datetime(TBL.date,'InputFormat','yyyy-MM-dd');
Y = TBL{:,2:end}';                % 20×T, NaN for missing
Y = 
% 3. Preallocate
n = size(TTT,1);  m = size(ZZ,1);  T = size(Y,2);
s = zeros(n,T);   P = zeros(n,n,T);
Ps= zeros(n,n,T);

% 4. Init diffuse prior
s(:,1)=0;  P(:,:,1)=eye(n)*1e6;

% 5. KF loop
for t=2:T
  % predict
  sp = TTT*s(:,t-1) + CCC;  
  Pp = TTT*P(:,:,t-1)*TTT' + RRR*Q*RRR';
  % update only if any obs
  y = Y(:,t);  idx = ~isnan(y);
  if any(idx)
    Zz = ZZ(:,idx)';  d = DD(idx);
    v = y(idx) - (Zz*sp + d);
    F = Zz*Pp*Zz';
    K = (Pp*Zz')/F;
    s(:,t) = sp + K*v;
    P(:,:,t) = Pp - K*Zz*Pp;
  else
    s(:,t) = sp;  P(:,:,t)=Pp;
  end
end

% 6. RTS smoother
ss = zeros(n,T);
ss(:,T)=s(:,T);  Ps(:,:,T)=P(:,:,T);
for t=T-1:-1:1
  Pp = TTT*P(:,:,t)*TTT' + RRR*Q*RRR';
  J  = P(:,:,t)*TTT'/Pp;
  ss(:,t)=s(:,t)+J*(ss(:,t+1)-TTT*s(:,t)-CCC);
  Ps(:,:,t)=P(:,:,t)+J*(Ps(:,:,t+1)-Pp)*J';
end

% 7. Extract smoothed obs_gdp
smGDP = ZZ(:,1)'*s + DD(1);

% 8. Plot & save
figure; plot(dates,smGDP,'LineWidth',1.5); hold on; plot(dates,Y(1,:));
datetick('x','yyyy-qq'), grid on
xlabel('Date'), ylabel('Smoothed obs\_gdp')

writetable(table(dates,smGDP','VariableNames',{'Date','smoothed_gdp'}),...
           'smoothed_gdp.csv');
disp('Done: smoothed_gdp.csv written.')
