%% kalman_smooth.m
clear; clc;
 % 2. Load data
    TBL = readtable('statespace/dataset.csv');
    dates = datetime(TBL.date,'InputFormat','yyyy-MM-dd');
    Y = TBL{:,2:end}';                % 20×T, NaN for missing
    Y_onlyFIN = nan(size(Y));
    Y_onlyFIN([6 9:11 14],:) = Y([6 9:11 14],:);
% 1. Load model
tic
for i =1:100
   warning off;
   if mod(i,10) ==0 
       disp(['Running posterior draw:' num2str(i)]); toc;
   end
    TTT = readmatrix(['statespace/' num2str(i) '_TTT.csv']);      % 91×91
    RRR = readmatrix(['statespace/' num2str(i) '_RRR.csv']);      % 91×29
    CCC = readmatrix(['statespace/' num2str(i) '_CCC.csv'])';      % 91×1
    ZZ  = readmatrix(['statespace/' num2str(i) '_ZZ.csv']);       % 20×91
    DD  = readmatrix(['statespace/' num2str(i) '_DD.csv'])';       % 20×1
    QQ   = readmatrix(['statespace/' num2str(i) '_QQ.csv']);% 29×29
    
    Q= eye(size(QQ))*10^-10;
    Q([2:5 10 22:23], [2:5 10 22:23])=QQ([2:5 10 22:23], [2:5 10 22:23]);
    
   
    % 3. Preallocate
    n = size(TTT,1);  m = size(ZZ,1);  T = size(Y,2);
    s = zeros(n,T);   P = zeros(n,n,T);
    Ps= zeros(n,n,T);
    
    % 4. Init diffuse prior
    s(:,1)=0;  P(:,:,1)=eye(n)*1e6;
    
    for timet = 5:T
        Y_rt = [Y(:,1:timet-4) Y_onlyFIN(:,timet-3:timet)]; 

        % 5. KF loop
        for t=2:timet
            % predict
            sp = TTT*s(:,t-1) + CCC;
            Pp = TTT*P(:,:,t-1)*TTT' + RRR*Q*RRR';
            % update only if any obs
            y = Y_rt(:,t);  idx = ~isnan(y);
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
        ss = zeros(n,timet+1);
        ss(:,timet+1)=zeros(size(s(:,t)));  Ps(:,:,timet+1)=P(:,:,timet);
        for t=timet:-1:1
            Pp = TTT*P(:,:,t)*TTT' + RRR*Q*RRR';
            J  = P(:,:,t)*TTT'/Pp;
            ss(:,t)=s(:,t)+J*(ss(:,t+1)-TTT*s(:,t)-CCC);
            Ps(:,:,t)=P(:,:,t)+J*(Ps(:,:,t+1)-Pp)*J';
        end
        smGDP(i,timet) = ZZ(:,1)'*ss(:,timet) + DD(1);

    end
    % 7. Extract smoothed obs_gdp
end
% 8. Plot & save
% [num_draws, T] = size(smGDP);
% smoothed_gdp_yoy = NaN(num_draws, T);  % initialize
% 
% % Convert to quarterly growth factors
% gdp_growth_factors = 1 + smGDP/100 ;
% 
% % Compute YoY growth
% for t = 5:T
%     smoothed_gdp_yoy(:, t) = (prod(gdp_growth_factors(:, t-3:t), 2) - 1)*100 ;
% end

quantiles = prctile(smGDP, [15 50 85], 1);

% Plot
figure;
hold on;
% Plot shaded confidence interval
fill([dates; flipud(dates)], ...
     [quantiles(1,:)'; flipud(quantiles(3,:)')], ...
     [0.8 0.8 1], 'EdgeColor', 'none');
plot(dates, quantiles(2,:), 'b-', 'LineWidth', 2); % Median
plot(dates, Y(1,:),'-k','LineWidth',3);
xlabel('Time');
ylabel('Smoothed gd[p\_obs]');
title('RTS Smoothed State with Parameter Uncertainty');
legend('70% CI', 'Median');
xlim([dates(160) dates(end)])
datetick('x', 'yyyy-QQ');

grid on;


writetable(table(dates,smGDP','VariableNames',{'Date','smoothed_gdp'}),...
    'smoothed_gdp.csv');
disp('Done: smoothed_gdp.csv written.')
