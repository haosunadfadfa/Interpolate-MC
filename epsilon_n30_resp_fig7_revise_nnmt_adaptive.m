% we leave o(b) term in calculation of bias 
close all;
clear;
warning('off');
L = 2;            % diameter of the area, in km (to simulate the source location)
% K = 10;            % number of sources   
sigmaa=1;
f = 5;            % frequency in kHz
af_dB = 0.11 * f^2 / (1 + f^2) + 44 * f^2 / (4100 + f^2) + 2.75 * 1e-4 * f^2 + 0.003;
% af_dB = 0.002 + 0.11 * f^2 / (1 + f^2) + 0.011 * f^2; % low frequency  model (< 500 Hz)
alpha = 1.5;
power = 1;
% Source signal model (for evaluation only)
std_ambient_noise = 0.02; %measurement noise
constant_noise= std_ambient_noise;
ii=30;delta=2;n=30;height=0.4;  %dimension of observation matrix
type_h = 'fixed';               % Type of bandwidth
%kernel ='uniform';
 kernel = 'epanechnikov';
%  kernel = 'gaussian';            % Type of kernel function
h_select = 'Leave-one-out CV';  % Bandwidth selection method
r=1;
C = 13;
% b=20;
% MM = [40 60 80 100 120 140 160 180 200];% 
% M1=MM(1);
% nn=50;
% hh2 = [0.9 0.9 0.8 0.7 0.6 0.5 0.5 0.5 0.4];
% hh1 = [0.3 0.2 0.2 0.2 0.2 0.2 0.2 0.1 0.1];
b=60;
% MM = [40 60 80 100 120 140 160 180 200];% 
MM = [40 100 160 220 280 340 400];% 
M_0 = 6;
M1=MM(1);
nn=50;
hh2 = [1.2 0.8 0.7 0.6 0.6 0.6 0.6];
hh1 = [0.3 0.2 0.2 0.2 0.1 0.1 0.1];
for index = 1 : length(MM)
    M = MM(index);
    h2=[hh2(index) hh2(index)];
    h1=[hh1(index) hh1(index)];

Gx=zeros(n+1,1);
Gy=zeros(n+1,1);
LH=L;
Gxinitial = (- LH / 2 + LH / (2 * n): LH / n: LH / 2) ; %grid center x
Gyinitial = (- LH / 2 + LH / (2 * n): LH / n: LH / 2) ;
Gxinitial=Gxinitial';
Gyinitial=Gyinitial';
Gx(1)=-1*L/2;
Gx(n+1)=1*L/2;
for i=2:n
    Gx(i)=Gxinitial(i-1)/2+Gxinitial(i)/2;
end
Gy(1)=-1*L/2;
Gy(n+1)=1*L/2;
for i=2:n
    Gy(i)=Gyinitial(i-1)/2+Gyinitial(i)/2;
end

bb=0.02;
Data = load('shad_200_1M400.mat');
sh = Data.Data.sh;
options.polytrend=2;

parfor ii=1:ii
% ======================== ENVIRONMENT SETUP ==============================
Z = Data.Data.location{ceil((M-40)/b)+1}{ii};
hZ = Data.Data.rss{ceil((M-40)/b)+1}{ii};
S =Data.Data.source{ceil((M-40)/b)+1}{ii};
K = length(S);E = zeros(n,n);V = zeros(n,n);

for hh=h1(1):bb:h2(1)
    h=[hh hh];
    count = 0;
    for i = 1:n
        for j=1:n
            [~,D] = knnsearch(Z,[Gxinitial(i) Gyinitial(j)],'K',M_0+1);
            if D(end)>h(1)
                continue
            end
            count = count+1;
        end
    end
    
    if count>(C*n*(log10(n))^2)
        h_start = hh;
        break
    end
end

%% global optimization
lpr = cell(round((h2(1)-h_start)/bb)+1,1);
error = cell(round((h2(1)-h_start)/bb)+1,1);
sr = cell(round((h2(1)-h_start)/bb)+1,1);
MSE_c= zeros(round((h2(1)-h_start)/bb)+1,1);
mse_b= zeros(round((h2(1)-h_start)/bb)+1,1);

for hh=h_start:bb:h2(1)
    h=[hh hh];
    LPR = zeros(n);
    quad=zeros(M,1);
    error_cpr = zeros(n);
    mse_c = 0;
    % count = 0;
    for i = 1:n
        for j=1:n
            [~,D] = knnsearch(Z,[Gxinitial(i) Gyinitial(j)],'K',M_0+1);
            if D(end)>h(1)
                continue
            end
            [LPR(i,j),~,KW,W] = LP(Z,hZ',1,type_h,kernel,h,[Gxinitial(i) Gyinitial(j)]);
            [~,beta] = LP(Z,hZ',2,type_h,kernel,h,[Gxinitial(i) Gyinitial(j)]);
            A = W' * KW * W;
            for iii=1 : M
                quad(iii)=(Z(iii,:)-[Gxinitial(i) Gyinitial(j)])*([beta(3) beta(4);beta(4) beta(5)])*(Z(iii,:)-[Gxinitial(i) Gyinitial(j)])';
            end
            E(i,j) = [1 0 0]*inv(A)*W'*KW*0.5*quad;
            V(i,j) = [1 0 0]*inv(A)*W'*(KW.*KW)*W*inv(A)*[1 0 0]'*std_ambient_noise^2;
           
            error_cpr(i,j) = delta*sqrt(abs(V(i,j)));
            LPR(i,j) = LPR(i,j)-E(i,j);
            mse_c = mse_c + E(i,j)^2 + V(i,j);
            
        end
    end
 
    lpr{round((hh-h_start)/bb)+1}=LPR;
    error{round((hh-h_start)/bb)+1}=error_cpr;
    % mse=mse/length(find(LPR));
    MSE_c(round((hh-h_start)/bb)+1) = mse_c/length(find(LPR));
end
[~,idx] = min(MSE_c);
H_0 = lpr{idx};
% sample(ii,:) = sr{idx};
err = error{idx};
h3 = (idx-1)*bb+h_start;
H_u = zeros(n);

for j=1:n
    for i=1:n
        hm=0;
        for k=1:K
            dmk = sqrt((norm([Gxinitial(i),Gyinitial(j)] - S(k, :)))^2+height^2);   % distance from sensor m to source k
            Amk = dmk ^ alpha * 10 ^ (- af_dB / 10)^dmk;
            Pmk = power/(Amk);
            hm = hm+Pmk;
         end
%         H_u(i,j) = hm;
        lx = floor((Gxinitial(i) + LH / 2) / (LH / nn)) + 1;
        ly = floor((Gyinitial(j) + LH / 2) / (LH / nn)) + 1;
        H_u(i,j) = hm*sh(lx,ly);
    end
end

[Hc6] = slumf_1st_mc_nn_bij...
    (Z, hZ,n,std_ambient_noise,Gxinitial,Gyinitial,type_h,kernel,delta,h1,h2,bb,h3);
mse_adaptive1(ii,:)=  (norm(Hc6-H_u,'fro')^2/(n^2));

end
ceil((M-M1)/b)+1

MSE_adaptive1(ceil((M-M1)/b)+1) = mean(mse_adaptive1)

end
% save('minor')
% save('result-n=30-s3.mat','MSE_adaptive1','MSE_lpr1','MSE_tps','MSE_kri')
% MSE_mc =[13.3068    10.5151    8.3663    7.9918    7.5735    6.4773 6.1138    5.9008  5.7811];  
% M=M1:b:M2;
figure 
% plot(MM, MSE_kri, '*-', 'LineWidth', 2, 'MarkerSize', 9)
% % hold on
% % plot(MM, MSE_mc, '^-','LineWidth', 2, 'MarkerSize', 9)
% hold on
% plot(MM, MSE_lpr1, 'x-','LineWidth', 2, 'MarkerSize', 9)
% hold on
% plot(MM, MSE_tps, 's-','LineWidth', 2, 'MarkerSize', 9)
% hold on
plot(MM, MSE_adaptive1, 'd-','LineWidth', 2, 'MarkerSize', 9)
% ylim([0,0.3])
set(gca, 'FontSize', 14);
xlabel('Number of Sensors M')
ylabel('Reconstruction MSE')
legend('k-LP','NNM-t')
set(gca, 'FontSize', 12);
BoxLineWidth = 2;
% Box
box on
set(gca, 'LineWidth', BoxLineWidth);
% Legend box
lgn = legend;
set(lgn, 'LineWidth', 1.5);
