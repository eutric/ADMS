clear
close all
% close all
clc
%load data
FRF_data = load("FRF_H1.mat");
f_vect = FRF_data.f;
FRFs_og = FRF_data.FRF_z;

f_res = length(f_vect)/f_vect(end);

f_max_interest = 8; % Hz
f_min_interest = 2;
n_f_max = ceil(f_max_interest*f_res); % number of elements of the first 8 Hz
n_f_min = ceil(f_min_interest*f_res); % number of elements of the first 2 Hz

f_vect = f_vect(n_f_min:n_f_max); % Frequency vector of interest
om_vect=f_vect*2*pi;
FRFs_og = FRFs_og(n_f_min:n_f_max,:)';  % FRFs of interest (0-8 Hz)
[m,n]=size(FRFs_og);

% look for FRFs in antimodal positions
%compute max value of each FRF to select only the stronger ones
max_v=[];
for i=1:m
    max_v(i)=max(FRFs_og(i,:));
end
[maxes,maxes_i]=sort(max_v,'descend');
%select first 50 FRF
FRFi=maxes_i(1:100);
FRFs=FRFs_og(FRFi,:);
absavg=mean(abs(FRFs),1);
phavg=mean(angle(FRFs),1);
figure
subplot(2,1,1)
semilogy(f_vect, abs(FRFs), 'cyan', LineWidth=.01)
grid on
xlim([2,8])
hold on
semilogy(f_vect, absavg, 'magenta', LineWidth=2)

subplot(2,1,2)
plot(f_vect, angle(FRFs), 'cyan',LineWidth=.01)
xlim([2,8])
grid on


%% look fer peaks to identify the resonant frequencies

[pks,pks_i]=findpeaks(absavg,"MinPeakWidth",5,"Threshold",2.05*1e-5);%using this criterias only the relevant peaks are found
subplot(2,1,1)
plot(f_vect(pks_i),pks,'bo',LineWidth=2)
res_f=f_vect(pks_i);
for i = 1:length(res_f)
    rf=res_f(i);%resonant frequency
    rfi=pks_i(i);%index of resonsnt freq in f_vect
    i1=find(absavg(1:rfi)<absavg(rfi)/sqrt(2),1,'last');
    subplot(2,1,1)
    plot(f_vect(i1),absavg(i1),'o')
    i2=find(absavg(rfi:end)<absavg(rfi)/sqrt(2),1,'first')+rfi;
    subplot(2,1,1)
    plot(f_vect(i2),absavg(i2),'o')
    % compute f1 and f2
    om1(i)=om_vect(i1);
    om2(i)=om_vect(i2);
    om0(i)=om_vect(rfi);
end
h_vect=(om2.^2-om1.^2)./4./om0.^2;
%for each FRF the mode shape are scaled of a factor, fixed the input
%position in k, the mode shape  for FRF_kj relative to a i mode
for i =1:length(res_f)
    for j=1:m
        unnnormed_mode(i,j)=-imag(2*h_vect(i).*(om_vect(pks_i(i))).^2*FRFs_og(j,pks_i(i)));
        
    end
end








