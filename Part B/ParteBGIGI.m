clear
close all
clc
%load data
FRF_data = load("FRF_H1.mat");
f_vect = FRF_data.f;
FRFs = FRF_data.FRF_z;

f_res = length(f_vect)/f_vect(end);

f_max_interest = 8; % Hz
f_min_interest = 2;
n_f_max = ceil(f_max_interest*f_res); % number of elements of the first 8 Hz
n_f_min = ceil(f_min_interest*f_res); % number of elements of the first 2 Hz

f_vect = f_vect(n_f_min:n_f_max); % Frequency vector of interest

FRFs = FRFs(n_f_min:n_f_max,:)';  % FRFs of interest (0-8 Hz)
[m,n]=size(FRFs);

% look for FRFs in antimodal positions
%compute max value of each FRF to select only the stronger ones
max_v=[];
for i=1:m
    max_v(i)=max(FRFs(i,:));
end
[maxes,maxes_i]=sort(max_v,'descend');
%select first 50 FRF
FRFi=maxes_i(1:100);
FRFs=FRFs(FRFi,:);
absavg= mean(abs(FRFs),1);
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
hold on
plot(f_vect, phavg, 'magenta', LineWidth=2)



