function val_mod=interpolate_data(val,sampling_rate)
% INTERPOLATE_DATA - Performs linear interpolation on the NaN removed Data
% to make the data uniformly sampled using given sampling rate
% 
% Usage:  val_mod=interpolate_data(val,sampling_rate)
%
% Input arguments:   
%     val: A matrix of size 2 by 'PBU-Time Series length'. 
%          First row contains Time axis value, Second row contains Signal Derivative value
%     sampling_rate: The time axis value difference between two consecutive data point    
% Output value:
%     val_mod:  A matrix of size 2 by 'PBU-Time Series length'. Time axis values has been linearly 
%               interpolated using given sampling rate. First row contains Time axis value,
%               Second row contains Signal Derivative value
% Other functions required:
%      None
%
% Farhan Asif Chowdhury, 2018
% aumyfarhan@gmail.com

time1=val(:,1);
derv1=val(:,2);

[c,ia,ic]=unique(time1,'stable');

time=time1(ia);
derv=derv1(ia);

time_interp=time(1):sampling_rate:time(end);
derv_interp= interp1(time,derv,time_interp);

val_mod=[];

val_mod(1:length(time_interp),1)=time_interp;
val_mod(1:length(derv_interp),2)=derv_interp;
