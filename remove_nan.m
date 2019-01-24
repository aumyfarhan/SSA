function val_mod=remove_nan(val)
% REMOVE_NAN - Removes the data points where signal derivative is NaN
%
% Usage:  val_mod=remove_nan(val)
%
% Input arguments:   
%     val: A matrix of size 2 by 'PBU-Time Series length'. First row contains Time axis value
%          Second row contains Signal Derivative value 
% Output value:
%     val_mod: Datapoints that has signal derivative value = -999, are removed from time and signal derivative.A matrix of size 2 by 'PBU-Time Series length'.
%              First row contains Time axis value, Second row contains Signal Derivative value             
% Other functions required:
%      None
%
% Farhan Asif Chowdhury, 2018
% aumyfarhan@gmail.com

new_time=val(:,1);
new_derv=val(:,2);

derv_ind=find(new_derv~=-999);
time=new_time(derv_ind);
derv=new_derv(derv_ind);

val_mod=[];

val_mod(1:length(time),1)=time;
val_mod(1:length(derv),2)=derv;