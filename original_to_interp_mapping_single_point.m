function [ind]=original_to_interp_mapping_single_point(point_time,time)
% INTERP_TO_ORIGINAL_MAPPING_SINGLE_POINT - Matches/maps a original time axis value
%                                           with the nearest/closest (based on time axis value)
%                                           interpolated time axis value and returns the corresponding index 
% 
% Usage:  [ind]=interp_to_original_mapping_single_point(point_time,time)
%
% Input arguments:   
%     point_time: A time value in the original data points (scaler)
%     time: Interpolated time axis values (vctor)      
% Output value:
%     ind:Index of the matched/mapped interpolated time value.      
% Other functions required:
%      None
%
% Farhan Asif Chowdhury, 2018
% aumyfarhan@gmail.com

ww=find((time-point_time)>0);

if (~isempty(ww))
    start_ind=max(min(ww),0);
else
    start_ind=1;
end

LL=length(time);

kk=find((time-point_time)<=0);

if (~isempty(kk))
    end_ind=min(max(kk),LL);

else
    end_ind=LL;
end


s_diff=abs(point_time-time(start_ind));
e_diff=abs(point_time-time(end_ind));

if (s_diff>e_diff)
    ind=end_ind;
else
    ind=start_ind;
end

