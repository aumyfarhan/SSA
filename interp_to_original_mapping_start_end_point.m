function [start_ind,end_ind]=interp_to_original_mapping_start_end_point...
    (start_time,end_time,time)
% INTERP_TO_ORIGINAL_MAPPING_START_END_POINT - Maps interpolated time
%                                              values to the original
%                                              datapoints and returns the
%                                              corresponding indexes
% 
% Usage:  [start_ind,end_ind]=interp_to_original_mapping_start_end_point...
%    (start_time,end_time,time)
%
% Input arguments:   
%     start_time: starting time of a segment from the interpolated values.
%     end_time:   ending time of a segment from the interpolated values.
%     time: original time values
% Output value:
%     start_ind: Index of the matched/mapped original time value of start_time.
%     end_ind: Index of the matched/mapped original time value of end_time.  
% Other functions required:
%     None
%
% Farhan Asif Chowdhury, 2018
% aumyfarhan@gmail.com

ww=find((time-start_time)>0);
start_ind=min(ww);

kk=find((end_time-time)>0);
end_ind=max(kk);

