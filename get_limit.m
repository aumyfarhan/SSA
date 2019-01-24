function [xlv,ylv]=get_limit(time,derv)
% GET_LIMIT - Calculates the minimum, maximum of time axis value
%             and signal derivative value for plot axis limit
% 
% Usage:  [xlv,ylv]=get_limit(time,derv)
%
% Input arguments:   
%     time: Time Axis Values
%     derv: Signal derivative values     
% Output value:
%     xlv: A 1 by 2 vector. [minimum x-axis(time) value, maximum x-axis(time) value]
%     ylv: A 1 by 2 vector. [minimum y-axis(signal derivative) value, maximum y-axis(signal derivative) value]
% Other functions required:
%      None
%
% Farhan Asif Chowdhury, 2018
% aumyfarhan@gmail.com

if (min(time)>-3.2 && max(time)+1<2.2)
    xlv=[-3.2,2.2];
else
    xlv=[floor(min(time))-0.2,max(floor(min(time))+5,ceil(max(time)))+0.2];
end


ylv=[min(-0.2,floor(min(derv)))-0.2,max(3,ceil(max(derv)))+0.2];
    