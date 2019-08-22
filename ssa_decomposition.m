function ssa_rc = ssa_decomposition(X, window, ssa_no)
% SSA_DECOMPOSITION - Performs SSA decomposition on the NaN removed, interpolated
%                     Signal derivative and returns the first five
%                     Reconstructed Components(RC) in a matrix
% 
% Usage:  ssa_rc=ssa_decomposition(X, window, ssa_no)
%
% Input arguments:   
%     X: Signal derivative values.
%     window: Lag window for building covariance matrix.
% Output value:
%     ssa_rc: The first five Reconstructed Components(RC) in a matrix ssa_rc
%             of size 'ssa_no' by 'PBU Time Series length'. First row of ssa_rc contains first RC,
%             second row contains second RC and so on.
%     ssa_no: Number of SSA reconstructed components to return (from begining).
% Other functions required:
%      None
%
% Farhan Asif Chowdhury, 2018
% aumyfarhan@gmail.com

N = length(X);   % length of time series

if ~exist('ssa_no','var')
    ssa_no = 5;
end
    
% if the window size is greater than the time series length
% then reduce window length

if (window>=N)  
    window=5*floor((N-5)/5);
end

trajectory_matrix=zeros(N-window+1,window);

for i=1:window
  trajectory_matrix(:,i) = X((1:N-window+1)+i-1);
end

covariance_matrix=trajectory_matrix'*trajectory_matrix / (N-window+1);

[eigen_vector,eigen_value] = eig(covariance_matrix);

eigen_value = diag(eigen_value);               
[eigen_value,ind]=sort(eigen_value,'descend'); 
eigen_vector = eigen_vector(:,ind);                  

principal_component = trajectory_matrix*eigen_vector;

reconstructed_component=zeros(N,window);

for i=1:window
  temp=principal_component(:,i)*eigen_vector(:,i)'; 
  temp=temp(end:-1:1,:);
  for n=1:N 
    reconstructed_component(n,i)=mean( diag(temp,-(N-window+1)+n) );
  end
end

ssa_rc = reconstructed_component(:,1:ssa_no)';
