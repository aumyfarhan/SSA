function ssa_rc=ssa_decomposition(X,window)
% SSA_DECOMPOSITION - Performs SSA decomposition on the NaN removed, interpolated
%                     Signal derivative and returns the first five
%                     Reconstructed Components(RC) in a matrix
% 
% Usage:  ssa_rc=ssa_decomposition(X,window)
%
% Input arguments:   
%     X: Signal derivative values.
%     window: Lag window for building covariance matrix.
% Output value:
%     ssa_rc: The first five Reconstructed Components(RC) in a matrix ssa_rc
%             of size 5 by 'PBU Time Series length'. First row of ssa_rc contains first RC,
%             second row contains second RC and so on.       
% Other functions required:
%      None
%
% Farhan Asif Chowdhury, 2018
% aumyfarhan@gmail.com

N = length(X);   % length of time series


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

ssa_rc=[];
ssa_rc(1,:)=reconstructed_component(:,1);
ssa_rc(2,:)=reconstructed_component(:,2);
ssa_rc(3,:)=reconstructed_component(:,3);
ssa_rc(4,:)=reconstructed_component(:,4);
ssa_rc(5,:)=reconstructed_component(:,5);
