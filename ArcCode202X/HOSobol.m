function [X] = HOSobol(m,s,d)
% Higher order Sobol sequence
% Create a higher order Sobol sequence.
% 2^m number of points
% s dimension of final point set
% d interlacing factor
% X Output Sobol sequence

N = pow2(m); % Number of points;
P = sobolset(d*s); % Get Sobol sequence;
sobolpoints = net(P,N); % Get net from Sobol sequence with N points;

% Create binary representation of digits;

W = sobolpoints* N;
Z = transpose(W);
Y = zeros(s,N);
for j = 1:s,
for i = 1:m,
for k = 1:d
Y(j,:) = bitset( Y(j,:),(m*d+1) - k - (i-1)*d,bitget( Z((j-1)*d+k,:),(m+1) - i));
end;
end;
end;
Y = Y * pow2(-m*d);

X=transpose(Y); % X is matrix of higher order Sobol points,
% where the number of columns equals the dimension
% and the number of rows equals the number of points;

end