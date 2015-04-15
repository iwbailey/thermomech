% test arrays

% stiffness matrix for row 1
a = [1,2,3,4,3,2,1; 
    2,3,4,5,4,3,2;
    4,5,6,7,6,5,4];

% slip def
b = [2,1,5,4;
    4,2,7,7;
    3,6,3,4];
[m,n] = size(b);

% desired answer
result1 = zeros( 1, size(b,2) );
for i=1:n,
    j1 = i;
    j2 = i+n-1;
    result1(i) = sum(sum(b.*a(:,j1:j2)));
end
fprintf('Result1:\n'); disp(result1)

% Get via conv2
tmp = conv2(a,rot90(b,2));
result2 = tmp(m,n:2*n-1);
fprintf('Result2:\n'); disp(result2)

% Get via fft2
b2 = [rot90(b,2), zeros(m,n-1)];
tmp = ifft2( fft2(a).*fft2(b2) );
result3 = tmp(m,n:2*n-1);
fprintf('Result3:\n'); disp(result3)