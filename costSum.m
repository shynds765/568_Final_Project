function J = costSum(x,Q,u,R,fp,tf)
% Author: Brady Beck (Created: 02/05/2023, Updated: 02/06/2023)
% Function Description: Used to find total cost of maneuver 
% Inputs: x (deltaX) = X_ref - X_sc(1:6), nx6
%         u - control vector, nx3
%         Q - matrix, 6x6
%         R - matrix, 3x3 
%         fp - penalty factor
%         tf - final time
% Output: J - total cost, vector (1) includes penalty, (2) does not  
    Q = diag(Q);
    R = diag(R);
    cost = 0;
    costp = 0;
    for i = 1:size(x,1)
        xtest = transpose(x(i,:));
        utest = u(:,i);
        costp = costp + (transpose(xtest)*Q*xtest + transpose(utest)*R*utest + fp(i));
        cost = cost + (transpose(xtest)*Q*xtest + transpose(utest)*R*utest);  
    end
    
    J = [(costp)*tf; cost*tf];
return 