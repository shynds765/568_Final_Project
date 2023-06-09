function J = costSum(x,Q,u,R,fp,t)
% Author: Brady Beck (Created: 02/05/2023, Updated: 02/06/2023)
% Function Description: Used to find total cost of maneuver 
% Inputs: x (deltaX) = X_ref - X_sc(1:6), nx6
%         u - control vector, nx3
%         Q - matrix, 6x6
%         R - matrix, 3x3 
%         fp - penalty factor
%         t -  time vector
% Output: J - total cost, vector (1) includes penalty, (2) does not  
    Q = diag(Q);
    R = diag(R);
    xtest = transpose(x(1,:));
    utest = u(:,1);
    dt = t(1) - 0;
    costp = (0.5*(transpose(xtest)*Q*xtest + transpose(utest)*R*utest) + fp(1))*dt;
    cost = (transpose(xtest)*Q*xtest + transpose(utest)*R*utest)*dt;
    for i = 2:size(x,1)
        xtest = transpose(x(i,:));
        utest = u(:,i);
        dt = t(i) - t(i-1);
        costp = costp + (0.5*(transpose(xtest)*Q*xtest + transpose(utest)*R*utest) + fp(i))*dt;
        cost = cost + (transpose(xtest)*Q*xtest + transpose(utest)*R*utest)*dt;  
    end
    
    J = [costp; cost];
return 