function [m, w_new] = GMM2(params,data,w)

gamma = params(1);
sigma = params(2);
lambda = params(3);
my = params(4);
delta = params(5);
delta_t = 1/12;
N = length(data);
f = zeros(N, 5);
s=0;
%epsilon = yN - x - gamma * (theta - x);

 for i = 1: length(data)
f(i, 1) = (data(i)-(gamma+my*lambda)*delta_t);
f(i, 2) = data(i)^2-delta_t*(sigma^2+(my^2+delta^2)*lambda);
f(i, 3) = data(i)^3-delta_t*(lambda*my*(my^2)+3*delta^2);
f(i, 4) = data(i)^4-delta_t*(my^4*lambda+6*my^2*gamma^2*lambda+3*delta^4*lambda)+3*delta_t^2*(sigma^2+(my^2+delta^2)*lambda)^2;
%f(:, 5) = f(:,1).^5;

        
        s = s + f(i,:)'*f(i,:);
       
 end
    g = [mean(f(:,1)); mean(f(:,2)); mean(f(:,3)); mean(f(:,4))];
    m = g'*w* g;
  w_new = 1/length(data)*s;
    w_new = inv(w_new);
    
 

end