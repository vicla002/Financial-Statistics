function [m, w_new] = GMM_CKLS(params,data,w)

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
for i = 1:length(data)
f(:, 1) = (data-((gamma*delta_t+my)-sqrt(lambda*delta_t)));
f(:, 2) = f(var(data)-n
f(:, 3) = f(:, 1).^3;
f(:, 4) = f(:, 1).^4;
f(:, 5) = f(:,1).^5;
end
    for i = 1: length(data)
        
        s = s + f(i,:)'*f(i,:);
       
    end
    g = [mean(f(:,1)); mean(f(:,2)); mean(f(:,3)); mean(f(:,4)); mean(f(:,5))];
    m = g'*w* g;
  w_new = 1/length(data)*s;
    w_new = inv(w_new);
    
 

end
