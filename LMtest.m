

f = @(x) x(1).^2 + 6*x(1).*x(2) + x(2).^2;

grad_f = @(x) [2*x(1) + 6*x(2); 6*x(1) + 2*x(2)];

x = [-5;7];

for k = 1:10
   p = -inv([2 6;6 2]+6*eye(2))*grad_f(x);
   alpha = 1;
   while (f(x+alpha*p) - f(x)) > 0.8*alpha*dot(grad_f(x),p)
       alpha = alpha/2;
   end
   x = x + alpha*p
end