function [obj, cs, pi] = func_foc_costs(x,param,cost)

p1 = x(1);
p2 = x(2);

q = param(1);
alpha = param(2);
beta = param(3);
sigma = param(4);

fees1 = param(5);
fees2 = 1-param(6);

c1 = cost(1);
c2 = cost(2);

c1 = c1 + fees1;
c2 = c2 + fees1;

aa1 = exp((q-beta*p1)/sigma);
aa2 = exp((-beta*p2)/sigma);
a1 = aa1 / (aa1+aa2);
a2 = aa2 / (aa1+aa2);

ap = -beta/sigma * a1*a2;
app1 = beta^2/sigma^2 * a1*a2*(a2-a1);
app2 = beta^2/sigma^2 * a1*a2*(a1-a2);

% figure out the case for seller 1
cond1 = (p1 < q/alpha);

if cond1==0
    foc1 = fees2*a1*(1+q-alpha*p1)+(fees2*p1-c1)*ap*(1+q-alpha*p1)-alpha*(fees2*p1-c1)*a1;
end

if cond1==1
    foc1 = fees2*a1+(fees2*p1-c1)*ap;
end

foc2 = a2*(1-alpha*p2)+(fees2*p2-c2)*ap*(1-alpha*p2)-alpha*(fees2*p2-c2)*a2;

obj = foc1^2+foc2^2;

if cond1==0
    soc1 = fees2*ap*(1+q-alpha*p1)-fees2*alpha*a1+fees2*ap*(1+q-alpha*p1)+(fees2*p1-c1)*app1*(1+q-alpha*p1) ...
        -alpha*(fees2*p1-c1)*ap-fees2*alpha*a1-alpha*(fees2*p1-c1)*ap;
end

if cond1==1
    soc1 = 2*fees2*ap+(fees2*p1-c1)*app1;
end

soc2 = fees2*ap*(1-alpha*p2)-fees2*alpha*a2+fees2*ap*(1-alpha*p2)+(fees2*p2-c2)*app2*(1-alpha*p2) ...
    -alpha*(fees2*p2-c2)*ap-fees2*alpha*a2-alpha*(fees2*p2-c2)*ap;

% rule out unrealistic prices
if fees2*p1<c1 || fees2*p2<c2 || p1>(1+q+alpha*c1)/(2*alpha) || p2>(1+alpha*c2)/(2*alpha) || soc1>=0 || soc2>=0
  obj = 1000;
end

d1 = 1+q-alpha*p1;
d2 = 1-alpha*p2;
    
cs = a1*(cond1*(q/alpha-p1 + 1/(2*alpha))+(1-cond1)*(((1+q)/alpha-p1)*d1/2))+a2*(1/alpha-p2)*d2/2;
pi = fees1*(a1*(cond1+(1-cond1)*d1)+a2*d2)+(1-fees2)*(a1*(cond1*p1+(1-cond1)*p1*d1)+a2*p2*d2);