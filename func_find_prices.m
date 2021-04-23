function [peqm,cond1out] = func_find_prices(pnum,param,cost)


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

pmin = [c1/fees2;c2/fees2];
pmax = [max((1+q+alpha*c1)/(2*alpha),(1+q)/alpha);1/alpha];

pstep = (pmax-pmin)/pnum;

p1vec = (pmin(1):pstep(1):pmax(1))';
p2vec = (pmin(2):pstep(2):pmax(2))';

aa1 = exp((q-beta*p1vec)/sigma)*ones(1,length(p2vec));
aa2 = ones(length(p1vec),1)*exp((-beta*p2vec')/sigma);
a1 = aa1 ./ (aa1+aa2);
a2 = aa2 ./ (aa1+aa2);

cond1 = (p1vec < q/alpha);

d1mat = a1.*(cond1*ones(1,length(p2vec))) + a1.*(((1-cond1).*(1+q-alpha*p1vec))*ones(1,length(p2vec)));
d2mat = a2.*(ones(length(p1vec),1)*(1-alpha*p2vec'));

pi1mat = (fees2*p1vec-c1)*ones(1,length(p2vec)).*d1mat;
pi2mat = (ones(length(p1vec),1)*(fees2*p2vec'-c2)).*d2mat;

maxmat1 = ones(length(p1vec),1)*max(pi1mat,[],1);
maxmat2 = max(pi2mat,[],2)*ones(1,length(p2vec));

mateqm = (maxmat1==pi1mat)+(maxmat2==pi2mat);
[eqmrow,eqmcol] = find(mateqm==2,1,'first');

peqm = [p1vec(eqmrow),p2vec(eqmcol)];

cond1out = cond1(eqmrow);