funcprot(0)
clear
T=5
dt=1/252

//Heston parameters:
theta=0.16^2; //long run vol sigma^2
sigma=0.25;//vol of vol
k=round(sigma^2/theta); //mean reversion speed

r=0; //riskless rate
s0=100;//initial stock price

N=T/dt//Time steps
m=15000;//Number of paths

grand("setsd",1)//set seed for simulation

dB=grand((N-1),m, 'nor',0,1);
dW=grand((N-1),m, 'nor',0,1);
//for correlated Brownian motion
dW1=grand((N-1),m, 'nor',0,1);

function f=performance(heston_par,kappa,trading_boundaries)
ETF=trading_boundaries
v0=theta
v=zeros(N,m);
//Simulation Part: Generate Heston paths (v,S) stoch_vol=1
//Simulation Part: Black-Scholes, if stoch_vol=0

//Creating approximating sample from stationary distribution of the stoch. variance
v(1,:)=v0;
a=k*theta
for i=2:N
vnew=((sigma*sqrt(dt)*dW(i-1,:)+sqrt( sigma^2*dt*(dW(i-1,:)).^2+4*(v(i-1,:)+(a-sigma^2/2)*dt)*(1+k*dt))) /(2*(1+k*dt))).^2;
v(i,:)=vnew;//stochastic variance
end

//Simulating daily Heston data for 5 years (m paths)
v=zeros(N,m);
v(1,:)=vnew;

for i=2:N

vnew=((sigma*sqrt(dt)*dW1(i-1,:)+sqrt( sigma^2*dt*(dW1(i-1,:)).^2+4*(v(i-1,:)+(a-sigma^2/2)*dt)*(1+k*dt))) /(2*(1+k*dt))).^2;
v(i,:)=vnew;//stochastic variance

end

v=v*(heston_par)+(v*0+theta)*(1-heston_par);

rho=-0.9;
//Heston Stock Price
s1=exp(cumsum((kappa-1/2)*v(1:(N-1),:)*dt+sqrt(dt)*sqrt(v(1:(N-1),:)).*(rho*dW1+sqrt(1-rho^2)*dB),'r'));
s=s0*[ones(1,m);s1];
mut=kappa*v(1:N,:);

   

epsi=0.001

   
M=length(ETF(:,1));

//Trading Part
   
for l=1:M
l 
etf=ETF(l,:);
pim=etf(2);
pip=etf(3);

pi0=(pim+pip)/2;

w0=100;

phi0=w0*pi0/s0;
cash0=w0*(1-pi0);


w=zeros(N,m);

phi=zeros(N,m);
prw=zeros(N,m);
dphi=zeros(N,m);

dcash=zeros(N,m);
cash=zeros(N,m);


cash(1,:)=cash0;
phi(1,:)=phi0;
w(1,:)=w0;

//wealth right before next trade
for i=2:N
    disp([l,i])

wnew=phi(i-1,:) .*s(i,:)+cash(i-1,:);//marking to market
//prop wealth before trade
pinew=phi(i-1,:) .*s(i,:) ./wnew;

pitarget=min(max(pinew,pim),pip);

sale=pinew>pip;

buy=pinew<pim;

//buy
dphi1=pim .*wnew ./s(i,:)-phi(i-1,:);
//sell    
dphi2=(pip .* wnew-phi(i-1,:) .*s(i,:)) ./ (s(i,:) .*(1-pip*0.001));

dphinew=dphi1 .*buy+dphi2 .*sale;



phi(i,:)=phi(i-1,:)+dphinew;
dphi(i,:)=dphinew;

w_after=wnew+epsi*dphinew .* sale .*s(i,:);
w(i,:)=w_after;

cash_after=cash(i-1,:)-dphinew .*s(i,:)+epsi .*dphinew .*sale .*s(i,:);
cash(i,:)=cash_after;

retS(i,:)=(s(i,:)-s(i-1,:)) ./s(i-1,:);
retW(i,:)=(w(i,:)-w(i-1,:)) ./w(i-1,:);

end
prw=phi .* s ./ w;

tc=epsi.* dphi .*(dphi<0)  ./phi .*prw; 

ATC=sum(tc,'r')/T;
mATC=mean(ATC);

mATC=mean(retW-3*retS-mut .*(prw-3)*dt)/dt
ETF(l,7)=mATC;

err=abs(mATC-ETF(l,5))/ETF(l,5);
TRE=sqrt(mean((retW-3*retS-mut .*(prw-3)*dt).^2)/dt);
ETF(l,8)=TRE
   

end
f=ETF
endfunction


ETF0=fscanfMat('Bull30bp10kappa0');
//Since trading for 5 years along 15000 paths takes long, we only compute a few trading boundaries.
//For the plots below, we will interpolate the results, to avoid unnatural kinks.
//Warning: If one changes the range of gammas in the files, then one needs also to change the below selection
//so to arrive at the same tracking error range of Figure 6.

ETF0=ETF0([980, 990,998,1000,1005,1010,1025,1050,1100,1200,1300,1400,1600,1700],:)
pim=ETF0(:,2);
pip=ETF0(:,3);
ExR=-theta/2*pim .*pip .*(pip-1)^2 ./(pip-pim) ./(1/0.001-pip);
beta=pip .*pim  .*log(pip ./pim) ./(pip-pim);
TrE=sqrt(theta)*sqrt(pip .*pim+3*(3-2*beta));
ETF0(:,5)=ExR;
ETF0(:,6)=TrE;


ETF0BS=performance(0,0,ETF0)
ETF0BS(:,9)=ETF0BS(:,5)-ETF0BS(:,1)/2 .*ETF0BS(:,6).^2
ETF0BS(:,10)=ETF0BS(:,7)-ETF0BS(:,1)/2 .*ETF0BS(:,8).^2

ETF0Heston=performance(1,0,ETF0)
ETF0Heston(:,9)=ETF0Heston(:,5)-ETF0Heston(:,1)/2 .*ETF0Heston(:,6).^2
ETF0Heston(:,10)=ETF0Heston(:,7)-ETF0Heston(:,1)/2 .*ETF0Heston(:,8).^2


ETF4=fscanfMat('Bull30bp10kappa1.5625');
ETF4=ETF4([980, 990,998,1000,1005,1010,1025,1050,1100,1200,1300,1400,1600,1700],:)

ETF4BS=performance(0,0.04/0.16^2,ETF4)
ETF4BS(:,9)=ETF4BS(:,5)-ETF4BS(:,1)/2 .*ETF4BS(:,6).^2
ETF4BS(:,10)=ETF4BS(:,7)-ETF4BS(:,1)/2 .*ETF4BS(:,8).^2

ETF4Heston=performance(1,0.04/0.16^2,ETF4)
ETF4Heston(:,9)=ETF4Heston(:,5)-ETF4Heston(:,1)/2 .*ETF4Heston(:,6).^2
ETF4Heston(:,10)=ETF4Heston(:,7)-ETF4Heston(:,1)/2 .*ETF4Heston(:,8).^2


   
ETF8=fscanfMat('Bull30bp10kappa3.125');
ETF8=ETF8([980, 990,998,1000,1005,1010,1025,1050,1100,1200,1300,1400,1600,1700],:)

ETF8BS=performance(0,0.08/0.16^2,ETF8)
ETF8BS(:,9)=ETF8BS(:,5)-ETF8BS(:,1)/2 .*ETF8BS(:,6).^2
ETF8BS(:,10)=ETF8BS(:,7)-ETF8BS(:,1)/2 .*ETF8BS(:,8).^2

ETF8Heston=performance(1,0.08/0.16^2,ETF8)
ETF8Heston(:,9)=ETF8Heston(:,5)-ETF8Heston(:,1)/2 .*ETF8Heston(:,6).^2
ETF8Heston(:,10)=ETF8Heston(:,7)-ETF8Heston(:,1)/2 .*ETF8Heston(:,8).^2
  

fprintfMat('ETF0BS',ETF0BS, "%5.10f")

fprintfMat('ETF4BS',ETF4BS, "%5.10f")

fprintfMat('ETF8BS',ETF8BS, "%5.10f")

fprintfMat('ETF0Heston',ETF0Heston, "%5.10f")

fprintfMat('ETF4Heston',ETF4Heston, "%5.10f")

fprintfMat('ETF8Heston',ETF8Heston, "%5.10f")

////////////////////////////////////////////

ETF0BS=fscanfMat('ETF0BS');

ETF4BS=fscanfMat('ETF4BS');

ETF8BS=fscanfMat('ETF8BS');

ETF0Heston=fscanfMat('ETF0Heston');

ETF4Heston=fscanfMat('ETF4Heston');

ETF8Heston=fscanfMat('ETF8Heston');





x=ETF0BS(:,6);
y=-ETF0BS(:,5);

x=flipdim(x,1)
y=flipdim(y,1)

xx=linspace(min(x),max(x),1000)

d=splin(x,y,"natural")
yy=interp(xx,x,y,d)
plot(xx,yy,'black')



x=ETF0BS(:,8);
y=-ETF0BS(:,7);

x=flipdim(x,1)
y=flipdim(y,1)

xx=linspace(min(x),max(x),1000)

d=splin(x,y,"natural")
yy=interp(xx,x,y,d)
plot(xx,yy,'b')




x=ETF0Heston(:,8);
y=-ETF0Heston(:,7);

x=flipdim(x,1)
y=flipdim(y,1)

xx=linspace(min(x),max(x),1000)
   
d=splin(x,y,"natural")
yy=interp(xx,x,y,d)
plot(xx,yy,'g')
   


 x=ETF4BS(:,6);
y=-ETF4BS(:,5);

x=flipdim(x,1)
y=flipdim(y,1)

xx=linspace(min(x),max(x),1000)

d=splin(x,y,"natural")
yy=interp(xx,x,y,d)
plot(xx,yy,'black--')



x=ETF4BS(:,8);
y=-ETF4BS(:,7);

x=flipdim(x,1)
y=flipdim(y,1)

xx=linspace(min(x),max(x),1000)

d=splin(x,y,"natural")
yy=interp(xx,x,y,d)
plot(xx,yy,'b--')




x=ETF4Heston(:,8);
y=-ETF4Heston(:,7);

x=flipdim(x,1)
y=flipdim(y,1)

xx=linspace(min(x),max(x),1000)
   
d=splin(x,y,"natural")
yy=interp(xx,x,y,d)
plot(xx,yy,'g--')





x=ETF8BS(:,6);
y=-ETF8BS(:,5);

x=flipdim(x,1)
y=flipdim(y,1)

xx=linspace(min(x),max(x),1000)

d=splin(x,y,"natural")
yy=interp(xx,x,y,d)
plot(xx,yy,'black-.')



   
x=ETF8BS(:,8);
y=-ETF8BS(:,7);

x=flipdim(x,1)
y=flipdim(y,1)

xx=linspace(min(x),max(x),1000)

d=splin(x,y,"natural")
yy=interp(xx,x,y,d)
plot(xx,yy,'b-.')




x=ETF8Heston(:,8);
y=-ETF8Heston(:,7);

x=flipdim(x,1)
y=flipdim(y,1)

xx=linspace(min(x),max(x),1000)
   
d=splin(x,y,"natural")
yy=interp(xx,x,y,d)
plot(xx,yy,'g-.')








mtlb_axis([0.01,0.03,0.05/100, 0.2/100])

a=gca();
a.log_flags = "ln" ; 


a=get("current_axes")
old_margins=a.margins
a.margins=old_margins/ins/2


new_locations=a.y_ticks.locations
new_labels=string(new_locations*100)+"%"
a.y_ticks = tlist(["ticks", "locations", "labels"], new_locations, new_labels);

new_locations=a.x_ticks.locations
new_labels=string(new_locations*100)+"%"
a.x_ticks = tlist(["ticks", "locations", "labels"], new_locations, new_labels);

xs2eps(0,'FigQF6.eps')
























