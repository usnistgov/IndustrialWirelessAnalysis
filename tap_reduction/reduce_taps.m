function [r_t,r_h,r_ph] =reduce_taps(h,r_paths)
% By Mohamed Hany

if size(h,1)<=r_paths
    r_h=h;
    r_t=repmat(1:size(h,1),1,size(h,2));
    return;
end

num_channels=1;
channel_energy = sum(abs(h).^2);
h_len = length(h(:,1));
t = [0:(h_len-1)] * 4; % for use in computing excess & RMS delays
excess_delay = zeros(1,num_channels);
RMS_delay = zeros(1,num_channels);
num_sig_paths = zeros(1,num_channels);
num_sig_e_paths = zeros(1,num_channels);
for k=1:num_channels
% determine excess delay and RMS delay
sq_h = abs(h(:,k)).^2 / channel_energy(k);
t_norm = t - 0; % remove the randomized arrival time of first cluster
excess_delay(k) = t_norm * sq_h;
RMS_delay(k) = sqrt( ((t_norm-excess_delay(k)).^2) * sq_h );
end

for count=1:size(h,2)
    count
loop_h=conj(h(:,count)');
loop_t=0:size(h,1)-1;
while  length(loop_h)>r_paths
    [min_gain,min_pos]=min(abs(loop_h).^2);
    
  
    A=abs(loop_t-loop_t(min_pos));
    [nb_t_abs,nb_pos]=min(A(A>0));
    if nb_pos==min_pos
        nb_pos=nb_pos+1;
    end
    loop_h(nb_pos)=loop_h(nb_pos)+loop_h(min_pos);
    loop_h(min_pos)=1000;
    l1=length(loop_h);
    loop_h=loop_h(abs(loop_h)<1000);
    l2=length(loop_h);
    loop_t(min_pos)=-1;
    loop_t=loop_t(loop_t>-1);
    length(loop_h);
    length(loop_t);
end
fun=@(x)sum((x-abs(loop_h).^2).^2);
% x0=abs(loop_h).^2;
x0=0.3*ones(1,r_paths);
Aeq(1,:)=ones(1,r_paths);
Aeq(2,:)=loop_t;
Aeq(3,:)=loop_t.^2;
Beq(1,:)=channel_energy(count);
Beq(2,:)=excess_delay(count)*channel_energy(count);
Beq(3,:)=channel_energy(count)*(excess_delay(count)^2+RMS_delay(count)^2);
lb=zeros(1,r_paths);
options=optimoptions('fmincon','TolX',1e-10,'TolCon',1e-3,'MaxFunEvals',25000,'MaxIter',5000);
% x = fmincon(fun,x0,[],[],Aeq,Beq,lb,[],[],options);
x=abs(loop_h).^2; % alternative to optimization problem to keep loop
Da=exp(1j*2*pi*[1:size(h,1)]'*[1:size(h,1)]/size(h,1));

Dh=exp(1j*2*pi*[1:size(h,1)]'*[1:r_paths]/size(h,1));

r_ph(:,count)=angle(Dh'*Da*h(:,count));
r_h(:,count)=sqrt(x);
r_t(:,count)=loop_t;
end

return;