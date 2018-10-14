function[formant]=formant_tracker(input,dimension,step,window_size,p,threshold)
% %Formant_tracker based on LPC Analysis
% format:[formants]=formant_tracker(input,dimension,step,window_size,p,threshold)
% where,
% input:input speech waveform;
% dimension:input size, in number of samples.
% step:step size , in number of samples;
% window_size: frame/window size, in number of samples;
% p: number of LP coeffcients;
% threshold: the threshold value for pole magnitude.Poles having higher magnitude are candidates for being formant
% frequencies.
% 
% Author:Narsimh Kamath.
% NITK, India.Thanks!
num=1;
dimension-window_size;
    
for i=1:step:dimension-window_size
    
    frame=input(num*step-window_size/2:num*step+window_size/2);
     
 

for k=1:p+1
    autocorrelation(k)=0;
    for l=1:window_size-k
    autocorrelation(k)=autocorrelation(k)+frame(l)*frame(l+k-1);
    end 
    autocorrelation(k)=autocorrelation(k)/window_size;
end

t=window_size/2;
for k=20:t
    acf(k)=0;
    for l=1:window_size-k
    acf(k)=acf(k)+frame(l)*frame(l+k-1);
    end 
    acf(k)=acf(k)/window_size;
end
[rr,ii]=max(acf(10:t));
tp=10000/ii;

alpha=zeros(p,p);
a=zeros(p,1);
E=zeros(p,1);
E(1)=autocorrelation(1);
for i=1:p
    sum=0;
    for j=1:i-1
        if(i)>1
        sum=sum+alpha(j,(i-1))*autocorrelation(abs(i-j)+1);
        end
    end
    
        k(i)=(autocorrelation(i+1)-sum)/E(i);
        alpha(i,i)=k(i);
        if i>1
            for j=1:i-1
                alpha(j,i)=alpha(j,i-1)-k(i)*alpha(i-j,i-1);
            end
        end
        E(i+1)=(1-k(i)^2)*E(i);
        var=E(i+1);
    end
    for j=1:p
        a(j)=alpha(j,p);
    end
    f=[1 -a'];
  gain=0;
  p;
    for i =1:p+1
        gain=gain+f(i)*autocorrelation(i);
    end
   gain;
    gain=gain^0.5;
   

root1=roots([1,-a']);

    mag_root=abs(root1);
    arg_root=angle(root1);
    k=1;
    for j=1:p
        if mag_root(j)>threshold
            if arg_root(j)>0  &arg_root(j)<pi
            formant(num,k)=arg_root(j)/pi*5000;
            
            num;
            k;
            k=k+1;
            end
        end
    end
   
    
    num=num+1;
end
 s=size(formant);
    for i=1:length(formant)
        for j=1:s(2)
            if formant(i,j)==0
                formant(i,j)=NaN;
            end
        end
    end
   
    plot(formant,'.b');
    return;
end