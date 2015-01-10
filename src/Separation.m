clc;clear all;close all; 
I1=imread ('input1.jpg');          
I2=imread ('input2.jpg'); 
I3=imread ('input3.jpg'); 
figure()
subplot(3,3,1),imshow(I1),title('Input1'),
subplot(3,3,2),imshow(I2),title('Input2'), 
subplot(3,3,3),imshow(I3),title('Input3'), 
[x1,y1,z1]=size(I1); 
[x2,y2,z2]=size(I2); 
[x3,y3,z3]=size(I3); 
s1=[reshape(I1,[1,x1*y1*z1])]; 
s2=[reshape(I2,[1,x2*y2*z2])]; 
s3=[reshape(I3,[1,x3*y3*z3])]; 
S_all=[s1;s2;s3];                           
S=double(S_all);                          
Sweight=rand(size(S_all,1));              
MixedS=Sweight*S;                           
ms1=reshape(MixedS(1,:),[x1,y1,z1]); 
ms2=reshape(MixedS(2,:),[x2,y2,z2]); 
ms3=reshape(MixedS(3,:),[x3,y3,z3]); 
I1_mixed=uint8(round(ms1)); 
I2_mixed=uint8(round(ms2)); 
I3_mixed=uint8(round(ms3)); 
subplot(3,3,4),imshow(I1_mixed),title('Mixed1'), 
subplot(3,3,5),imshow(I2_mixed),title('Mixed2'), 
subplot(3,3,6),imshow(I3_mixed),title('Mixed3'), 
MixedS_bak=MixedS;                         
MixedS_mean=zeros(3,1); 
for i=1:3 
    MixedS_mean(i)=mean(MixedS(i,:)); 
end                                        
 
for i=1:3 
    for j=1:size(MixedS,2) 
        MixedS(i,j)=MixedS(i,j)-MixedS_mean(i); 
    end 
end 
MixedS_cov=cov(MixedS');                 
[E,D]=eig(MixedS_cov);                   
Q=inv(sqrt(D))*(E)';                      
MixedS_white=Q*MixedS;                    
IsI=cov(MixedS_white');                        
X=MixedS_white;                            
[VariableNum,SampleNum]=size(X); 
numofIC=VariableNum;                       
B=zeros(numofIC,VariableNum);               
for r=1:numofIC                           
    i=1;maxIterationsNum=150;             
    b=2*(rand(numofIC,1)-.5);                   
    b=b/norm(b);                         
    while i<=maxIterationsNum+1 
       if i == maxIterationsNum           
            fprintf('No convergence¡£', r,maxIterationsNum); 
            break; 
        end 
        bOld=b;                           
        u=1; 
        t=X'*b; 
        g=t.^3; 
        dg=3*t.^2; 
        b=((1-u)*t'*g*b+u*X*g)/SampleNum-mean(dg)*b;     
                                           
        b=b-B*B'*b;                       
        b=b/norm(b);  
        if abs(abs(b'*bOld)-1)<1e-10      
             B(:,r)=b;                  
             break; 
         end 
        i=i+1;         
    end 
end 
ICAedS=B'*Q*MixedS_bak;                    
ICAedS_bak=ICAedS; 
ICAedS=abs(55*ICAedS); 
is1=reshape(ICAedS(1,:),[x1,y1,z1]); 
is2=reshape(ICAedS(2,:),[x2,y2,z2]); 
is3=reshape(ICAedS(3,:),[x3,y3,z3]); 
I1_icaed =uint8 (round(is1)); 
I2_icaed =uint8 (round(is2)); 
I3_icaed =uint8 (round(is3)); 
subplot(3,3,7),imshow(I1_icaed),title('Restored1'); 
subplot(3,3,8),imshow(I2_icaed),title('Restored2');
subplot(3,3,9),imshow(I3_icaed),title('Restored3');
