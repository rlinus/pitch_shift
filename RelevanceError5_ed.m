% check compensation as a function of sensory variance
n=2000; % number of trials
Vsall=2:1:100; % Sensory variance
Vi=10; % Variance internal model
delp=-30; % pitch shfit
p0=2e-7; % the constant k in the document

pzero=200; % pitch target
pzerom=200; % pitch motor output

% some variables used for plotting
pestall=zeros(size(Vsall));
peststd=pestall;
postall=pestall;
pestA=zeros(1,n); postA=pestA; pzeromA=pestA;
pbias=0;pbiasA=pestA; % pitch bias (epsilon in document)
pbiasall=pestall;
for i=1:length(Vsall)
    Vs=Vsall(i);
    aiA=randn(n,1)*sqrt(Vi);
    afA=randn(n,1)*sqrt(Vs)+aiA+delp;
    Vi=20;
    for j=1:n
        ai=aiA(j)+pzerom;
        af=afA(j)+pzerom;
        p=exp(-(ai-af).^2/2/(Vi+Vs))/sqrt(2*pi*(Vs+Vi));
        post=p./(p+p0);
        pest=post.*(ai*Vs+af*Vi)/(Vs+Vi) +(1-post).*ai; %!! pest=ai+post.*(af-ai)*Vi/(Vs+Vi);
        %was1=post>0.5;
        if 1%was1 % one source
            pbias=pbias-.1*(pest-pzero);
             pzerom=pzero+pbias;
           % pzerom=pzerom+(pest-ai)/100; % bring internal mean closer to visual mean
         %   Vi=Vi/1.001;
        else
        %    Vi=Vi*1.001; % two sources; widen internal variance (more uncertainty)
        end
        pbiasA(j)=pbias;
        pzeromA(j)=pzerom;
        pestA(j)=pest+pzero-ai;
        postA(j)=post;
    end
    
    %  pest=post.*(ai*Vs+af*Vi)/(Vs+Vi)+(1-post).*(af.*hh+ai.*(1-hh));
    pestall(i)=mean(pestA(end/2:end));
    peststd(i)=std(pestA(end/2:end));
    postall(i)=mean(postA(end/2:end));
    pbiasall(i)=mean(pbiasA(end/2:end));
end
figure(4);clf;plot(pzeromA); title('Motor');
figure(3); clf; plot(pestA); title('Pest');
figure(1);clf;plot(Vsall,100*abs(pbiasall/delp));%plot(pestall-pzero);
ylabel('Compensation (Percent)');
xlabel('Sensory Variance');
figure(2);clf;plot(peststd)
figure(5);clf;plot(pestall); title('adaptation of ppd')
