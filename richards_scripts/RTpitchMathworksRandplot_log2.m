Alim=3; % look for Alim consecutive samples of voiced pitch
WinSize_ms=5*ceil(WindowLen*1000/Fs);
StepSize_ms=5*ceil(AnalysisLen*1000/Fs);
verbose=0;
f0_low=95; % take only pitch values above that level
f0_high=140; % below that level
max_pitch_change=5; % disregard pitch changes (from one point to the next) by more than 10.

f0s=cell(length(levs)+1,1+do_var);
h0=zeros(Fs/10,cs);
for i=1:size(f0s,1)
    f0s{i,1}=h0;
    if do_var
        f0s{i,2}=h0;
    end
end
counts=f0s;
f0sneg=f0s; countsneg=counts;
if ~computed_pitch
    f0_times=cell(1,cs); f0_values=f0_times; f0_candidates=f0_times;
end
for i=1:cs
    y0=Ys0{i}'; % what i said
    y=Ys{i}'; % what i heard
    % [voiced, pitch_plot] = f_VOICED(y, Fs, WinSize_ms/1000);
    if computed_pitch
        f0_time=f0_times{i}; f0_value=f0_values{i}; f0_candidate=f0_candidates{i};
    else
        [f0_time,f0_value,SHR,f0_candidate]=shrp(y,Fs,[50 550],WinSize_ms,StepSize_ms,.7,1250,0,0); % compute pitch
        f0_times{i}=f0_time; f0_values{i}=f0_value; f0_candidates{i}=f0_candidate;
    end
    flower_nonzero=find(f0_candidate(:,1)>f0_low);
    f0_value(flower_nonzero)=f0_candidate(flower_nonzero,1);
    
    cand=abs(diff(f0_value))<max_pitch_change & f0_value(2:end)>f0_low & f0_value(2:end)<f0_high; % Rich's crude way of extracting voiced region
    
    if verbose
        figure(1); clf; subplot(211); plot(y0); axis tight; subplot(212); plot(f0_value);
        hold on; %plot(f0_candidates(:,1),'r');
        f0nan=nan(size(f0_value)); f0nan(cand)=f0_value(cand);
        plot(f0nan,'r');
        plot([1 length(f0_value)],f0_high*[1 1],'g'); plot([1 length(f0_value)],f0_low*[1 1],'g');
        axis tight;   pause
    end
    %   [f0_time,f0_value,SHR,f0_candidates]=shrp(y,Fs,[50 550],WinSize_ms,StepSize_ms,.2,1250,0,0); % compute pitch
    
    % f0_value2=func_pitch(y,Fs);
    %    figure(1);clf;plot(y,'r'); hold on;plot(voiced,'b')
    %pause
    %    cand=voiced;
    a=regionprops(cand,'Area','Pixellist');
    
    t_on=0;     % ton=1 means relative to onset of new pitch level, ton=0 means relative to onset of first voiced region
    for j=1:length(a)
        if a(j).Area>Alim
            if t_on==0
                t_on=a(j).PixelList(1,2);
            end
            hh=1+a(j).PixelList(1,2):a(j).PixelList(end,2);
            f0s{Levs(i),Vars(i)+1}(hh-t_on+1,i)=f0s{Levs(i),Vars(i)+1}(hh-t_on+1,i)+f0_value(hh);
            counts{Levs(i),Vars(i)+1}(hh-t_on+1,i)=counts{Levs(i),Vars(i)+1}(hh-t_on+1,i)+1;
            
            l=length(f0_value);
            f0sneg{Levs(i+1),Vars(i)+1}(l-(hh-t_on+1),i)=f0sneg{Levs(i+1),Vars(i)+1}(l-(hh-t_on+1),i)+f0_value(hh);
            countsneg{Levs(i+1),Vars(i)+1}(l-(hh-t_on+1),i)=countsneg{Levs(i+1),Vars(i)+1}(l-(hh-t_on+1),i)+1;
            %           figure(1);clf;plot(f0s(hh));pause
        end
    end
    %             figure(1);clf;subplot(311); plot(y); axis tight; subplot(312);plot(f0_value); hold on;plot(f0_value(cand),'r'); axis tight;
    %              subplot(313); plot(cand);axis tight;
    %              pause
end
f0s(lev0,2)=f0s(lev0,1); counts(lev0,2)=counts(lev0,1);

computed_pitch=1;
%%

cols={'y','r','m','g','c','b','k'};
Stds=cell(size(f0s,1),1+do_var); Means=Stds; Dts=Means;
for j=1:1+do_var
    tstart=1; % start time for outlier computation
    tlim=75; % end time for outlier computation
    nout=0; % number of outliers to remove (defined by std of f0)
    f0s2var=cell(length(levs)+1,1+do_var);
    %  figure(j);clf;
    figure(100+j);clf;
    figure(200+j);clf;
    leg_str=cell(1,length(levs)+1);
    Var_tit={'Var 0 ','Var 1'};
    
    for i=size(f0s,1):-1:1
        hh=std(f0s{i,j}(tstart:tlim,:),[],1);
        [hhsort,hhisort]=sort(hh,'descend');
        
        which_el=hhisort(nout+1:end);
        f0s2{i,j}=sum(f0s{i,j}(:,which_el),2);
        counts2{i,j}=sum(counts{i,j}(:,which_el),2);
        f0s2var{i,j}=sum((f0s{i,j}(:,which_el)).^2,2);
        
        f0sneg2{i,j}=sum(f0sneg{i,j}(:,which_el),2);
        countsneg2{i,j}=sum(countsneg{i,j}(:,which_el),2);
        
        % remove outliers
        
        % Now plot the heard pitch (heard by the experimenter)
        %         figure(j);
        %
        %         subplot(211); hold on
        h1=find(counts2{i,j});
        %         h1neg=find(countsneg2{1});
        %         plot(StepSize_ms*(1-length(h1neg):length(h1))', [f0sneg2{i,j}(h1neg)./(countsneg2{i,j}(h1neg)+eps); f0s2{i,j}(h1)./(counts2{i,j}(h1)+eps)],cols{i});
        %
        %         subplot(212); hold on
        %         plot(StepSize_ms*(1-length(h1neg):length(h1))',[countsneg2{i,j}(h1neg) ; counts2{i,j}(h1)],cols{i});
        %
        
        figure(100+j);
        subplot(211); hold on;
        tmax=100;
        lw=3;
        
        h1b=h1(h1<tmax)';
        t1=StepSize_ms*(h1b);
        n=counts2{i,j}(h1b)+eps;
        m1=f0s2{i,j}(h1b)./n;
        std1=sqrt(f0s2var{i,j}(h1b)./n-m1.^2);
        
        Means{i,j}=m1; Stds{i,j}=std1; Dts{i,j}=t1;
        patch([t1 t1(end:-1:1)],[m1+std1;m1(end:-1:1)-std1(end:-1:1)],cols{i},'FaceAlpha',.5);
        hold on;
        plot(t1,m1,cols{i},'linewidth',lw);
        subplot(212); hold on;
        plot(t1,n,cols{i},'linewidth',lw);       
        leg_str{size(f0s,1)+1-i}=['Shift ' num2str(i-lev0)];
      end
    
    figure(j);
    subplot(211);
    legend(leg_str,'Location','best');
    title(Var_tit{j},'fontsize',14);
    
    figure(100+j);
    subplot(211);
    title(Var_tit{j},'fontsize',14);
    subplot(212);
    legend(leg_str,'Location','best');
    
    
    %
    figure(200+j);
  for i=size(f0s,1):-1:1
      if i==lev0
          continue
      end
        l=min(length(Means{i,j}),length(Means{lev0,j}));
         subplot(311); hold on;
        % plot(Dts{i,j}(1:l),Means{i,j}(1:l)-Means{lev0,j}(1:l),cols{i},'linewidth',lw);
        plot(Dts{i,j}(1:l),1200*log2(Means{i,j}(1:l)./Means{lev0,j}(1:l)),cols{i},'linewidth',lw);
     %    1200*log(SynthesisLens/AnalysisLen);
         subplot(312); hold on;
         plot(Dts{i,j},Stds{i,j},cols{i},'linewidth',lw);    
         subplot(313); hold on;
%         plot(Dts{i,j}(1:l),(Means{i,j}(1:l)-Means{lev0,j}(1:l))/,cols{i},'linewidth',lw);         
        plot(Dts{i,j}(1:l),100*1200*log2(Means{i,j}(1:l)./Means{lev0,j}(1:l))/cents_pitches(i),cols{i},'linewidth',lw);
  end
    legend(leg_str(levs),'Location','best');
      subplot(311); ylabel('Relative pitch (cents)'); xlabel('Time');
      subplot(312); ylabel('Pitch Std (Hz)'); xlabel('Time');
      subplot(313); ylabel('Percent compensation'); xlabel('Time');
end

if do_var
    nn=size(f0s,1);
    for i=1:nn
        if i==lev0
            continue
        end
        figure(1000+i); clf;
        subplot(211);hold on;
        plot(Dts{i,1},Means{i,1},'b','linewidth',lw);
        plot(Dts{i,2},Means{i,2},'r','linewidth',lw);
        
        legend(Var_tit,'Location','Best');
        ylabel('Mean'); title(['Shift ' num2str(i-lev0)],'fontsize',14);
        subplot(212); hold on;
        plot(Dts{i,1},Stds{i,1},'b','linewidth',lw); plot(Dts{i,2},Stds{i,2},'r','linewidth',lw);
        ylabel('Std');
    end
end

