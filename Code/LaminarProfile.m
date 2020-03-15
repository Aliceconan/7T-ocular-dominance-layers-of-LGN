clear;
subj={'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10'};

EPI2d_people=[1 2 3 4 5 8 9];EPI3d_people=[1 2 4 5 6 7 10];SSFP3d_people=[1 2 3 4 5];
sequence_people={EPI2d_people,SSFP3d_people,EPI3d_people};

%%best LREye pattern slices
S01_left=157:162;S01_right=155:160;
S02_left=151:156;S02_right=151:156;
S03_left=152:157;S03_right=138:143;
S04_left=152:157;S04_right=153:158;
S05_left=164:169;S05_right=165:170;
S06_left=159:164;S06_right=157:162;
S07_left=163:168;S07_right=157:162;
S08_left=152:157;S08_right=149:154;
S09_left=152:157;S09_right=156:161;
S10_left=152:157;S10_right=156:161;

left_coord={};right_coord={};angle={};
for people=1:10
    left_coord{people}=eval([subj{people} '_left']);
    right_coord{people}=eval([subj{people} '_right']);
    angle{people}=[45,45,45,45,45,45;-45,-45,-45,-45,-45,-45];
    
end


sequence_name={'EPI2d','SSFP3d','EPI3d'};
conds={'LE-RE','M-P'};
cond=1;
interp_num=16;column=3;
for sequence=1:3
    layer_profile1={}; 
    
    for people=sequence_people{sequence}
        
        LEFT=[];RIGHT=[];
        
        % load the 1D file that contains the beta value of LE, RE and LE-RE
        % of LGN. The first 3 columns are coordinate values, the 4th column
        % is LE beta, the 5th column is RE beta, the 6th column is LE-RE beta          
        slice = load(['./Doc/' subj{people} '/Beta_L_R_L-R_' sequence_name{sequence} '.1D']);   
        slice=slice(:,[1,2,3,column+3]);
        
        leftdata=slice(slice(:,1)>130,:);
        rightdata=slice(slice(:,1)<130,:);
        
       
        for coronal=1:6
            
            % left LGN's data
            line=0;slice_data=[];data=[];
            for slices=left_coord{people}(coronal)
                stats_slice=leftdata((leftdata(:,2)==slices),:);
                data=[data;stats_slice];
            end
            COM = round(mean(data(:,1:3)));y=COM(2);
            data(:,1:3) = -(data(:,1:3) - ones(size(data,1),1)*COM)+10;
            for coord1=min(data(:,1)):max(data(:,1))
                for coord3=min(data(:,3)):max(data(:,3))
                    if any((data(:,1)==coord1).*(data(:,3)==coord3))~=0
                        line=line+1;
                        slice_data(line,1)=coord1;
                        slice_data(line,2)=y;
                        slice_data(line,3)=coord3;
                        slice_data(line,4)=mean(data(logical((data(:,1)==coord1).*(data(:,3)==coord3)),4));   %when x=x1,z=z1, y is slices, so this step is slice_mean;
                    else
                    end
                end
            end
            slice_image = zeros(20,20);
            for vi = 1:size(slice_data,1)
                x = slice_data(vi,1);
                z = slice_data(vi,3);
                slice_image(z,x) = slice_data(vi,4);
            end
            subplot_k = 0;
            for rot_angle = angle{people}(1,coronal)
                rot_slice_image = imrotate(slice_image,rot_angle,'bicubic');rot_slice_image = rot_slice_image(end/2-9:end/2+10,end/2-9:end/2+10);
                subplot_k = subplot_k+1;
                if rot_angle==angle{people}(1,coronal)
                    a=sum(rot_slice_image)./sum(rot_slice_image~=0);a=a(~isnan(a));
                    left=fliplr(a);
                    xleft=1:1:size(left,2);xinterpleft=1:(size(left,2)-1)/(interp_num-1):size(left,2);
                    left=interp1(xleft,left,xinterpleft,'linear');
                    LEFT=[LEFT;left];
                end
            end
            
            
            % right LGN's data
            line=0;slice_data=[];data=[];
            for slices=right_coord{people}(coronal)
                stats_slice=rightdata((rightdata(:,2)==slices),:);
                data=[data;stats_slice];
            end
            COM = round(mean(data(:,1:3)));y=COM(2);
            data(:,1:3) = -(data(:,1:3) - ones(size(data,1),1)*COM)+10;
            for coord1=min(data(:,1)):max(data(:,1))
                for coord3=min(data(:,3)):max(data(:,3))
                    if any((data(:,1)==coord1).*(data(:,3)==coord3))~=0
                        line=line+1;
                        slice_data(line,1)=coord1;
                        slice_data(line,2)=y;
                        slice_data(line,3)=coord3;
                        slice_data(line,4)=mean(data(logical((data(:,1)==coord1).*(data(:,3)==coord3)),4));
                    else
                    end
                end
            end
            slice_image = zeros(20,20);
            for vi = 1:size(slice_data,1)
                x = slice_data(vi,1);
                z = slice_data(vi,3);
                slice_image(z,x) = slice_data(vi,4);
            end
            subplot_k = 0;
            for rot_angle = angle{people}(2,coronal)
                rot_slice_image = imrotate(slice_image,rot_angle,'bicubic');rot_slice_image = rot_slice_image(end/2-9:end/2+10,end/2-9:end/2+10);
                subplot_k = subplot_k+1;
                if rot_angle==angle{people}(2,coronal)
                    a=sum(rot_slice_image)./sum(rot_slice_image~=0);a=a(~isnan(a));
                    right=a;
                    xright=1:1:size(right,2);xinterpright=1:(size(right,2)-1)/(interp_num-1):size(right,2);
                    right=interp1(xright,right,xinterpright,'linear');
                    RIGHT=[RIGHT;right];
                end
            end
        end
        layer_profile1{people,1}=mean(LEFT,1);
        layer_profile1{people,2}=mean(RIGHT,1);
        
    end
    
    
    linewidth=0.1;fontsize=15;
    
    
    %% subjects layer profile
    people_num=size(sequence_people{sequence},2);    
    figure;k=0;
    for i=sequence_people{sequence}
        k=k+1;
        subplot(1,people_num,k)
        x1=-floor(size(layer_profile1{i,1},2)/2);x2=-floor(size(layer_profile1{i,2},2)/2);
        plot([x1:x1+size(layer_profile1{i,1},2)-1],zscore(layer_profile1{i,1}'),'r-',[x2:x2+size(layer_profile1{i,2},2)-1],zscore(layer_profile1{i,2}'),'b-'),xlim([-10 10]);
        pbaspect([1.5,1,1]);
    end
    
    
    %% subjects layer profile mean
    interp_num=16;LEFT=[];RIGHT=[];
    for people = sequence_people{sequence}
        left=layer_profile1{people,1};
        right=layer_profile1{people,2};
        LEFT=[LEFT;left];
        RIGHT=[RIGHT;right];
    end
    left_mean=mean(LEFT);left_SEM=std(LEFT)/sqrt(size(LEFT,1));
    right_mean=mean(RIGHT);right_SEM=std(RIGHT)/sqrt(size(RIGHT,1));
    figure;
    plot(1:16,left_mean,'r',1:16,right_mean,'b');
    hold on;errorbar(1:16,left_mean,left_SEM,'r.');
    hold on;errorbar(1:16,right_mean,right_SEM,'b.');
    title([sequence_name{sequence} ' ' conds{cond} ' layer profile'],'FontSize',fontsize);
    legend('Left LGN','Right LGN');legend('boxoff');set(legend,'FontSize',12);
    xlabel('Medial to Lateral','FontSize',fontsize);ylabel('signal change%','FontSize',fontsize);
    pbaspect([1.5,1,1]);
    
    left_for_clac_corr{sequence}=left_mean;
    right_for_clac_corr{sequence}=right_mean;

end


figure;
for sequence=1:3
    corr_seq=corrcoef(left_for_clac_corr{sequence},right_for_clac_corr{sequence});
    lLGNrLGN_beta{sequence}=[left_for_clac_corr{sequence},right_for_clac_corr{sequence}];
    corr_lLGNrLGN(sequence)=corr_seq(1,2);
    xmax=max(left_for_clac_corr{sequence});xmin=min(left_for_clac_corr{sequence});
    ymax=max(right_for_clac_corr{sequence});ymin=min(right_for_clac_corr{sequence});
    subplot(1,3,sequence)
    plot(left_for_clac_corr{sequence},right_for_clac_corr{sequence},'.','MarkerFaceColor','b');h=lsline;set(h,'LineWidth',1,'Color',[0 0 0])
    xlabel('left LGN','FontSize',15);ylabel('right LGN','FontSize',15);
    title(['r=' num2str(corr_lLGNrLGN(sequence))],'FontSize',15);axis([xmin-0.1 xmax+0.1 ymin-0.1 ymax+0.1]);     
    axis square;
    
end

figure;
seq_group=[1 2;2 3;1 3];
for k=1:3
    corr_seq=corrcoef(lLGNrLGN_beta{seq_group(k,1)},lLGNrLGN_beta{seq_group(k,2)});
    xmax=max(lLGNrLGN_beta{seq_group(k,1)});xmin=min(lLGNrLGN_beta{seq_group(k,1)});
    ymax=max(lLGNrLGN_beta{seq_group(k,2)});ymin=min(lLGNrLGN_beta{seq_group(k,2)});
    subplot(1,3,k)
    plot(lLGNrLGN_beta{seq_group(k,1)},lLGNrLGN_beta{seq_group(k,2)},'.','MarkerFaceColor','b');h=lsline;set(h,'LineWidth',1,'Color',[0 0 0])
    xlabel(sequence_name{seq_group(k,1)},'FontSize',15);ylabel(sequence_name{seq_group(k,2)},'FontSize',15);
    title(['r=' num2str(corr_seq(1,2))],'FontSize',15);axis([xmin-0.1 xmax+0.1 ymin-0.1 ymax+0.1]);     
    axis square;
end

