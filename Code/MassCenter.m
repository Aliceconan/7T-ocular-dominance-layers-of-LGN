clear;
subj={'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10'};

%Mass centers distribution of Inner-Middle-Outer clusters, ROIs from 2dEPI.LREye

COM_Inner=[];COM_Middle=[];COM_Outer=[];dotAll=[];
sort_number=14;  %sort_number means the voxels in each cluster
for people=[1 2 3 4 5]
    LGN_All=load(['./Doc/' subj{people} '/LGN1Inner2Middle3Outer4.1D']);  %Middle refers to Ipsilateral cluster, Outer refers to Contralateral cluster
    a=load(['./Doc/' subj{people} '/Tvalue_LE-RE_2dEPI_3dEPI_3dSSFP.1D']);
    LGN_All(:,5)=a(:,4);  %the 4th column of a is the Tvalue of LE-RE of 2dEPI data

    left=LGN_All(LGN_All(:,1)>130,:);
    COM_left=ones(size(left,1),1)*mean(left(:,1:3));left(:,1:3)=left(:,1:3)-COM_left;
    left_Inner=left(left(:,4)==2,1:3);
    left_Middle=left(left(:,4)==3,[1 2 3 5]);LM_sort=sort(left_Middle(:,4),'descend');thred=LM_sort(sort_number);left_Middle=left_Middle(left_Middle(:,4)>=thred,1:3);
    left_Outer=left(left(:,4)==4,[1 2 3 5]);LO_sort=sort(left_Outer(:,4),'ascend');thred=LO_sort(sort_number);left_Outer=left_Outer(left_Outer(:,4)<=thred,1:3);
    COM_Inner=[COM_Inner;mean(left_Inner)];COM_Middle=[COM_Middle;mean(left_Middle)];COM_Outer=[COM_Outer;mean(left_Outer)];
    
    right=LGN_All(LGN_All(:,1)<130,:);
    COM_right=ones(size(right,1),1)*mean(right(:,1:3));right(:,1:3)=right(:,1:3)-COM_right;right(:,1)=-right(:,1);
    right_Inner=right(right(:,4)==2,1:3);
    right_Middle=right(right(:,4)==3,[1 2 3 5]);RM_sort=sort(right_Middle(:,4),'ascend');thred=RM_sort(sort_number);right_Middle=right_Middle(right_Middle(:,4)<=thred,1:3);
    right_Outer=right(right(:,4)==4,[1 2 3 5]);RO_sort=sort(right_Outer(:,4),'descend');thred=RO_sort(sort_number);right_Outer=right_Outer(right_Outer(:,4)>=thred,1:3);
    COM_Inner=[COM_Inner;mean(right_Inner)];COM_Middle=[COM_Middle;mean(right_Middle)];COM_Outer=[COM_Outer;mean(right_Outer)];
    
    
end

%coronal:xz;  sagittal:yz;  axial:xy : first letter is x axis, second letter is y axis
slices={'coronal','sagittal','axial'};
x_orient={'Medial-to-Lateral','Anterior-to-Posterior','Medial-to-Lateral'};
y_orient={'Inferior-to-Superior','Inferior-to-Superior','Posterior-to-Anterior'};


COM_Middle=0.6*COM_Middle;COM_Outer=0.6*COM_Outer;
Middle_mean=mean(COM_Middle);Middle_SEM=std(COM_Middle)/sqrt(size(COM_Middle,1));
Outer_mean=mean(COM_Outer);Outer_SEM=std(COM_Outer)/sqrt(size(COM_Outer,1));
errorbarwidth=2;dotSize=100;fontsize=20;dotfontsize=10;
coord=[1 3;2 3;1 2];
%simulation mass center:
Ipsi_simu=[-0.1266   0.2602     -0.2152];Contra_simu=[1.2926    -0.1109     0.5483];
% Ipsi:        -0.1266   0.2602     -0.2152    
% Contra:       1.2926   -0.1109     0.5483  
% axis: Medial(-) to lateral(+)   Anterior (-) to Posterior(+)    Ventral(-) to Dorsal(+)

figure;
for slice=1:3   %1 is coronal; 2 is sagittal; 3 is axial
    subplot(1,3,slice);
    
    for LGN=1:10
        hold on;text(COM_Middle(LGN,coord(slice,1)),COM_Middle(LGN,coord(slice,2)),num2str(LGN),'FontSize',dotfontsize,'Color','b');
        hold on;text(COM_Outer(LGN,coord(slice,1)),COM_Outer(LGN,coord(slice,2)),num2str(LGN),'FontSize',dotfontsize,'Color','r');
    end
    
    
    
    hold on;dotAll(1)=scatter(Middle_mean(:,coord(slice,1)),Middle_mean(:,coord(slice,2)),dotSize,'b','filled');
    hold on;plot([Middle_mean(:,coord(slice,1))-Middle_SEM(coord(slice,1)) Middle_mean(:,coord(slice,1))+Middle_SEM(coord(slice,1))], ...
        [Middle_mean(coord(slice,2)) Middle_mean(coord(slice,2))],'b',[Middle_mean(:,coord(slice,1)) Middle_mean(:,coord(slice,1))], ...
        [Middle_mean(:,coord(slice,2))-Middle_SEM(coord(slice,2)) Middle_mean(:,coord(slice,2))+Middle_SEM(coord(slice,2))],'b','LineWidth',errorbarwidth);
    
    hold on;dotAll(2)=scatter(Outer_mean(:,coord(slice,1)),Outer_mean(:,coord(slice,2)),dotSize,'r','filled');
    hold on;plot([Outer_mean(:,coord(slice,1))-Outer_SEM(coord(slice,1)) Outer_mean(:,coord(slice,1))+Outer_SEM(coord(slice,1))], ...
        [Outer_mean(coord(slice,2)) Outer_mean(coord(slice,2))],'r',[Outer_mean(:,coord(slice,1)) Outer_mean(:,coord(slice,1))], ...
        [Outer_mean(:,coord(slice,2))-Outer_SEM(coord(slice,2)) Outer_mean(:,coord(slice,2))+Outer_SEM(coord(slice,2))],'r','LineWidth',errorbarwidth);
    
    hold on;scatter(Ipsi_simu(:,coord(slice,1)),Ipsi_simu(:,coord(slice,2)),dotSize,'b');
    hold on;scatter(Contra_simu(:,coord(slice,1)),Contra_simu(:,coord(slice,2)),dotSize,'r');
    
    legend(dotAll,'Middle','Lateral');
    title([slices(slice)],'FontSize',fontsize);
    xlabel([x_orient{slice} ' '],'FontSize',fontsize);ylabel([y_orient{slice} ' '],'FontSize',fontsize);
    if slice==3
        set(gca,'ydir','reverse');
    end
    xlim([min(min([COM_Middle(:,coord(slice,1));COM_Outer(:,coord(slice,1));Ipsi_simu(:,coord(slice,1));Contra_simu(:,coord(slice,1))]))-0.5, ...
        max(max([COM_Middle(:,coord(slice,1));COM_Outer(:,coord(slice,1));Ipsi_simu(:,coord(slice,1));Contra_simu(:,coord(slice,1))]))+0.5]);
    ylim([min(min([COM_Middle(:,coord(slice,2));COM_Outer(:,coord(slice,2));Ipsi_simu(:,coord(slice,2));Contra_simu(:,coord(slice,2))]))-0.5, ...
        max(max([COM_Middle(:,coord(slice,2));COM_Outer(:,coord(slice,2));Ipsi_simu(:,coord(slice,2));Contra_simu(:,coord(slice,2))]))+0.5]);
end
