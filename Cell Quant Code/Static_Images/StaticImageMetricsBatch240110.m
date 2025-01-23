%Code for obtaining nuclei and myotube stats from stained images
   
clearvars -except MCmean0

[Metrics, Auto]=Metrics_choose;
Nuc_count=Metrics(1);
Wid=Metrics(2);
Uni=Metrics(3);
Ach=Metrics(4);
Stri_im=Metrics(5);
FI=Metrics(6);
%input image variables
prompt = {'Pixel width (microns):','Filename prefix:','Start Position:','Number of Positions:','Enhancement factor:','image length:','Dilation factor:','% Myo-marker threshold:','Mean nuclei pixel size:','Min myotube length:','Myotube image size:'};
dlgtitle = 'Input';
dims = [1 35];
%Pericentrin threshold will change with fluorescence intensity
definput = {'0.454','Day3','1','1','20','1000','3','0.15','300','250','250'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

%pixel size
ps=str2double(answer(1));
Expt=cell2mat(answer(2));
pos=str2double(answer(3));
n_pos=str2double(answer(4));
EF=str2double(answer(5));
L2=str2double(answer(6));
DilF=str2double(answer(7));
Pericentrin_Thresh=str2double(answer(8));
MCmean0=str2double(answer(9));
Min_MT_Length=str2double(answer(10));
L1=str2double(answer(11));
L2=round(L2/ps,0);
L_half=round(L2/2,0);
count_no=0;
AA=0;

for p=pos:pos+(n_pos-1)
count_no=count_no+1;
p
%if p<10
%    PlateLetter=sprintf('%s0',row);
%else
%    PlateLetter=sprintf('%s',row);
%end
%Read images, binarize images and crop to length
filename1=sprintf('%s-Pos%d-MTDNA.tif',Expt,p);
filename2=sprintf('%s-Pos%d-DNABW.tif',Expt,p);
filename3=sprintf('%s-Pos%d-DNA.tif',Expt,p);
filename4=sprintf('%s-Pos%d-MTDNABW.tif',Expt,p);
filename5=sprintf('%s-Pos%d-ActinBW.tif',Expt,p);
filename6=sprintf('%s-Pos%d-Actin.tif',Expt,p);
if Ach==1
    filename7=sprintf('%s-Pos%d-BTX.tif',Expt,p);
end

Myotube=imread(filename5);
%Determine centrepoint of image;
Image_dimensions=size(Myotube);
X_C=round(Image_dimensions(1)/2,0);
Y_C=round(Image_dimensions(2)/2,0);

if Stri_im==1
    Actin=imread(filename6);
    Actin=Actin(X_C-L_half:X_C+L_half,Y_C-L_half:Y_C+L_half);
    Actin=mat2gray(2*Actin);
end
%Define myotube regions from actin staining and remove small regions

Myotube=imbinarize(Myotube);
Myotube=Myotube(X_C-L_half:X_C+L_half,Y_C-L_half:Y_C+L_half);
Myotube=bwareaopen(Myotube,2000); %Check
%Read Pericentrin stained regions for observation
Pericentrin=imread(filename1)*5*EF;
Pericentrin=Pericentrin(X_C-L_half:X_C+L_half,Y_C-L_half:Y_C+L_half);
%Read DAPI stained nuclei regions
NucBW=imread(filename2);
NucBW=NucBW(X_C-L_half:X_C+L_half,Y_C-L_half:Y_C+L_half);
NucBW=imbinarize(NucBW);
seMB = strel('disk',DilF);
%Dilate Nuclei
NucDil=imdilate(NucBW,seMB);
%Read unbinarized DAPI nuclei staining for observation
Nuclei=imread(filename3)*EF;
Nuclei=Nuclei(X_C-L_half:X_C+L_half,Y_C-L_half:Y_C+L_half);
%Read pericentrin stained regions
Pericentrin_Nuclei=imread(filename4);
Pericentrin_Nuclei=Pericentrin_Nuclei(X_C-L_half:X_C+L_half,Y_C-L_half:Y_C+L_half);
Pericentrin_Nuclei=imbinarize(Pericentrin_Nuclei);

Pericentrin_Nuclei=imdilate(Pericentrin_Nuclei,seMB);
%Calculate overlapping regions of Pericentrin and DAPI (pericentrin inside
%nuclei)
Pericentrin_Nuclei=Pericentrin_Nuclei & NucDil;

%figure;imshow(MBNuc);
%Centroid/area stats for DAPI stained nuclei regions
statsNuc = regionprops('table',NucBW,'Centroid','area');
   
%Define average nuclei size to use as a threshold for excluding smaller
%regions.
emc = exist('MCmean0');
if emc == 1
   MCmean=MCmean0;
else
   MCmean=mean(statsNuc.Area);
end
%remove areas which are smaller than 1/2 of the mean nucleus area.
Nuc=bwareaopen(NucBW,ceil(1*MCmean/2));
% Define small nuclei as 'Dead cells'
Dead_cells=bwareaopen(NucBW,ceil(1*MCmean/30));
Dead_cells=Dead_cells & ~Nuc;
%figure; imshow(imoverlay(Nuclei,Dead_cells))
statsDeadCell = regionprops('table',Dead_cells,'Centroid','area');
Dead_cell_count(count_no)=length(statsDeadCell.Area);
clear statsMC;
clear statsDeadCell;
%statsNuc=regionprops('table',Nuc,'Centroid','area');
%Choose whether to dilate the nuclei image again? 
%NucDil=imdilate(Nuc,seMB);
%NucDil=Nuc;
%Search through nuclei regions and count the number of pericentrin positive
%pixels
StatsNucPx=regionprops('table',Nuc,'PixelList','Centroid','area');
pxcount=0;
%test=0;
Mx_size=size(Myotube);
NucMT(1:Mx_size(1),1:Mx_size(2))=0;
NucMB(1:Mx_size(1),1:Mx_size(2))=0;
MT_Mature=0;
for px1=1:length(StatsNucPx.PixelList)
    Pixels1=cell2mat(StatsNucPx.PixelList(px1));
    for px2=1:length(Pixels1)
        if logical(Pericentrin_Nuclei(Pixels1(px2,2),Pixels1(px2,1)))
            %test(Pixels1(px2,2),Pixels1(px2,1))=1;
            pxcount=pxcount+1;
        end
    end
    %Calculate proportion of cell with pericentrin staining
    Proportion_pericentrin(px1)=pxcount/length(Pixels1);
    pxcount=0;
        %Count cells with > 50% pericentrin staining
        if Proportion_pericentrin(px1)> 0.5
          MT_Mature=MT_Mature+1;
        end
        if Proportion_pericentrin(px1)< Pericentrin_Thresh
            for px2=1:length(Pixels1)
                NucMB(Pixels1(px2,2),Pixels1(px2,1))=1;
            end
        else 
            for px2=1:length(Pixels1)
                NucMT(Pixels1(px2,2),Pixels1(px2,1))=1;
            end
        end
end
NucMB_BW=imbinarize(NucMB);
[MBSingNuc,MBMultiNuc,MBMultiEst]=nuc_count(NucMB_BW,MCmean);
NucMT_BW=imbinarize(NucMT);
[MTSingNuc,MTMultiNuc,MTMultiEst]=nuc_count(NucMT_BW,MCmean);
AddMBx=0;
AddMBy=0;
AddMB=0;
RemMBx=0;
RemMBy=0;
RemMB=0;
AddMTx=0;
AddMTy=0;
AddMT=0;
RemMTx=0;
RemMTy=0;
RemMT=0;
%MCSub=0;
%MCTot=0;

if Auto==0
    %Show overlay of DAPI and Pericentrin
    CC=imfuse(NucMB, NucMT,'falsecolor','Scaling','joint','ColorChannels',[1 0 2]);  
    figure(1);imshow(imoverlay(CC,Pericentrin_Nuclei,'yellow'))

%Function to manually add/remove myoblasts missed by segmentation 

    title('Left click->swap myoblast (red) for myonucleus (blue), Right click -> swap myonucleus for myoblast, spacebar-> finish')

    [RemMBx,RemMBy,RemMTx,RemMTy]=ManLab(Pericentrin);
    
if sum(RemMTx)>0
    AddMB=length(RemMTx);
    RemMT=AddMB;
else
    AddMB=0;
    RemMT=0;
end

if sum(RemMBx)>0
    RemMB=length(RemMBx);
    AddMT=RemMB;
else
    RemMB=0;
    AddMT=0;
end

end

if Wid==1
%Calculate Myotube density as a fraction of the sample area

MTsize=size(Myotube);
s1=MTsize(1);
s2=MTsize(2);

[MTDens,MTD]=MTDensity(Myotube,s1,s2);
Proportional_Area(count_no)=MTD;

%remove holes from myotube image
Myotube=imdilate(Myotube,seMB);
Filled=imfill(Myotube,'holes');
holes = Filled & ~Myotube;
bigholes = bwareaopen(holes, 200);
smallholes = holes & ~bigholes;
I = Myotube | smallholes;
%Calculate a reference Myotube width by dividing density by total length.
%I=Myotube(:,:);
%figure;imshow(I);
BWskel = bwskel(I,'MinBranchLength',Min_MT_Length);
%figure;imshow(BWskel);
TotLength=sum(sum(BWskel));
ReferenceWidth=MTD/TotLength;
ReferenceWidthPx=ReferenceWidth*ps;
ReferenceWidthDil=round(ReferenceWidth/2,0);
Reference_Width(count_no)=ReferenceWidthPx;         
end

%Calculate mean intensity of BTX per myonucleus
%Find myonuclei by subtracting nuclei from myotubes
Myonuclei=I & NucBW;
Myonuclei=bwareaopen(Myonuclei,round(0.25*MCmean));
%figure;imshow(Myonuclei);
%Remove groups of myonuclei
LgMyonuclei=bwareaopen(Myonuclei,2*MCmean);
Myonuclei= Myonuclei & ~LgMyonuclei;
%figure;imshow(Myonuclei);
if Ach==1
    BTX=imread(filename7);
    BTX=BTX(X_C-L_half:X_C+L_half,Y_C-L_half:Y_C+L_half);
%convert to 8 bit
%BTX=uint8(BTX);
    StatsBTXNucPx=regionprops(Myonuclei,'PixelIdxList','Centroid');
    for k=1:numel(StatsBTXNucPx)
    idxBTX=StatsBTXNucPx(k).PixelIdxList;
    MeanInt(k)=mean(BTX(idxBTX));
    end
    matrixname=sprintf('Mean_BTX_Intensity-Pos%d',p);
    writematrix(MeanInt,matrixname);
end
%if Uni==1
%Calculate distance between myonuclei
%1. Split skeleton at branch points
BWBranch = bwmorph(BWskel,'branchpoints');
BWBranch = bwmorph(BWBranch,'thicken',1);
BWsplit = BWskel&~BWBranch;        
%2. remove smaller myotube sections
BWsplitLG=bwareaopen(BWsplit, Min_MT_Length);
%3. get coordinates for regions
statsBWskel=regionprops(BWsplitLG,'PixelIdxList','PixelList');
MT_count(count_no)=numel(statsBWskel);
%4. Create image with a single myotube
cellno=0;
MT_Count=0;
for kk=1:numel(statsBWskel)
if MT_Count <=300
    MT_Count =MT_Count + 1;
    SingleMTBW=zeros(size(BWsplitLG));
    idx1=statsBWskel(kk).PixelIdxList;
    SingleMTBW(idx1)=1;
%5. Widen single myotube
    seMT = strel('disk',ReferenceWidthDil);
    SingleMTBW=imdilate(SingleMTBW,seMT);
    SingleMTimagestats=regionprops(SingleMTBW,'PixelIdxList','PixelList','BoundingBox');

    if Stri_im==1
        %'Cut out' actin stained myotube section for striation analysis
        SingleActin=Actin.*SingleMTBW;
        Bounds=floor(SingleMTimagestats(1).BoundingBox);
        Bounds(Bounds==0)=1;

        SingleActin=SingleActin(Bounds(2):Bounds(2)+Bounds(4),Bounds(1):Bounds(1)+Bounds(3));
        ImW=Bounds(3);
        ImH=Bounds(4);
    %Give cut images same dimensions
        SingleActin=trim_pad_sample(ImW, ImH, L1, SingleActin);
    %MaskSample=trim_pad_sample(ImW, ImH, L1, MaskSample);
    %ActinSample=trim_pad_sample(ImW, ImH, L1, ActinSample);
    %Store the coordinate of the centre of the image
        Centre_x_coord=Bounds(1)+(ImW/2);
        Centre_y_coord=Bounds(2)+(ImH/2);

        AA=AA+1;
        SingleActin_Image(:,:,AA)=SingleActin(1:L1,1:L1);
    %MaskSample_Image(:,:,AA)=MaskSample(1:L1,1:L1);
    %ActinSample_Image(:,:,AA)=ActinSample(1:L1,1:L1);
        Centre_coords_Image(AA,:)=[Centre_x_coord,Centre_y_coord,p,count_no];
    %figure;imshow(SingleActin)
    end
    if Uni==1
%6. Check for myonuclei within myotube
%Find regions with both myotube and myonuclei
    SingleMTBWnuc=SingleMTBW & Myonuclei;
%define centroids of these regions
    tempstatsnuc=regionprops(SingleMTBWnuc,'centroid');
%7. Find distances between centroids of myonuclei
    if numel(tempstatsnuc)>1
        for m=1:numel(tempstatsnuc)
            A(m,:)=tempstatsnuc(m).Centroid;
        end
        B=pdist(A);
        D=squareform(B);
        cellno=cellno+1;
        A=[];
    %find distance between nuclei in order
        [Dist_mean(cellno), Dist_SD(cellno)]= MeanDistanceFunction(D);
    %figure;imshow(SingleMTBWnuc);
    end
    end
end
end
if Uni==1
    matrixname_Dist=sprintf('Nuc_Dist_Mean-%s%d',Expt,p);
    matrixname_DistSD=sprintf('Nuc_Dist_SD-%s%d',Expt,p);
    %writematrix(Dist_mean,matrixname_Dist);
    %writematrix(Dist_SD,matrixname_DistSD);
    csvwrite(matrixname_Dist,Dist_mean);
    csvwrite(matrixname_DistSD,Dist_SD);
end

%Sum of single nuclei, estimated multi-nuclei and manually added myonuclei 
%minus manually removed nuclei and nuclei labelled as myocytes.

MBTot=length(MBSingNuc)+AddMB-RemMB+MBMultiEst;
MBTot=round(MBTot,0);

MTEst2=length(MTSingNuc)+AddMT-RemMT+MTMultiEst;
MTEst2=round(MTEst2,0);

Myoblast_estimate(count_no)=MBTot;

MyonucleiEstimate(count_no)=MTEst2;

MatureMyonucleiCount(count_no)=MT_Mature;

if Wid==1
Myonuclei_per_um(count_no)=MTEst2/TotLength;
end

if FI==1
%Fusion index
FIndex(count_no)=MTEst2/(MBTot+MTEst2);
end
clear Actin C MBNuc Nuc NucBW NucDil Pericentrin Nuclei Pericentrin_Nuclei I Myotube BWskel ...
    holes bigholes smallholes Filled BTX Myonuclei LgMyonuclei MeanInt BWBranch BWsplit BWsplitLG...
    NucMB NucMB_BW NucMT NucMT_BW SingleMTBW SingleMTBWnuc Dead_cells
end
function [MTDens,MTD]=MTDensity(Myotubes,s1,s2)
MTD=sum(sum(Myotubes));
TotA=s1*s2;
MTDens=MTD/TotA;
end


function [xt1,yt1,xt2,yt2]=ManLab(Iin)
[sizex]=size(Iin(:,1,:,1));
[sizey]=size(Iin(1,:,:,1));
for row = 1 : 500 : sizex(1)%3154
 line([1, sizey(2)], [row, row]);
end
for column = 1 : 500 : sizey(2)%3954
 line([column , column ], [1, sizex(1)]);
end
xt1 = 0;
yt1 = 0;
xt2 = 0;
yt2 = 0;
x = 0;
y = 0;
button = 1;
i = 0;
j = 0;
while button <=3
    [xg,yg,button] = ginput(1);
    
    if button==1
        j = j+1;
        xt1(j) = xg; yt1(j)= yg;
        hold on
        figure(1);plot(xt1(j),yt1(j),'b.','MarkerSize',10)
        a = [j]'; b = num2str(a); c = cellstr(b);
        drawnow
    end
    
    if button==3
        i = i+1;
        xt2(i) = xg; yt2(i)= yg;
        hold on
        figure(1);plot(xt2(i),yt2(i),'r.','MarkerSize',10)       
        a = [j]'; b = num2str(a); c = cellstr(b);
        drawnow
    end
end
end

function [Dist_mean, Dist_SD]=MeanDistanceFunction(D)
%function for finding the nearest neighbour distances of myonuclei within a
%single myotube section.
maximum = max(max(D));
[x,y]=find(D==maximum);
i=y(1);
j=0;
while i~=x(1)
    j=j+1;
    MxLine=D(i,:);
    Dist(j)=min(MxLine(MxLine>0));
    D(i,:)=0;
    D(:,i)=0;
    i=find(MxLine==Dist(j));
end
Dist_mean=mean(Dist);
Dist_SD=std(Dist);
end

function [Metrics,Auto]=Metrics_choose 
%List to choose metrics to analyse
list = {'Nuclei count','Myotube area & width','Nuclei spatial distribution','Acetylcholine receptors','Striation images','Fusion Index'};
[indx] = listdlg('ListString',list,'PromptString','Choose metrics');
for i=1:6
    if ismember((i),indx)
        Metrics(i)=1;
    else
        Metrics(i)=0;
    end
end
if Metrics(1)==1
    list2 = {'Automated','Semi-Automated'};
                    [indx2] = listdlg('ListString',list2,'PromptString','Automation of nuclei count');
                    if ismember(1,indx2)
                        Auto=1;
                    else
                        Auto=0;
                    end
else
    Auto=0;
end
end

function SingleActin=trim_pad_sample(ImW,ImH,L1,SingleActin)
if ImW > L1
    Diff=ceil((ImW-L1)/2);
    SingleActin=SingleActin(:,Diff:ImW-Diff);
end
if ImH > L1
    Diff=ceil((ImH-L1)/2);
    SingleActin=SingleActin(Diff:ImH-Diff,:);
end 
ImSize=size(SingleActin);
if ImSize(1)<L1
    Diff=ceil((L1-ImSize(1))/2);
    SingleActin=padarray(SingleActin,[Diff 0],0,'both');
end
if ImSize(2)<L1
    Diff=ceil((L1-ImSize(2))/2);
    SingleActin=padarray(SingleActin,[0 Diff],0,'both');
end
end

function [SingNuc,MultiNuc,MultiEst]=nuc_count(Nuc,MCmean)
StatsNucPx=regionprops('table',Nuc,'PixelList','Centroid','area');
SingNuc=[];
MultiNuc=[];
SingNucLg=[];
%nuclei size correction based on mean
NucArea=StatsNucPx.Area;
i=0;
j=0;
MultiNuc=0;
for a=1:length(NucArea)
    if NucArea(a)<= 2*MCmean
        i=i+1;
        SingNuc(i)=NucArea(a);
    else 
        j=j+1;
        MultiNuc(j)=NucArea(a);
    end
end

i=0;
for a=1:length(SingNuc)
    if SingNuc(a)>= 100
        i=i+1;
        SingNucLg(i)=SingNuc(a);
    end
end
%Estimate of area of multiple nuclei by dividing total area by mean of
%larger individual nuclei
MultiEst=sum(MultiNuc)/(mean(SingNucLg));
end