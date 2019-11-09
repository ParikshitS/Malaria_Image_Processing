function varargout = home(varargin)
% HOME MATLAB code for home.fig
%      HOME, by itself, creates a new HOME or raises the existing
%      singleton*.
%
%      H = HOME returns the handle to a new HOME or the handle to
%      the existing singleton*.
%
%      HOME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HOME.M with the given input arguments.
%
%      HOME('Property','Value',...) creates a new HOME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before home_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to home_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help home

% Last Modified by GUIDE v2.5 21-Apr-2017 09:17:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @home_OpeningFcn, ...
                   'gui_OutputFcn',  @home_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before home is made visible.
function home_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to home (see VARARGIN)

% Choose default command line output for home
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes home wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = home_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = msgbox({'Doctors use thick and thin blood smears to determine whether you have malaria','Thick blood smears are most useful for detecting the presence of parasites, because they examine a larger sample of blood.','Thin Films are used for Percentage of red blood cells infected (parasite density), not the number of parasites'});
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set(ch, 'FontName', 'comic');
set(ch, 'FontSize', 10);
pos = get( h, 'Position' );
set( h, 'Position', pos+50 );

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
chos=0;
possibility=7;

while chos~=possibility
    chos=menu('AUTOMATED MALARIA PARASITE DETECTION ','LOAD IMAGE','PREPROCESS','EXTRACTION','SEGMENTATION','MORPHOLOGY OPERATION','TOTAL CELLS','TOTAL INFECTED');

if chos==1,
[f,p]=uigetfile('.jpg','Input test RBC image');
I=imread([p,f]);
figure,imshow(I),title('input image');
end


if chos==2,
%% RGB2gray image 
I2=rgb2gray(I);
subplot(3,3,1);imshow(I2);title('RGB TO GRAYSCALE');
%% Preprocessing 


I3=medfilt2(I2,[3 3]);
subplot(3,3,2);imshow(I3);title('FILTERED IMAGE');

y1=histeq(I3);
subplot(3,3,3);imshow(y1);title('HISTOGRAMED IMAGE');

end
if chos==3,    %Extraction
%% Extracting the blue plane 
bPlane = I(:,:,3)  - 0.5*(I(:,:,1)) - 0.5*(I(:,:,2));
subplot(3,3,1);
imshow(bPlane);title('EXTRACTION-1');

%% Extract out purple cells

BW = bPlane > 35;

subplot(3,3,2);
imshow(BW);
title('EXTRACTION-2');

%% Remove noise 100 pixels or less
BW2 = bwareaopen(BW, 100);
subplot(3,3,3)
himage7=imshow(BW2);
title('NOISE REMOVAL');

%% contrast 
I5=imadjust(I3);
subplot(3,3,4);imshow(I5);
title('INTENSITY ADJUSTMENT');
f=graythresh(I2);

I6=im2bw(I5,f);

subplot(3,3,5);imshow(I6);
title('BINARY IMAGE');

end
if chos==4,% Segmentation


chos1=0;
possib_seg=2;
while chos1~=possib_seg
   chos1=menu('Choose segmentation method', 'Sobel', 'Prewitt');
    
I7=bwareaopen(I6,50);
%subplot(3,3,1);imshow(I7)
%title('AREA OPENING');
if chos1==1,
 [~, threshold] = edge(I7, 'sobel');
fudgeFactor = .5;
BWs = edge(I7,'sobel', threshold * fudgeFactor);
subplot(1,1,1);imshow(BWs);
title('Sobel Segmentation');
break

else if chos1==2,
   
    BWs = edge(I7,'Prewitt')
    subplot(1,1,1);imshow(BWs);
    title('Prewitt Segmentation');
    end 
end 

end

end

if chos==5,
%%morphology

imshow(BWs), title('GRADIENT MASK');
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);

BWsdil = imdilate(BWs, [se90 se0]);
subplot(3,3,1);
imshow(BWsdil);
title('DILATION');

BWdfill = imfill(BWsdil, 'holes');
label=bwlabel(BWdfill);
RBCcount=max(max(label));
% Rmin = 30;
% Rmax = 65;
% [centersBright, radiiBright] = imfindcircles(BWdfill,[Rmin Rmax],'ObjectPolarity','bright');
% %viscircles(centersBright, radiiBright,'Color','b');
subplot(3,3,2);imshow(BWdfill);
title('HOLES FILLING')


I2=rgb2gray(I);
I_eq = adapthisteq(I2);
subplot(3,3,3);imshow(I_eq);
title('ADAPTIVE HISTOGRAM');
bw = im2bw(I_eq, graythresh(I_eq));
%bw_comp=imcomplement(bw)
 %%figure, imshow(bw)
% %  bw2 = imfill(bw,'holes');
% % bw3 = imopen(bw2, ones(5,5));
% % bw4 = bwareaopen(bw3, 40);
bw4_hbreak = bwmorph(bw,'hbreak');
bw4_diag = bwmorph(bw4_hbreak,'diag');
bw4_skel = bwmorph(bw4_diag,'skel');
bw4_perim = bwperim(bw);
bw_comp=imcomplement(bw4_perim);
overlay1 = imoverlay(I_eq, bw4_perim, [.3 1 .3]);
% label=bwlabel(bw_comp)
% totalRBC=max(max(label))
overlay1 = imoverlay(mat2gray(I_eq), bw4_perim, [.3 1 .3]);
subplot(3,3,4);imshow(overlay1);
title('overlay of original image');
mask_em = imextendedmax(I_eq, 30);
subplot(3,3,5);
imshow(mask_em);title('masked image');

%%title('MASKED IMAGE');
%mask_em = imclose(mask_em, ones(5,5));
%mask_em = imfill(mask_em, 'holes');
%mask_em = bwareaopen(mask_em, 40);
overlay2 = imoverlay(I_eq, bw4_perim | mask_em, [.3 1 .3]);
subplot(3,3,6);imshow(overlay2);
title('overlay of masked image');
end
if chos==6,

figure,imshow(overlay1)
title('TOTAL NO OF CELLS')

I8=imfill(I6,'holes');
%%figure,imshow(I8)
end
if chos==7,

L=bwlabel(BW2);
%% Superimpose onto original image
figure, imshow(I), hold on
himage = imshow(BW2);

% % % % BW=imbinarize(I8);
% BW_temp1 = im2double(I8);
% BW_temp=graythresh(BW_temp1);
% 
% BW=im2bw(BW_temp);
set(himage, 'AlphaData', 0.5);
%k=imfill(BW,'holes');
totalinfectedcells=max(max(L));
title('TOTAL NO Of Malaria Parasites Found:')

%totalcells=max(max(centersBright))
%title('TOTAL NO OF RBC1 CELLS')
% % msgbox( sprintf('The total number of cells are %f', RBCcount));
% % if totalinfectedcells==0
% %     msgbox( sprintf('The total infected cells are 0'));
% % else
% % %msgbox( sprintf('The total number of parasites are %f', totalinfectedcells));
% % 
% % end
parasite_percent=totalinfectedcells/RBCcount*100;
set(handles.text9,'String',RBCcount);
set(handles.text10,'String',totalinfectedcells);
set(handles.text13,'String',parasite_percent);
end
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
choss=0;
possibilitys=5;
while choss~=possibilitys
    choss=menu('AUTOMATED MALARIA PARASITE DETECTION ','LOAD IMAGE','PREPROCESS','WBC Detection','SEGMENTATION','Infected');
   if choss==1,
[f,p]=uigetfile('.jpg','Input test RBC image');
I=imread([p,f]);
figure,imshow(I),title('input image');
   end
if choss==2,
I_bw = im2bw(I, 0.5)
I_comp=imcomplement(I_bw)
I_WBC = bwareaopen(I_comp, 400);
true_parasite=imsubtract(I_comp,I_WBC)
figure,imshow(I_comp),title('Preprocess image');
end
if choss==3,
  figure,imshow(I_WBC),title('WBC image');  
end
if choss==4,
     cform = makecform('srgb2lab');
lab_he = applycform(I,cform);

ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

nColors = 3;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);
                                  
pixel_labels = reshape(cluster_idx,nrows,ncols);
imshow(pixel_labels,[]), title('image labeled by cluster index');

segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 3]);

for k = 1:nColors
    color = I;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
end

subplot(2,2,1);imshow(segmented_images{1});title('objects in cluster 1');


subplot(2,2,2);imshow(segmented_images{2});title('objects in cluster 2');

subplot(2,2,3);imshow(segmented_images{3});title('objects in cluster 3');

mean_cluster_value = mean(cluster_center,2);
[tmp, idx] = sort(mean_cluster_value);
blue_cluster_num = idx(1);

L = lab_he(:,:,1);
blue_idx = find(pixel_labels == blue_cluster_num);
L_blue = L(blue_idx);
%%is_light_blue = imbinarize(L_blue);
level = graythresh(L_blue)
is_light_blue = im2bw(L_blue, level)

nuclei_labels = repmat(uint8(0),[nrows ncols]);
nuclei_labels(blue_idx(is_light_blue==false)) = 1;
nuclei_labels = repmat(nuclei_labels,[1 1 3]);
blue_nuclei = I;
blue_nuclei(nuclei_labels ~= 1) = 0;
%subplot(1,1,1);imshow(blue_nuclei);title('blue nuclei');
 
end
if choss==5,
figure,imshow(true_parasite),title('Infected  image');
L1= bwlabel(true_parasite)
true_parasite_count= max(max(L1))
%figure, imshow(I), hold on
%himage1 = imshow(true_parasite);
%msgbox( sprintf('The total parasites are %f', true_parasite_count));
set(handles.text11,'String',true_parasite_count);
end
end
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;
