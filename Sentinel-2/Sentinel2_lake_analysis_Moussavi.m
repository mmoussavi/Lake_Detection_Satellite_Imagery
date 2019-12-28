%%% This code requires python gdal library
%%%set python environment: contains gdal 
setenv('PATH', '/Users/...../Anaconda-3/anaconda3/envs/py37/bin:/Users/mmoussavi/google-cloud-sdk/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin');
% choose directory where Sentinel-2 image folders are stored
pathFolder = uigetdir('/..../Sentinel-2/Data/');
directory=horzcat(' ',pathFolder, '/'); 
d = dir(pathFolder);
isub = [d(:).isdir]; %# returns logical vector
Scenes = {d(isub).name}';
Scenes(ismember(Scenes,{'.','..'})) = [];

addpath ('.......') %% to be filled out by user (directory where codes and dependencies are stored)
home = pwd;
addpath(pwd);

%%% lake depth calculation  parameters
g_B4 =  0.83; % Band 4 %g values calculated according to Williamson et al (2018)
date_number = repmat(datetime(0,0,0), (size(Scenes)));


%%%Reflectance of optically-deep water for S2 coastal scenes around
%%%Antarctica -- Compiled in Google Earth Engine -- fore more information
%%%on this dataset, please refer to "Antarctic Supraglacial Lake Detection Using
%%%Landsat 8 and Sentinel-2 Imagery: towards Continental Generation of Lake
%%%Volumes" published in Remote Sensing.
Rinf_Table=tabularTextDatastore('Sentinel2_Antarctica_Rinf_2014_2019.csv','TreatAsMissing','NA','MissingValue',0, 'Delimiter', ',');
while hasdata(Rinf_Table)
    Rinf_data = read(Rinf_Table);
end

Rinf_B4_values=Rinf_data.MIN_B4_5perc;
Rinf_Tile_Numbers=string(Rinf_data.TILE);
Rinf_UTMzones=double(extractBetween(Rinf_Tile_Numbers, 1, 2));
Rinf_ID_values=string(Rinf_data.SCENE_ID);
Rinf_dates=extractBetween(Rinf_ID_values, 20, 27);
 

volume = zeros(size(Scenes));
cloud_percent = zeros(size(Scenes));
pixel = 10;%% meters

for scene = 1:size(Scenes,1)

    %change to appropriate directory
    expression = strcat('cd ', ' ', directory,Scenes{scene},'/GRANULE/');
    eval(expression);
    folder=ls;
    IMG_Dir=horzcat(' ',folder); 
    expression = strcat('cd ', ' ', IMG_Dir,'/IMG_DATA');
    eval(expression);
    
    
    %%%%%%%% find Rinf values for B4 and B8 using Sentinel-2 scenes (closest
    %%%%%%%% in time and space) that contain open ocean water pixels (very dark pixels)
    
    %%%%find image date
    SceneID=string(Scenes{scene});
    Date=datetime(extractBetween(SceneID,12, 19), 'InputFormat',  'yyyyMMdd');

    %%% find closest dates
    NumDays = daysact(datetime(Rinf_dates, 'InputFormat',  'yyyyMMdd'), Date);
    date_threshold=prctile(abs(NumDays), 5);
    index_close_dates=find(abs(NumDays) <= date_threshold);
    close_dates=Rinf_dates(index_close_dates);

    %%%% find closest UTM zones
    image_Tile=extractBetween(SceneID,40, 44);
    image_utmzone=double(extractBetween(SceneID,40, 41));
    
    distance=Rinf_UTMzones-image_utmzone;
    
    index_close_pairs=find(abs(distance) <= prctile(abs(distance), 5));
    close_pairs=Rinf_Tile_Numbers(index_close_pairs);
    %%%%%% find closest Tiles in space and time
    index_Rinf_candidates_space_Time=intersect(index_close_dates,  index_close_pairs);
    date_candidates=Rinf_dates(index_Rinf_candidates_space_Time);
    Tile_candidates=Rinf_Tile_Numbers(index_Rinf_candidates_space_Time);
    B4_min_candidates=Rinf_B4_values(index_Rinf_candidates_space_Time);
    index_1000_B4=find((B4_min_candidates >0 & B4_min_candidates <1000 ));
    final_values_B4_Rinf=(B4_min_candidates(index_1000_B4 ))/10000.;
    final_dates_B4_Rinf=date_candidates(index_1000_B4);
    final_Tiles_B4_Rinf=Tile_candidates(index_1000_B4);
    
    B4_Rinf_final=median(final_values_B4_Rinf);
    
    %%%%% If no values are found, then the long-term Rinf values will be
    %%%%% assigned to the Rinf parameter --- the code was tested many times and
    %%%%% this scenario never happened.
    if isfile('S2_B4.tif')
    %%%tiff file exists
       B4_tiffinfo=geotiffinfo('S2_B4.tif'); 
       B4_info= B4_tiffinfo.SpatialRef;
    else 
        %%%% convert one of the .jp2 files to get the geographic info
        commandstr=('gdal_translate *B04.jp2 S2_B4.tif');
        [status4, commandOut4]=system(commandstr); 
        B4_tiffinfo=geotiffinfo('S2_B4.tif'); 
        B4_info= B4_tiffinfo.SpatialRef;
    end

  
    B22=dir('*B02.jp2'); B2= imread(B22.name);
    B33=dir('*B03.jp2');B3= imread(B33.name);
    B44=dir('*B04.jp2');B4= imread(B44.name);
    B1010=dir('*B10.jp2');B10= imread(B1010.name);
    B1111=dir('*B11.jp2');B11= imread(B1111.name);

    %%%resample the original image (20 m resolution) to 10 m resolution to
    %%%%match Band 3, this is necessary step for NDSI computation later on

    B11= imresize(B11, 2, 'nearest');

    %%%% resample B10 from 60 m to 10 m

    B10 = imresize(B10, 6, 'nearest');

    %%%% change data type
    B2=double(B2);
    B3=double(B3);
    B4=double(B4); 
    B10=double(B10);
    B11 = double(B11);

    %%%%%%%%% Calculation of NDWI
    NDWI = (B2-B4)./(B2+B4); % This is NDWI
    NDWI = NDWI*10000;
    NDWI((NDWI==NDWI(1))) = 0;
    NDWI(1,1) = 2^16-1;
    NDWI = uint16(NDWI);
    NDWI(1,1) = 0;

    %%%%%%%Calculation of  NDSI
    NDSI = (B3-B11)./(B3+B11); % This is NDSI
    NDSI = NDSI*10000;
    NDSI((NDSI==NDSI(1))) = 0;
    NDSI(1,1) = 2^16-1;
    NDSI = uint16(NDSI);
    NDSI(1,1) = 0;


    %%%%%%%% Rock Masking
    rock= zeros(size(B4));
    rock= uint8(rock);
    rock(NDSI <9000 & B2>0 & B3<4000 &  B2 < 4000 )=1; %%% version 3
    rock_mask=rock;


    %%%% Cloud mask
    cloud_mask = zeros(size(B4));
    cloud_mask = uint8(cloud_mask);
    cloud_mask(B10 > 30 & B11 > 1100 & B2 > 6000 & B2<9700 )=1;%%% version 2
    cloud_mask=standardizeMissing(cloud_mask, 0);
    cloud_mask(rock_mask ==1)=0;


    %%%%%% lake_mask
    lake_mask = zeros(size(B4));
    lake_mask = uint8(lake_mask);
    lake_mask(NDWI>1800 & (B3-B4)>800 & (B3-B4)< 4000  &(B2-B3)>700 & B11 > 10 & B2> 4000)=1; %%%%% version 3 (B11 > 10 gets rid of some cloud shadows)
    lake_mask=standardizeMissing(lake_mask, 0);
    lake_mask(cloud_mask==1)=0;
    lake_mask(rock_mask==1)=0;

    CC = bwconncomp(lake_mask,26);
    min_pixel=45;%%45x (10x10)m2)=4500 m2 :this is equivalent to 5 landsat pixels (30x30)m2: 5x900=4500 m2
    min_width=6;%%% this is equivalent to two landsat pixels
    
    for lake = 1:CC.NumObjects

        if size(CC.PixelIdxList{lake},1) < min_pixel 
            lake_mask(CC.PixelIdxList{lake}) = 0;

        else %all same row or all same column
            [r,c] = ind2sub(size(lake_mask),CC.PixelIdxList{lake});
            r = unique(r);
            c = unique(c);
            if size(r,1) == min_width || size(c,1) == min_width
                lake_mask(CC.PixelIdxList{lake}) = 0;
            end
        end
    end           

    %%%write all masks to a geotiff file
    All_Masks = zeros(size(B4));
    All_Masks = uint8(All_Masks);
    All_Masks(lake_mask==1)=1;
    All_Masks(rock_mask==1)=2;
    All_Masks(cloud_mask==1)=3;

        %write All_Masks
        expression = strcat('geotiffwrite(''',Scenes{scene},'_All_Masks.TIF'', All_Masks, B4_info, ''GeoKeyDirectoryTag'',B4_tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);');
        eval(expression)

   %%%%%%%%%%%% Create Browse images of masks superimposed on RGB composite
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RGBimage = cat(3, imadjust(uint16(B4)),  imadjust(uint16(B3)),imadjust(uint16(B2)));
    RGB_resized=imresize(RGBimage, 0.05, 'nearest');           
    B = imoverlay(RGB_resized,imresize(rock_mask, 0.05, 'nearest'),[0.12 0.32 0.8]); 
    C= imoverlay(B,imresize(cloud_mask, 0.05, 'nearest'), [0.4 .58 .5]);
    D=imoverlay(C,imresize(lake_mask, 0.05, 'nearest'), [0.8 0.2 0.1 ]);

    name_mask=strcat(Scenes{scene},'_masks_overlaid_2.png' );
    name_orig=strcat(Scenes{scene},'_original.png' );
    imwrite(((RGB_resized)), name_orig)
    imwrite((D), name_mask) 

   %%%%% Depth Calculation using Band 4: Red Band
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mask=lake_mask;
    image=B4; 

    %read in cloud mask
    %calculate cloud cover percentage
    id=find(cloud_mask == 1);
    id0=find(image > 0);
    siz_id=size(id);
    siz_id0=size(id0);
    cloud_cover_percent=100*siz_id(1)/siz_id0(1);
    cloud_percent(scene) = cloud_cover_percent; 
 
    %%% lake mask 
    mask = mask(:,:,1);
    mask_thick = bwmorph(mask,'thicken',3);
    mask = bwmorph(mask_thick,'thin',3);
    mask_thick = double(mask_thick);
    mask = double(mask);
       
    %format files appropriately
    image = double(image);
    image = image/10000;
    z = zeros(size(mask));

    %calculating lake-specific albedo
    CC = bwconncomp(mask_thick);
    num_lakes = CC.NumObjects;
    Ad_lakes = zeros(num_lakes,1);

    lake_boundaries = bwboundaries(mask_thick,8,'noholes');

    for lake=1:num_lakes
        bound = lake_boundaries{lake};
        bound = sub2ind(size(mask),bound(:,1),bound(:,2));
        Ad_lakes(lake) = mean(image(bound));
    end

    CC = bwconncomp(mask);
    props = regionprops(CC,mask,{'PixelIdxList'}); 
    Ad = zeros(size(mask));
    
    for lake = 1:num_lakes
        Ad(props(lake).PixelIdxList) = Ad_lakes(lake);
    end

    %filtering values which would cause a negative log
    temp = image-B4_Rinf_final;
    index = find(temp>0);

    temp2 = Ad-B4_Rinf_final;
    index2 = find(temp2>0);

    index = intersect(index,index2);

    %calculating lake depth
    z(index) = (log(Ad(index)-B4_Rinf_final)-log(temp(index)))/(g_B4);
    z_masked = mask.*z*1000; %can do this because the LakeMask is a 1/0 binary, written in mm
    z_masked(z_masked < 0) = 0;

    %write lake depths out
    z_masked(1,1) = 2^16-1;
    z_masked = uint16(z_masked);
    z_masked(1,1) = 0;
    
    expression = strcat('geotiffwrite(''',Scenes{scene},'_B4_depth.tif'', z_masked, B4_info, ''GeoKeyDirectoryTag'',B4_tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);');
    eval(expression)

     B4_Depth=z_masked;


    %%%% Calculating total volume per scene
    s = Scenes{scene};
    year = s(18:21);
    month=s(22:23);
    dayy = s(24:25);
    date_vol=s(18:25);
    
    date_number(scene)=datetime(str2double(year),str2double(month),str2double(dayy));
    
    image=z_masked;
    
    vol = sum(image(image > 0)) * pixel * pixel / 1000; %depth in mm, 10 m pixels
    volume(scene) = vol;    
    
   %%%%% generating lake depth and area distribution histograms
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Depths=double(image(image > 0)); 
   
   avg_Depth=mean(Depths(:))/1000; % in meters
   max_depth=ceil(max(Depths(:))/1000);% in meters
   
   %%%%% Lake depth distribution
   if isempty(Depths) ==0
   figure('rend','painters','pos',[10 10 1100 600],'visible','off')
   subplot(1,2,1);
   number_of_bins = 100;
   % Fit a non-parametric kernel smoothing distribution
   h=histfit(Depths/1000,number_of_bins,'kernel');
   h(1).FaceColor = [.85 .85 0.85];
   h(2).LineWidth=4;h(2).Color=[.8 .3 0.4];
   xlabel('Depth (m)', 'FontWeight', 'bold', 'FontSize', 25);ylabel('Pixel Count','FontWeight', 'bold','FontSize', 25);
   ax=gca; ax.FontSize=22;ax.XLim=[0,max_depth];ax.LineWidth=2;ax.GridLineStyle=':';
   e1=strcat('Total Volume=  ', '  ', ' ', num2str(fix(vol/1e6)),'x10^6 m^3');
   t=text (avg_Depth+0.4,max(h(1).YData)*0.9 , e1);t.FontSize=22;t.FontWeight='bold';
   e2=strcat('Mean Depth=  ', '  ', num2str(round(mean(Depths)/1000, 2)),'  m');
   t2=text(avg_Depth+0.4, max(h(1).YData)*0.8, e2);t2.FontSize=22;t2.FontWeight='bold';
   
   %%%%% Lake area distribution   
       CC = bwconncomp(image);
       area_pixel=zeros(CC.NumObjects);
       for lake=1:CC.NumObjects
           area_pixel(lake)=size(CC.PixelIdxList{lake},1);
       end
   subplot(1,2,2) ;
   Area_1=area_pixel*100*1e-6; %% this calculates area in square kilometers, each sentinel pixel is 10 m.
   idx_outliers_removed=find(Area_1 > 0.0005 );%%% removing lakes that are less than 5 pixels in area.
   Area=Area_1(idx_outliers_removed);
   max_Area=ceil(max(Area(:)));
   %%% find percentage of lakes that are less than 1 square kilometer in area
   id_1=find(Area < 1); 
   percentage_less_1=numel(id_1)*100/numel(Area); 
   text_less_1=strcat( num2str(fix(percentage_less_1)), '% of lakes < 1 km^2');
   text_total=strcat('Total lake area= ', num2str(fix(sum(Area))), ' km^2');
   cdf_area=cdfplot(Area);
   hold on;
   x=1;
   y =(percentage_less_1/100); 
   Horz_line=line([0,fix(max(Area))],[y,y]);Horz_line.LineWidth=4; Horz_line.LineStyle=':';hold on;
   vert_line=line([x,x], [0,y] );vert_line.LineWidth=4; vert_line.LineStyle=':';
   ax=gca; ax.FontSize=22;ax.XLim=[0,max_Area];ax.LineWidth=2;cdf_area.LineWidth=7;ax.GridLineStyle=':';cdf_area.Color=[0.64 0.08 0.18] ;
   xlabel('Lake Area (km^2)', 'FontWeight', 'bold', 'FontSize', 25);ylabel('Cumulative Distribution','FontWeight', 'bold','FontSize', 25);
      t_less_1=text(1.2, 0.85, text_less_1);t_less_1.FontSize=22; t_less_1.FontWeight='Bold';
      t_total=text(1.2, 0.75, text_total);t_total.FontSize=22;t_total.FontWeight='Bold';
   exp_txt2=strcat(Scenes{scene}, '_lake_Depth_area_Distribution.png');
   saveas(gcf,exp_txt2)
   end

    
end
cd ..
filename_output_volume='output_volume.mat';
save (filename_output_volume,'Scenes',  'date_number','volume', 'cloud_percent');

expression = strcat('cd ''',home,'''');
eval(expression)

