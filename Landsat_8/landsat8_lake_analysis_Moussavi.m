function landsat8_lake_analysis_Moussavi(pathFolder)

%%% pathfolder: direcotry in which your Landsat 8 image folders are stored.
directory=horzcat(' ',pathFolder, '/'); 
d = dir(pathFolder);
isub = [d(:).isdir]; %# returns logical vector
Scenes = {d(isub).name}';
Scenes(ismember(Scenes,{'.','..'})) = [];
addpath ('/Users/mmoussavi/Documents/Antarctic_Lakes_king_Mapped/code/Main_Codes')
home = pwd;
addpath(pwd);



%%% lake depth calculation model parameters: Attenuation Coefficient
g_B4 =  0.8; % Band 4 -- (regressed g values), Pope et al., 2015
g_B8 =  0.36; % Band 4 -- (regressed g values), Pope et al., 2015

date_number = repmat(datetime(0,0,0), (size(Scenes)));       
volume = zeros(size(Scenes));
cloud_percent = zeros(size(Scenes));
pixel = 30; % meters


%%%Reflectance of optically-deep water for L8 coastal scenes around
%%%Antarctica -- Compiled in Google Earth Engine -- fore more information
%%%on this dataset, please refer to "Antarctic Supraglacial Lake Detection Using
%%%Landsat 8 and Sentinel-2 Imagery: towards Continental Generation of Lake
%%%Volumes" published in Remote Sensing.

Rinf_Table=tabularTextDatastore('/Users/mmoussavi/Google Drive/something/Antarctica_Landsat8_Rinf_2014_2019.csv','TreatAsMissing','NA','MissingValue',0);
while hasdata(Rinf_Table)
    Rinf_data = read(Rinf_Table);
end
    Rinf_dates=Rinf_data.DATE;
    Rinf_B4_values=Rinf_data.MIN_B4_5perc;
    Rinf_B8_values=Rinf_data.MIN_B8_5perc;
    Rinf_paths=Rinf_data.PATH;
    Rinf_rows=Rinf_data.ROW;

  
for scene = 1:size(Scenes,1)
    %change to appropriate directory
    expression = strcat('cd ', ' ', directory,Scenes{scene},'/');
    eval(expression);
    
    %read metadata
    expression = strcat('mtl =ls8_meta_parser(''',Scenes{scene},'_MTL.txt'');');
    eval(expression);
    
    %%%%%%%% find Rinf values for B4 and B8 using landsat scenes (closest
    %%%%%%%% in time and space) that contain open ocean water pixels (very dark pixels)
    
    Date=mtl.PRODUCT_METADATA.DATE_ACQUIRED;
    %%% find closest dates
    NumDays = daysact(datetime(Rinf_dates), Date);
    date_threshold=prctile(abs(NumDays), 5);
    index_close_dates=find(abs(NumDays) <= date_threshold);

    %%%% find closest path-row pairs
    image_path=mtl.PRODUCT_METADATA.WRS_PATH;
    image_row=mtl.PRODUCT_METADATA.WRS_ROW;
    
    distance=sqrt((Rinf_paths-image_path).^2+(Rinf_rows-image_row).^2);
    index_close_pairs=find(abs(distance) <= prctile(distance, 5));
    
    %%%%%% find close path-row pairs in space and time
    index_Rinf_candidates_space_Time=intersect(index_close_dates,  index_close_pairs);
    
    B4_min_candidates_dates=Rinf_B4_values(index_Rinf_candidates_space_Time);
    B8_min_candidates_dates=Rinf_B8_values(index_Rinf_candidates_space_Time);
    
    
    B4_Rinf_final=median(B4_min_candidates_dates(B4_min_candidates_dates >0 & B4_min_candidates_dates < 0.1 ));
    B8_Rinf_final=median(B8_min_candidates_dates(B8_min_candidates_dates >0 & B8_min_candidates_dates < 0.1));
    
    
    %%%%% If no values are found, then the long-term Rinf values will be
    %%%%% assigned to the Rinf parameter --- the code was tested many times and
    %%%%% this scenario never happened.
    if isnan(B4_Rinf_final)
        B4_Rinf_final=0.04;
    end
    if isnan(B8_Rinf_final)
        B8_Rinf_final=0.06;  
    end
    
    %%%%% Load Landsat8 Bands needed for Analysis
    
    %read in B2: Blue Band   
    expression = strcat('[B2, B2_info] = geotiffread(''',Scenes{scene},'_B2.TIF'');');
    eval(expression);
    expression = strcat('B2_tiffinfo = geotiffinfo(''',Scenes{scene},'_B2.TIF'');');
    eval(expression);
    
    B2 = double(B2);
    B2 = ((B2*mtl.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_2)+mtl.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_2)/sin(mtl.IMAGE_ATTRIBUTES.SUN_ELEVATION * pi / 180); %DN to TOA reflectance conversion
    B2 = B2*10000;
    B2((B2==B2(1))) = 0;
    B2(1,1) = 2^16-1;
    B2 = uint16(B2);
    B2(1,1) = 0;

    %read in B3: Green Band
    expression = strcat('[B3, B3_info] = geotiffread(''',Scenes{scene},'_B3.TIF'');');
    eval(expression);
    expression = strcat('B3_tiffinfo = geotiffinfo(''',Scenes{scene},'_B3.TIF'');');
    eval(expression);

    B3 = double(B3);
    B3 = ((B3*mtl.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_3)+mtl.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_3)/sin(mtl.IMAGE_ATTRIBUTES.SUN_ELEVATION * pi / 180); %DN to TOA reflectance conversion
    B3 = B3*10000;
    B3((B3==B3(1))) = 0;
    B3(1,1) = 2^16-1;
    B3 = uint16(B3);
    B3(1,1) = 0;
        
    %read in B4: Red Band
    expression = strcat('[B4, B4_info] = geotiffread(''',Scenes{scene},'_B4.TIF'');');
    eval(expression);
    expression = strcat('B4_tiffinfo = geotiffinfo(''',Scenes{scene},'_B4.TIF'');');
    eval(expression);

    B4 = double(B4);
    B4 = ((B4*mtl.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_4)+mtl.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_4)/sin(mtl.IMAGE_ATTRIBUTES.SUN_ELEVATION * pi / 180); %DN to TOA reflectance conversion
    B4 = B4*10000;
    B4((B4==B4(1))) = 0;
    B4(1,1) = 2^16-1;
    B4 = uint16(B4);
    B4(1,1) = 0;

    %read in B6; SWIR Band
    expression = strcat('[B6, B6_info] = geotiffread(''',Scenes{scene},'_B6.TIF'');');
    eval(expression);
    expression = strcat('B6_tiffinfo = geotiffinfo(''',Scenes{scene},'_B6.TIF'');');
    eval(expression);
    
    B6 = double(B6);
    B6 = ((B6*mtl.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_6)+mtl.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_6)/sin(mtl.IMAGE_ATTRIBUTES.SUN_ELEVATION * pi / 180); %DN to TOA reflectance conversion
    B6 = B6*10000;
    B6((B6==B6(1))) = 0;
    B6(1,1) = 2^16-1;
    B6 = uint16(B6);
    B6(1,1) = 0;

    %read in B8; Pnachromatic Band
    expression = strcat('[B8, ~] = geotiffread(''',Scenes{scene},'_B8.TIF'');');
    eval(expression);
    expression = strcat('B8_tiffinfo = geotiffinfo(''',Scenes{scene},'_B8.TIF'');');
    eval(expression);
  
    B8 = double(B8);
    B8=imresize(B8, 0.5, 'bilinear');
    B8 = ((B8*mtl.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_8)+mtl.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_8)/sin(mtl.IMAGE_ATTRIBUTES.SUN_ELEVATION * pi / 180); %DN to TOA reflectance conversion
    B8 = B8*10000;
    B8((B8==B8(1))) = 0;
    B8(1,1) = 2^16-1;
    B8 = uint16(B8);
    B8(1,1) = 0;
    
    % read in B10: Thermal Band
    expression = strcat('[B10, B10_info] = geotiffread(''',Scenes{scene},'_B10.TIF'');');
    eval(expression);
    expression = strcat('B10_tiffinfo = geotiffinfo(''',Scenes{scene},'_B10.TIF'');');
    eval(expression);

    B10 = double(B10);
    B10=(B10*mtl.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_10)+mtl.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_10;
    B10=mtl.TIRS_THERMAL_CONSTANTS.K2_CONSTANT_BAND_10./log((mtl.TIRS_THERMAL_CONSTANTS.K1_CONSTANT_BAND_10./B10)+1);
    B10 = B10*10;
    B10((B10==B10(1))) = 0;
    B10(1,1) = 2^16-1;
    B10 = uint16(B10);
    B10(1,1) = 0;
        

    %%%%%%% Calculation of NDWI
    B2 = double(B2);
    B4 = double(B4);
    NDWI = (B2-B4)./(B2+B4); % This is NDWI
    NDWI = NDWI*10000;
    NDWI((NDWI==NDWI(1))) = 0;
    NDWI(1,1) = 2^16-1;
    NDWI = uint16(NDWI);
    NDWI(1,1) = 0;  
    
    %%%%%% Calculation of NDSI
    B3 = double(B3);
    B6 = double(B6);
    NDSI = (B3-B6)./(B3+B6); % This is NDSI
    NDSI = NDSI*10000;
    NDSI((NDSI==NDSI(1))) = 0;
    NDSI(1,1) = 2^16-1;
    NDSI = uint16(NDSI);
    NDSI(1,1) = 0;

    %%%%%% Calculation of TIRS/Blue
    B10 = double(B10);
    B2 = double(B2);
    TIRS_Blue= B10./B2; % This is TIRS/Blue 
    TIRS_Blue = TIRS_Blue*10000;
    TIRS_Blue((TIRS_Blue==TIRS_Blue(1))) = 0;
    TIRS_Blue(1,1) = 2^16-1;
    TIRS_Blue = uint16(TIRS_Blue);
    TIRS_Blue(1,1) = 0;
 
    %%%%%%% Rock outcrop/Seawater Masking
    rock_mask = zeros(size(B4));
    rock_mask =uint8(rock_mask);
    rock_mask(TIRS_Blue > 6500 & B4 > 0 & B2 <3500)=1; %%% B4 filtering --> to avoid including image borders
    rock_mask=standardizeMissing(rock_mask, 0);
   
    %%%%%%% Cloud Masking
    cloud_mask = zeros(size(B4));
    cloud_mask = uint8(cloud_mask);
    cloud_mask(B6 > 1000 & NDSI < 8000 & B2 >6000  & B2 < 9500) = 1;
    cloud_mask=standardizeMissing(cloud_mask, 0);
    cloud_mask(rock_mask ==1)=0;
    
    %%%%% Lake Masking
    lake_mask = zeros(size(B4));
    lake_mask = uint8(lake_mask);
    lake_mask(NDWI > 1900  &(B3-B4)> 700 & (B2-B3)>1100) = 1; 
    lake_mask(rock_mask==1)=0;
    lake_mask(cloud_mask==1)=0;
    lake_mask=standardizeMissing(lake_mask, 0);

    %Remove "small" lakes(less than 5 pixels in area or narrower than 2 pixels)
    CC = bwconncomp(lake_mask,4);
    min_pixel=5;
    min_width=2;
    
    for lake = 1:CC.NumObjects  
        if size(CC.PixelIdxList{lake},1) < min_pixel 
            lake_mask(CC.PixelIdxList{lake}) = 0;    
        else 
            [r,c] = ind2sub(size(lake_mask),CC.PixelIdxList{lake});
            r = unique(r);
            c = unique(c);
            if size(r,1) == min_width || size(c,1) == min_width
                lake_mask(CC.PixelIdxList{lake}) = 0;
            end
        end
    end           
          
   %%%%% Depth Calculation using Band 4: Red Band
 
   mask=lake_mask;
   image=B4; 

   %read in cloud mask to calculate cloud cover percentage
   id=find(cloud_mask == 1);
   id0=find(image > 0);
   siz_id=size(id);
   siz_id0=size(id0);
   cloud_cover_percent=100*siz_id(1)/siz_id0(1);
   cloud_percent(scene) = cloud_cover_percent; 
 
   %%% lake mask 
   mask = mask(:,:,1);
   mask_thick = bwmorph(mask,'thicken',1);
   mask = bwmorph(mask_thick,'thin',1);
   mask_thick = double(mask_thick);
   mask = double(mask);
       
   %format files appropriately
   image = double(image);
   image = image/10000;
   
   %calculating lake-specific albedo
   z = zeros(size(mask));
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
    
   B4_Depth=z_masked;
%    expression = strcat('geotiffwrite(''',Scenes{scene},'_B4_depth.tif'', z_masked, B3_info, ''GeoKeyDirectoryTag'',B3_tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);');
%    eval(expression)
   
   %%%%% Depth Calculation using Band 4: Red Band
   mask=lake_mask;
   image=B8; 

   %%% lake mask 
   mask = mask(:,:,1);
   mask_thick = bwmorph(mask,'thicken',1);
   mask = bwmorph(mask_thick,'thin',1);
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
   temp = image-B8_Rinf_final;
   index = find(temp>0);

   temp2 = Ad-B8_Rinf_final;
   index2 = find(temp2>0);
   index = intersect(index,index2);

   %calculating lake depth
   z(index) = (log(Ad(index)-B8_Rinf_final)-log(temp(index)))/(g_B8);
   z_masked = mask.*z*1000; %can do this because the LakeMask is a 1/0 binary, written in mm
   z_masked(z_masked < 0) = 0;

    %write lake depths out
    z_masked(1,1) = 2^16-1;
    z_masked = uint16(z_masked);
    z_masked(1,1) = 0;
    
   B8_Depth=z_masked;
%    expression = strcat('geotiffwrite(''',Scenes{scene},'_B8_depth.tif'', B8_Depth, B3_info, ''GeoKeyDirectoryTag'',B3_tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);');
%    eval(expression)
    
   %%%%% calculating average depths from B4 and B8
   Depth_Average=(B4_Depth+B8_Depth)./2;
   Depth_Average = uint16(Depth_Average);

   expression = strcat('geotiffwrite(''',Scenes{scene},'_AverageB4B8_depth.tif'', Depth_Average, B3_info, ''GeoKeyDirectoryTag'',B3_tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);');
   eval(expression)
  
   %%%write union mask 
   All_Masks = zeros(size(B4));
   All_Masks(lake_mask==1)=1; % lakes
   All_Masks(rock_mask==1)=2; % rocks and seawater
   All_Masks(cloud_mask==1)=3; % clouds
   %write All_Masks
   expression = strcat('geotiffwrite(''',Scenes{scene},'_All_Masks.TIF'', All_Masks, B4_info, ''GeoKeyDirectoryTag'',B4_tiffinfo.GeoTIFFTags.GeoKeyDirectoryTag);');
   eval(expression)


   %%%% Calculating volume
   s = Scenes{scene};
   year = s(18:21);
   month=s(22:23);
   dayy = s(24:25);
   date_vol=s(18:25);

   date_number(scene)=datetime(str2double(year),str2double(month),str2double(dayy));  
   image=Depth_Average;
   vol = sum(image(image > 0)) * pixel * pixel / 1000; %depth in mm, 30 m pixels
   volume(scene) = vol;    
    
   %%%%% generating lake depth and area distribution

   Depths=double(image(image > 0)); 
   avg_Depth=mean(Depths(:))/1000; % in meters
   max_depth=ceil(max(Depths(:))/1000);% in meters 
   %%% 

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
   Area_1=area_pixel*900*1e-6; %% this calculates area in square kilometers, each landsat pixel is 30 m.
   idx_outliers_removed=(Area_1 > 0.0005 );%%% removing lakes that are less than 5 pixels in area.
   Area=Area_1(idx_outliers_removed);
   max_Area=ceil(max(Area(:)));
   
   %%% find percentage of lakes that are less than 1 square kilometer in
   %%% area
   id_1=find(Area < 1); 
   percentage_less_1=numel(id_1)*100/numel(Area); 
   text_less_1=strcat( num2str(fix(percentage_less_1)), '% of lakes < 1 km^2');
   text_total=strcat('Total Area= ', num2str( round(sum(Area),2)), ' km^2');
   cdf_area=cdfplot(Area); 
   hold on;
   x=1;
   y =(percentage_less_1/100); 
   Horz_line=line([0,fix(max(Area))],[y,y]);Horz_line.LineWidth=4; Horz_line.LineStyle=':';hold on;
   vert_line=line([x,x], [0,y] );vert_line.LineWidth=4; vert_line.LineStyle=':';
   ax=gca; ax.FontSize=22;ax.XLim=[0,max_Area];ax.LineWidth=2;cdf_area.LineWidth=7;ax.GridLineStyle=':';cdf_area.Color=[0.64 0.08 0.18] ;
   xlabel('Lake Area (km^2)', 'FontWeight', 'bold', 'FontSize', 25);ylabel('Cumulative Distribution','FontWeight', 'bold','FontSize', 25);
   t_less_1=text(0.5, 0.85, text_less_1);t_less_1.FontSize=22; t_less_1.FontWeight='Bold';
   t_total=text(0.5, 0.75, text_total);t_total.FontSize=22;t_total.FontWeight='Bold';
   exp_txt2=strcat(Scenes{scene}, '_Lake_Depth_Area_Distribution.png');
   saveas(gcf,exp_txt2)
   
   end
   
   %%%% Create Browse images of masks superimposed on RGB composite
   RGBimage = cat(3, imadjust(uint16(B4)),  imadjust(uint16(B3)),imadjust(uint16(B2)));
   RGB_resized=imresize(RGBimage, 0.2, 'nearest');
            
   B = imoverlay(RGB_resized,imresize(rock_mask, 0.2, 'nearest'),[0.12 0.32 0.8]); 
   C= imoverlay(B,imresize(cloud_mask, 0.2, 'nearest'), [0.4 .58 .5]);
   D=imoverlay(C,imresize(lake_mask, 0.2, 'nearest'), [0.8 0.2 0.1 ]);
           
   name_mask=strcat(Scenes{scene},'_masks_overlaid.png' );
   name_orig=strcat(Scenes{scene},'_original.png' );
   imwrite(((RGB_resized)), name_orig)
   imwrite((D), name_mask) 
   
   %change to appropriate directory
   expression = strcat('cd ', ' ', directory,Scenes{scene},'/');
   eval(expression);
 
end
  cd ..
  filename_output_volume='output_volume.mat';
  save (filename_output_volume,'Scenes',  'date_number','volume', 'cloud_percent');

  expression = strcat('cd ''',home,'''');
  eval(expression)
end

