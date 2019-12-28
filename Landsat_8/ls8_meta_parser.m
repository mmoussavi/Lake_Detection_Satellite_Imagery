function MTL_list=ls8_meta_parser(MTL_filename)
%ls8_meta_parser Metadata parser for Landsat 8 data
% Usage meta=ls8_meta_parser(finename)
% meta = variable to contain metadata
% filname = filename of external metadata file (read in)
%
% Neil Arnold, Scott Polar Research Institute, 2015
%
% Based on original version for LS 7 from Mathworks file exchange Mar 2011 by Seongsu Jeong
% Geomatics and Remote Sensing Laboratory (GRSLAB), Yonsei University
%

if nargin()==1    
    nMTL=size(MTL_filename,1);    
    
else
    error('Incorrect input : Enter metadata file name');
    
end

for cnt=1:nMTL
    % File open and line string input
    fin=fopen(MTL_filename(cnt,:),'r');
    str_in=fgetl(fin);

    while ~strcmp(str_in,'END')
        % String parsing and refinement
        % input line refinement and processing
        
        [field, value]=strtok(str_in,'=');
        
        %field name refinement
        %remove unnecessary space character from the field description
        if sum(field(1)==' ') %space character detector
            field=strtok(field,' ');
        end
        
        %Value refinement
        value=strtok(value,'=');
        
        %remove unnecessary space character from the field description
        while value(1,1)==' '
            value=value(1,2:size(value,2));  
        end

        %disp(strcat(field,': ',value));
        
        if sum(value(1)=='"' && value(size(value,2))=='"') %If the value is a string wrapped by large quotation mark (")
            value=value(1,2:size(value,2)-1);
        elseif isempty(findstr(field,'TIME')) && isempty(findstr(field,'DATE'))
            value=str2num(value);
        end
        
        % Field detection and assignment routine
        switch field
            
            % GROUP=METADATA_FILE_INFO
            case 'ORIGIN'
                MTL.METADATA_FILE_INFO.ORIGIN=value;                
            case 'REQUEST_ID'
                MTL.METADATA_FILE_INFO.REQUEST_ID=value;
            case 'LANDSAT_SCENE_ID'
                MTL.METADATA_FILE_INFO.LANDSAT_SCENE_ID=value;
            case 'STATION_ID'
                MTL.METADATA_FILE_INFO.STATION_ID=value;
            case 'FILE_DATE'
                MTL.METADATA_FILE_INFO.FILE_DATE=value;   
            case 'PROCESSING_SOFTWARE_VERSION'
                MTL.METADATA_FILE_INFO.PROCESSING_SOFTWARE_VERSION=value;
 
            % GROUP=PRODUCT_METADATA
            case 'DATA_TYPE'        
                MTL.PRODUCT_METADATA.DATA_TYPE=value;
            case 'ELEVATION_SOURCE'        
                MTL.PRODUCT_METADATA.ELEVATION_SOURCE=value;
            case 'OUTPUT_FORMAT'        
                MTL.PRODUCT_METADATA.OUTPUT_FORMAT=value;
            case 'SPACECRAFT_ID'        
                MTL.PRODUCT_METADATA.SPACECRAFT_ID=value;
            case 'SENSOR_ID'        
                MTL.PRODUCT_METADATA.SENSOR_ID=value;
            case 'WRS_PATH'        
                MTL.PRODUCT_METADATA.WRS_PATH=value;
            case 'WRS_ROW'        
                MTL.PRODUCT_METADATA.WRS_ROW=value;
            case 'NADIR_OFFNADIR'        
                MTL.PRODUCT_METADATA.NADIR_OFFNADIR=value;
            case 'TARGET_WRS_PATH'        
                MTL.PRODUCT_METADATA.TARGET_WRS_PATH=value;
            case 'TARGET_WRS_ROW'        
                MTL.PRODUCT_METADATA.TARGET_WRS_ROW=value;
            case 'DATE_ACQUIRED'        
                MTL.PRODUCT_METADATA.DATE_ACQUIRED=value;
            case 'SCENE_CENTER_TIME'        
                MTL.PRODUCT_METADATA.SCENE_CENTER_TIME=value;
                %
            case 'CORNER_UL_LAT_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_UL_LAT_PRODUCT=value;
            case 'CORNER_UL_LON_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_UR_LON_PRODUCT=value;
            case 'CORNER_UR_LAT_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_UR_LAT_PRODUCT=value;
            case 'CORNER_UR_LON_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_UR_LON_PRODUCT=value;
            case 'CORNER_LL_LAT_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_LL_LAT_PRODUCT=value;
            case 'CORNER_LL_LON_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_LL_LON_PRODUCT=value;
            case 'CORNER_LR_LAT_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_LR_LAT_PRODUCT=value;
            case 'CORNER_LR_LON_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_LR_LON_PRODUCT=value;
                 %
            case 'CORNER_UL_PROJECTION_X_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_UL_PROJECTION_X_PRODUCT=value;
            case 'CORNER_UL_PROJECTION_Y_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_UR_PROJECTION_Y_PRODUCT=value;
            case 'CORNER_UR_PROJECTION_X_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_UR_PROJECTION_X_PRODUCT=value;
            case 'CORNER_UR_PROJECTION_Y_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_UR_PROJECTION_Y_PRODUCT=value;
            case 'CORNER_LL_PROJECTION_X_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_LL_PROJECTION_X_PRODUCT=value;
            case 'CORNER_LL_PROJECTION_Y_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_LL_PROJECTION_Y_PRODUCT=value;
            case 'CORNER_LR_PROJECTION_X_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_LR_PROJECTION_X_PRODUCT=value;
            case 'CORNER_LR_PROJECTION_Y_PRODUCT'        
                 MTL.PRODUCT_METADATA.CORNER_LR_PROJECTION_Y_PRODUCT=value;
                 %
            case 'PANCHROMATIC_LINES'        
                MTL.PRODUCT_METADATA.PANCHROMATIC_LINES=value;     
            case 'PANCHROMATIC_SAMPLES'        
                MTL.PRODUCT_METADATA.PANCHROMATIC_SAMPLES=value;     
            case 'REFLECTIVE_LINES'        
                MTL.PRODUCT_METADATA.REFLECTIVE_LINES=value;     
            case 'REFLECTIVE_SAMPLES'        
                MTL.PRODUCT_METADATA.REFLECTIVE_SAMPLES=value;     
            case 'THERMAL_LINES'        
                MTL.PRODUCT_METADATA.THERMAL_LINES=value;     
            case 'THERMAL_SAMPLES'        
                MTL.PRODUCT_METADATA.THERMAL_SAMPLES=value;
                %
            case 'FILE_NAME_BAND1'        
                MTL.PRODUCT_METADATA.FILE_NAME_BAND1=value;     
            case 'FILE_NAME_BAND2'        
                MTL.PRODUCT_METADATA.FILE_NAME_BAND2=value;     
            case 'FILE_NAME_BAND3'        
                MTL.PRODUCT_METADATA.FILE_NAME_BAND3=value;     
            case 'FILE_NAME_BAND4'        
                MTL.PRODUCT_METADATA.FILE_NAME_BAND4=value;     
            case 'FILE_NAME_BAND5'        
                MTL.PRODUCT_METADATA.FILE_NAME_BAND5=value;     
            case 'FILE_NAME_BAND6'        
                MTL.PRODUCT_METADATA.FILE_NAME_BAND6=value;
            case 'FILE_NAME_BAND7'        
                MTL.PRODUCT_METADATA.FILE_NAME_BAND7=value;
            case 'FILE_NAME_BAND8'
                MTL.PRODUCT_METADATA.FILE_NAME_BAND8=value;
            case 'FILE_NAME_BAND9'
                MTL.PRODUCT_METADATA.FILE_NAME_BAND9=value;
            case 'FILE_NAME_BAND10'
                MTL.PRODUCT_METADATA.FILE_NAME_BAND10=value;
            case 'FILE_NAME_BAND11'
                MTL.PRODUCT_METADATA.FILE_NAME_BAND11=value;
                %
            case 'METADATA_FILE_NAME'
                MTL.PRODUCT_METADATA.METADATA_FILE_NAME=value;
            case 'BPF_NAME_OLI'
                MTL.PRODUCT_METADATA.BPF_NAME_OLI=value;
            case 'BPF_NAME_TIRS'
                MTL.PRODUCT_METADATA.BPF_NAME_TIRS=value;
            case 'CPF_FILE_NAME'
                MTL.PRODUCT_METADATA.CPF_FILE_NAME=value;
            case 'RLUT_FILE_NAME'
                MTL.PRODUCT_METADATA.CPF_FILE_NAME=value;
                
            %GROUP = IMAGE_ATTRIBUTES
            case 'CLOUD_COVER'
                MTL.IMAGE_ATTRIBUTES.CLOUD_COVER=value;
            case 'IMAGE_QUALITY_OLI'
                MTL.IMAGE_ATTRIBUTES.IMAGE_QUALITY_OLI=value;
            case 'IMAGE_QUALITY_TIRS'
                MTL.IMAGE_ATTRIBUTES.IMAGE_QUALITY_TIRS=value;
            case 'ROLL_ANGLE'
                MTL.IMAGE_ATTRIBUTES.ROLL_ANGLE=value;
            case 'SUN_AZIMUTH'
                MTL.IMAGE_ATTRIBUTES.SUN_AZIMUTH=value;
            case 'SUN_ELEVATION'
                MTL.IMAGE_ATTRIBUTES.SUN_ELEVATION=value;
            case 'EARTH_SUN_DISTANCE'
                MTL.IMAGE_ATTRIBUTES.EARTH_SUN_DISTANCE=value;
            case 'GROUND_CONTROL_POINTS_MODEL'
                MTL.IMAGE_ATTRIBUTES.GROUND_CONTROL_POINTS_MODEL=value;
            case 'GEOMETRIC_RMSE_MODEL'
                MTL.IMAGE_ATTRIBUTES.GEOMETRIC_RMSE_MODEL=value;
            case 'GEOMETRIC_RMSE_MODEL_Y'
                MTL.IMAGE_ATTRIBUTES.GEOMETRIC_RMSE_MODEL_Y=value;
            case 'GEOMETRIC_RMSE_MODEL_X'
                MTL.IMAGE_ATTRIBUTES.GEOMETRIC_RMSE_MODEL_X=value;
            case 'GROUND_CONTROL_POINTS_VERIFY'
                MTL.IMAGE_ATTRIBUTES.GROUND_CONTROL_POINTS_VERIFY=value;
            case 'GEOMETRIC_RMSE_VERIFY'
                MTL.IMAGE_ATTRIBUTES.GEOMETRIC_RMSE_VERIFY_X=value;
   
            %GROUP = MIN_MAX_RADIANCE
            case 'RADIANCE_MAXIMUM_BAND_1'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_1=value;        
            case 'RADIANCE_MINIMUM_BAND_1'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_1=value;        
            case 'RADIANCE_MAXIMUM_BAND_2'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_2=value;        
            case 'RADIANCE_MINIMUM_BAND_2'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_2=value;        
            case 'RADIANCE_MAXIMUM_BAND_3'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_3=value;        
            case 'RADIANCE_MINIMUM_BAND_3'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_3=value;        
            case 'RADIANCE_MAXIMUM_BAND_4'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_4=value;        
            case 'RADIANCE_MINIMUM_BAND_4'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_4=value;        
            case 'RADIANCE_MAXIMUM_BAND_5'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_5=value;        
            case 'RADIANCE_MINIMUM_BAND_5'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_5=value;        
            case 'RADIANCE_MAXIMUM_BAND_6'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_6=value;        
            case 'RADIANCE_MINIMUM_BAND_6'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_6=value;        
            case 'RADIANCE_MAXIMUM_BAND_7'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_7=value;        
            case 'RADIANCE_MINIMUM_BAND_7'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_7=value;        
            case 'RADIANCE_MAXIMUM_BAND_8'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_8=value;        
            case 'RADIANCE_MINIMUM_BAND_8'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_8=value;        
            case 'RADIANCE_MAXIMUM_BAND_9'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_9=value;        
            case 'RADIANCE_MINIMUM_BAND_9'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_9=value;        
            case 'RADIANCE_MAXIMUM_BAND_10'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_10=value;        
            case 'RADIANCE_MINIMUM_BAND_10'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_10=value;        
            case 'RADIANCE_MAXIMUM_BAND_11'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MAXIMUM_BAND_11=value;        
            case 'RADIANCE_MINIMUM_BAND_11'        
                MTL.MIN_MAX_RADIANCE.RADIANCE_MINIMUM_BAND_11=value;         
                
            %GROUP = MIN_MAX_REFLECTANCE
            case 'REFLECTANCE_MAXIMUM_BAND_1'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_1=value;        
            case 'REFLECTANCE_MINIMUM_BAND_1'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_1=value;        
            case 'REFLECTANCE_MAXIMUM_BAND_2'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_2=value;        
            case 'REFLECTANCE_MINIMUM_BAND_2'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_2=value;        
            case 'REFLECTANCE_MAXIMUM_BAND_3'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_3=value;        
            case 'REFLECTANCE_MINIMUM_BAND_3'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_3=value;        
            case 'REFLECTANCE_MAXIMUM_BAND_4'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_4=value;        
            case 'REFLECTANCE_MINIMUM_BAND_4'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_4=value;        
            case 'REFLECTANCE_MAXIMUM_BAND_5'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_5=value;        
            case 'REFLECTANCE_MINIMUM_BAND_5'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_5=value;        
            case 'REFLECTANCE_MAXIMUM_BAND_6'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_6=value;        
            case 'REFLECTANCE_MINIMUM_BAND_6'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_6=value;        
            case 'REFLECTANCE_MAXIMUM_BAND_7'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_7=value;        
            case 'REFLECTANCE_MINIMUM_BAND_7'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_7=value;        
            case 'REFLECTANCE_MAXIMUM_BAND_8'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_8=value;        
            case 'REFLECTANCE_MINIMUM_BAND_8'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_8=value;        
            case 'REFLECTANCE_MAXIMUM_BAND_9'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_9=value;        
            case 'REFLECTANCE_MINIMUM_BAND_9'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_9=value;        
            case 'REFLECTANCE_MAXIMUM_BAND_10'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_10=value;        
            case 'REFLECTANCE_MINIMUM_BAND_10'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_10=value;        
            case 'REFLECTANCE_MAXIMUM_BAND_11'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MAXIMUM_BAND_11=value;        
            case 'REFLECTANCE_MINIMUM_BAND_11'        
                MTL.MIN_MAX_REFLECTANCE.REFLECTANCE_MINIMUM_BAND_11=value;        

            %GROUP = MIN_MAX_PIXEL_VALUE
            case 'QUANTIZE_CAL_MAX_BAND_1'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_1=value;        
            case 'QUANTIZE_CAL_MIN_BAND_1'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_1=value;        
            case 'QUANTIZE_CAL_MAX_BAND_2'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_2=value;        
            case 'QUANTIZE_CAL_MIN_BAND_2'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_2=value;        
            case 'QUANTIZE_CAL_MAX_BAND_3'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_3=value;        
            case 'QUANTIZE_CAL_MIN_BAND_3'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_3=value;        
            case 'QUANTIZE_CAL_MAX_BAND_4'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_4=value;        
            case 'QUANTIZE_CAL_MIN_BAND_4'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_4=value;        
            case 'QUANTIZE_CAL_MAX_BAND_5'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_5=value;        
            case 'QUANTIZE_CAL_MIN_BAND_5'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_5=value;        
            case 'QUANTIZE_CAL_MAX_BAND_6'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_6=value;        
            case 'QUANTIZE_CAL_MIN_BAND_6'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_6=value;        
            case 'QUANTIZE_CAL_MAX_BAND_7'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_7=value;        
            case 'QUANTIZE_CAL_MIN_BAND_7'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_7=value;        
            case 'QUANTIZE_CAL_MAX_BAND_8'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_8=value;        
            case 'QUANTIZE_CAL_MIN_BAND_8'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_8=value;        
            case 'QUANTIZE_CAL_MAX_BAND_9'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_9=value;        
            case 'QUANTIZE_CAL_MIN_BAND_9'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_9=value;        
            case 'QUANTIZE_CAL_MAX_BAND_10'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_10=value;        
            case 'QUANTIZE_CAL_MIN_BAND_10'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_10=value;        
            case 'QUANTIZE_CAL_MAX_BAND_11'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MAX_BAND_11=value;        
            case 'QUANTIZE_CAL_MIN_BAND_11'        
                MTL.MIN_MAX_PIXEL_VALUE.QUANTIZE_CAL_MIN_BAND_11=value;        

                %GROUP = RADIOMETRIC_RESCALING        
            case 'RADIANCE_MULT_BAND_1'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_1=value;
            case 'RADIANCE_ADD_BAND_1'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_1=value;
            case 'RADIANCE_MULT_BAND_2'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_2=value;
            case 'RADIANCE_ADD_BAND_2'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_2=value;
            case 'RADIANCE_MULT_BAND_3'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_3=value;
            case 'RADIANCE_ADD_BAND_3'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_3=value;
            case 'RADIANCE_MULT_BAND_4'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_4=value;
            case 'RADIANCE_ADD_BAND_4'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_4=value;
            case 'RADIANCE_MULT_BAND_5'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_5=value;
            case 'RADIANCE_ADD_BAND_5'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_5=value;
            case 'RADIANCE_MULT_BAND_6'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_6=value;
            case 'RADIANCE_ADD_BAND_6'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_6=value;
            case 'RADIANCE_MULT_BAND_7'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_7=value;
            case 'RADIANCE_ADD_BAND_7'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_7=value;
            case 'RADIANCE_MULT_BAND_8'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_8=value;
            case 'RADIANCE_ADD_BAND_8'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_8=value;
            case 'RADIANCE_MULT_BAND_9'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_9=value;
            case 'RADIANCE_ADD_BAND_9'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_9=value;
            case 'RADIANCE_MULT_BAND_10'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_10=value;
            case 'RADIANCE_ADD_BAND_10'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_10=value;
            case 'RADIANCE_MULT_BAND_11'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_11=value;
            case 'RADIANCE_ADD_BAND_11'
                MTL.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_11=value;
            case 'REFLECTANCE_MULT_BAND_1'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_1=value;
            case 'REFLECTANCE_ADD_BAND_1'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_1=value;
            case 'REFLECTANCE_MULT_BAND_2'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_2=value;
            case 'REFLECTANCE_ADD_BAND_2'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_2=value;
            case 'REFLECTANCE_MULT_BAND_3'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_3=value;
            case 'REFLECTANCE_ADD_BAND_3'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_3=value;
            case 'REFLECTANCE_MULT_BAND_4'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_4=value;
            case 'REFLECTANCE_ADD_BAND_4'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_4=value;
            case 'REFLECTANCE_MULT_BAND_5'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_5=value;
            case 'REFLECTANCE_ADD_BAND_5'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_5=value;
            case 'REFLECTANCE_MULT_BAND_6'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_6=value;
            case 'REFLECTANCE_ADD_BAND_6'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_6=value;
            case 'REFLECTANCE_MULT_BAND_7'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_7=value;
            case 'REFLECTANCE_ADD_BAND_7'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_7=value;
            case 'REFLECTANCE_MULT_BAND_8'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_8=value;
            case 'REFLECTANCE_ADD_BAND_8'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_8=value;
            case 'REFLECTANCE_MULT_BAND_9'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_9=value;
            case 'REFLECTANCE_ADD_BAND_9'
                MTL.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_9=value;
                
                %GROUP = TIRS_THERMAL_CONSTANTS
            case 'K1_CONSTANT_BAND_10'
                MTL.TIRS_THERMAL_CONSTANTS.K1_CONSTANT_BAND_10=value;
            case 'K2_CONSTANT_BAND_10'
                MTL.TIRS_THERMAL_CONSTANTS.K2_CONSTANT_BAND_10=value;
            case 'K1_CONSTANT_BAND_11'
                MTL.TIRS_THERMAL_CONSTANTS.K1_CONSTANT_BAND_11=value;
            case 'K2_CONSTANT_BAND_11'
                MTL.TIRS_THERMAL_CONSTANTS.K2_CONSTANT_BAND_11=value;

            %GROUP = PROJECTION_PARAMETERS
            case 'MAP_PROJECTION'
                MTL.PROJECTION_PARAMETERS.MAP_PROJECTION=value;
            case 'DATUM'
                MTL.PROJECTION_PARAMETERS.DATUM=value;
            case 'ELLIPSOID'
                MTL.PROJECTION_PARAMETERS.ELLIPSOID=value;
            case 'UTM_ZONE'
                MTL.PROJECTION_PARAMETERS.UTM_ZONE=value;                
            case 'GRID_CELL_SIZE_PANCHROMATIC'
                MTL.PROJECTION_PARAMETERS.GRID_CELL_SIZE_PANCHROMATIC=value;
            case 'GRID_CELL_SIZE_REFLECTIVE'
                MTL.PROJECTION_PARAMETERS.GRID_CELL_SIZE_REFLECTIVE=value;
            case 'GRID_CELL_SIZE_THERMAL'
                MTL.PROJECTION_PARAMETERS.GRID_CELL_SIZE_THERMAL=value;
            case 'ORIENTATION'
                MTL.PROJECTION_PARAMETERS.ORIENTATION=value;
            case 'RESAMPLING_OPTION'
                MTL.PROJECTION_PARAMETERS.RESAMPLING_OPTION=value;
                    
        end
        
        str_in=fgetl(fin);
    end
    
    MTL_list(cnt,1)=MTL;
    fclose(fin);
    
end