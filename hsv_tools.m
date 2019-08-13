function hsv_tools(arg)

switch arg
    case 1
        x_y_z_data_curser
    case 2
        hsvp_gen
    case 3
        dist_along_frame
    otherwise
        disp('No match found')
end


function x_y_z_data_curser

    %fg=findall(gcf,'Type','figure')
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@curser3)
    
    
function hsvp_gen

    % Get hsv_analyser data
    h=guidata(findall(0,'Tag','hsv_analyzer'));
    b=h.b;
    
    fn = [b.dir b.file_name(1:end-4) 'txt'];
       
    
    if ~exist(fn,'file')
        % If the file not exist, just create it
        fID = fopen(fn , 'a+');
        fprintf(fID,'This file contains high speed video data positions (hsvp)\n\n');
        fprintf(fID,'cam_x:\t %0.2f\n',b.cam_x);
        fprintf(fID,'cam_y:\t %0.2f\n',b.cam_y);
        fprintf(fID,'cam_z:\t %0.2f\n',b.cam_z);
        fprintf(fID,'ref_x:\t %0.2f\n',b.ref_x);
        fprintf(fID,'ref_y:\t %0.2f\n',b.ref_y);
        fprintf(fID,'ref_z:\t %0.2f\n',b.ref_z);
        fprintf(fID,'scr_x:\t %i\n',b.scr_x);
        fprintf(fID,'scr_y:\t %i\n',b.scr_y);
        fprintf(fID,'align_time:\t %0.7f\n',b.align_time);
        fprintf(fID,'align_frame:\t %i\n',b.align_frame);
        fprintf(fID,'frame_rate:\t %i\n\n',h.ph.frame_rate);
        
        fprintf(fID,'*****************************************\n');
        fprintf(fID,'frameNu\t\tN_x\t\tN_y\t\tIntensity\n');
        fprintf(fID,'*****************************************\n\n');
    else
        % if the file exist, just open it for append
        fID = fopen(fn);
        % read the header lines
        hdr=textscan(fID,'%s %f',11,'HeaderLines',2);
        %clc
        %if hdr{2}(10)~= b.align_frame
        %    tt1 = hdr{2}(10)
        %    tt2 = b.align_frame
        %    disp('yes')
        %end
        %return
        if hdr{2}(1)~= b.cam_x || hdr{2}(2)~= b.cam_y || hdr{2}(3)~= b.cam_z || ...
            hdr{2}(4)~= b.ref_x || hdr{2}(5)~= b.ref_y || hdr{2}(6)~= b.ref_z || ...
            hdr{2}(7)~= b.scr_x || hdr{2}(8)~= b.scr_y || hdr{2}(9)~= round(b.align_time*1e7)/1e7 || ...
            hdr{2}(10)~= b.align_frame || hdr{2}(11)~= round(h.ph.frame_rate)
         
            opt = questdlg('Header values are different from the current values',...
                'Header Error','Ok','Modify header','Cancel','Ok');
            
            switch opt
                case 'Modify header'
                    errordlg('Sorry I am working on it. Do it manually for now.')
                otherwise
                    % do nothing
            end            
        end
           
        fclose(fID);
        fID = fopen(fn , 'a+');
    end
    
    % handles to the previos points
    hl = nan;
    
    xy = [];
    n = 0;
    % Loop, picking up the points.
    disp('Left mouse button picks points.')
    disp('Right mouse button picks last point.')
    but = 1;
    
    while but == 1
                
        [xi,yi,but] = ginput(1);
        
        
        hp1 = plot(xi,yi,'ro','MarkerSize',2,'MarkerFaceColor','r');
        n = n+1;
        xy(:,n) = [xi;yi];
        
        % Get the frame number
        h=guidata(findall(0,'Tag','hsv_analyzer'));
        b=h.b;
       
        % Get brightest pixel
        matlabIm = get_image(h,b.current_frame);
        matlabIm = matlabIm';
        
        midX = round(xi);
        midY = round(yi);
        maxI = -inf;
        maxXind = nan;
        maxYind = nan;
        
        for i1 = midX - 1:midX + 1
            for j1 = midY - 1:midY + 1;
                try 
                    tempMaxI = matlabIm(i1,j1);
                    if maxI < tempMaxI
                        maxI = tempMaxI;
                        maxXind = i1;
                        maxYind = j1;
                    end
                end
            end
        end
        
        hp2 = plot(maxXind,maxYind,'+');
        
        
        str = 'Add this to the file?';
        str = sprintf('%s\n\nOriginal:\n    X:%0.1f   \n    Y:%0.1f \n',str,xi,yi);
        [x,y,z] = scr2ldar(xi,yi,b);
        str = sprintf('%s\nx:%0.1fm   y:%0.1fm   z:%0.1fm',str,x,y,z);
        
        str = sprintf('%s\n\nModified:\n    X:%0.1f   \n    Y:%0.1f \n',str,maxXind,maxYind);
        [x,y,z] = scr2ldar(maxXind,maxYind,b);
        str = sprintf('%s\nx:%0.1fm   y:%0.1fm   z:%0.1fm',str,x,y,z);
        
      
        
        opt = questdlg(str,'HSVP file gen','Add Modified','Add the original','No','Add Modified');
        
        delete(hp1)
        delete(hp2)
        
        switch opt
            case 'Add Modified'
                fprintf(fID,'%i\t\t%0.1f\t\t%0.1f\t\t%i\n',b.current_frame,maxXind,maxYind,maxI);
                plot(maxXind,maxYind,'g+','MarkerSize',5)
                
                if isfield(h,'hsvpd')
                    [L1 ~]  = size(h.hsvpd);
                else
                    L1 = 0;
                end
                
                h.hsvpd(L1+1,1) = b.current_frame;
                h.hsvpd(L1+1,2) = maxXind;
                h.hsvpd(L1+1,3) = maxYind;
                h.hsvpd = sortrows(h.hsvpd,1);
                
                guidata(findall(0,'Tag','hsv_analyzer'), h)
            case 'Add the original'
                fprintf(fID,'%i\t\t%0.1f\t\t%0.1f\t\t%i\n',b.current_frame,xi,yi,matlabIm(midX,midY));
                plot(xi,yi,'g+','MarkerSize',5)
                
                if isfield(h,'hsvpd')
                    [L1 ~]  = size(h.hsvpd);
                else
                    L1 = 0;
                end
                
                h.hsvpd(L1+1,1) = b.current_frame;
                h.hsvpd(L1+1,2) = xi;
                h.hsvpd(L1+1,3) = yi;
                h.hsvpd = sortrows(h.hsvpd,1);
                
            otherwise                
                % Do nothing
        end
    end
    
    fclose(fID);
    
    
function [x y z] = scr2ldar(N_x,N_y,b)

f = 8.0e-3;  %Focal length
p = 20.0e-6; %Fixel Size


r = sqrt((((b.ref_x-b.cam_x)).^2+(b.ref_y-b.cam_y ).^2 ).*(1+(p^2.*(b.scr_x-N_x ).^2)./f^2 ));

theta = atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x ))+atan(p*(b.scr_x-N_x)/f);

x = b.cam_x + r.*cos(theta);
y = b.cam_y + r.*sin(theta);
z = b.cam_z + r.*((b.scr_y-N_y).*p)./f;

function dist_along_frame
    %try
        h=guidata(findall(0,'Tag','hsv_analyzer'));
        b=h.b;

        [Nx,Ny]=ginput(1);
        [x1 y1 z1] = scr2ldar(Nx,Ny,b);
        pause(0.25)
        [Nx,Ny]=ginput(1);
        [x2 y2 z2] = scr2ldar(Nx,Ny,b);
        dist = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
        msg = sprintf('Distance = %0.1f m',dist);
        msgbox(msg,'Distance')
%     catch E
%         msgbox('An error occured' ,'Ditance')
%         rethrow(E)
%     end
    
function matlabIm = get_image(handles,imageNu)

handles.ph.imgRange.First = imageNu;

%Read the cine image into the buffer
[~, unshiftedIm, imgHeader] = PhGetCineImage(handles.ph.cineHandle,...
                                                handles.ph.imgRange, ...
                                                handles.ph.imgSizeInBytes);
pImCount = libpointer('int32Ptr',1);
                                                                                      

% Transform 1D image pixels to 1D/3D image pixels to be used with MATLAB
[unshiftedIm] = ExtractImageMatrixFromImageBuffer(unshiftedIm, imgHeader);

bps = GetEffectiveBitsFromIH(imgHeader);
[matlabIm, ~] = ConstructMatlabImage(unshiftedIm, imgHeader.biWidth, imgHeader.biHeight, 1, bps);