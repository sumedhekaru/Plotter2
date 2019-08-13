function data = hsvp_extract(fn)

fID = fopen(fn,'r');
% Header
hdr=textscan(fID,'%s %f',11,'HeaderLines',2);

% data
tmp = textscan(fID,'%f %f %f %f','HeaderLines',5);
fclose(fID);

tmp = cell2mat(tmp);
tmp = sortrows(tmp,1);

data.frameN = tmp(:,1);
data.N_x    = tmp(:,2);
data.N_y    = tmp(:,3);
data.I      = tmp(:,4);

data.cam_x  = hdr{2}(1);
data.cam_y  = hdr{2}(2);
data.cam_z  = hdr{2}(3);
data.ref_x  = hdr{2}(4);
data.ref_y  = hdr{2}(5);
data.ref_z  = hdr{2}(6);
data.scr_x  = hdr{2}(7);
data.scr_y  = hdr{2}(8);
data.align_time  = hdr{2}(9);
data.align_frame = hdr{2}(10);
data.frame_rate  = hdr{2}(11);
