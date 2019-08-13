function eps2pdf

%dn = 'E:\Sumedhe\Documents\PBFA_paper\LocatingIBP4-JGR - REV\';
dn = 'C:\Users\Sumedhe\Desktop\IBP_Modeling_JGR\Version2\pics\';
%dn = uigetdir;

fn = dir([dn '*.eps']);
L = length(fn);


wbh = waitbar(0,'Converting');

for i = 3:L
   
   msg = sprintf('Converting ... %0.1f%%',(i-2)/(L-2)*100);
   
   try
        waitbar((i-2)/(L-2),wbh,msg)
   catch
       return
   end       
   
   name =  fn(i).name;
   cmnd = sprintf('epstopdf %s%s',dn,name)
   dos(cmnd);
end

delete(wbh)