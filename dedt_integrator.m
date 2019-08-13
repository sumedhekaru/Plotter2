function dedt_integrator

h=guidata(findall(0,'Tag','plotter2'));
g = h.g;
settings = h.sen_set;

% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

% finding the file extention number
ext = 3;
i = 57;

% Finding the stattion ID
sid=settings.sen_IDs{ceil(i/3)};

% finding the file name
filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
    bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);

hfilename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
    bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);


[t y ch_freq_str] = FA_Extract1(filename,hfilename,g.t1,g.t2,0,settings,i);

%[t,y,filtered] = ch_high_pass(t,y,3000);
y = y + 3.23799;

figure
hold all
plot(t,y)

inds = find(isnan(y));

y([1:2]) = [];
t([1:2]) = [];

intY = cumtrapz(y);

%[t,intY,filtered] = ch_high_pass(t,intY,3000);

plot(t,-intY)

tools2fig
xlim([g.t1 g.t2])


ext = 1;

% finding the file name
filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
    bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);

hfilename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
    bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);


[t y ch_freq_str] = FA_Extract1(filename,hfilename,g.t1,g.t2,0,settings,i);
plot(t,y*10)