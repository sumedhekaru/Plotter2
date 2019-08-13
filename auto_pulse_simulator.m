function auto_pulse_simulator

a=struct( ...
             'maxA' , {500000}       , ...
               't1' , {4.0000e-006} , ...
               't2' , {2.1500e-005} , ...
                'v' , {0.1e8}   , ...
               'H1' , {5146.2}        , ...
               'H2' , {6000}        , ...
                'x' , {-9582.6}, ...
                'y' , {6013.7} , ...
               'T1' , {0}           , ...
               'T2' , {2.5000e-004} , ...
            'lamda' , {200}         , ...
           't_step' , {1.0000e-006} , ...
                'N' , {1}           , ...
        't_reflect' , {5.0000e-006} , ...
    'coeff_reflect' , {-0.5000}     , ...
          't_shift' , {0} , ...
         'is_saved' , {0}           , ...
          'save_fn' , {''}          , ...
        'file_name' , {'C:\Users\sumedhe\Desktop\Plotter2\pulse_parameter_test.mat'}, ...
         'plots_on' , {[0 0 0 0]}   , ...
       'addi_plots' , {[0 0 0 0]}   , ...
            's_num' , {1}             ...
        );
    
    % Open plotter 2 data    
    try; h=guidata(findall(0,'Tag','plotter2'));
    catch; disp('Run plotter2 first!'); return;
    end
    
    settings = h.sen_set;
    g = h.g;
    g.chgraphs = zeros(1,60);
    g.lpgraphs = zeros(1,60);
        
    % Detection parameters (real)
    para.sn = 1;
    para.ppV = 12.97;
    para.ppt = 77636.33377517;
    para.npV = -27.73;
    para.npt = 77636.33376917;
    para.begt = 77636.33376517;
    para.endt = 77636.33378617;
    
    g.chgraphs(3) = 1;
    plot_all4(g)
    hold all
    
    
    data = pulse_simulator(a);
    ft = find_features(data);
    data.t = data.t + para.begt - ft.begt;
    ft = find_features(data);
    
    
    
    
    
    
%     figure
%     hold all
    plot(data.t,data.E_tot)
%     plot(data.t,data.E_rad)
%     plot(data.t,data.E_ind)
%     plot(data.t,data.E_stat)

    plot([ft.ppt ft.npt],[ft.ppV ft.npV],'ro','markerfacecolor','r')
    plot([ft.begt ft.endt],[0 0],'go','markerfacecolor','g')
    
    
    
    
    
    
function ft = find_features(data)

% Positive peak
[ft.ppV ind] = max(data.E_tot);
ft.ppt       = data.t(ind);
% Negative peak
[ft.npV ind] = min(data.E_tot);
ft.npt       = data.t(ind);

% Begining time
ind          = find(data.E_tot<ft.npV*.0001);
ft.begt      = data.t(ind(1));

% End time
ind         = find(data.E_stat > min(data.E_stat)*.999);
ft.endt     = data.t(ind(end));
