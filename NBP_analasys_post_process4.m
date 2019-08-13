function NBP_analasys_post_process4
% This is to analyze 

indx = [523
524
525
526
529
534
538
542
543
545
548
553
560
569
571
574
580
581
582
584
585
591
595
597
599
602
606
609
620

];

for i = 1:length(indx)
    load_data2plotter(indx(i))
    PBFA_auto6
end


