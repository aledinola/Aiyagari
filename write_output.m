function [] = write_output(ED,agg,mom,par)
% Write output (parameters, equilib prices and model moments) on the screen
% and (to be yet implemented) on a file.

if ~isstruct(par)
    error('Input argument par must be a struct')
end
if ~isstruct(agg)
    error('Input argument agg must be a struct')
end
if ~isstruct(mom)
    error('Input argument mom must be a struct')
end
% if ~isstruct(mom)
%     error('mom must be a struct')
% end

fprintf("************************************************ \n")
    fprintf(" \n")
%     fprintf("MODEL PARAMETERS \n")
%     % ECONOMIC PARAMETERS
%     fprintf("beta       = %f \n",par.beta)
%     fprintf("delta      = %f \n",par.delta)
%     fprintf("alpha      = %f \n",par.alpha)
%     fprintf("vi         = %f \n",par.vi)
%     fprintf("lambda     = %f \n",par.lambda)
%     fprintf("gamma      = %f \n",par.gamma)
%     fprintf("le         = %f \n",par.le)
%     fprintf("th_bar     = %f \n",shock.th_bar)
%     fprintf("th_std     = %f \n",shock.th_std)
%     fprintf("p11        = %f \n",shock.p11)
%     fprintf("p21        = %f \n",shock.p21)
%     fprintf("p22        = %f \n",shock.p22)
%     fprintf("p44        = %f \n",shock.p44)
%     fprintf("e6         = %f \n",shock.e6)
%     fprintf("pi6        = %f \n",shock.pi6)
%     fprintf("pi63       = %f \n",shock.pi63)
%     fprintf("lambda_hsv = %f \n",par.lambda_hsv)
%     fprintf("tau_hsv    = %f \n",par.tau_hsv)
%     % COMPUTATIONAL PARAMETERS
%     fprintf("Size grid assets = %d \n",par.na)
%     fprintf("Size grid theta  = %d \n",par.nz)
%     fprintf("a_max            = %f \n",par.amax)
%     fprintf("a_space          = %f \n",par.curv)
%     fprintf("Golden method    = %d \n",par.do_golden)
%     fprintf("St.dev. taste shock = %f \n",par.sig_e)
%     fprintf("do_GE      = %d \n",par.do_GE)
%     fprintf("do_howard  = %d \n",par.do_howard)
%     fprintf("n_howard   = %d \n",par.n_howard)
%     fprintf("monotoni   = %d \n",par.monotonicity)
%     fprintf("do_super   = %d \n",par.do_super)
%     fprintf(" \n")
    
    fprintf("MODEL MOMENTS \n")
    fprintf("Excess demand       = %f \n",ED)
    fprintf("Residual mkt clear  = %f \n",agg.res)
    fprintf("Residual1 mkt clear = %f \n",agg.res1)
    %fprintf("Share of workers    = %f \n",mom.share_w)
    %fprintf("Share of entrepr    = %f \n",mom.share_e)
    fprintf("Interest rate       = %f \n",par.r)
    fprintf("Wage                = %f \n",par.wage)
    fprintf("K/Y                 = %f \n",agg.Ks/agg.YY)
    fprintf("C/Y                 = %f \n",agg.CC/agg.YY)
    fprintf("I/Y                 = %f \n",agg.II/agg.YY)
    fprintf("K/L                 = %f \n",agg.Ks/agg.Ls)
    %fprintf("share_entre_inc     = %f \n",mom.share_entre_inc)
    fprintf("YY (total output)   = %f \n",agg.YY)
    fprintf("KK (total capital)  = %f \n",agg.Ks)
    fprintf("CC (total consump)  = %f \n",agg.CC)
    
    fprintf('CV labor product    = %f \n',mom.cv_z)
    fprintf('CV consumption      = %f \n',mom.cv_c)
    fprintf('CV wealth           = %f \n',mom.cv_a)

    %fprintf("Tot. Taxes          = %f \n",mom.TT)
    %fprintf("Tot. Taxes/YY       = %f \n",mom.TT/mom.YY)
    %fprintf("Gini wealth         = %f \n",mom.gini_a)
    %fprintf("Gini income         = %f \n",mom.gini_inc)
    %fprintf("Gini consump.       = %f \n",mom.gini_c)
    
%     fprintf("share_inc_top1      = %f \n",mom.share_inc_top1)
%     fprintf("Med wealth E-to-W   = %f \n",mom.medwealth_entre_work)
%     fprintf("Share entre hiring  = %f \n",mom.frac_hiring)
%     
%     fprintf("%% entre inc 0-p20    = %f \n",mom.share_entre_incbins(1))
%     fprintf("%% entre inc p20-p40  = %f \n",mom.share_entre_incbins(2))
%     fprintf("%% entre inc p40-p60  = %f \n",mom.share_entre_incbins(3))
%     fprintf("%% entre inc p60-p80  = %f \n",mom.share_entre_incbins(4))
%     fprintf("%% entre inc p80-p90  = %f \n",mom.share_entre_incbins(5))
%     fprintf("%% entre inc p90-p99  = %f \n",mom.share_entre_incbins(6))
%     fprintf("%% entre inc >p99     = %f \n",mom.share_entre_incbins(7))
%     
%     fprintf("%% inc 0-p20    = %f \n",mom.mass_inc_bins(1))
%     fprintf("%% inc p20-p40  = %f \n",mom.mass_inc_bins(2))
%     fprintf("%% inc p40-p60  = %f \n",mom.mass_inc_bins(3))
%     fprintf("%% inc p60-p80  = %f \n",mom.mass_inc_bins(4))
%     fprintf("%% inc p80-p90  = %f \n",mom.mass_inc_bins(5))
%     fprintf("%% inc p90-p99  = %f \n",mom.mass_inc_bins(6))
%     fprintf("%% inc >p99     = %f \n",mom.mass_inc_bins(7))
%     
%     fprintf("Time solution (sec) = %f \n",time)


end %end function

