$title SDP-DICE Model
* Fully replicable to DICE2016R
* Used for Chang and Yoo (2022) "How Learning Helps Mitigate the Worst When the Downside of Climate Change is Extreme"
* Uncertainty lies in Equilibrium Climate Sensitivity
* Assumes Learning

$offlisting
$eolcom #

* set bs and threads to number of cores
$if not set bs $set bs 1
$if not set threads $set threads 1
option threadsAsync=%threads%;

*       Number of VFI iterations
$if not set niter       $set niter        106

*       Set grid points uniquely for each state variable
$if not set gridn_k     $set gridn_k 6
$if not set gridn_m     $set gridn_m 4
$if not set gridn_p     $set gridn_p 5
$if not set gridn_u     $set gridn_u 4
$if not set gridn_l     $set gridn_l 2
$if not set gridn_ot    $set gridn_ot 3
$if not set gridn_bm    $set gridn_bm 6
$if not set gridn_bs    $set gridn_bs 2

$eval gridn (max(%gridn_k%,%gridn_m%,%gridn_u%,%gridn_l%,%gridn_ot%,%gridn_p%,%gridn_bm%,%gridn_bs%))

* Parameter definition for DICE 2016
$eval t_last (%niter%)
$eval t_last_1 (%t_last%-1)
$eval t_last_5 (%t_last%-5)

*        Set Nodes and Smoyak sparse grid
sets    t               Time horizon                    /1*%t_last%/
        t_first(t)      Initial time period             /1/
        t_last(t)       Final time period               /%t_last%/
        t_exc_last(t)   Time excluding the last         /1*%t_last_1%/,
        t_exc_last5(t)   Time excluding the last        /1*%t_last_5%/,
        ct(t)           Current time period;

ct(t) = yes$t_last(t);

* Model parameters

Parameters
*       Preferences
        elasmu   Elasticity of marginal utility of consumption   /1.45 /
        prstp    Initial rate of social time preference per year /.015 /
        dk       Depreciation rate on capital (per year)         /.100 /
        gama     Capital elasticity in production function       /.300 /
        disc     discount factor in DICE 2016                              ,
*       Climate sensitivity parameters
        t2xco2   Equilibrium temp impact (oC per doubling CO2)   /3.1 /
        lambda   Reference system's climate sensitivity (C-Wm2)  /1.2 /
        var_e    Variance of temperature shock                   /0.011/
*       Roe-Baker Distribution for the feedback factor N(0.6,0.1^2)
        lev_bs_init Variance of the PDF                           /0.01/
        lev_s_init Meadian of the PDF                             /3/
        true_sens True Climate Sensitivity
*       Climate damage parameters
        a2        Damage quadratic term                          /0.010038/
        a3        Damage exponent                                /2.00/
        expcost2  Exponent of control cost function              /2.6/
*       Climate model parameters
        fco22x   Forcings of equilibrium CO2 doubling (Wm-2)      /3.6813 /
        c1       Climate equation coefficient for upper level     /0.1005 /
        c3       Transfer coefficient upper to lower stratum      /0.088  /
        c4       Transfer coefficient for lower level             /0.025  /
*       Flow paramaters
        b12      Carbon cycle transition matrix                   /.12   /
        b23      Carbon cycle transition matrix                   /0.007 /
        b11      Carbon cycle transition matrix
        b21      Carbon cycle transition matrix
        b22      Carbon cycle transition matrix
        b32      Carbon cycle transition matrix
        b33      Carbon cycle transition matrix
        mate     Equilibrium concentration atmosphere  (GtC)      /588  /
        mue      Equilibrium concentration in upper strata (GtC)  /360  /
        mle      Equilibrium concentration in lower strata (GtC)  /1720 /
*       Climate model parameters
        fex0     2015 forcings of non-CO2 GHG (Wm-2)              /0.5  /
        fex1     2100 forcings of non-CO2 GHG (Wm-2)              /1.0  /
        eland0   Carbon emissions from land 2015 (GtCO2 per year) /2.6  /
        deland   Decline rate of land emissions (per period)      /.115 /
*       Abatement cost
        pback    Cost of backstop 2010$ per tCO2 2015             /550  /
        gback    Initial cost decline backstop cost per period    /.025 /;
*       Discounting factor
        disc=1/(1+prstp)**5 ;
*       Parameters for long-run consistency of carbon cycle
        b11 = 1 - b12;
        b21 = b12*mate/mue;
        b22 = 1 - b21 - b23;
        b32 = b23*mue/mle;
        b33 = 1 - b32;

*        Exogenous variables (originated from DICE 2016R)
Parameters
         tstep          Years per Period                                 /5      /
         pop0           Initial world population 2015 (millions)         /7403   /
         popasym        Asymptotic level of labor                        /11500  /
         popadj         Growth rate to calibrate to 2050 pop projection  /0.134  /
         q0             Initial world gross output 2015 (trill 2010 USD) /105.5  /
         e0             Industrial emissions 2015 (GtCO2 per year)       /35.85  /
         miu0           Initial emissions control rate for base case 2015/.03    /
         a0             Initial level of total factor productivity       /5.115  /
         ga0            Initial growth rate for TFP per 5 years          /0.076  /
         dela           Decline rate of TFP per 5 years                  /0.005  /
         gsigma1        Initial growth of sigma (per year)               /-0.0152/
         sig0           Carbon intensity 2010 (kgCO2 per output 2005 USD 2010)
         dsig           Decline rate of decarbonization (per period)     /-0.001 /
         ga(t)          Growth rate of productivity from
         gsig(t)        Change in sigma (cumulative improvement of energy efficiency)
         lbar_t(t)      Level of population at a current period
         abar_t(t)      Level of total factor productivity at a current period
         sbar_t(t)      CO2-equivalent-emissions output ratio at a current period
         forcoth(t)     Exogenous forcing for other greenhouse gases
         etree(t)       Emissions from deforestation
         cost1(t)       Adjusted cost for backstop
         to_basis(t)    Approximation of the law of ocean temperature;

sig0 = e0/(q0*(1-miu0));
lbar_t(t)$(t.val eq 1) = pop0;
loop(t, lbar_t(t+1)=lbar_t(t););
loop(t, lbar_t(t+1)=lbar_t(t)*(popasym/lbar_t(t))**popadj ;);
lbar_t(t) = lbar_t(t)/1000;
ga(t) = ga0*exp(-dela*5*((t.val-1)));
abar_t(t)$(t.val eq 1) = a0;  loop(t, abar_t(t+1)=abar_t(t)/((1-ga(t))););
gsig(t)$(t.val eq 1)=gsigma1; loop(t,gsig(t+1)=gsig(t)*((1+dsig)**tstep) ;);
sbar_t(t)$(t.val eq 1)=sig0;  loop(t,sbar_t(t+1)=(sbar_t(t)*exp(gsig(t)*tstep)););
forcoth(t) = fex0 + (1/17)*(fex1-fex0)*((t.val -1))$(t.val lt 18) + (fex1-fex0)$(t.val ge 18);
etree(t)   = eland0*(1-deland)**((t.val-1));
cost1(t)   = pback*(1-gback)**((t.val-1))*SBAR_T(t)/expcost2/1000;

* Data for Smolyak grid

*        Set Nodes and Smoyak sparse grid
sets    sp      State variable set indices      /CAP,MAT,PHI,MU,ML,OT,BMIU,BSIG/,
        ia      Universal node index            /1*%gridn%/,
        ik(ia)  Capital node index              /1*%gridn_k%/,
        im(ia)  MAT node index                  /1*%gridn_m%/,
        iu(ia)  mu node index                   /1*%gridn_u%/,
        il(ia)  ml node index                   /1*%gridn_l%/,
        ip(ia)  phit node index                 /1*%gridn_p%/,
        iot(ia) OT node index                   /1*%gridn_ot%/,
        ibm(ia) BMIU node index                  /1*%gridn_bm%/,
        ibs(ia) BSIG node index                  /1*%gridn_bs%/,
        ic(ia)  Order of Chebyshev polynomial   /1*%gridn%/,
        iter    Dynamic programming iterations  /1*%niter%/;

alias(ik,jk,kk)
alias(im,jm,km)
alias(iu,ju,ku)
alias(il,jl,kl)
alias(ip,jp,kp)
alias(iot,jot,kot)
alias(ibm,jbm,kbm)
alias(ibs,jbs,kbs)
alias(universal,*);

$setglobal isn        ik,im,iu,il,ip,iot,ibm,ibs
$setglobal jsn        jk,jm,ju,jl,jp,jot,jbm,jbs
$setglobal ksn        kk,km,ku,kl,kp,kot,kbm,kbs

set     sgrid(%isn%)      grid points for Smolyak sparse grid;
sgrid(%isn%) = yes;
alias(sgrid,ssgrid);

*-------------------------------------------------------------------------------------------
*       simplicial complete polynomial
*-------------------------------------------------------------------------------------------

set     cgrid(%isn%)      grid points for Smolyak sparse grid;

*       compute alpha/n for each dimension.
parameter       alphan(%isn%);
alphan(%isn%) = (ik.val-1)/(card(ik)-1) +
                (im.val-1)/(card(im)-1) +
                (iu.val-1)/(card(iu)-1) +
                (il.val-1)/(card(il)-1) +
*                (ip.val-1)/(card(ip)-1) ;
                (ip.val-1)/(card(ip)-1) +
                (iot.val-1)/(card(iot)-1) +
                (ibm.val-1)/(card(ibm)-1) +
                (ibs.val-1)/(card(ibs)-1);

cgrid(%isn%) = yes$(alphan(%isn%) le 1);
alias(cgrid,ccgrid)

parameter ncgrid, nsgrid;
ncgrid = sum(cgrid(%isn%),1);
nsgrid = sum(sgrid(%isn%),1);

*       check number of grid points
display ncgrid, nsgrid;

*       Specify lower and upper bounds for each state var
Parameters
        lo(sp)          Lowerbound on state variable,
        up(sp)          Upperbound on state variable,
        arg             Argument of cosine weighting function,
        x               Node value for the state variable on the unit interval,
        x,xx            Node value for the state variable on the unit interval,
        x_k             Node value for the state variable on the unit interval,
        x_m             Node value for the state variable on the unit interval,
        x_p             Node value for the state variable on the unit interval,
        x_u             Node value for the state variable on the unit interval,
        x_l             Node value for the state variable on the unit interval,
        x_ot            Node value for the state variable on the unit interval,
        x_bm            Node value for the state variable on the unit interval,
        x_bs            Node value for the state variable on the unit interval,
        lev_k           Level value at node for grid point calculation,
        lev_kt          Level value of augmented capital (TFP and labor),
        lev_m           Level value at node for grid point calculation,
        lev_p           Level value at node for grid point calculation,
        lev_mu          Level value at node for grid point calculation,
        lev_ml          Level value at node for grid point calculation,
        lev_ot          Level value at node for grid point calculation,
        lev_bm          Level value at node for grid point calculation,
        lev_s           Level value at node for grid point calculation,
        lev_bs          Level value at node for grid point calculation,
        sigma           Normalized standard deviation of phi /0.1/;

*       Define basis for chebyshev polynomial expansion
arg(ik) = ((2*ord(ik)-1)*pi)/(2*card(ik));
xx(ik) = cos(arg(ik))*.85;
x_k(ik) = xx(ik);

arg(im) = ((2*ord(im)-1)*pi)/(2*card(im));
xx(im) = cos(arg(im))*.7;
x_m(im) = xx(im);

arg(iu) = ((2*ord(iu)-1)*pi)/(2*card(iu));
xx(iu) = cos(arg(iu))*1;
x_u(iu) = xx(iu);

arg(il) = ((2*ord(il)-1)*pi)/(2*card(il));
xx(il) = cos(arg(il))*1;
x_l(il) = xx(il);

arg(ip) = ((2*ord(ip)-1)*pi)/(2*card(ip));
xx(ip) = cos(arg(ip))*.7;
x_p(ip) = xx(ip);

arg(iot) = ((2*ord(iot)-1)*pi)/(2*card(iot));
xx(iot) = cos(arg(iot))*1;
x_ot(iot) = xx(iot);

arg(ibm) = ((2*ord(ibm)-1)*pi)/(2*card(ibm));
xx(ibm) = cos(arg(ibm))*1;
x_bm(ibm) = xx(ibm);

arg(ibs) = ((2*ord(ibs)-1)*pi)/(2*card(ibs));
xx(ibs) = cos(arg(ibs))*1;
x_bs(ibs) = xx(ibs);

*        Boundary
up("CAP") = 6;
up("PHI") = 9.9 ;
up("OT")  = up("PHI");
up("MAT") = 1300/mate;
up("MU")  = up("MAT")*(1-b11)/b21;
up("ML")  = up("MU")*b23/(1-b33);
up("BMIU")= 12;
up("BSIG")= 0.013;

lo("CAP") = .5;
lo("PHI") = (2*(0.75)-up("PHI")-up("PHI")*(smin(ip,x_p(ip))+(1e-15)))/(1-(smin(ip,x_p(ip))+(1e-15)));
lo("OT")  = (2*(0.0068)-up("OT")-up("OT")*(smin(iot,x_ot(iot))+(1e-15)))/(1-(smin(iot,x_ot(iot))+(1e-15)));
lo("MAT") = (2*(851/mate)-up("MAT")-up("MAT")*(smin(im,x_m(im))+(1e-15)))/(1-(smin(im,x_m(im))+(1e-15)));
lo("MU")  = (2*(460/mate)-up("MU")-up("MU")*(smin(iu,x_u(iu))+(1e-15)))/(1-(smin(iu,x_u(iu))+(1e-15)));
lo("ML")  = (2*(1740/mate)-up("ML")-up("ML")*(smin(il,x_l(il))+(1e-15)))/(1-(smin(il,x_l(il))+(1e-15)));
lo("BMIU")= (2*(2.8)-up("BMIU")-up("BMIU")*(smin(ibm,x_bm(ibm))+(1e-15)))/(1-(smin(ibm,x_bm(ibm))+(1e-15)));
lo("BSIG")= (2*(1e-8)-up("BSIG")-up("BSIG")*(smin(ibs,x_bs(ibs))+(1e-15)))/(1-(smin(ibs,x_bs(ibs))+(1e-15)));

*        Grid points
lev_k(ik)   = (lo("cap")+up("cap")+(up("cap")-lo("cap"))*x_k(ik))/2;
lev_m(im)   = (lo("mat")+up("mat")+(up("mat")-lo("mat"))*x_m(im))/2;
lev_p(ip)   = (lo("phi")+up("phi")+(up("phi")-lo("phi"))*x_p(ip))/2;
lev_mu(iu)  = (lo("mu")+up("mu")+(up("mu")-lo("mu"))*x_u(iu))/2;
lev_ml(il)  = (lo("ml")+up("ml")+(up("ml")-lo("ml"))*x_l(il))/2;
lev_kt(t,ik)= lev_k(ik)*lbar_t(t)*abar_t(t)**(1/(1-gama));
lev_ot(iot) = (lo("OT")+up("OT")+(up("OT")-lo("OT"))*x_ot(iot))/2;
lev_s(ibm)  = (lo("BMIU")+up("BMIU")+(up("BMIU")-lo("BMIU"))*x_bm(ibm))/2;
lev_bs(ibs) = (lo("BSIG")+up("BSIG")+(up("BSIG")-lo("BSIG"))*x_bs(ibs))/2;
lev_bm(ibm) = 1-lambda/lev_s(ibm) ;

display lev_k, lev_m, lev_mu, lev_ml, lev_p, lev_ot, lev_s, lev_bm, lev_bs, lev_kt;

parameters par_lev;
par_lev('lev_k',ia)    = lev_k(ia);
par_lev('lev_m',ia)    = lev_m(ia);
par_lev('lev_p',ia)    = lev_p(ia);
par_lev('lev_mu',ia)   = lev_mu(ia);
par_lev('lev_ml',ia)   = lev_ml(ia);
par_lev('lev_ot',ia)   = lev_ot(ia);
par_lev('lev_s',ia)    = lev_s(ia);
par_lev('lev_bs',ia)   = lev_bs(ia);
par_lev('lev_mlev',ia) = lev_m(ia)*588;
par_lev('lev_mulev',ia)= lev_mu(ia)*588;
par_lev('lev_mllev',ia)= lev_ml(ia)*588;

* Chebyshev polynomial terms
$eval ptcard ( round(card(ia)/2) )
$eval ptcard_k ( round(card(ik)/2) )
$eval ptcard_m ( round(card(im)/2) )
$eval ptcard_u ( round(card(iu)/2) )
$eval ptcard_l ( round(card(il)/2) )
$eval ptcard_p ( round(card(ip)/2) )
$eval ptcard_ot ( round(card(iot)/2) )
$eval ptcard_bm ( round(card(ibm)/2) )
$eval ptcard_bs ( round(card(ibs)/2) )

set     pt              Chebyshev polynomial terms /1 * %ptcard%/,
        pt_k(pt)        Chebyshev polynomial terms /1 * %ptcard_k%/,
        pt_m(pt)        Chebyshev polynomial terms /1 * %ptcard_m%/,
        pt_u(pt)        Chebyshev polynomial terms /1 * %ptcard_u%/,
        pt_l(pt)        Chebyshev polynomial terms /1 * %ptcard_l%/,
        pt_p(pt)        Chebyshev polynomial terms /1 * %ptcard_p%/,
        pt_ot(pt)        Chebyshev polynomial terms /1 * %ptcard_ot%/,
        pt_bm(pt)        Chebyshev polynomial terms /1 * %ptcard_bm%/,
        pt_bs(pt)        Chebyshev polynomial terms /1 * %ptcard_bs%/,

        csp_k(ia,pt)    Index for Chebyshev polynomial coefficients and exponents,
        csp_m(ia,pt)    Index for Chebyshev polynomial coefficients and exponents,
        csp_u(ia,pt)    Index for Chebyshev polynomial coefficients and exponents,
        csp_l(ia,pt)    Index for Chebyshev polynomial coefficients and exponents,
        csp_p(ia,pt)    Index for Chebyshev polynomial coefficients and exponents,
        csp_ot(ia,pt)   Index for Chebyshev polynomial coefficients and exponents,
        csp_bm(ia,pt)   Index for Chebyshev polynomial coefficients and exponents,
        csp_bs(ia,pt)   Index for Chebyshev polynomial coefficients and exponents;

parameters
        ce_k(ia,pt)     Exponents of state variables for Chebyshev Polynomial terms,
        ce_m(ia,pt)     Exponents of state variables for Chebyshev Polynomial terms,
        ce_u(ia,pt)     Exponents of state variables for Chebyshev Polynomial terms,
        ce_l(ia,pt)     Exponents of state variables for Chebyshev Polynomial terms,
        ce_p(ia,pt)     Exponents of state variables for Chebyshev Polynomial terms,
        ce_ot(ia,pt)     Exponents of state variables for Chebyshev Polynomial terms,
        ce_bm(ia,pt)     Exponents of state variables for Chebyshev Polynomial terms,
        ce_bs(ia,pt)     Exponents of state variables for Chebyshev Polynomial terms,
        cc_k(*,*)       Coefficients of state variables for Chebyshev Polynomial terms,
        cc_m(*,*)       Coefficients of state variables for Chebyshev Polynomial terms,
        cc_u(*,*)       Coefficients of state variables for Chebyshev Polynomial terms,
        cc_l(*,*)       Coefficients of state variables for Chebyshev Polynomial terms,
        cc_p(*,*)       Coefficients of state variables for Chebyshev Polynomial terms,
        cc_ot(*,*)       Coefficients of state variables for Chebyshev Polynomial terms,
        cc_bm(*,*)       Coefficients of state variables for Chebyshev Polynomial terms,
        cc_bs(*,*)       Coefficients of state variables for Chebyshev Polynomial terms;

*       Assign polynomial coefficients to express each chebyshev polynomial basis
*       e.g. cc("1","1") = coefficient for 1st polynomial term of 1st order Chebyshev basis

table  cc(*,*)  Coefficients of state variables for Chebyshev Polynomial terms

                1       2       3       4       5       6
        1       1
        2       1
        3       2       1
        4       4       3
        5       8       8       1
        6       16      20      5
        7       32      48      18      1
        8       64      112     56      7
        9       128     256     160     32      1
        10      256     576     432     120     9
        11      512     1280    1120    400     50      1;

parameters
        bar_k(ik,ia)            Polynomial terms used for least squares estimation,
        bar_m(im,ia)            Polynomial terms used for least squares estimation,
        bar_mu(iu,ia)           Polynomial terms used for least squares estimation,
        bar_ml(il,ia)           Polynomial terms used for least squares estimation,
        bar_p(ip,ia)            Polynomial terms used for least squares estimation,
        bar_ot(iot,ia)            Polynomial terms used for least squares estimation,
        bar_bm(ibm,ia)            Polynomial terms used for least squares estimation,
        bar_bs(ibs,ia)            Polynomial terms used for least squares estimation;

parameter       coeff(t,%isn%)       Estimated Coefficients;

*-------------------------------------------------------------------------------------------
*       capital
*-------------------------------------------------------------------------------------------

csp_k(ik,pt_k) = yes$cc(ik,pt_k);

*       Exponents are assigned
ce_k(csp_k(ik,pt_k)) = ik.val - 1 - (pt_k.val-1)*2;

*       Signs of coefficients alternate
cc_k(csp_k(ik,pt_k)) = cc(ik,pt_k) * power(-1,(pt_k.val-1));

*       Apply above algorithm to compute chebyshev polynomial terms
bar_k(ik,jk) = sum(pt_k$csp_k(jk,pt_k),cc_k(jk,pt_k)*power(x_k(ik),ce_k(jk,pt_k)));

*-------------------------------------------------------------------------------------------
*       mat
*-------------------------------------------------------------------------------------------

csp_m(im,pt_m) = yes$cc(im,pt_m);

*       Exponents are assigned
ce_m(csp_m(im,pt_m)) = im.val - 1 - (pt_m.val-1)*2;

*       Signs of coefficients alternate
cc_m(csp_m(im,pt_m)) = cc(im,pt_m) * power(-1,(pt_m.val-1));

*       Apply above algorithm to compute chebyshev polynomial terms
bar_m(im,jm) = sum(pt_m$csp_m(jm,pt_m),cc_m(jm,pt_m)*power(x_m(im),ce_m(jm,pt_m)));

*-------------------------------------------------------------------------------------------
*       mu
*-------------------------------------------------------------------------------------------

csp_u(iu,pt_u) = yes$cc(iu,pt_u);

*       Exponents are assigned
ce_u(csp_u(iu,pt_u)) = iu.val - 1 - (pt_u.val-1)*2;

*       Signs of coefficients alternate
cc_u(csp_u(iu,pt_u)) = cc(iu,pt_u) * power(-1,(pt_u.val-1));

*       Apply above algorithm to compute chebyshev polynomial terms
bar_mu(iu,ju) = sum(pt_u$csp_u(ju,pt_u),cc_u(ju,pt_u)*power(x_u(iu),ce_u(ju,pt_u)));

*-------------------------------------------------------------------------------------------
*       ml
*-------------------------------------------------------------------------------------------

csp_l(il,pt_l) = yes$cc(il,pt_l);

*       Exponents are assigned
ce_l(csp_l(il,pt_l)) = il.val - 1 - (pt_l.val-1)*2;

*       Signs of coefficients alternate
cc_l(csp_l(il,pt_l)) = cc(il,pt_l) * power(-1,(pt_l.val-1));

*       Apply above algorithm to compute chebyshev polynomial terms
bar_ml(il,jl) = sum(pt_l$csp_l(jl,pt_l),cc_l(jl,pt_l)*power(x_l(il),ce_l(jl,pt_l)));

*-------------------------------------------------------------------------------------------
*       phi
*-------------------------------------------------------------------------------------------

csp_p(ip,pt_p) = yes$cc(ip,pt_p);

*       Exponents are assigned
ce_p(csp_p(ip,pt_p)) = ip.val-1 - (pt_p.val-1)*2;

*       Signs of coefficients alternate
cc_p(csp_p(ip,pt_p)) = cc(ip,pt_p) * power(-1,(pt_p.val-1));

*       Apply above algorithm to compute chebyshev polynomial terms
bar_p(ip,jp) = sum(pt_p$csp_p(jp,pt_p),cc_p(jp,pt_p)*power(x_p(ip),ce_p(jp,pt_p)));

*-------------------------------------------------------------------------------------------
*       OT
*-------------------------------------------------------------------------------------------

csp_ot(iot,pt_ot) = yes$cc(iot,pt_ot);

*       Exponents are assigned
ce_ot(csp_ot(iot,pt_ot)) = iot.val-1 - (pt_ot.val-1)*2;

*       Signs of coefficients alternate
cc_ot(csp_ot(iot,pt_ot)) = cc(iot,pt_ot) * power(-1,(pt_ot.val-1));

*       Apply above algorithm to compute chebyshev polynomial terms
bar_ot(iot,jot) = sum(pt_ot$csp_ot(jot,pt_ot),cc_ot(jot,pt_ot)*power(x_ot(iot),ce_ot(jot,pt_ot)));

*-------------------------------------------------------------------------------------------
*       BMIU
*-------------------------------------------------------------------------------------------

csp_bm(ibm,pt_bm) = yes$cc(ibm,pt_bm);

*       Exponents are assigned
ce_bm(csp_bm(ibm,pt_bm)) = ibm.val-1 - (pt_bm.val-1)*2;

*       Signs of coefficients alternate
cc_bm(csp_bm(ibm,pt_bm)) = cc(ibm,pt_bm) * power(-1,(pt_bm.val-1));

*       Apply above algorithm to compute chebyshev polynomial terms
bar_bm(ibm,jbm) = sum(pt_bm$csp_bm(jbm,pt_bm),cc_bm(jbm,pt_bm)*power(x_bm(ibm),ce_bm(jbm,pt_bm)));

*-------------------------------------------------------------------------------------------
*       BSIG
*-------------------------------------------------------------------------------------------

csp_bs(ibs,pt_bs) = yes$cc(ibs,pt_bs);

*       Exponents are assigned
ce_bs(csp_bs(ibs,pt_bs)) = ibs.val-1 - (pt_bs.val-1)*2;

*       Signs of coefficients alternate
cc_bs(csp_bs(ibs,pt_bs)) = cc(ibs,pt_bs) * power(-1,(pt_bs.val-1));

*       Apply above algorithm to compute chebyshev polynomial terms
bar_bs(ibs,jbs) = sum(pt_bs$csp_bs(jbs,pt_bs),cc_bs(jbs,pt_bs)*power(x_bs(ibs),ce_bs(jbs,pt_bs)));

*-------------------------------------------------------------------------------------------
*       Macros used in the program
*-------------------------------------------------------------------------------------------

*       Normalized value of state variables that appear in the equations
$macro K_N(%isn%) ((EF_K(t+1,%isn%)-(lo("cap")+up("cap"))/2)/((up("cap")-lo("cap"))/2))

$macro MAT_N(%isn%) ((MAT(t+1,%isn%)-(lo("mat")+up("mat"))/2)/((up("mat")-lo("mat"))/2))

$macro PHIT_N(%isn%) ((PHIT(t+1,%isn%,ig,ih)-(lo("phi")+up("phi"))/2)/((up("phi")-lo("phi"))/2))

$macro K_CS(%isn%,jk) (sum(pt_k$csp_k(jk,pt_k),cc_k(jk,pt_k) * power(K_N(%isn%),ce_k(jk,pt_k))))
$macro MAT_CS(%isn%,jm) (sum(pt_m$csp_m(jm,pt_m),cc_m(jm,pt_m) * power(MAT_N(%isn%),ce_m(jm,pt_m))))
$macro PHIT_CS(%isn%,jp) (sum(pt_p$csp_p(jp,pt_p),cc_p(jp,pt_p) * power(PHIT_N(%isn%),ce_p(jp,pt_p))))
$macro MU_CS(%isn%,ju) (sum(pt_u$csp_u(ju,pt_u),cc_u(ju,pt_u) * power(mu_n(%isn%),ce_u(ju,pt_u))))
$macro ML_CS(%isn%,jl) (sum(pt_l$csp_l(jl,pt_l),cc_l(jl,pt_l) * power(ml_n(%isn%),ce_l(jl,pt_l))))
$macro OT_CS(%isn%,jot) (sum(pt_ot$csp_ot(jot,pt_ot),cc_ot(jot,pt_ot) * power(ot_n(%isn%),ce_ot(jot,pt_ot))))
$macro BM_CS(%isn%,jbm) (sum(pt_bm$csp_bm(jbm,pt_bm),cc_bm(jbm,pt_bm) * power(s_n(%isn%,ig,ih),ce_bm(jbm,pt_bm))))
$macro BS_CS(%isn%,jbs) (sum(pt_bs$csp_bs(jbs,pt_bs),cc_bs(jbs,pt_bs) * power(bsig_n(%isn%),ce_bs(jbs,pt_bs))))

$macro PVL (sum(ccgrid(%jsn%), A(t,%jsn%) * \
                                bar_k(ik,jk) * \
                                bar_m(im,jm) * \
                                bar_mu(iu,ju) * \
                                bar_ml(il,jl) * \
                                bar_p(ip,jp) * \
                                bar_ot(iot,jot) * \
                                bar_bm(ibm,jbm) * \
                                bar_bs(ibs,jbs)))

$macro PV(A) (sum(ccgrid(%jsn%), A(t+1,%jsn%) * \
                                 K_CS(%isn%,jk) * \
                                 MAT_CS(%isn%,jm) * \
                                 MU_CS(%isn%,ju) * \
                                 ML_CS(%isn%,jl) * \
                                 PHIT_CS(%isn%,jp) * \
                                 OT_CS(%isn%,jot)  * \
                                 BM_CS(%isn%,jbm)  * \
                                 BS_CS(%isn%,jbs) ))

*-------------------------------------------------------------------------------------------
* Gauss-Hermite Numerical Integration for the feedback factor and the temperature shock (Double expectation)
*-------------------------------------------------------------------------------------------
$if not set nGH       $set nGH        2    #choose among 2, 3, 5

set     ig               Gaussian-Hermite grid points /1*%nGH%/
alias (ig,ih);

parameter ghdata(ig,*)    Tabulated Gaussian-Hermite approximation points;

$ifi %nGH%=='2' ghdata('1','zeta') = -0.70710678119; ghdata('1','omega') = 0.88622692545; ghdata('2','zeta') = (-1)*ghdata('1','zeta'); ghdata('2','omega') = ghdata('1','omega');
$ifi %nGH%=='3' ghdata('1','zeta') = -1.224744871; ghdata('1','omega') = 0.295408975; ghdata('2','zeta') = 0; ghdata('2','omega') = 1.181635901; ghdata('3','zeta') = (-1)*ghdata('1','zeta'); ghdata('3','omega') = ghdata('1','omega');
$ifi %nGH%=='5' ghdata('1','zeta') = -2.02018287 ; ghdata('1','omega') = 0.019953242; ghdata('2','zeta') = -0.95857247; ghdata('2','omega') = 0.393619323; ghdata('3','zeta') = 0; ghdata('3','omega') =  0.94530872;
$ifi %nGH%=='5' ghdata('4','zeta') = (-1)*ghdata('2','zeta'); ghdata('4','omega') = ghdata('2','omega'); ghdata('5','zeta') = (-1)*ghdata('1','zeta'); ghdata('5','omega') = ghdata('1','omega');

parameters
        omega(ig)                Normalized GH weights
        stoc_bm(%isn%,ig)        Quadrature points for the feedback factor ~N(lev_bm lev_bs)
        stoc_bm_init(ig)         Quadrature points for the feedback factor ~N(lev_bm lev_bs)
        stoc_ts(ig)              Quadrature points for the temperature shock: ~N(0 var_e) ;

*       Take expectation of function value
omega(ig) = ghdata(ig,"omega")/sqrt(3.141592);
stoc_bm(sgrid(%isn%),ig) = sqrt(2*lev_bs(ibs))*ghdata(ig,"zeta") + (min(lev_bm(ibm),0.99));
stoc_bm_init(ig) = sqrt(2*lev_bs_init)*ghdata(ig,"zeta") + (1-lambda/lev_s_init) ;
stoc_ts(ih) = sqrt(2*var_e)*ghdata(ih,"zeta");

parameter lo_bm, up_bm, check_snodes(%isn%,ig),check_snodes_init(ig);
up_bm = 1-lambda/up("bmiu");
lo_bm = 1-lambda/max(0.000001,lo("bmiu"));
check_snodes(sgrid(%isn%),ig) = yes$(stoc_bm(%isn%,ig) lt lo_bm or stoc_bm(%isn%,ig) gt up_bm);
check_snodes_init(ig) = yes$(stoc_bm_init(ig) lt lo_bm or stoc_bm_init(ig) gt up_bm);

parameter        exp_f    Mean of the feedback factor
                 median_s Mean of the climate sensitivity;
median_s = sum(ig, omega(ig) * lambda/(1-(stoc_bm_init(ig))));
t2xco2 = median_s ;
exp_f = 1-lambda/t2xco2;

*-------------------------------------------------------------------------------------------
* Equations
*-------------------------------------------------------------------------------------------
nonnegative variables
        C(t,%isn%)              Consumption,
        EF_K(t,%isn%)           Effective Capital,
        MAT(t,%isn%)            Carbon concentration (proportion to pre-industrial level 588 GtC),
        PHIT(t,%isn%,ig,ih)     Atmospheric Temperature (degree C),
        MIU(t,%isn%)            Carbon control rate (0<=MIU<=1);

variables
        OBJ                Objective,
        OBJL               Least Squares Objective,
        TU(t,%isn%)        Nodal approximations of total utility,
        A(t,%isn%)         Estimated Coefficients,
        I(t,%isn%)         Investment,
        E(t,%isn%)         Emissions,
        ESTGRAD_K(t,%isn%) Gradient of value function wrt capital,
        TRGRAD_K(t,%isn%)  Gradient of value function wrt capital,
        ESTGRAD_M(t,%isn%) Gradient of value function wrt carbon stock
        ESTSCC(t,%isn%)    Estimated SCC
        ESTSCC_valfn(t,%isn%)    Estimated SCC from the value function;

parameters
        mu(%isn%)         Measure of carbon concentration increase in shallow oceans,
        mu_n(%isn%)       normalized mu,
        ml(%isn%)         Measure of carbon concentration increase in deeper oceans,
        ml_n(%isn%)       normalized ml,
        ot(%isn%)         Measure of carbon concentration increase in Ocean temperature,
        ot_n(%isn%)       normalized ot,
        ygross(t,%isn%)   Gross world product including abatement and damages (tril 2005 USD per year),
        damfrac(ip)       Damages as fraction of gross output,
        s_n(%isn%,ig,ih)  normalized bmiu,
        bsig_n(%isn%)     normalized bsig,
        bmiu(%isn%,ig,ih) Mean of climate sensitivity (belief)
        bsig(%isn%)       Variance of climate sensitivity (a measure of belief) ;

*-------------------------------------------------------------------------------------------
*       Macros to reduce problem size
*-------------------------------------------------------------------------------------------
$macro K(%isn%)        (EF_K(t+1,%isn%)*lbar_t(t+1)*abar_t(t+1)**(1/(1-gama)))
$macro ABTFRC(t,%isn%) (cost1(t)*MIU(t,%isn%)**expcost2)
$macro FORCE(%isn%)    (fco22x * (log(MAT(t+1,%isn%))/log(2)) + forcoth(t+1))
$macro Y(t,%isn%)      (ygross(t,%isn%)*(1-damfrac(ip)-ABTFRC(t,%isn%)))
$macro margin          (c1*fco22x/lambda*lev_p(ip))

*-------------------------------------------------------------------------------------------
*       Parameters
*-------------------------------------------------------------------------------------------
ygross(t,sgrid(%isn%)) = abar_t(t)*lbar_t(t)**(1-gama)*lev_kt(t,ik)**gama;
damfrac(ip) = a2*lev_p(ip)**a3;
mu(sgrid(%isn%)) = lev_m(im)*b12 + lev_mu(iu)*b22;
mu_n(sgrid)  = (mu(sgrid)-(lo("MU")+up("MU"))/2)/((up("MU")-lo("MU"))/2);
ml(sgrid(%isn%)) = lev_ml(il)*b33 + lev_mu(iu)*b23;
ml_n(sgrid)  = (ml(sgrid)-(lo("ML")+up("ML"))/2)/((up("ML")-lo("ML"))/2);
ot(sgrid(%isn%)) = lev_ot(iot)+c4*(lev_p(ip)-lev_ot(iot));
ot_n(sgrid)  = (ot(sgrid)-(lo("OT")+up("OT"))/2)/((up("OT")-lo("OT"))/2);
bmiu(sgrid(%isn%),ig,ih) = max( lo_bm, min( up_bm ,((margin*stoc_bm(%isn%,ig)+stoc_ts(ih))*margin*lev_bs(ibs)+lev_bm(ibm)*var_e)/(margin**2*lev_bs(ibs)+var_e))); # Assume Bayesian Learning
bsig(sgrid(%isn%))    = lev_bs(ibs)*var_e/((margin**2)*lev_bs(ibs)+var_e) ;  # Assume Bayesian Learning
s_n(sgrid,ig,ih) = min(1, ((lambda/(1-bmiu(sgrid,ig,ih))) -(lo("BMIU")+up("BMIU"))/2)/((up("BMIU")-lo("BMIU"))/2));
bsig_n(sgrid)= (bsig(sgrid)-(lo("BSIG")+up("BSIG"))/2)/((up("BSIG")-lo("BSIG"))/2);

*-------------------------------------------------------------------------------------------
*       Bellman Equations
*-------------------------------------------------------------------------------------------
equations
        objdef,phitndef,matdef,cdef,idef,edef,tu_nonneg;

cdef(ct(t),sgrid(%isn%))..      C(t,%isn%) =e= Y(t,%isn%) - I(t,%isn%);

idef(ct(t),sgrid(%isn%))..      I(t,%isn%) =e= (K(%isn%) - (1-dk)**5*lev_kt(t,ik))/5;

edef(ct(t),sgrid(%isn%))..      E(t,%isn%) =e= (1-MIU(t,%isn%))*sbar_t(t)*ygross(t,%isn%) + etree(t);

phitndef(t+1,sgrid(%isn%),ig,ih)$ct(t)..
                                PHIT(t+1,%isn%,ig,ih) =e=
                                   lev_p(ip) +
                                   c1 * (FORCE(%isn%) - fco22x/lambda*lev_p(ip) +
                                          fco22x/lambda*stoc_bm(%isn%,ig)*lev_p(ip) -
                                          c3*(lev_p(ip)-lev_ot(iot))) +
                                   stoc_ts(ih);

matdef(t+1,sgrid(%isn%))$ct(t)..
                                MAT(t+1,%isn%) =e=
                                lev_m(im)*b11 + lev_mu(iu)*b21 + (E(t,%isn%)*5/3.666)/mate;

objdef..        OBJ =e= sum(sgrid(%isn%),

                         sum(ct(t)$t_exc_last(t),

                            (1-disc)*(C(t,%isn%))**(1-elasmu)*lbar_t(t)**elasmu/(1-elasmu)+
                            disc*sum((ig,ih), omega(ig)*omega(ih)*PV(coeff))) +

                         sum(ct(t)$t_last(t),
                            (1-disc)*(C(t,%isn%))**(1-elasmu)*lbar_t(t)**elasmu/(1-elasmu)));

tu_nonneg(ct(t),sgrid(%isn%)).. (1-disc)*(C(t,%isn%))**(1-elasmu)*lbar_t(t)**elasmu/(1-elasmu)+
                            disc*sum((ig,ih), omega(ig)*omega(ih)*PV(coeff)) =l= 0;

*-------------------------------------------------------------------------------------------
*       LSQR Problem
*-------------------------------------------------------------------------------------------
equations
        lsqrdef         lsqr objective,
        lsqrdef_gt      lsqr objective with gradient targeting,
        gradkdef        gradient expression wrt capital,
        gradmdef        gradient expression wrt carbon stock,
        sccdef          social cost of carbon calculation,
        sccdef_val      social cost of carbon calculation from the estimated value function ;

lsqrdef..       OBJL =e= sum((ct(t),%isn%),sqr(PVL - TU(t,%isn%)));

gradkdef(ct(t),sgrid(%isn%))..
                        ESTGRAD_K(t,%isn%) =e=
                        1/((up("cap") - lo("cap"))/2) * (
                        sum(cgrid(%jsn%), A(t,%jsn%) *
                          sum(km,bar_m(im,km)$(km.val eq jm.val))  *
                          sum(ku,bar_mu(iu,ku)$(ku.val eq ju.val)) *
                          sum(kl,bar_ml(il,kl)$(kl.val eq jl.val)) *
                          sum(kp,bar_p(ip,kp)$(kp.val eq jp.val))  *
                          sum(kot,bar_ot(iot,kot)$(kot.val eq jot.val)) *
                          sum(kbm,bar_bm(ibm,kbm)$(kbm.val eq jbm.val)) *
                          sum(kbs,bar_bs(ibs,kbs)$(kbs.val eq jbs.val)) *

                          sum(kk$(kk.val eq jk.val),
                                sum(pt_k$csp_k(kk,pt_k),
                                  (cc_k(kk,pt_k)*ce_k(kk,pt_k))$(ce_k(kk,pt_k) ge 1) *
                                  (power(x_k(ik),ce_k(kk,pt_k)-1)$(x_k(ik) and ce_k(kk,pt_k) ge 1) +
                                   1$(not x_k(ik) and ce_k(kk,pt_k) eq 1))
                                )
                          )
                        ));

gradmdef(ct(t),sgrid(%isn%))..
                        ESTGRAD_M(t,%isn%) =e=
                        1/((up("mat") - lo("mat"))/2) * (
                        sum(cgrid(%jsn%), A(t,%jsn%) *
                          sum(kk,bar_k(ik,kk)$(kk.val eq jk.val))  *
                          sum(ku,bar_mu(iu,ku)$(ku.val eq ju.val)) *
                          sum(kl,bar_ml(il,kl)$(kl.val eq jl.val)) *
                          sum(kp,bar_p(ip,kp)$(kp.val eq jp.val))  *
                          sum(kot,bar_ot(iot,kot)$(kot.val eq jot.val)) *
                          sum(kbm,bar_bm(ibm,kbm)$(kbm.val eq jbm.val)) *
                          sum(kbs,bar_bs(ibs,kbs)$(kbs.val eq jbs.val)) *

                          sum(km$(km.val eq jm.val),
                                sum(pt_m$csp_m(km,pt_m),
                                  (cc_m(km,pt_m)*ce_m(km,pt_m))$(ce_m(km,pt_m) ge 1) *
                                  (power(x_m(im),ce_m(km,pt_m)-1)$(x_m(im) and ce_m(km,pt_m) ge 1) +
                                   1$(not x_m(im) and ce_m(km,pt_m) eq 1))
                                )
                          )
                        ));


sccdef_val(ct(t),%isn%)$(sgrid(%isn%) and t.val ne card(t))..
                         ESTSCC_valfn(t,%isn%) =e= -1000*(ESTGRAD_M(t,%isn%)/3.666/mate)/(ESTGRAD_K(t,%isn%)/lbar_t(t)*abar_t(t)**(1/(gama-1)));

sccdef(ct(t),%isn%)$(sgrid(%isn%) and t.val ne card(t))..
                         ESTSCC(t,%isn%) =e= -1000*edef.m(t,%isn%)/(cdef.m(t,%isn%)+1e-10);

model   bellman_reduced /objdef,phitndef,matdef,cdef,idef,edef,tu_nonneg/;
model   lsqr /lsqrdef,gradkdef,gradmdef/;
option  nlp=conopt3;

*-------------------------------------------------------------------------------------------
* Traditional Value Iteration Solve
*-------------------------------------------------------------------------------------------
bellman_reduced.holdfixed = 1;
bellman_reduced.dictfile  = 0;
lsqr.holdfixed = 1;
lsqr.dictfile  = 0;

* Initial guess and boundary conditions
$ondotl
MIU.L(t,sgrid(%isn%))      = .5;
EF_K.L(t,sgrid(%isn%))   = lev_k(ik);
C.L(t,sgrid(%isn%)) = Y(t,%isn%)*0.8;
I.L(t,sgrid(%isn%)) = Y(t,%isn%)*0.2;
E.L(t,sgrid(%isn%)) = (1-MIU(t,%isn%))*sbar_t(t)*ygross(t,%isn%) + etree(t);
MAT.L(t,sgrid(%isn%))    = lev_m(im);
PHIT.L(t,sgrid(%isn%),ig,ih) = lev_p(ip);

ESTGRAD_K.l(t,sgrid(%isn%))=1e-4;
ESTGRAD_M.l(t,sgrid(%isn%))=-(1e-4);
ESTSCC.l(t,%isn%)=500;

MIU.UP(t,sgrid(%isn%))     = 1;
EF_K.UP(t,sgrid(%isn%)) = up("cap")*0.99;
MAT.UP(t+1,sgrid(%isn%))   = up("mat");
PHIT.UP(t+1,sgrid(%isn%),ig,ih) = up("Phi");
C.UP(t_exc_last(t),sgrid(%isn%))     = ygross(t,%isn%);
I.UP(t_exc_last(t),sgrid(%isn%))     = ygross(t,%isn%);
C.LO(t_exc_last,sgrid(%isn%))     = 10;
I.LO(t_exc_last,sgrid(%isn%)) = 0;
MIU.LO(t,sgrid(%isn%))     = 1e-13;
EF_K.LO(t_exc_last,sgrid(%isn%)) = max(0,lo("CAP"));
MAT.LO(t,sgrid(%isn%))   = lo("MAT")*1.01;
PHIT.LO(t,sgrid(%isn%),ig,ih)  = 1e-6;

ESTGRAD_M.up(t_exc_last,sgrid(%isn%)) =0;
ESTGRAD_K.lo(t_exc_last,sgrid(%isn%)) =0;
$offdotl

*        Pre-check or Post-check on the boundary
parameter itlog,itlog2,errorlog,errorlog2,errorlog3,prechk,c_chk,i_chk,tulog;
prechk("mu") = sum(sgrid, mu(sgrid)$(mu(sgrid) gt up("mu") or mu(sgrid) lt lo("mu")));
prechk("ml") = sum(sgrid, ml(sgrid)$(ml(sgrid) gt up("ml") or ml(sgrid) lt lo("ml")));
prechk("ot") = sum(sgrid, ot(sgrid)$(ot(sgrid) gt up("ot") or ot(sgrid) lt lo("ot")));
prechk("bsig") = sum(sgrid, bsig(sgrid)$(bsig(sgrid) gt up("bsig") or bsig(sgrid) lt lo("bsig")));
abort$(sum(universal, prechk(universal)) gt 0) 'Check the domain of state variables';

*-------------------------------------------------------------------------------------------
*       Parameters for gradient targeting
*-------------------------------------------------------------------------------------------
parameter       grad_e, grad_c;
$macro dlev_kt_dlev_k(t)     (lbar_t(t)*abar_t(t)**(1/(1-gama)))
$macro dygross_dlev_kt(t,ik) (abar_t(t)*lbar_t(t)**(1-gama)*gama*lev_kt(t,ik)**(gama-1))

*-------------------------------------------------------------------------------------------
*       Define current time period to terminal time period
*-------------------------------------------------------------------------------------------
A.L(t,cgrid) = 0;
coeff(t,cgrid) = A.L(t,cgrid);

*-------------------------------------------------------------------------------------------
*       VFI Loop
*-------------------------------------------------------------------------------------------
parameter       randu(%isn%)    Assign random values to each isn between 0 and 1;
randu(sgrid) = uniform(0,1);

set             b               Bellman submodels to be solved in parallel (number of cores) /1*%bs%/,
                isnmap(%isn%,b) Association of states and submodels
                isa(%isn%)      Grid index to keep track of parallel solve;

parameter       gap(b)          Gap according to the number of cores;
gap(b) = 1/card(b);
gap(b) = gap(b)*b.val;

isnmap(sgrid,b) = yes$(randu(sgrid) gt (gap(b)-1/card(b)) and randu(sgrid) le gap(b));

option isnmap:0:0:1;
display isnmap;

parameter       handle(b)       Pointer for grid solution;

set             submit(b)       List of regions to submit,
                done(b)         List of regions which are completed;

parameter       nlplog(b,iter,*)        Bellmand solution log;

option RESLIM = 50000000;
option solprint = on;
file ktitle; ktitle.lw=0; ktitle.nd = 0;
lsqr.solvelink = 2;
bellman_reduced.solvelink = 2;

singleton set ctn(t), ctnn(t);
loop(iter,
        ct(t) = yes$(t.val eq card(t)-iter.val+1);
        ctn(t) = yes$(t.val eq card(t)-iter.val+2);
        ctnn(t) = yes$(t.val eq card(t)-iter.val+3);

         MIU.L(ct(t),sgrid(%isn%)) = MIU.L(ctn,%isn%);
         EF_K.L(ct(t),sgrid(%isn%)) = EF_K.L(ctn,%isn%);
         C.L(ct(t),sgrid(%isn%))$(iter.val > 1) = C.L(ctn,%isn%);
         I.L(ct(t),sgrid(%isn%)) = I.L(ctn,%isn%);

         MAT.L(ctn(t),sgrid(%isn%)) = MAT.L(ctnn,%isn%);
         PHIT.L(ctn(t),sgrid(%isn%),ig,ih) = PHIT.L(ctnn,%isn%,ig,ih);

*       Solve in parallel
*       Bellman max problem
$ifthen.threads %threads%==1

        bellman_reduced.solvelink = 2;
        solve bellman_reduced using nlp max OBJ;

$else.threads
        bellman_reduced.solvelink = %solvelink.AsyncThreads%;

        done(b) = no;
        handle(b) = 0;
        repeat

        submit(b) = no;

        loop(b$(not (done(b) or handle(b))),
                submit(b) = yes$(card(submit)+card(handle)<%threads%););

        loop(submit(b),
                isa(%isn%) = isnmap(%isn%,b);
                solve bellman_reduced using nlp max OBJ;
                handle(b) = bellman_reduced.handle; );

        put_utility 'title'/ktitle 'NLP solution.  Finished: ',card(done),', remaining: ',
                     (card(b)-card(done)),', running: ',card(handle),'.';

        display$ReadyCollect(handle) 'Waiting for next instance to collect.';

        loop(b$handlecollect(handle(b)),
            nlplog(b,iter,"modelstat") = bellman_reduced.modelstat;
            nlplog(b,iter,"solvestat") = bellman_reduced.solvestat;
            nlplog(b,iter,"OBJ.L")     = OBJ.L;
            display$handledelete(handle(b)) 'trouble deleting handles' ;
            handle(b) = 0;
            done(b) = yes;);

        until (card(done)=card(b));
        abort$(card(done)<>card(b)) 'Bellman jobs did not return:', handle;

$endif.threads

$ondotl
*       Fix TU to solve LSQR problem
        TU.FX(ct(t),%isn%) =
                (1-disc)*(C(t,%isn%))**(1-elasmu)*lbar_t(t)**elasmu/(1-elasmu)+
                disc*sum((ig,ih), omega(ig)*omega(ih)*PV(coeff));
$offdotl

*       LSQR problem

        A.L(ct(t),cgrid) = A.L(ctn,cgrid);

        solve lsqr using nlp minimizing OBJL;

coeff(ct(t),cgrid) = A.L(t,cgrid);

*-------------------------------------------------------------------------------------------
*       Reporting
*-------------------------------------------------------------------------------------------
$ondotl
        c_chk(t,sgrid(%isn%)) = yes$(C(t,%isn%) lt 1); c_chk(t_last,sgrid(%isn%)) = no;
        abort$(sum((t,sgrid(%isn%)), c_chk(t,%isn%)) gt 0) 'Check the value function: Consumption(C) should be non-negative';
        i_chk(t,sgrid(%isn%)) = yes$(I(t,%isn%) lt 0); i_chk(t_last,sgrid(%isn%)) = no;
        abort$(sum((t,sgrid(%isn%)), i_chk(t,%isn%)) gt 0) 'Check the value function: Investment(I) should be non-negative';

* Iteration Log
        tulog(ct(t))                        = TU(t,"1","1","1","1","1","1","1","1");
        itlog("K_N",ct(t),sgrid(%isn%))     = K_N(%isn%);
        itlog("Y",ct(t),sgrid(%isn%))       = Y(t,%isn%);
        itlog("I",ct(t),sgrid(%isn%))       = I(t,%isn%);
        itlog("C",ct(t),sgrid(%isn%))       = C(t,%isn%);
        itlog("K",ct(t),sgrid(%isn%))       = K(%isn%);
        itlog("EF_K",ct(t),sgrid)           = EF_K(t+1,sgrid);
        itlog("MAT",ct(t),sgrid)            = MAT(t+1,sgrid)*mate;
        itlog("lev_mat",ct(t),sgrid(%isn%)) = lev_m(im)*mate ;
        itlog("MIU",ct(t),sgrid)            = MIU(t,sgrid);
        itlog("Ygross",ct(t),sgrid(%isn%))  = YGROSS(t,%isn%);
        itlog("E",ct(t),sgrid(%isn%))       = E(t,%isn%);
        itlog("lev_kt",ct(t),sgrid(%isn%))  = lev_kt(t,ik) ;
        itlog("delta_K",ct(t),sgrid)        = itlog("K",t,sgrid)-itlog("lev_kt",t,sgrid);
        itlog("delta_M",ct(t),sgrid(%isn%)) = itlog("MAT",t,sgrid)-lev_m(im)*mate;
        itlog("Inv_share",ct(t),sgrid)      = itlog("I",t,sgrid)/itlog("Y",t,sgrid);
        itlog2("PHIT",ct(t),sgrid,ig,ih)    = PHIT(t+1,sgrid,ig,ih);
        itlog("SCC",ct(t),sgrid(%isn%))     = ESTSCC(t,%isn%);
        itlog("CPRICE",ct(t),sgrid(%isn%))  = ((pback*(1-gback)**(t.val-1))*MIU.l(t,%isn%)**(expcost2-1));
$offdotl

errorlog("error_K",ct(t_exc_last5),sgrid)        = yes$(itlog("EF_K",t_exc_last5,sgrid) ge up("CAP") or itlog("EF_K",t_exc_last5,sgrid) le lo("CAP"));
errorlog("error_MAT",ct(t_exc_last5),sgrid)      = yes$(itlog("MAT",t_exc_last5,sgrid)/mate ge up("MAT") or itlog("MAT",t_exc_last5,sgrid)/mate le lo("MAT"));
errorlog2("error_PHIT",ct(t),sgrid,ig,ih)        = yes$(itlog2("PHIT",t,sgrid,ig,ih) ge up("PHI") or itlog2("PHIT",t,sgrid,ig,ih) le lo("PHI"));
errorlog3(ct(t)) = yes$((bellman_reduced.solvestat ne 1) or (bellman_reduced.modelstat ne 2));
abort$(sum((universal,t(t_exc_last5),sgrid), errorlog(universal,t,sgrid) + sum((ig,ih),errorlog2(universal,t,sgrid,ig,ih)) + errorlog3(t)) gt 0) 'Check the optimal solution';

option sysout   = off;
option solprint = off;

execute_unload 'coeff.gdx', coeff;
execute_unload 'itlog.gdx', itlog;
execute_unload 'itlog2.gdx', itlog2;
execute_unload 'err.gdx', errorlog;
execute_unload 'err2.gdx', errorlog2;
);

*-------------------------------------------------------------------------------------------
*       Simulation
*-------------------------------------------------------------------------------------------
* Set the true value for ECS
true_sens = 8;

$eval te_last (%niter%-6)
Sets    te(t)       Simulation time periods                   /1*%te_last%/
        sim_iter    Simulation iteration for stochasticity    /1*1000/ ;

parameter
        k0              Initial capital level                           /223/
        mat0            Initial carbon concentration                    /851/
        mu0             Initial carbon concentration in upper ocean     /460/
        ml0             Initial carbon concentration in lower ocean     /1740/
        phit0           Initial atmospheric temperature                 /.85/
        ocn0            Initial ocean temperature                       /.0068/

        PRE_EFK         Current effective capital level
        PRE_MAT         Current carbon concentration
        PRE_MU          Current carbon concentration in upper ocean
        PRE_ML          Current carbon concentration in lower ocean
        PRE_PHIT        Current atmospheric temperature
        PRE_OCN         Current ocean temperature
        PRE_BMIU        Prior belief on the Feedback factor
        PRE_S           Prior belief on the Feedback factor
        PRE_BSIG        Prior degree of belief
        PRE_T           Current Time period                             /1/
        PST_T           Next Time period
        coef_t(%isn%)   Coefficients of the current value function
        coef_tn(%isn%)  Coefficients of the next value function
        temp_shock      Random Temperature Shock                        /0/

        grad_k          Gradient of value function wrt capital
        grad_m          Gradient of value function wrt mat
        log_sim         Log for simulations
        log_avg         Log for average of simulations                   ;

nonnegative variables
         Sim_C           Consumption,
         emission        Emissions,
         RI              Investment,
         Sim_MIU         Carbon control level,
         PST_EFK         Effective Capital for next period,
         PST_OCN         Next period ocean change level,
         PST_MAT         Next period carbon concentration level,
         PST_PHIT(ig,ih) Possible temperature change in next period,
         OBS_PST_AT      Realized temperature change in next period,
         PST_MU          Next period carbon concentration level,
         PST_ML          Next period carbon concentration level,
         PST_S(ig,ih)    Posterior belief on the climate sensitivity,
         PST_BSIG        Posterior degree of belief (variance);
variable
         PST_BMIU(ig,ih) Posterior belief on the feedback factor,
         sim_OBJ         Objective;

$macro K_sim_CS(A)    (sum(pt_k$csp_K(jk,pt_K),cc_k(jk,pt_k) * power(A,ce_k(jk,pt_k))))
$macro MAT_sim_CS(A)  (sum(pt_m$csp_m(jm,pt_m),cc_m(jm,pt_m) * power(A,ce_m(jm,pt_m))))
$macro MU_sim_CS(A)   (sum(pt_u$csp_u(ju,pt_u),cc_u(ju,pt_u) * power(A,ce_u(ju,pt_u))))
$macro ML_sim_CS(A)   (sum(pt_l$csp_l(jl,pt_l),cc_l(jl,pt_l) * power(A,ce_l(jl,pt_l))))
$macro PHIT_sim_CS(A) (sum(pt_p$csp_p(jp,pt_p),cc_p(jp,pt_p) * power(A,ce_p(jp,pt_p))))
$macro OT_sim_CS(A)   (sum(pt_ot$csp_ot(jot,pt_ot),cc_ot(jot,pt_ot) * power(A,ce_ot(jot,pt_ot))))
$macro BMIU_sim_CS(A) (sum(pt_bm$csp_bm(jbm,pt_bm),cc_bm(jbm,pt_bm) * power(A,ce_bm(jbm,pt_bm))))
$macro BSIG_sim_CS(A) (sum(pt_bs$csp_bs(jbs,pt_bs),cc_bs(jbs,pt_bs) * power(A,ce_bs(jbs,pt_bs))))

$macro PRE_ST_PV  ( sum(ccgrid(%jsn%), coef_t(%jsn%) *\
                         K_sim_CS(PRE_K_N) * MAT_sim_CS(PRE_MAT_N) * MU_sim_CS(PRE_MU_N) * ML_sim_CS(PRE_ML_N) * PHIT_sim_CS(PRE_PHIT_N) * OT_sim_CS(PRE_OT_N) * BMIU_sim_CS(PRE_BMIU_N) * BSIG_sim_CS(PRE_BSIG_N)))

$macro PST_ST_PV  ( sum(ccgrid(%jsn%), coef_tn(%jsn%) *\
                         K_sim_CS(PST_K_N) * MAT_sim_CS(PST_MAT_N) * MU_sim_CS(PST_MU_N) * ML_sim_CS(PST_ML_N) * PHIT_sim_CS(PST_PHIT_N) * OT_sim_CS(PST_OT_N) * BMIU_sim_CS(PST_BMIU_N) * BSIG_sim_CS(PST_BSIG_N)))

$macro PRE_K_N    ((PRE_EFK-(lo("cap")+up("cap"))/2)/((up("cap")-lo("cap"))/2))
$macro PRE_MAT_N  ((PRE_MAT-(lo("mat")+up("mat"))/2)/((up("mat")-lo("mat"))/2))
$macro PRE_MU_N   ((PRE_MU-(lo("mu")+up("mu"))/2)/((up("mu")-lo("mu"))/2))
$macro PRE_ML_N   ((PRE_ML-(lo("ml")+up("ml"))/2)/((up("ml")-lo("ml"))/2))
$macro PRE_PHIT_N ((PRE_PHIT-(lo("phi")+up("phi"))/2)/((up("phi")-lo("phi"))/2))
$macro PRE_OT_N   ((PRE_OCN-(lo("OT")+up("OT"))/2)/((up("OT")-lo("OT"))/2))
$macro PRE_BMIU_N  ((PRE_S-(lo("BMIU")+up("BMIU"))/2)/((up("BMIU")-lo("BMIU"))/2))
$macro PRE_BSIG_N  ((PRE_BSIG-(lo("BSIG")+up("BSIG"))/2)/((up("BSIG")-lo("BSIG"))/2))

$macro PST_K_N    ((PST_EFK-(lo("cap")+up("cap"))/2)/((up("cap")-lo("cap"))/2))
$macro PST_MAT_N  ((PST_MAT-(lo("mat")+up("mat"))/2)/((up("mat")-lo("mat"))/2))
$macro PST_MU_N   ((PST_MU-(lo("mu")+up("mu"))/2)/((up("mu")-lo("mu"))/2))
$macro PST_ML_N   ((PST_ML-(lo("ml")+up("ml"))/2)/((up("ml")-lo("ml"))/2))
$macro PST_PHIT_N ((PST_PHIT(ig,ih)-(lo("phi")+up("phi"))/2)/((up("phi")-lo("phi"))/2))
$macro PST_OT_N   ((PST_OCN-(lo("OT")+up("OT"))/2)/((up("OT")-lo("OT"))/2))
$macro PST_BMIU_N  ((PST_S(ig,ih)-(lo("BMIU")+up("BMIU"))/2)/((up("BMIU")-lo("BMIU"))/2))
$macro PST_BSIG_N  ((PST_BSIG-(lo("BSIG")+up("BSIG"))/2)/((up("BSIG")-lo("BSIG"))/2))

$macro LBAR_TE   (sum(t,lbar_t(t)$(t.val eq PRE_T)))
$macro LBAR_TEN  (sum(t,lbar_t(t)$(t.val eq PST_T)))
$macro ABAR_TE   (sum(t,abar_t(t)$(t.val eq PRE_T)))
$macro ABAR_TEN  (sum(t,abar_t(t)$(t.val eq PST_T)))
$macro SBAR_TE   (sum(t,sbar_t(t)$(t.val eq PRE_T)))
$macro SBAR_TEN  (sum(t,sbar_t(t)$(t.val eq PST_T)))
$macro forcoth_TE (sum(t,forcoth(t)$(t.val eq PST_T)))
$macro etree_TE (eland0*(1-deland)**(PRE_T-1))
$macro cost1_TE  (pback*(1-gback)**(PRE_T-1)*SBAR_TE/expcost2/1000)

$macro PRE_rk       (PRE_EFK*ABAR_TE**(1/(1-gama)) *LBAR_TE)
$macro PST_rk       (PST_EFK*ABAR_TEN**(1/(1-gama))*LBAR_TEN)
$macro sim_ygross   (ABAR_TE*PRE_RK**gama*LBAR_TE**(1-gama))
$macro sim_y        (sim_Ygross*(1-(a2*PRE_PHIT**a3)) - sim_Ygross * COST1_TE * (Sim_MIU**expcost2))
$macro PST_force    (fco22x * (log(PST_MAT)/log(2)) + forcoth_TE)
$macro det_temp_sim (PRE_PHIT + c1 * (PST_force - fco22x/lambda*PRE_PHIT - c3*(PRE_PHIT-PRE_OCN)) )
$macro margin_sim   (c1*fco22x/lambda*PRE_PHIT)
$macro unc_bm       (max(min(sqrt(2*PRE_BSIG)*ghdata(ig,"zeta")+PRE_BMIU,up_bm),lo_bm))

$macro mymax(t1,t2,d) [0.5 * [t1 + t2 + sqrt(sqr(t1-t2) + sqr(d))] ]     #numerical max function
$macro mymin(t1,t2,d) [-0.5 * [-t1 - t2 + sqrt(sqr(-t1+t2) + sqr(d))] ]  #numerical min function

equations       sim_objective,ceq,eeq,ieq,otneq,matneq,muneq,mlneq,bmiuneq,bsigneq,obs_tempeq,phiteq,bseq;

sim_objective..     sim_OBJ =E= ((1-disc)*(Sim_C/LBAR_TE)**(1-elasmu)/(1-elasmu)*LBAR_TE + disc * sum((ig,ih), omega(ig)* omega(ih) * PST_ST_PV));
ceq..    Sim_c + RI =E= sim_y;
eeq..    emission =E= (1-Sim_MIU) * SBAR_TE * sim_ygross + etree_TE ;
ieq..    RI*5 =E= PST_rk - (1-dk)**5*PRE_rk ;
otneq..  PST_OCN =E= PRE_OCN + c4*(PRE_PHIT-PRE_OCN);
muneq..  PST_MU =E= PRE_MAT*b12 + PRE_MU*b22 ;
mlneq..  PST_ML =E= PRE_ML*b33 + PRE_MU*b23 ;
matneq.. PST_MAT =E= emission*5/3.666/mate + PRE_MAT*b11 + PRE_MU*b21 ;
bmiuneq(ig,ih).. PST_BMIU(ig,ih) =E= ((PST_PHIT(ig,ih)-det_temp_sim)*margin_sim*PRE_BSIG+PRE_BMIU*var_e)/((margin_sim**2)*PRE_BSIG+var_e); # Assume Bayesian Learning
bsigneq..        PST_BSIG =E= PRE_BSIG*var_e/((margin_sim**2)*PRE_BSIG+var_e); # Assume Bayesian Learning
bseq(ig,ih)..    PST_S(ig,ih) =E= mymax(lambda/(1-lo_bm), mymin(lambda/(1-up_bm), lambda/(1-PST_BMIU(ig,ih)), 1e-8), 1e-8);
phiteq(ig,ih)..  PST_PHIT(ig,ih) =E= PRE_PHIT + c1*(PST_force - fco22x/lambda*PRE_PHIT + fco22x/lambda*unc_bm*PRE_PHIT - c3*(PRE_PHIT-PRE_OCN)) + stoc_ts(ih) ;
Obs_tempeq..     OBS_PST_AT =E= PRE_PHIT + c1*(PST_force-fco22x/true_sens*PRE_PHIT-c3*(PRE_PHIT-PRE_OCN)) + temp_shock ;  #Realization of temperature change

model   sim_lsqr /sim_objective, ceq, eeq, ieq, otneq, matneq, muneq, mlneq, bmiuneq, bsigneq, phiteq, bseq, obs_tempeq  /;
option nlp=conopt4;
sim_lsqr.optfile = 1;
sim_lsqr.solprint = 2;
execseed = 1 + gmillisec(jnow);

loop(sim_iter,
*       Initial values and boundary
PRE_T     = 1;
PST_T     = 2;
PRE_EFK   = k0/LBAR_TE*ABAR_TE**(1/(gama-1));
PRE_OCN   = ocn0;
PRE_PHIT  = phit0;
PRE_MAT   = mat0/mate;
PRE_MU    = mu0/mate;
PRE_ML    = ml0/mate;
PRE_BSIG  = lev_bs_init;
PRE_S     = lev_s_init;
PRE_BMIU  = 1-lambda/PRE_S;

PST_EFK.L = PRE_EFK*1.01;
PST_MAT.L = PRE_MAT*1.01;
PST_MU.L  = PRE_MU*1.01;
PST_ML.L  = PRE_ML*1.01;
PST_PHIT.L(ig,ih) = PRE_PHIT*1.01;
OBS_PST_AT.L      = PRE_PHIT*1.01;
Sim_MIU.L = 0.1;

PST_EFK.UP = up("cap");
PST_MAT.UP = up("mat");
PST_MU.UP  = up("mu");
PST_ML.UP  = up("ml");
PST_PHIT.UP(ig,ih) = up("phi");
OBS_PST_AT.UP = up("phi");
PST_OCN.UP = up("ot");
Sim_MIU.UP = 1;

PST_EFK.LO = lo("cap");
PST_MAT.LO = lo("mat");
PST_MU.LO  = lo("mu");
PST_ML.LO  = lo("ml");
PST_OCN.LO = lo("ot");
PST_PHIT.LO(ig,ih) = lo("phi");
OBS_PST_AT.LO = lo("phi");
Sim_MIU.LO = 0;
Sim_C.LO   = 1e-6;
RI.LO      = 1e-6;

loop(te,
         PST_T = PRE_T + 1 ;
         coef_t(%isn%) = coeff(te,%isn%) ;
         coef_tn(%isn%) = coeff(te+1,%isn%) ;
         temp_shock = normal(0,sqrt(var_e));

         solve sim_lsqr using nlp maximizing sim_OBJ;
         abort$(sim_lsqr.solvestat>1) "Lsqr terminates with abnormal solver status";
         abort$(sim_lsqr.modelstat>2) "Lsqr terminates with abnormal model status";

$ondotl
         log_sim(sim_iter,te,"Y(net)") = Sim_Y;
         log_sim(sim_iter,te,"Ygross") = Sim_Ygross ;
         log_sim(sim_iter,te,"V1")  = PRE_ST_PV;
         log_sim(sim_iter,te,"V2")  = (1-disc)*(Sim_C/LBAR_TE)**(1-elasmu)/(1-elasmu)*LBAR_TE + disc * sum((ig,ih), omega(ig)*omega(ih)* PST_ST_PV);
         log_sim(sim_iter,te,"Present Utility") = (1-disc)*(Sim_C/LBAR_TE)**(1-elasmu)/(1-elasmu)*LBAR_TE ;
         log_sim(sim_iter,te,"OBJL")= abs(log_sim(sim_iter,te,"V1") - log_sim(sim_iter,te,"V2")) ;
         log_sim(sim_iter,te,"C")   = Sim_C ;
         log_sim(sim_iter,te,"I")   = RI ;
         log_sim(sim_iter,te,"MIU") = Sim_MIU;
         log_sim(sim_iter,te,"CPRICE") = pback*(1-gback)**(PRE_T-1)*Sim_MIU**(expcost2-1) ;
         log_sim(sim_iter,te,"Emit")= Emission ;
         log_sim(sim_iter,te,"Dam_Frac") = a2*PRE_PHIT**a3;
         log_sim(sim_iter,te,"Damage") = sim_Ygross*a2*PRE_PHIT**a3 ;
         log_sim(sim_iter,te,"Cost")= sim_Ygross * COST1_TE * (Sim_MIU**expcost2);
         log_sim(sim_iter,te,"L")   = LBAR_TE;
         log_sim(sim_iter,te,"TFP") = ABAR_TE;
         log_sim(sim_iter,te,"Sig") = SBAR_TE;
         log_sim(sim_iter,te,"forcing_next") = PST_force;
         log_sim(sim_iter,te,"forcoth") = forcoth_TE;
         log_sim(sim_iter,te,"K")   = PRE_RK ;
         log_sim(sim_iter,te,"Year")= PRE_T*5+2010 ;
         log_sim(sim_iter,te,"OCN") = PRE_OCN ;
         log_sim(sim_iter,te,"MAT") = PRE_MAT*mate ;
         log_sim(sim_iter,te,"MU")  = PRE_MU*mate ;
         log_sim(sim_iter,te,"ML")  = PRE_ML*mate ;
         log_sim(sim_iter,te,"BMIU")= PRE_BMIU;
         log_sim(sim_iter,te,"BSIG")= PRE_BSIG;
         log_sim(sim_iter,te,"True_Sensitivity") = true_sens;
         log_sim(sim_iter,te,"ECS") = PRE_S;
         log_sim(sim_iter,te,"Temp_Shock") = temp_shock;
         log_sim(sim_iter,te,"Actual AT")  = PRE_PHIT ;
         log_sim(sim_iter,te+1,"Anticipated AT") = sum((ig,ih), omega(ig)*omega(ih)*PST_PHIT.L(ig,ih));
         log_sim(sim_iter,te,"Effective_K")= PRE_EFK;
         log_sim(sim_iter,te,"bnd_err_UP") = Yes$(PRE_MAT gt smax(im,lev_m(im)) or PRE_EFK gt smax(ik,lev_k(ik)) or PRE_PHIT gt smax(ip,lev_p(ip)) or PRE_MU gt smax(iu,lev_mu(iu)) or PRE_ML gt smax(il,lev_ml(il)) or PRE_BMIU gt smax(ibm,lev_bm(ibm)) or PRE_BSIG gt smax(ibs,lev_bs(ibs)));
         log_sim(sim_iter,te,"bnd_err_LO") = Yes$(PRE_MAT lt smin(im,lev_m(im)) or PRE_EFK lt smin(ik,lev_k(ik)) or PRE_PHIT lt smin(ip,lev_p(ip)) or PRE_MU lt smin(iu,lev_mu(iu)) or PRE_ML lt smin(il,lev_ml(il)) or PRE_BMIU lt smin(ibm,lev_bm(ibm)) or PRE_BSIG lt smin(ibs,lev_bs(ibs)));

*        Calculation of SCC
         grad_k =   1/((up("cap") - lo("cap"))/2) * (
                        sum(cgrid(%jsn%), coef_t(%jsn%) *
                          sum(km,MAT_sim_CS(PRE_MAT_N)$(km.val eq jm.val)) *
                          sum(ku,MU_sim_CS(PRE_MU_N)$(ku.val eq ju.val)) *
                          sum(kl,ML_sim_CS(PRE_ML_N)$(kl.val eq jl.val)) *
                          sum(kp,PHIT_sim_CS(PRE_PHIT_N)$(kp.val eq jp.val))  *
                          sum(kot,OT_sim_CS(PRE_OT_N)$(kot.val eq jot.val))  *
                          sum(kbm,BMIU_sim_CS(PRE_BMIU_N)$(kbm.val eq jbm.val))  *
                          sum(kbs,BSIG_sim_CS(PRE_BSIG_N)$(kbs.val eq jbs.val)) *
                          sum(kk$(kk.val eq jk.val),
                                sum(pt_k$csp_k(kk,pt_k),
                                  (cc_k(kk,pt_k)*ce_k(kk,pt_k))$(ce_k(kk,pt_k) ge 1) *
                                  (power(PRE_K_N,ce_k(kk,pt_k)-1)$(ce_k(kk,pt_k) ge 1) +
                                   1$(not ce_k(kk,pt_k) eq 1))))));

         grad_m =   1/((up("mat") - lo("mat"))/2) * (
                        sum(cgrid(%jsn%), coef_t(%jsn%) *
                          sum(kk,K_sim_CS(PRE_K_N)$(kk.val eq jk.val)) *
                          sum(ku,MU_sim_CS(PRE_MU_N)$(ku.val eq ju.val)) *
                          sum(kl,ML_sim_CS(PRE_ML_N)$(kl.val eq jl.val)) *
                          sum(kp,PHIT_sim_CS(PRE_PHIT_N)$(kp.val eq jp.val))  *
                          sum(kot,OT_sim_CS(PRE_OT_N)$(kot.val eq jot.val))  *
                          sum(kbm,BMIU_sim_CS(PRE_BMIU_N)$(kbm.val eq jbm.val))  *
                          sum(kbs,BSIG_sim_CS(PRE_BSIG_N)$(kbs.val eq jbs.val)) *
                          sum(km$(km.val eq jm.val),
                                sum(pt_m$csp_m(km,pt_m),
                                  (cc_m(km,pt_m)*ce_m(km,pt_m))$(ce_m(km,pt_m) ge 1) *
                                  (power(PRE_MAT_N,ce_m(km,pt_m)-1)$(ce_m(km,pt_m) ge 1) +
                                   1$(not ce_m(km,pt_m) eq 1))))));

         grad_k=grad_k/LBAR_TE*ABAR_TE**(1/(gama-1));
         grad_m=grad_m/3.666/mate;

         log_sim(sim_iter,te,"SCC") = -1000 * eeq.m/ceq.m;
         log_sim(sim_iter,te,"SCC_val") = -1000 * grad_m/grad_k;
         log_sim(sim_iter,te,"gradk") = grad_k;
         log_sim(sim_iter,te,"gradm") = grad_m;

*        Update variables for the next period
         PRE_BMIU = (max(lo_bm,min(up_bm,((OBS_PST_AT.L-det_temp_sim)*margin_sim*PRE_BSIG+PRE_BMIU*var_e)/((margin_sim**2)*PRE_BSIG+var_e)))); # Assume Bayesian Learning
         PRE_S    = lambda/(1-PRE_BMIU);
         PRE_EFK  = PST_EFK.L;
         PRE_OCN  = PST_OCN.L;
         PRE_MAT  = PST_MAT.L;
         PRE_PHIT = OBS_PST_AT.L;
         PRE_T    = PST_T;
         PRE_MU   = PST_MU.L + PRE_ML*b32;
         PRE_ML   = PST_ML.L;
         PRE_BSIG = PST_BSIG.L;

$offdotl
);
);
log_avg(te,universal)=sum(sim_iter, log_sim(sim_iter,te,universal))/card(sim_iter);

execute_unload 'log_sim.gdx', log_sim;
