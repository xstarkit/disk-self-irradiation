import graph;



void read_spectrum(string filename, real[] E, real[] I)
{
    file f = input(filename).line();
    real[][] raw = f.dimension(0,0);
    close(f);

    // if row[0] contains only one number, we have module 15 and we skip first 4 lines
    // if row[0] contains two numbers, we have module 15i and we skip first 2 lines
    int row_offset = (raw[0].length == 1) ? 3 : 2;

    real[][] tmp = transpose(raw[row_offset:]);

    E.delete();
    E.append(tmp[0]*4.135667e-18); // 4.135667e-18 = Hz->keV conversion factor (planck_h/(1e3*electronvolt))
    I.delete();
    I.append(tmp[1]);
}




// --- settings ----------------

//real mass = 10;
//real spin = 0.65;
//real alpha = 0.100;
//real mdot[] = {0.1, 0.3, 0.6};
//pen  colors[] = {black,deepblue,heavyred}; if (settings.gray) colors = new pen[] {black,black,black};
//real inc  = 66;
real xmin = 0.01;
real xmax = 30.0;
real ymin = -infinity;
real ymax = infinity;
real size_factor = 1.0;
real size_aspect = 3/2;
defaultpen(currentpen+linewidth(1pt));
usersetting();


pen  pen_original = black+linetype("3 3");
pen  pen_retrad   = black+linetype("0 3");
pen  pen_final    = black+solid;


//------ plot graph --------------

picture plot_spectrum() {
    picture pic;
    size(pic, 80mm*size_factor*size_aspect, 80mm*size_factor, IgnoreAspect);
    scale(pic, Log, Log);

    real[] E;
    real[] I;

    read_spectrum("../a000/l010/0013/stage1.15", E, I);
    draw(pic, graph(pic, E, I), pen_original, "original");

    read_spectrum("../a000/l010/0013/stage1.15i", E, I);
    draw(pic, graph(pic, E, I), pen_retrad, "retrad");

    read_spectrum("../a000/l010/0013/stage2.15", E, I);
    draw(pic, graph(pic, E, I), pen_final, "final");
    real ymax = 10^(round(log10(max(I)))+0.5);
    real ymin = ymax*1e-5;

    // label flux curves
//    label(pic, rotate(-62)*"$L=0.1$", (log10(9.5),log10(5*ymin)), W, Helvetica()+fontsize(9pt), NoFill);
//    label(pic, rotate(-62)*"$L=0.3$", (log10(13.5),log10(5*ymin)), W, Helvetica()+fontsize(9pt), NoFill);
//    label(pic, rotate(-62)*"$L=0.6$", (log10(21),log10(5*ymin)), W, Helvetica()+fontsize(9pt), NoFill);

    xlimits(pic, xmin, xmax);
    ylimits(pic, ymin, ymax);
    crop(pic);

    xaxis(pic, "energy [keV]", BottomTop, LeftTicks(DefaultFormat), above=true);
    yaxis(pic, "flux [erg/cm$^2$/keV/s]", LeftRight, RightTicks, above=true);
/*
    // attach legend
    legendlinelength = 40;
    legendmargin = 6;
    Legend l1 = Legend(" NT", black+linewidth(1pt)+dashed+Helvetica()+fontsize(9pt));
    Legend l2 = Legend(" SD", black+linewidth(1pt)+solid+Helvetica()+fontsize(9pt));
    pic.legend.push(l1);
    pic.legend.push(l2);
    attach(pic, legend(pic, black+linewidth(.5pt)),point(pic,NE),15S+15W,Fill(gray(0.98)));

    // labels for spin & mass
    label(pic, "$M\!=\!"+string(mass)+" M_\odot,\;a\!=\!"+string(spin)+"$", point(pic,NE), 25S+13.6W, fontsize(8pt), NoFill);
    label(pic, "$\alpha="+string(alpha)+",\;{\rm inc}\!=\!"+string(inc)+"^\circ$", point(pic,NE), 30S+13W, fontsize(8pt), NoFill);
    //label(pic, "$h_f=\rm{BHSPEC}$", point(pic,NE), 35S+14W, fontsize(8pt), NoFill);
    label(pic, "$h_f=1.6$", point(pic,NE), 35S+11.5W, fontsize(8pt), NoFill);
*/

    return pic;
}


picture P1 = plot_spectrum();
add(P1.fit(),(0,0),N);

