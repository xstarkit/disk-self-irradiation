import graph;
import palette;



bounds shade_image(picture pic=currentpicture, pair[] z, real[] f,
             range range=Full, pen[] palette)
{
  if(z.length != f.length)
    abort("z and f arrays have different lengths");

  real m=min(f);
  real M=max(f);
  bounds bounds=range(pic,m,M);
  real rmin=pic.scale.z.T(bounds.min);
  real rmax=pic.scale.z.T(bounds.max);

  palette=adjust(pic,m,M,rmin,rmax,palette);

  // Crop data to allowed range and scale
  if(range != Full || pic.scale.z.scale.T != identity ||
     pic.scale.z.postscale.T != identity) {
    scalefcn T=pic.scale.z.T;
    real m=bounds.min;
    real M=bounds.max;
    f=map(new real(real x) {return T(min(max(x,m),M));},f);
  }

  int[] edges={0,0,1};
  int N=palette.length-1;

  int[][] trn=triangulate(z);
  real step=rmax == rmin ? 0.0 : N/(rmax-rmin);
  for(int i=0; i < trn.length; ++i) {
    int[] trni=trn[i];
    int i0=trni[0], i1=trni[1], i2=trni[2];
    pen color(int i) {return palette[round((f[i]-rmin)*step)];}
    gouraudshade(pic,z[i0]--z[i1]--z[i2]--cycle,
                 new pen[] {color(i0),color(i1),color(i2)},edges);
  }
  return bounds; // Return bounds used for color space
}



string filename = "../a000/l010/0050/stage1.15i_ang";

file f = input(filename).line();
real[][] data = f.dimension(0,0);
close(f);

pair[] xy;
real[] f;


for (int i=0; i<data.length; ++i) {
    real theta = data[i][0]/180*pi;  // read first value (value of theta angle)
    data[i].delete(0);               // remove first value from array
    pair[] theta_xy = sequence(
        new pair(int j) {
            real phi = (j+0.5) * 2*pi/data[i].length;
            real r   = theta/(pi/2);
            return (500*r*cos(phi), 500*r*sin(phi));
        }, 
        data[i].length
    );
    real[] theta_f = sequence(
        new real(int j) {
            return (data[i][j]<=0.0) ? data[i][j] : log10(data[i][j]);
        }, 
        data[i].length
    );

    xy.append(theta_xy);
    f.append(theta_f);
}


//write(data.length);
//write(xy.length);


write(xy.length);

while (xy.length > 1000) {
    int[][] trn=triangulate(xy);

    // step 1 - find triangle with least difference
    int[] toremove;
    for (int i=0; i<trn.length; ++i) {
        int i0=trn[i][0], i1=trn[i][1], i2=trn[i][2];
        real diff = max(f[i0], f[i1], f[i2]);
        if (diff == 0.0) {
            pair p12 = 0.5(xy[i1]+xy[i2]);
            pair t = 1/3*xy[i0] + 2/3*p12;
            real val = f[i0];
            //xy.push(t); f.push(val);
            toremove.push(i0);
            toremove.push(i1);
            toremove.push(i2);
        }
    }
    write(xy.length);

    toremove = sort(toremove);
    int lastremoved = -1;
    for (int i=toremove.length-1; i>=0; --i) {
        if (lastremoved == toremove[i]) continue;
        xy.delete(toremove[i]);
        f.delete(toremove[i]);
        lastremoved = toremove[i];
    }
    write(xy.length);

//abort("iter 0");
//continue;
break;

    // step 1 - find triangle with least difference
    real mindiff = 1e200;
    int imindiff = -1;
    int[][] trn=triangulate(xy);
    for(int i=0; i<trn.length; ++i) {
        int i0=trn[i][0], i1=trn[i][1], i2=trn[i][2];
        real diff = max(f[i0], f[i1], f[i2]);
        if (diff < mindiff) {
            mindiff = diff;
            imindiff = i;
        }
    }

    // step 2 - remove one point from the array
/*    
    int[] trni=trn[imindiff];
    int i0=trni[0], i1=trni[1], i2=trni[2];
    
    int iremove;
    if ((length(xy[i0])>=length(xy[i1])) && (length(xy[i0])>=length(xy[i2]))) {
        iremove = (abs(f[i0]-f[i1]) < abs(f[i0]-f[i2])) ? 1 : 2;
    } else
    if ((length(xy[i1])>=length(xy[i0])) && (length(xy[i1])>=length(xy[i2]))) {
        iremove = (abs(f[i1]-f[i0]) < abs(f[i1]-f[i2])) ? 0 : 2;
    } else
    if ((length(xy[i2])>=length(xy[i0])) && (length(xy[i2])>=length(xy[i1]))) {
        iremove = (abs(f[i2]-f[i0]) < abs(f[i2]-f[i1])) ? 0 : 1;
    }
*/
    int[] trni=trn[imindiff];
    int i0=trni[0], i1=trni[1], i2=trni[2];
    write(imindiff, mindiff);

    int iremove;
    real d01 = length(xy[i0]-xy[i1]);
    real d02 = length(xy[i0]-xy[i2]);
    real d12 = length(xy[i1]-xy[i2]);
    if ((d01<=d02) && (d01<=d12)) {
        iremove = (d02<d12) ? 0 : 1;
    } else
    if ((d02<=d12) && (d02<=d01)) {
        iremove = (d01<d12) ? 0 : 2;
    } else
    if ((d12<=d01) && (d12<=d02)) {
        iremove = (d01<d02) ? 1 : 2;
    }

    xy.delete(iremove);
    f.delete(iremove);
write(xy.length);
} 

//write(xy.length);






bounds range = shade_image(xy, f, Full, Rainbow());
//palette(bar,"$A$",range,(0,0),(0.5cm,8cm),Right,Palette,
//        PaletteTicks("$%+#.1f$"));



