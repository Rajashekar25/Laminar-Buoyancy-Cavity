
#include <iostream>
# include <cmath>

using namespace std;
const int N=102;
//const double Re=100.0;
const double Gr = 1e4;

void vupwind( double v[N][N], double vp[N][N],double fe[N][N],double fw[N][N],double fn[N][N],double fs[N][N], double theta[N][N], float dx, float dy, float volp, float dt, int Nreal)
{
     int i,j, gs=0;
     float Rp, ap, rms;
     double ffe,ffw,ffn,ffs,ff,tf,ve,vw,vn,vs;

        do
        {
        rms=0;
   for( i=1;i<N-1;++i)
    {
    for(j=1;j<N-1;++j)
        {
Rp = (volp/dt)*(v[i][j]-vp[i][j])+(volp/sqrt(Gr))*((vp[i+1][j]-2*vp[i][j]+vp[i-1][j])/(dx*dx)+(vp[i][j+1]-2*vp[i][j]+vp[i][j-1])/(dy*dy))+theta[i][j]*volp;

if (fe[i][j]>=0 )
   {
     ffe=fe[i][j]; ve=vp[i][j];
   }
 else
     {
    ffe=0; ve=vp[i+1][j];
     }
 if (fw[i][j]>=0 )
    {
       ffw=fw[i][j]; vw=vp[i][j];
    }
 else
    {
        ffw=0;  vw=vp[i-1][j];
    }
 if (fn[i][j]>=0 )
   {
      ffn=fn[i][j]; vn=vp[i][j];
   }
 else
    {
        ffn=0; vn=vp[i][j+1];

    }
 if (fs[i][j]>=0 )
    {
      ffs=fs[i][j]; vs=vp[i][j];
    }
 else
  {
      ffs=0; vs=vp[i][j-1];
  }
    tf=fe[i][j]*ve + fw[i][j]*vw + fn[i][j]*vn + fs[i][j]*vs;
     ff=ffe+ffw+ffn+ffs;
    Rp = Rp - tf;
ap = (volp/dt) + (2.0*volp/sqrt(Gr))*(1.0/(dx*dx)+1.0/(dy*dy)) + ff;
vp[i][j]=vp[i][j] + Rp/ap;
rms=rms+Rp*Rp;
}
    }

//boundary conditions ---- fictitious cells

    for (j=1;j<N-1;++j)    //left wall
        vp[0][j] = -vp[1][j];

    for (j=1;j<N-1;++j)      // right wall
        vp[N-1][j] = -vp[N-2][j];

    for (j=1;j<N-1;++j)      // bottom wall
        vp[j][0] = -vp[j][1];

    for (j=1;j<N-1;++j)    // top wall
        vp[j][N-1] = -vp[j][N-2];

    rms=sqrt(rms/Nreal);
   gs=gs+1;
        }
    while (rms >1e-6);

}
