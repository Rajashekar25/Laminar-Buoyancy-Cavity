#include <iostream>
# include <cmath>

using namespace std;
const int N=102;
const double Gr = 1e4;
const float pr=1;

// solves energy equation in semi implicit form

void Energy( double thetanew[N][N], double theta[N][N], double unew[N][N], double vnew[N][N], double fe[N][N],double fw[N][N],double fn[N][N],double fs[N][N], float dx, float dy, float volp, float dt,int Nreal)
{
   int i,j, gs=0;
     float Rp, ap, rms;
     double ffe,ffw,ffn,ffs,ff,tf,te,tw,tn,ts;

        do
        {
        rms=0;
   for( i=1;i<N-1;++i)
    {
    for(j=1;j<N-1;++j)
        {
Rp = (volp/dt)*(theta[i][j]-thetanew[i][j])+(volp/(pr*sqrt(Gr)))*((thetanew[i+1][j]-2*thetanew[i][j]+thetanew[i-1][j])/(dx*dx)+(thetanew[i][j+1]-2*thetanew[i][j]+thetanew[i][j-1])/(dy*dy)) ;

if (fe[i][j]>=0 )
   {
     ffe=fe[i][j]; te=thetanew[i][j];
   }
 else
     {
    ffe=0; te=thetanew[i+1][j];
     }
 if (fw[i][j]>=0 )
    {
       ffw=fw[i][j]; tw=thetanew[i][j];
    }
 else
    {
        ffw=0;  tw=thetanew[i-1][j];
    }
 if (fn[i][j]>=0 )
   {
      ffn=fn[i][j]; tn=thetanew[i][j];
   }
 else
    {
        ffn=0; tn=thetanew[i][j+1];

    }
 if (fs[i][j]>=0 )
    {
      ffs=fs[i][j]; ts=thetanew[i][j];
    }
 else
  {
      ffs=0; ts=thetanew[i][j-1];
  }
     tf = fe[i][j]*te + fw[i][j]*tw + fn[i][j]*tn + fs[i][j]*ts;
     ff = ffe+ffw+ffn+ffs;
    Rp = Rp - tf;
ap = (volp/dt) + (2.0*volp/(pr*sqrt(Gr)))*(1.0/(dx*dx)+1.0/(dy*dy)) + ff;
thetanew[i][j] = thetanew[i][j] + Rp/ap;
rms=rms+Rp*Rp;
        }
    }

//boundary conditions ---- fictitious cells

    for (j=1;j<N-1;++j)    //left wall
        thetanew[0][j] = 1.0 -thetanew[1][j];


    for (j=1;j<N-1;++j)      // right wall
        thetanew[N-1][j] = -1.0 -thetanew[N-2][j];


    for (j=1;j<N-1;++j)      // bottom wall
        thetanew[j][0] = thetanew[j][1];


    for (j=1;j<N-1;++j)    // top wall
        thetanew[j][N-1] = thetanew[j][N-2];


   rms=sqrt(rms/Nreal);
   gs=gs+1;
        }
    while (rms >1e-6);
}

