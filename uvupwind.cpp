#include <iostream>
# include <cmath>

using namespace std;
const int N=102;
const double Re=100.0;
const double Gr = 1e4;

void uvupwind( double u[N][N],  double up[N][N], double fe[N][N],double fw[N][N],double fn[N][N],double fs[N][N], float dx, float dy, float volp, float dt, int Nreal)
{
   int i,j, gs=0;
     float Rp, ap, rms;
     double ffe,ffw,ffn,ffs,ff,tf,ue,uw,un,us;

        do
        {
        rms=0;
   for( i=1;i<N-1;++i)
    {
    for(j=1;j<N-1;++j)
        {
Rp = (volp/dt)*(u[i][j]-up[i][j])+(volp/sqrt(Gr))*((up[i+1][j]-2*up[i][j]+up[i-1][j])/(dx*dx)+(up[i][j+1]-2*up[i][j]+up[i][j-1])/(dy*dy));

if (fe[i][j]>=0 )
   {
     ffe=fe[i][j]; ue=up[i][j];
   }
 else
     {
    ffe=0; ue=up[i+1][j];
     }
 if (fw[i][j]>=0 )
    {
       ffw=fw[i][j]; uw=up[i][j];
    }
 else
    {
        ffw=0;  uw=up[i-1][j];
    }
 if (fn[i][j]>=0 )
   {
      ffn=fn[i][j]; un=up[i][j];
   }
 else
    {
        ffn=0; un=up[i][j+1];
    }
 if (fs[i][j]>=0 )
    {
      ffs=fs[i][j]; us=up[i][j];
    }
 else
  {
      ffs=0; us=up[i][j-1];
  }
     tf=fe[i][j]*ue + fw[i][j]*uw + fn[i][j]*un + fs[i][j]*us;
     ff=ffe+ffw+ffn+ffs;
    Rp = Rp - tf;
ap = (volp/dt) + (2.0*volp/sqrt(Gr))*(1.0/(dx*dx)+1.0/(dy*dy)) + ff;
up[i][j]=up[i][j] + Rp/ap;
rms=rms+Rp*Rp;
}
    }

//boundary conditions ---- fictitious cells

    for (j=1;j<N-1;++j)    //left wall
        up[0][j] = -up[1][j];


    for (j=1;j<N-1;++j)      // right wall
        up[N-1][j] = -up[N-2][j];


    for (j=1;j<N-1;++j)      // bottom wall
        up[j][0] = -up[j][1];


    for (j=1;j<N-1;++j)    // top wall
        up[j][N-1] = -up[j][N-2];

 rms=sqrt(rms/Nreal);
   gs=gs+1;
        }
    while (rms >1e-6);
}


