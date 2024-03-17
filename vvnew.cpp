#include <iostream>
# include <cmath>

using namespace std;
const int N=102;
const double Gr = 1e4;

// corrected velocities

void vvnew( double vnew[N][N], double v[N][N], double p[N][N], double fe[N][N],double fw[N][N],double fn[N][N],double fs[N][N],double theta[N][N], float dx, float dy, float volp, float dt, int Nreal)
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
Rp = (volp/dt)*(v[i][j]-vnew[i][j])+(volp/sqrt(Gr))*((vnew[i+1][j]-2*vnew[i][j]+vnew[i-1][j])/(dx*dx)+(vnew[i][j+1]-2*vnew[i][j]+vnew[i][j-1])/(dy*dy)) - (p[i][j+1]*dx-p[i][j-1]*dx)/2.0 + theta[i][j]*volp;

if (fe[i][j]>=0 )
   {
     ffe=fe[i][j]; ve=vnew[i][j];
   }
 else
     {
    ffe=0; ve=vnew[i+1][j];
     }
 if (fw[i][j]>=0 )
    {
       ffw=fw[i][j]; vw=vnew[i][j];
    }
 else
    {
        ffw=0;  vw=vnew[i-1][j];
    }
 if (fn[i][j]>=0 )
   {
      ffn=fn[i][j]; vn=vnew[i][j];
   }
 else
    {
        ffn=0; vn=vnew[i][j+1];

    }
 if (fs[i][j]>=0 )
    {
      ffs=fs[i][j]; vs=vnew[i][j];
    }
 else
  {
      ffs=0; vs=vnew[i][j-1];
  }
     tf=fe[i][j]*ve + fw[i][j]*vw + fn[i][j]*vn + fs[i][j]*vs;
     ff=ffe+ffw+ffn+ffs;
    Rp = Rp - tf;
ap = (volp/dt) + (2.0*volp/sqrt(Gr))*(1.0/(dx*dx)+1.0/(dy*dy)) + ff;
vnew[i][j]=vnew[i][j] + Rp/ap;
rms=rms+Rp*Rp;

        }
    }

    for (j=1;j<N-1;++j)    //left wall
        vnew[0][j] = -vnew[1][j];


    for (j=1;j<N-1;++j)      // right wall
        vnew[N-1][j] = -vnew[N-2][j];


    for (j=1;j<N-1;++j)      // bottom wall
        vnew[j][0] = -vnew[j][1];


    for (j=1;j<N-1;++j)    // top wall
        vnew[j][N-1] = -vnew[j][N-2];


    rms=sqrt(rms/Nreal);
   gs=gs+1;
        }
    while (rms >1e-6);

}
