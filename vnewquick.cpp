
#include <iostream>
# include <cmath>

using namespace std;
const int N=102;
const double Gr = 1e4;

void vnewquick( double vnew[N][N],  double v[N][N], double p[N][N], double fe[N][N],double fw[N][N],double fn[N][N],double fs[N][N], double theta[N][N], float dx, float dy, float volp, float dt, int Nreal)
{
   int i,j,k, gs=0;
    int arr[]={1,N-2};
     float Rp, ap, rms;
     double ffe,ffw,ffn,ffs,ff,tf,ve,vw,vn,vs;

        do
        {
        rms=0;
        //interior cells quick scheme
   for( i=2;i<N-2;++i)
    {
    for(j=2;j<N-2;++j)
        {
Rp = (volp/dt)*(v[i][j]-vnew[i][j])+(volp/sqrt(Gr))*((vnew[i+1][j]-2*vnew[i][j]+vnew[i-1][j])/(dx*dx)+(vnew[i][j+1]-2*vnew[i][j]+vnew[i][j-1])/(dy*dy)) - (p[i][j+1]*dx-p[i][j-1]*dx)/2.0 + theta[i][j]*volp;

if (fe[i][j]>=0 )
   {
     ffe=0.75*fe[i][j]; ve = (0.75)*vnew[i][j]+(0.375)*vnew[i+1][j]-(0.125)*vnew[i-1][j];
   }
 else
     {
    ffe=0.375*fe[i][j];  ve = (0.75)*vnew[i+1][j]+(0.375)*vnew[i][j]-(0.125)*vnew[i+2][j];
     }
 if (fw[i][j]>=0 )
    {
         ffw=0.75*fw[i][j]; vw = (0.75)*vnew[i][j]+(0.375)*vnew[i-1][j]-(0.125)*vnew[i+1][j];
    }
 else
    {
         ffw=0.375*fw[i][j];  vw = (0.75)*vnew[i-1][j]+(0.375)*vnew[i][j]-(0.125)*vnew[i-2][j];
    }
 if (fn[i][j]>=0 )
   {
      ffn=0.75*fn[i][j]; vn = (0.75)*vnew[i][j]+(0.375)*vnew[i][j+1]-(0.125)*vnew[i][j-1];
   }
 else
    {
       ffn=0.375*fn[i][j]; vn = (0.75)*vnew[i][j+1]+(0.375)*vnew[i][j]-(0.125)*vnew[i][j+2];
    }
 if (fs[i][j]>=0 )
    {
       ffs=0.75*fs[i][j];  vs = (0.75)*vnew[i][j]+(0.375)*vnew[i][j-1]-(0.125)*vnew[i][j+1];
    }
 else
  {
      ffs=0.375*fs[i][j];  vs = (0.75)*vnew[i][j-1]+(0.375)*vnew[i][j]-(0.125)*vnew[i][j-2];
  }
     tf=fe[i][j]*ve + fw[i][j]*vw + fn[i][j]*vn + fs[i][j]*vs;
     ff=ffe+ffw+ffn+ffs;
    Rp = Rp - tf;
ap = (volp/dt) + (2.0*volp/sqrt(Gr))*(1.0/(dx*dx)+1.0/(dy*dy)) + ff;
vnew[i][j]=vnew[i][j] + Rp/ap;
rms=rms+Rp*Rp;
}
    }
//top and bottom cells// upwind scheme
    for (k=0;k<2;++k){
    j=arr[k];
for (i=1;i<N-1;++i)
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



//left and right cells  // upwind scheme
    for (k=0;k<2;++k){
    i=arr[k];
for (j=2;j<N-2;++j){
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

//boundary conditions ---- fictitious cells

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




