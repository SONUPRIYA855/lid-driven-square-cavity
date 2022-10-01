#include<stdio.h>
#include<math.h>
int main()
{
FILE *fp,*U,*V;
fp=fopen("stream.dat","w");
U=fopen("uv_velocity.dat","w");
double SUM,s,error,wold[129][129],siold[129][129],wnew[129][129],sinew[129][129],u[129][129],
v[129][129],
sqrdelx,sqrdely,delx=0.0078125,dely=0.0078125,x;
int i,j,k,count=0;
sqrdelx=delx*delx;
sqrdely=sqrdelx;
for(j=0;j<129;j++)
{
for(i=0;i<129;i++)
{
siold[i][j]=0;
wold[i][j]=0;
}
}
for(i=1;i<128;i++)
{
wold[i][128]=0-(2*128);
 }
for(j=0;j<129;j++)
{
for(i=0;i<129;i++)
{
sinew[i][j]=siold[i][j];
wnew[i][j]=wold[i][j];
}
}
//update sin
do
{
//start :
count++;

for(k=0;k<10;k++)
{
for(j=1;j<128;j++)
{
for(i=1;i<128;i++)
{
 sinew[i][j]=(sinew[i+1][j]+sinew[i-1][j]+sinew[i][j+1]+sinew[i][j-1]+(sqrdelx*wnew[i][j]))/4.0;
}
}
}
//update boundary condition
for(j=1;j<128;j++)
{
wnew[0][j]=(0-(2*sinew[1][j]))/sqrdelx;
wnew[128][j]=(0-(2*sinew[127][j]))/sqrdelx;
}
for(i=1;i<128;i++)
{
wnew[i][0]=(0-(2*sinew[i][1]))/sqrdely;
wnew[i][128]=(0-(2*(sinew[i][127]+dely)))/sqrdely;
}
//update vorticity
for(k=0;k<2;k++)
{
for(j=1;j<128;j++)
{
for(i=1;i<128;i++)
{
 wnew[i][j]=(wnew[i+1][j]+wnew[i-1][j]+wnew[i][j+1]+wnew[i][j-1]-(100.0*(wnew[i+1][j]-wnew[i-1][j])*(sinew[i][j+1]-sinew[i][j-1]))+(100.0*(wnew[i][j+1]-wnew[i][j-1])*(sinew[i+1][j]-sinew[i-1][j])))/4.0;
}
}
}
SUM=0;
s=0;
for(j=0;j<129;j++)
{
for(i=0;i<129;i++)
{
SUM=SUM+fabs(wnew[i][j]-wold[i][j]);

s=s+fabs(wnew[i][j]);
}
}
error=SUM/s;
printf("error=%lf\n",error);
for(j=0;j<129;j++)
{
for(i=0;i<129;i++)
{
siold[i][j]=sinew[i][j];
wold[i][j]=wnew[i][j];
}
}
}
while(error>0.000001);
for(j=0;j<129;j++)
{
for(i=0;i<129;i++)
{
fprintf(fp,"%d\t%d\t%lf\n",i+1,j+1,sinew[i][j]);
}
}
 for(j=1;j<128;j++)
{
for(i=1;i<128;i++)
{
u[i][j]=(sinew[i][j+1]-sinew[i][j-1])/(2*dely);
v[i][j]=(0-1)*(sinew[i+1][j]-sinew[i-1][j])/(2*delx);
}
}
for(j=0;j<129;j++)
{
for(i=0;i<129;i++)
{
fprintf(U,"%d\t%d\t%lf\t%lf\n",i+1,j+1,u[i][j],v[i][j]);
 }
 }
 //counting no. of iterations
 printf("\n%d",count);
 fclose(fp);
return 0;
}
