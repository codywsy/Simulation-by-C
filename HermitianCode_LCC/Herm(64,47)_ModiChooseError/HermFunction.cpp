//HermFunction.cpp
#include "FiniteFieldBasisCal.h"
#include "FiniteFieldBasisGF(16).h"




extern int root[];

void findpoints(int *point_temp)
{
	int i, j, u, a1, a2, a3, x, y;

	//Initialisation
	for(i=0;i<2;i++)
		for(j=0;j<n;j++)	//w^3
			point_temp[i][j]=0;

	//find points over x^3-y^2-y=0
	u=0;
	for(i=0;i<q;i++)	//number of elements in GF(q)
	{
		x=root[i];
		for(j=0;j<q;j++)	//number of elements in GF(q)
		{
			y=root[j];
			a1=power(x, w+1);
			a2=power(y, w);
			a3=y;
			if(add(a3, add(a1, a2))==0)
			{
				point[0][u]=x;
				point[1][u]=y;
				u++;
			}
		}
	}

}


void tgorder(int *tg_order_temp)
{
	int i,j,u;

	j=0;
	for(i=0;i<tg_size;i++)
		for(u=i;u>=0;u--)
			if(u<=w)
			{
				tg_order_temp[j][0]=u;	// 0<deg_x<=w
				tg_order_temp[j][1]=i-u;	// 0<deg_y
				j++;
			}

}