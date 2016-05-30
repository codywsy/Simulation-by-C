#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "main.h"
#include "FiniteFieldBasisGF(4).h"
#include "FiniteFieldBasisCal.h"
#include "BackInter.h"
#include "NewInter.h"

//from Herm(8,4)_GS_v1.0.cpp
extern int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];

void NewInter(int g[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int interpoint[][n-eta])
{
	int i, j, u, v, z, delta_temp, lod_temp;
	int delta[init_polyNum], lod[init_polyNum], J[init_polyNum], lod_min, j_min, flag_lod_min;
	int f[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int poly[init_polyNum][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int g1[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize+1];
	int g2[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];


	//Initialization
	for(j=0; j<init_polyNum; ++j)
		for(u=0; u<New_interpoly_Zsize; ++u)
			for(v=0; v<New_interpoly_Ysize; ++v)
				for(z=0; z<New_interpoly_Xsize; ++z)
					poly[j][u][v][z] = 0;

	for(u=0; u<(lm+1); ++u)
		for(v=0; v<w; ++v)
			poly[w*u+v][u][v][0] = 1;

	//Interpolation
	for(i=0; i<n-eta; ++i)
	{
		//Calculate each poly's leading order
		for(j=0; j<init_polyNum; ++j)
		{
			lod_temp = 0;
			lod[j] = 0;
			for(u=0; u<New_interpoly_Zsize; ++u)
				for(v=0; v<New_interpoly_Ysize; ++v)
					for(z=0; z<New_interpoly_Xsize; ++z)
					{
						if(poly[j][u][v][z] != 0)
						{
							lod_temp = MonoOrderConvert(u,v,z);
							if(lod_temp > lod[j])
								lod[j] = lod_temp;
						}
					}
		}

		//Calculate the hasse derivative mapping of each polynomials
		j_min = -1;
		for (flag_lod_min = 1, j=0; j < init_polyNum; ++j)
		{

			delta[j] = 0;
			J[j] = 0;
			for(u=0; u<New_interpoly_Zsize; ++u)
			{
				delta_temp = 0;

				for(v=0; v<New_interpoly_Ysize; ++v)
					for(z=0; z<New_interpoly_Xsize; ++z)
						if(poly[j][u][v][z] != 0)
						{
							int temp = mul( power(interpoint[0][i],z), power(interpoint[1][i],v) );
							temp = mul( temp, poly[j][u][v][z] );
							//Hasse derivative mapping
							delta_temp = add( delta_temp, temp );
						}
				
				if(u==0)
					delta[j] = delta_temp;
				else if(u>0)
				{
					delta_temp = mul( delta_temp, power(interpoint[2][i],u) );
					delta[j] = add( delta[j], delta_temp );
				}

			}

			if(delta[j] != 0 && flag_lod_min == 1)
			{
				flag_lod_min = 0;
				lod_min = lod[j];
				j_min = j;
			}
			
			if (delta[j] != 0)
				J[j] = 1;
		}	

		//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
		for(j=0; j<init_polyNum; ++j)
			if( lod[j]<lod_min && J[j]!=0)
			{
				lod_min = lod[j];
				j_min = j;
			}

		//Polynomial modification
		if(j_min != -1)
		{
			//f = g[j_min]
			for(u=0; u<New_interpoly_Zsize; ++u)
				for(v=0; v<New_interpoly_Ysize; ++v)
					for(z=0; z<New_interpoly_Xsize; ++z)
						f[u][v][z] = poly[j_min][u][v][z];

			//Modify nonzero polynomials
			for(j=0; j<init_polyNum; ++j)
			{
				int temp1, temp2;

				if (J[j] != 0)
				{
					if (j != j_min)
					{
						//delta*poly_k+delta_k*f
						for (u = 0; u < New_interpoly_Zsize; ++u)
							for (v = 0; v < New_interpoly_Ysize; ++v)
								for (z = 0; z < New_interpoly_Xsize; ++z)
								{
									if (poly[j][u][v][z] != 0)
										temp1 = mul(delta[j_min], poly[j][u][v][z]);
									else
										temp1 = 0;

									if (f[u][v][z] != 0)
										temp2 = mul(delta[j], f[u][v][z]);
									else
										temp2 = 0;

									if (temp1 != 0 || temp2 != 0)
										poly[j][u][v][z] = add(temp1, temp2);
									else
										poly[j][u][v][z] = 0;
								}
					}
					else if (j == j_min)
					{
						for (u = 0; u < New_interpoly_Zsize; ++u)
							for (v = 0; v < New_interpoly_Ysize; ++v)
							{
								for (z = 0; z < (New_interpoly_Xsize + 1); ++z)
									g1[u][v][z] = 0;
								for (z = 0; z < New_interpoly_Xsize; ++z)
									g2[u][v][z] = 0;
							}


						//g1 = x * f
						for (u = 0; u < New_interpoly_Zsize; ++u)
							for (v = 0; v < New_interpoly_Ysize; ++v)
								for (z = 0; z < New_interpoly_Xsize; ++z)
									if (poly[j][u][v][z] != 0)
										g1[u][v][z + 1] = poly[j][u][v][z];

						//g2 = xi * f
						for (u = 0; u < New_interpoly_Zsize; ++u)
							for (v = 0; v < New_interpoly_Ysize; ++v)
								for (z = 0; z < New_interpoly_Xsize; ++z)
									if (poly[j][u][v][z] != 0)
										g2[u][v][z] = mul(interpoint[0][i], poly[j][u][v][z]);

						//poly = g1 + g2
						for (u = 0; u < New_interpoly_Zsize; ++u)
							for (v = 0; v < New_interpoly_Ysize; ++v)
								for (z = 0; z < New_interpoly_Xsize; ++z)
									if (g1[u][v][z] != 0 || g2[u][v][z] != 0)
										poly[j][u][v][z] = add(g1[u][v][z], g2[u][v][z]);
									else
										poly[j][u][v][z] = 0;
					}
				}
			}
		}
	}

	//convert x^w+1=y^w+y
	for(j=0; j<init_polyNum; ++j)
	{
		//Initialization
		int poly_temp[New_interpoly_Zsize][interpoly_Ysize][New_interpoly_Xsize];
		for(u=0; u<New_interpoly_Zsize; ++u)
			for(v=0; v<interpoly_Ysize; ++v)
				if(v<New_interpoly_Ysize)
				{	
					for(z=0; z<New_interpoly_Xsize; ++z)
						poly_temp[u][v][z] = poly[j][u][v][z];
				}
				else
				{
					for(z=0; z<New_interpoly_Xsize; ++z)
						poly_temp[u][v][z] = 0;
				}

		//x^(w+1)=y^w+y
		for(u=0; u<New_interpoly_Zsize; u++)
		{
			int index_flag = New_interpoly_Xsize-1;	//the max deg_x of poly_temp
			while(index_flag > w)
			{
				int temp = index_flag - (w+1);	//deg_x - (w+1)
				for(v=0; v<New_interpoly_Ysize; ++v)
					if(poly_temp[u][v][index_flag] != 0)
					{
						poly_temp[u][v+1][temp] = add( poly_temp[u][v+1][temp], poly_temp[u][v][index_flag] );
						poly_temp[u][v+w][temp] = add( poly_temp[u][v+w][temp], poly_temp[u][v][index_flag] );
						poly_temp[u][v][index_flag] = 0;
					}
				index_flag = index_flag - 1;
			}
		}

		//return the polynomial group
		for (u = 0; u < interpoly_Zsize; ++u)
			for (v = 0; v < interpoly_Ysize; ++v)
				for (z = 0; z < interpoly_Xsize; ++z)
					g[j][u][v][z] = poly_temp[u][v][z];
	}




}


int MonoOrderConvert(int degz, int degy, int degx)
{
	int New_degy = degy;
	int New_degx = degx;

	if (degx > w)
	{
		New_degy = (degx / (w + 1)) * w;
		New_degx = degx % (w + 1);
	}

	return mono_order[degz][New_degy][New_degx];
}