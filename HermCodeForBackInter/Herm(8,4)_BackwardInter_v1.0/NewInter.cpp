#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "main.h"
#include "FiniteFieldBasisGF(4).h"
#include "FiniteFieldBasisCal.h"
#include "BackInter.h"
#include "NewInter.h"
#include "test_LM_validity.h"

//from Herm(8,4)_GS_v1.0.cpp
extern int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];
extern int x_ordered[2][n];
extern int large_vec[choose_num][n];
extern float test_set_ordered[2][n];
extern int Q_interpoly_BF[test_vec_num][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];	//chosen interpolated polynomial [z-deg+1][max(deg_y)+1][w+1]
extern int degree_test[test_vec_num];		//store the degree of the choosen polynomial
extern unsigned long int seq_num_Now;


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

	//debug:
	//Calculate each poly's leading order
	for (j = 0; j<init_polyNum; ++j)
	{
		lod_temp = 0;
		lod[j] = 0;
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
				{
					if (poly[j][u][v][z] != 0)
					{
						lod_temp = MonoOrderConvert(u, v, z);
						if (lod_temp > lod[j])
							lod[j] = lod_temp;
					}
				}
	}

	////debug: BF for Hermitian
	//int Bpoint[3] = { 3, 3, 2 };
	//int Fpoint[3] = { 3, 3, 3 };
	//BackInterpolation(poly, Bpoint);
	//InterOnce(poly, Fpoint);

	//Bpoint[0] = 0;	Bpoint[1] = 1;	Bpoint[2] = 0;
	//Fpoint[0] = 0;	Fpoint[1] = 1;	Fpoint[2] = 1;
	//BackInterpolation(poly, Bpoint);
	//InterOnce(poly, Fpoint);

	//Bpoint[0] = 3;	Bpoint[1] = 3;	Bpoint[2] = 3;
	//Fpoint[0] = 3;	Fpoint[1] = 3;	Fpoint[2] = 2;
	//BackInterpolation(poly, Bpoint);
	//InterOnce(poly, Fpoint);


	//convert x^w+1=y^w+y
	ConvertX2Y(poly, g, init_polyNum, w);




}


#ifndef _GS_Normal_

/************************************
Function:
	General Backward Interpolation with MonoOrderConvert()

Problem:
   Have something wrong with LM(Q) judgement
************************************/

void NewInter_1()
{
	int i, j, u, v, z;
	int delta_temp, delta[init_polyNum], lod_temp, lod_min, flag_lod_min, lod[init_polyNum], j_min, J[init_polyNum];
	int com_elem_interpoint[n][3];
	int poly[init_polyNum][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int f[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int g1[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize + 1];
	int g2[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int g[init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	int GrayOrder[test_vec_num];
	int uncom_elem_interpoint[eta][4];

	//Initialize
	for (i = 0; i<test_vec_num; i++)
		for (u = 0; u<interpoly_Zsize; u++)	//rs
			for (v = 0; v<interpoly_Ysize; v++)	//max(deg_y)+1
				for (z = 0; z<interpoly_Xsize; z++)	//w+1
					Q_interpoly_BF[i][u][v][z] = 0;

	//Interpolation for first path

	//set common element interpoint (xi,ri)
	for (i = 0; i<n; ++i)
	{
		com_elem_interpoint[i][0] = x_ordered[0][i];
		com_elem_interpoint[i][1] = x_ordered[1][i];
		com_elem_interpoint[i][2] = large_vec[0][(int)test_set_ordered[1][i]];;
	}

	//Initialization
	for (j = 0; j<init_polyNum; ++j)
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
					poly[j][u][v][z] = 0;

	for (u = 0; u<(lm + 1); ++u)
		for (v = 0; v<w; ++v)
			poly[w*u + v][u][v][0] = 1;

	//Interpolation
	for (i = 0; i<n; ++i)
	{
		//Calculate each poly's leading order
		for (j = 0; j<init_polyNum; ++j)
		{
			lod_temp = 0;
			lod[j] = 0;
			for (u = 0; u<New_interpoly_Zsize; ++u)
				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
					{
						if (poly[j][u][v][z] != 0)
						{
							lod_temp = MonoOrderConvert(u, v, z);
							if (lod_temp > lod[j])
								lod[j] = lod_temp;
						}
					}
		}

		//Calculate the hasse derivative mapping of each polynomials
		j_min = -1;
		for (flag_lod_min = 1, j = 0; j < init_polyNum; ++j)
		{

			delta[j] = 0;
			J[j] = 0;
			for (u = 0; u<New_interpoly_Zsize; ++u)
			{
				delta_temp = 0;

				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
						if (poly[j][u][v][z] != 0)
						{
							int temp = mul(power(com_elem_interpoint[i][0], z), power(com_elem_interpoint[i][1], v));
							temp = mul(temp, poly[j][u][v][z]);
							//Hasse derivative mapping
							delta_temp = add(delta_temp, temp);
						}

				if (u == 0)
					delta[j] = delta_temp;
				else if (u>0)
				{
					delta_temp = mul(delta_temp, power(com_elem_interpoint[i][2], u));
					delta[j] = add(delta[j], delta_temp);
				}

			}

			if (delta[j] != 0 && flag_lod_min == 1)
			{
				flag_lod_min = 0;
				lod_min = lod[j];
				j_min = j;
			}

			if (delta[j] != 0)
				J[j] = 1;
		}

		//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
		for (j = 0; j<init_polyNum; ++j)
			if (lod[j]<lod_min && J[j] != 0)
			{
				lod_min = lod[j];
				j_min = j;
			}

		//Polynomial modification
		if (j_min != -1)
		{
			//f = g[j_min]
			for (u = 0; u<New_interpoly_Zsize; ++u)
				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
						f[u][v][z] = poly[j_min][u][v][z];

			//Modify nonzero polynomials
			for (j = 0; j<init_polyNum; ++j)
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
										g2[u][v][z] = mul(com_elem_interpoint[i][0], poly[j][u][v][z]);

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

	for (j = 0; j<init_polyNum; ++j)
		for (u = 0; u<interpoly_Zsize; ++u)
			for (v = 0; v<interpoly_Ysize; ++v)
				for (z = 0; z<interpoly_Xsize; ++z)
					g[j][u][v][z] = 0;

	ConvertX2Y(poly, g, init_polyNum, w);


	j_min = ChooseTheMinPoly(g, init_polyNum);
	for (u = 0; u < interpoly_Zsize; ++u)
		for (v = 0; v < interpoly_Ysize; ++v)
			for (z = 0; z < interpoly_Xsize; ++z)
				Q_interpoly_BF[0][u][v][z] = g[j_min][u][v][z];
	//Interpolation for first path finish


	//Produce GrayCode used by the ordering of BFinterpolation
	for (i = 0; i < test_vec_num; ++i)
	{
		GrayOrder[i] = i ^ (i >> 1);
	}

	//produce interpoint group for Backward interpolation
	for (i = 0; i<eta; ++i)
	{
		uncom_elem_interpoint[i][0] = x_ordered[0][n - eta + i];
		uncom_elem_interpoint[i][1] = x_ordered[1][n - eta + i];
		uncom_elem_interpoint[i][2] = large_vec[0][(int)test_set_ordered[1][n - eta + i]];
		uncom_elem_interpoint[i][3] = large_vec[1][(int)test_set_ordered[1][n - eta + i]];
	}

	//Binary tree in Backward interpolation
	for (i = 1; i<test_vec_num; ++i)
	{
		//different location and value between GrayOrder[i-1] and GrayOrder[i]
		int result = GrayOrder[i - 1] ^ GrayOrder[i];
		int location = -1;
		for (int temp = result; temp > 0;)
		{
			temp = temp >> 1;
			++location;
		}
		location = eta - 1 - location;

		int Blocation = result & GrayOrder[i - 1];
		int Bpoint[3] = { uncom_elem_interpoint[location][0], uncom_elem_interpoint[location][1], uncom_elem_interpoint[location][2 + Blocation] };
		int Flocation = (result & GrayOrder[i]) ? 1 : 0;
		int Fpoint[3] = { uncom_elem_interpoint[location][0], uncom_elem_interpoint[location][1], uncom_elem_interpoint[location][2 + Flocation] };

		//Backward Interpolation
		BackInterpolation(poly, Bpoint);

		//Forward Interpolation
		InterOnce(poly, Fpoint);

		//find the min polynomial
		for (j = 0; j<init_polyNum; ++j)
			for (u = 0; u<interpoly_Zsize; ++u)
				for (v = 0; v<interpoly_Ysize; ++v)
					for (z = 0; z<interpoly_Xsize; ++z)
						g[j][u][v][z] = 0;

		ConvertX2Y(poly, g, init_polyNum, w);

		j_min = -1;
		j_min = ChooseTheMinPoly(g, init_polyNum);
		for (u = 0; u<interpoly_Zsize; ++u)
			for (v = 0; v<interpoly_Ysize; ++v)
				for (z = 0; z<interpoly_Xsize; ++z)
					Q_interpoly_BF[i][u][v][z] = g[j_min][u][v][z];
	}

}


/************************************
Function:
	Backward Interpolation(2.0) with New LM(Q) judgement type 
************************************/
void NewInter_2()
{
	int i, j, u, v, z;
	int delta_temp, delta[init_polyNum], lod_temp, lod_min, flag_lod_min, lod[init_polyNum], j_min, J[init_polyNum];
	int com_elem_interpoint[n][3];
	int poly_min_order[eta];
	int poly_min_order_all[n], location_all;
	int poly[init_polyNum][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int f[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int g1[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize + 1];
	int g2[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int g[init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	int GrayOrder[test_vec_num];
	int uncom_elem_interpoint[eta][4];
	
	//****debug***
	int LeadingMono[init_polyNum][3];
	//*********

	//Initialize
	for (i = 0; i<test_vec_num; i++)
		for (u = 0; u<interpoly_Zsize; u++)	//rs
			for (v = 0; v<interpoly_Ysize; v++)	//max(deg_y)+1
				for (z = 0; z<interpoly_Xsize; z++)	//w+1
					Q_interpoly_BF[i][u][v][z] = 0;

	//Interpolation for first path

	//set common element interpoint (xi,ri)
	for (i = 0; i<n; ++i)
	{
		com_elem_interpoint[i][0] = x_ordered[0][i];
		com_elem_interpoint[i][1] = x_ordered[1][i];
		com_elem_interpoint[i][2] = large_vec[0][(int)test_set_ordered[1][i]];
	}

	//debug - BF test
	//if (seq_num_Now == 13)
	//{
	//	com_elem_interpoint[n - 1][2] = large_vec[1][(int)test_set_ordered[1][i - 1]];
	//	//com_elem_interpoint[n - 2][2] = large_vec[1][(int)test_set_ordered[1][i - 2]];

	//	//int temp;
	//	//temp = com_elem_interpoint[n - 2][0];
	//	//com_elem_interpoint[n - 2][0] = com_elem_interpoint[n - 1][0];
	//	//com_elem_interpoint[n - 1][0] = temp;

	//	//temp = com_elem_interpoint[n - 2][1];
	//	//com_elem_interpoint[n - 2][1] = com_elem_interpoint[n - 1][1];
	//	//com_elem_interpoint[n - 1][1] = temp;

	//	//temp = com_elem_interpoint[n - 2][2];
	//	//com_elem_interpoint[n - 2][2] = com_elem_interpoint[n - 1][2];
	//	//com_elem_interpoint[n - 1][2] = temp;
	//	
	//}

	////debug - random order
	//com_elem_interpoint[0][0] = 0;	com_elem_interpoint[0][1] = 0;	com_elem_interpoint[0][2] = 1;
	//com_elem_interpoint[1][0] = 2;	com_elem_interpoint[1][1] = 2;	com_elem_interpoint[1][2] = 3;
	//com_elem_interpoint[2][0] = 1;	com_elem_interpoint[2][1] = 3;  com_elem_interpoint[2][2] = 3;
	//com_elem_interpoint[3][0] = 2;	com_elem_interpoint[3][1] = 3;	com_elem_interpoint[3][2] = 1;
	//com_elem_interpoint[4][0] = 3;	com_elem_interpoint[4][1] = 3;	com_elem_interpoint[4][2] = 3;
	//com_elem_interpoint[5][0] = 3;	com_elem_interpoint[5][1] = 2;	com_elem_interpoint[5][2] = 1;
	//com_elem_interpoint[6][0] = 1;	com_elem_interpoint[6][1] = 2;	com_elem_interpoint[6][2] = 0;
	//com_elem_interpoint[7][0] = 0;	com_elem_interpoint[7][1] = 1;	com_elem_interpoint[7][2] = 2;
	//
	////random num without repeat
	//int unorder_set[n];
	//for (j = 0; j < n; j++)
	//	unorder_set[j] = 0;

	
	//Initialization
	for (j = 0; j<init_polyNum; ++j)
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
					poly[j][u][v][z] = 0;

	for (u = 0; u<(lm + 1); ++u)
		for (v = 0; v<w; ++v)
			poly[w*u + v][u][v][0] = 1;

	//Interpolation
	for (i = 0; i<n; ++i)
	{
		//Calculate each poly's leading order
		//1st, poly converts to g
		for (j = 0; j<init_polyNum; ++j)
			for (u = 0; u<interpoly_Zsize; ++u)
				for (v = 0; v<interpoly_Ysize; ++v)
					for (z = 0; z<interpoly_Xsize; ++z)
						g[j][u][v][z] = 0;

		ConvertX2Y(poly, g, init_polyNum, w);

		//2nd, use g to calculate the lod[]
		for (j = 0; j<init_polyNum; ++j)
		{
			lod_temp = 0;
			lod[j] = 0;
			for (u = 0; u<New_interpoly_Zsize; ++u)
				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
					{
						if (g[j][u][v][z] != 0)
						{
							lod_temp = mono_order[u][v][z];
							if (lod_temp > lod[j])
								lod[j] = lod_temp;
						}
					}
		}

		//Calculate the hasse derivative mapping of each polynomials
		j_min = -1;
		for (flag_lod_min = 1, j = 0; j < init_polyNum; ++j)
		{

			delta[j] = 0;
			J[j] = 0;
			for (u = 0; u<New_interpoly_Zsize; ++u)
			{
				delta_temp = 0;

				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
						if (poly[j][u][v][z] != 0)
						{
							int temp = mul(power(com_elem_interpoint[i][0], z), power(com_elem_interpoint[i][1], v));
							temp = mul(temp, poly[j][u][v][z]);
							//Hasse derivative mapping
							delta_temp = add(delta_temp, temp);
						}

				if (u == 0)
					delta[j] = delta_temp;
				else if (u>0)
				{
					delta_temp = mul(delta_temp, power(com_elem_interpoint[i][2], u));
					delta[j] = add(delta[j], delta_temp);
				}

			}

			if (delta[j] != 0 && flag_lod_min == 1)
			{
				flag_lod_min = 0;
				lod_min = lod[j];
				j_min = j;
			}

			if (delta[j] != 0)
				J[j] = 1;
		}

		//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
		for (j = 0; j<init_polyNum; ++j)
			if (lod[j]<lod_min && J[j] != 0)
			{
				lod_min = lod[j];
				j_min = j;
			}

		//record the j_min seq for every poly
		if (i - n + eta >= 0)
			poly_min_order[i - n + eta] = j_min;
		
		poly_min_order_all[i] = j_min;


		//Polynomial modification
		if (j_min != -1)
		{
			//f = g[j_min]
			for (u = 0; u<New_interpoly_Zsize; ++u)
				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
						f[u][v][z] = poly[j_min][u][v][z];

			//Modify nonzero polynomials
			for (j = 0; j<init_polyNum; ++j)
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
										g2[u][v][z] = mul(com_elem_interpoint[i][0], poly[j][u][v][z]);

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

	for (j = 0; j<init_polyNum; ++j)
		for (u = 0; u<interpoly_Zsize; ++u)
			for (v = 0; v<interpoly_Ysize; ++v)
				for (z = 0; z<interpoly_Xsize; ++z)
					g[j][u][v][z] = 0;

	ConvertX2Y(poly, g, init_polyNum, w);

	j_min = ChooseTheMinPoly(g, init_polyNum);
	for (u = 0; u < interpoly_Zsize; ++u)
		for (v = 0; v < interpoly_Ysize; ++v)
			for (z = 0; z < interpoly_Xsize; ++z)
				Q_interpoly_BF[0][u][v][z] = g[j_min][u][v][z];
	//Interpolation for first path finish


	////debug
	//if (seq_num_Now == 13)
	//{
	//	int poly_min_order_temp[eta] = { 0, 2 };
	//	int location = 0;
	//	int point_BF[3] = { 1, 2, 0 };
	//	int point_FI[3] = { 1, 2, 1 };
	//	test_LM_validity(poly, point_BF, point_FI, poly_min_order_temp, location);
	//	
	//}



	//Produce GrayCode used by the ordering of BFinterpolation
	for (i = 0; i < test_vec_num; ++i)
	{
		GrayOrder[i] = i ^ (i >> 1);
	}

	//produce interpoint group for Backward interpolation
	for (i = 0; i<eta; ++i)
	{
		uncom_elem_interpoint[i][0] = x_ordered[0][n - eta + i];
		uncom_elem_interpoint[i][1] = x_ordered[1][n - eta + i];
		uncom_elem_interpoint[i][2] = large_vec[0][(int)test_set_ordered[1][n - eta + i]];
		uncom_elem_interpoint[i][3] = large_vec[1][(int)test_set_ordered[1][n - eta + i]];
	}

	//Binary tree in Backward interpolation
	for (i = 1; i<test_vec_num; ++i)
	{
		//different location and value between GrayOrder[i-1] and GrayOrder[i]
		int result = GrayOrder[i - 1] ^ GrayOrder[i];
		int location = -1;
		for (int temp = result; temp > 0;)
		{
			temp = temp >> 1;
			++location;
		}
		location = eta - 1 - location;

		int Blocation = result & GrayOrder[i - 1];
		int Bpoint[3] = { uncom_elem_interpoint[location][0], uncom_elem_interpoint[location][1], uncom_elem_interpoint[location][2 + Blocation] };
		int Flocation = (result & GrayOrder[i]) ? 1 : 0;
		int Fpoint[3] = { uncom_elem_interpoint[location][0], uncom_elem_interpoint[location][1], uncom_elem_interpoint[location][2 + Flocation] };

		//Calculate each poly's leading order
		for (j = 0; j<init_polyNum; ++j)
		{
			lod_temp = 0;
			lod[j] = 0;
			for (u = 0; u<New_interpoly_Zsize; ++u)
				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
					{
						if (poly[j][u][v][z] != 0)
						{
							lod_temp = MonoOrderConvert(u, v, z);
							if (lod_temp > lod[j])
								lod[j] = lod_temp;
						}
					}
		}


		//Backward Interpolation
		BackInterpolation_2(poly, Bpoint, poly_min_order, location, lod, poly_min_order_all, &location_all);

		//Calculate each poly's leading order
		for (j = 0; j<init_polyNum; ++j)
		{
			lod_temp = 0;
			lod[j] = 0;
			for (u = 0; u<New_interpoly_Zsize; ++u)
				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
					{
						if (poly[j][u][v][z] != 0)
						{
							lod_temp = MonoOrderConvert(u, v, z);
							if (lod_temp > lod[j])
								lod[j] = lod_temp;
						}
					}
		}


		//debug: print the lod after BF
#ifdef _PrintGroupLod_
		printf("\nAfter %dst BI, the lod =\t{%d,\t%d,\t%d,\t%d}-->", i, lod[0], lod[1], lod[2], lod[3]);
#endif

		//Forward Interpolation
		InterOnce_2(poly, Fpoint, poly_min_order, location, poly_min_order_all, location_all);

		//Calculate each poly's leading order
		for (j = 0; j<init_polyNum; ++j)
		{
			lod_temp = 0;
			lod[j] = 0;
			for (u = 0; u<New_interpoly_Zsize; ++u)
				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
					{
						if (poly[j][u][v][z] != 0)
						{
							lod_temp = MonoOrderConvert(u, v, z);
							if (lod_temp > lod[j])
								lod[j] = lod_temp;
						}
					}
		}

		//debug: print the lod after BF
#ifdef _PrintGroupLod_
		printf("After %dst FI, the lod =\t{%d,\t%d,\t%d,\t%d}", i, lod[0], lod[1], lod[2], lod[3]);
#endif

		//find the min polynomial
		for (j = 0; j<init_polyNum; ++j)
			for (u = 0; u<interpoly_Zsize; ++u)
				for (v = 0; v<interpoly_Ysize; ++v)
					for (z = 0; z<interpoly_Xsize; ++z)
						g[j][u][v][z] = 0;

		ConvertX2Y(poly, g, init_polyNum, w);

		j_min = -1;
		j_min = ChooseTheMinPoly(g, init_polyNum);
		for (u = 0; u<interpoly_Zsize; ++u)
			for (v = 0; v<interpoly_Ysize; ++v)
				for (z = 0; z<interpoly_Xsize; ++z)
					Q_interpoly_BF[GrayOrder[i]][u][v][z] = g[j_min][u][v][z];
	}
#ifdef _PrintGroupLod_
	printf("\n\n");
#endif

}

/************************************
Function:
	Backward Interpolation(3.0) with New delta[] judgement
************************************/
void NewInter_3()
{
	int i, j, u, v, z;
	int delta_temp, delta[init_polyNum], lod_temp, lod_min, flag_lod_min, lod[init_polyNum], j_min, J[init_polyNum];
	int com_elem_interpoint[n][3];
	int poly[init_polyNum][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int f[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int g1[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize + 1];
	int g2[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int g[init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	int GrayOrder[test_vec_num];
	int uncom_elem_interpoint[eta][4];

	//Initialize
	for (i = 0; i<test_vec_num; i++)
		for (u = 0; u<interpoly_Zsize; u++)	//rs
			for (v = 0; v<interpoly_Ysize; v++)	//max(deg_y)+1
				for (z = 0; z<interpoly_Xsize; z++)	//w+1
					Q_interpoly_BF[i][u][v][z] = 0;

	//Interpolation for first path

	//set common element interpoint (xi,ri)
	for (i = 0; i<n; ++i)
	{
		com_elem_interpoint[i][0] = x_ordered[0][i];
		com_elem_interpoint[i][1] = x_ordered[1][i];
		com_elem_interpoint[i][2] = large_vec[0][(int)test_set_ordered[1][i]];;
	}

	//Initialization
	for (j = 0; j<init_polyNum; ++j)
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
					poly[j][u][v][z] = 0;

	for (u = 0; u<(lm + 1); ++u)
		for (v = 0; v<w; ++v)
			poly[w*u + v][u][v][0] = 1;

	//Interpolation
	for (i = 0; i<n; ++i)
	{
		//Calculate each poly's leading order
		for (j = 0; j<init_polyNum; ++j)
		{
			lod_temp = 0;
			lod[j] = 0;
			for (u = 0; u<New_interpoly_Zsize; ++u)
				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
					{
						if (poly[j][u][v][z] != 0)
						{
							lod_temp = MonoOrderConvert(u, v, z);
							if (lod_temp > lod[j])
								lod[j] = lod_temp;
						}
					}
		}

		//Calculate the hasse derivative mapping of each polynomials
		//1st, poly converts to g
		for (j = 0; j<init_polyNum; ++j)
			for (u = 0; u<interpoly_Zsize; ++u)
				for (v = 0; v<interpoly_Ysize; ++v)
					for (z = 0; z<interpoly_Xsize; ++z)
						g[j][u][v][z] = 0;

		ConvertX2Y(poly, g, init_polyNum, w);

		//2nd, use g calculate delta[]
		j_min = -1;
		for (flag_lod_min = 1, j = 0; j < init_polyNum; ++j)
		{

			delta[j] = 0;
			J[j] = 0;
			for (u = 0; u<interpoly_Zsize; ++u)
			{
				delta_temp = 0;

				for (v = 0; v<interpoly_Ysize; ++v)
					for (z = 0; z<interpoly_Xsize; ++z)
						if (g[j][u][v][z] != 0)
						{
							int temp = mul(power(com_elem_interpoint[i][0], z), power(com_elem_interpoint[i][1], v));
							temp = mul(temp, g[j][u][v][z]);
							//Hasse derivative mapping
							delta_temp = add(delta_temp, temp);
						}

				if (u == 0)
					delta[j] = delta_temp;
				else if (u>0)
				{
					delta_temp = mul(delta_temp, power(com_elem_interpoint[i][2], u));
					delta[j] = add(delta[j], delta_temp);
				}

			}

			if (delta[j] != 0 && flag_lod_min == 1)
			{
				flag_lod_min = 0;
				lod_min = lod[j];
				j_min = j;
			}

			if (delta[j] != 0)
				J[j] = 1;
		}

		//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
		for (j = 0; j<init_polyNum; ++j)
			if (lod[j]<lod_min && J[j] != 0)
			{
				lod_min = lod[j];
				j_min = j;
			}

		//Polynomial modification
		if (j_min != -1)
		{
			//f = g[j_min]
			for (u = 0; u<New_interpoly_Zsize; ++u)
				for (v = 0; v<New_interpoly_Ysize; ++v)
					for (z = 0; z<New_interpoly_Xsize; ++z)
						f[u][v][z] = poly[j_min][u][v][z];

			//Modify nonzero polynomials
			for (j = 0; j<init_polyNum; ++j)
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
										g2[u][v][z] = mul(com_elem_interpoint[i][0], poly[j][u][v][z]);

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

	for (j = 0; j<init_polyNum; ++j)
		for (u = 0; u<interpoly_Zsize; ++u)
			for (v = 0; v<interpoly_Ysize; ++v)
				for (z = 0; z<interpoly_Xsize; ++z)
					g[j][u][v][z] = 0;

	ConvertX2Y(poly, g, init_polyNum, w);


	j_min = ChooseTheMinPoly(g, init_polyNum);
	for (u = 0; u < interpoly_Zsize; ++u)
		for (v = 0; v < interpoly_Ysize; ++v)
			for (z = 0; z < interpoly_Xsize; ++z)
				Q_interpoly_BF[0][u][v][z] = g[j_min][u][v][z];
	//Interpolation for first path finish


	//Produce GrayCode used by the ordering of BFinterpolation
	for (i = 0; i < test_vec_num; ++i)
	{
		GrayOrder[i] = i ^ (i >> 1);
	}

	//produce interpoint group for Backward interpolation
	for (i = 0; i<eta; ++i)
	{
		uncom_elem_interpoint[i][0] = x_ordered[0][n - eta + i];
		uncom_elem_interpoint[i][1] = x_ordered[1][n - eta + i];
		uncom_elem_interpoint[i][2] = large_vec[0][(int)test_set_ordered[1][n - eta + i]];
		uncom_elem_interpoint[i][3] = large_vec[1][(int)test_set_ordered[1][n - eta + i]];
	}

	//Binary tree in Backward interpolation
	for (i = 1; i<test_vec_num; ++i)
	{
		//different location and value between GrayOrder[i-1] and GrayOrder[i]
		int result = GrayOrder[i - 1] ^ GrayOrder[i];
		int location = -1;
		for (int temp = result; temp > 0;)
		{
			temp = temp >> 1;
			++location;
		}
		location = eta - 1 - location;

		int Blocation = result & GrayOrder[i - 1];
		int Bpoint[3] = { uncom_elem_interpoint[location][0], uncom_elem_interpoint[location][1], uncom_elem_interpoint[location][2 + Blocation] };
		int Flocation = (result & GrayOrder[i]) ? 1 : 0;
		int Fpoint[3] = { uncom_elem_interpoint[location][0], uncom_elem_interpoint[location][1], uncom_elem_interpoint[location][2 + Flocation] };

		//Backward Interpolation
		BackInterpolation(poly, Bpoint);

		//Forward Interpolation
		InterOnce(poly, Fpoint);

		//find the min polynomial
		for (j = 0; j<init_polyNum; ++j)
			for (u = 0; u<interpoly_Zsize; ++u)
				for (v = 0; v<interpoly_Ysize; ++v)
					for (z = 0; z<interpoly_Xsize; ++z)
						g[j][u][v][z] = 0;

		ConvertX2Y(poly, g, init_polyNum, w);

		j_min = -1;
		j_min = ChooseTheMinPoly(g, init_polyNum);
		for (u = 0; u<interpoly_Zsize; ++u)
			for (v = 0; v<interpoly_Ysize; ++v)
				for (z = 0; z<interpoly_Xsize; ++z)
					Q_interpoly_BF[i][u][v][z] = g[j_min][u][v][z];
	}

}

#endif

int MonoOrderConvert(int degz, int degy, int degx)
{
	int New_degy = degy;
	int New_degx = degx;

	if (degx > w)
	{
		New_degy = (degx / (w + 1)) * w + degy;
		New_degx = degx % (w + 1);
	}

	return mono_order[degz][New_degy][New_degx];
}

/************************************
Function:
	poly[] interpolate the interpoint[] once

************************************/
void InterOnce(int poly[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int interpoint[3])
{
	int j, j_min, u, v, z;
	int lod[init_polyNum], lod_temp, lod_min, delta[init_polyNum], delta_temp, flag_lod_min;
	int J[init_polyNum];
	int f[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int g1[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize + 1];
	int g2[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];

	//Calculate each poly's leading order
	for (j = 0; j<init_polyNum; ++j)
	{
		lod_temp = 0;
		lod[j] = 0;
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
				{
					if (poly[j][u][v][z] != 0)
					{
						lod_temp = MonoOrderConvert(u, v, z);
						if (lod_temp > lod[j])
							lod[j] = lod_temp;
					}
				}
	}

	//Calculate the hasse derivative mapping of each polynomials
	j_min = -1;
	for (flag_lod_min = 1, j = 0; j < init_polyNum; ++j)
	{

		delta[j] = 0;
		J[j] = 0;
		for (u = 0; u<New_interpoly_Zsize; ++u)
		{
			delta_temp = 0;

			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
					if (poly[j][u][v][z] != 0)
					{
						int temp = mul(power(interpoint[0], z), power(interpoint[1], v));
						temp = mul(temp, poly[j][u][v][z]);
						//Hasse derivative mapping
						delta_temp = add(delta_temp, temp);
					}

			if (u == 0)
				delta[j] = delta_temp;
			else if (u>0)
			{
				delta_temp = mul(delta_temp, power(interpoint[2], u));
				delta[j] = add(delta[j], delta_temp);
			}

		}

		if (delta[j] != 0 && flag_lod_min == 1)
		{
			flag_lod_min = 0;
			lod_min = lod[j];
			j_min = j;
		}

		if (delta[j] != 0)
			J[j] = 1;
	}

	//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
	for (j = 0; j<init_polyNum; ++j)
		if (lod[j]<lod_min && J[j] != 0)
		{
			lod_min = lod[j];
			j_min = j;
		}

	//Polynomial modification
	if (j_min != -1)
	{
		//f = g[j_min]
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
					f[u][v][z] = poly[j_min][u][v][z];

		//Modify nonzero polynomials
		for (j = 0; j<init_polyNum; ++j)
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
									g2[u][v][z] = mul(interpoint[0], poly[j][u][v][z]);

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

void InterOnce_2(int poly[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int interpoint[3], int poly_min_order[eta], int location, int poly_min_order_all[n], int location_all)
{
	int j, j_min, u, v, z;
	int lod[init_polyNum], lod_temp, lod_min, delta[init_polyNum], delta_temp, flag_lod_min;
	int J[init_polyNum];
	int f[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];
	int g1[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize + 1];
	int g2[New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];

	//Calculate each poly's leading order
	for (j = 0; j<init_polyNum; ++j)
	{
		lod_temp = 0;
		lod[j] = 0;
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
				{
					if (poly[j][u][v][z] != 0)
					{
						lod_temp = MonoOrderConvert(u, v, z);
						if (lod_temp > lod[j])
							lod[j] = lod_temp;
					}
				}
	}

	//Calculate the hasse derivative mapping of each polynomials
	j_min = -1;
	for (flag_lod_min = 1, j = 0; j < init_polyNum; ++j)
	{

		delta[j] = 0;
		J[j] = 0;
		for (u = 0; u<New_interpoly_Zsize; ++u)
		{
			delta_temp = 0;

			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
					if (poly[j][u][v][z] != 0)
					{
						int temp = mul(power(interpoint[0], z), power(interpoint[1], v));
						temp = mul(temp, poly[j][u][v][z]);
						//Hasse derivative mapping
						delta_temp = add(delta_temp, temp);
					}

			if (u == 0)
				delta[j] = delta_temp;
			else if (u>0)
			{
				delta_temp = mul(delta_temp, power(interpoint[2], u));
				delta[j] = add(delta[j], delta_temp);
			}

		}

		if (delta[j] != 0 && flag_lod_min == 1)
		{
			flag_lod_min = 0;
			lod_min = lod[j];
			j_min = j;
		}

		if (delta[j] != 0)
			J[j] = 1;
	}

	//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
	for (j = 0; j<init_polyNum; ++j)
		if (lod[j]<lod_min && J[j] != 0)
		{
			lod_min = lod[j];
			j_min = j;
		}

	//update poly_min_order
	poly_min_order[location] = j_min;
	poly_min_order_all[location_all] = j_min;

	//Polynomial modification
	if (j_min != -1)
	{
		//f = g[j_min]
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
					f[u][v][z] = poly[j_min][u][v][z];

		//Modify nonzero polynomials
		for (j = 0; j<init_polyNum; ++j)
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
									g2[u][v][z] = mul(interpoint[0], poly[j][u][v][z]);

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



/***********************************
Function:
	According x^(w+1) = y^w + y, polySour is converted to polyDes
***********************************/
void ConvertX2Y(int polySour[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int polyDes[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int PolyNum, int omga)
{
	int j, u, v, z;

	for (j = 0; j<PolyNum; ++j)
	{
		//Initialization
		int poly_temp[New_interpoly_Zsize][interpoly_Ysize][New_interpoly_Xsize];
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<interpoly_Ysize; ++v)
				if (v<New_interpoly_Ysize)
				{
					for (z = 0; z<New_interpoly_Xsize; ++z)
						poly_temp[u][v][z] = polySour[j][u][v][z];
				}
				else
				{
					for (z = 0; z<New_interpoly_Xsize; ++z)
						poly_temp[u][v][z] = 0;
				}

		//x^(w+1)=y^w+y
		for (u = 0; u<New_interpoly_Zsize; u++)
		{
			int index_flag = New_interpoly_Xsize - 1;	//the max deg_x of poly_temp
			while (index_flag > omga)
			{
				int temp = index_flag - (omga + 1);	//deg_x - (w+1)
				for (v = 0; v<New_interpoly_Ysize; ++v)
					if (poly_temp[u][v][index_flag] != 0)
					{
						poly_temp[u][v + 1][temp] = add(poly_temp[u][v + 1][temp], poly_temp[u][v][index_flag]);
						poly_temp[u][v + omga][temp] = add(poly_temp[u][v + omga][temp], poly_temp[u][v][index_flag]);
						poly_temp[u][v][index_flag] = 0;
					}
				index_flag = index_flag - 1;
			}
		}

		//return the polynomial group
		for (u = 0; u < interpoly_Zsize; ++u)
			for (v = 0; v < interpoly_Ysize; ++v)
				for (z = 0; z < interpoly_Xsize; ++z)
					polyDes[j][u][v][z] = poly_temp[u][v][z];
	}
}

int ChooseTheMinPoly(int polySour[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int PolyNum)
{
	int temp, index_temp;
	int *degree_temp = (int *)malloc(sizeof(int)*PolyNum);


	for (int j = 0; j < PolyNum; ++j)
	{
		temp = -1;
		degree_temp[j] = 0;
		for (int u = 0; u < interpoly_Zsize; ++u)
			for (int v = 0; v < interpoly_Ysize; ++v)
				for (int z = 0; z < interpoly_Xsize; ++z)
					if (polySour[j][u][v][z] != 0 && temp < mono_order[u][v][z])
						temp = mono_order[u][v][z];

		degree_temp[j] = temp;
	}

	index_temp = 0;
	temp = degree_temp[0];
	for (int j = 1; j < PolyNum; ++j)
		if (temp>degree_temp[j])
		{
			temp = degree_temp[j];
			index_temp = j;
		}

	free(degree_temp);

	return index_temp;

}