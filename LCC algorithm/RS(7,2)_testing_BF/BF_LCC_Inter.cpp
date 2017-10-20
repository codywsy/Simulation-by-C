#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "main.h"
#include "FiniteFieldBasisGF(8).h"
#include "FiniteFieldBasisCal.h"
#include "BackInter.h"
#include "BF_LCC_Inter.h"

extern int mono_order[monoTable_Ysize][monoTable_Xsize];
extern int x_ordered[n];
extern int large_vec[choose_num][n];
extern float test_set_ordered[2][n];
extern int Q_interpoly_BF[test_vec_num][interpoly_Ysize][interpoly_Xsize];	//chosen interpolated polynomial [z-deg+1][max(deg_y)+1][w+1]
extern int degree_test[test_vec_num];		//store the degree of the choosen polynomial
extern unsigned long int seq_num_Now;

void BF_LCC_inter()
{
	int i, j, u, v;
	int j_min, lod_temp, lod[init_polyNum], lod_min, J[init_polyNum], delta[init_polyNum], delta_temp, delta_temp1, flag_lod_min;
	int com_elem_interpoint[n][2];
	int poly_min_order[eta];
	int poly_min_order_all[n], location_all;
	int poly[init_polyNum][interpoly_Ysize][interpoly_Xsize];
	int f[interpoly_Ysize][interpoly_Xsize];
	int g1[interpoly_Ysize][interpoly_Xsize + 1];
	int g2[interpoly_Ysize][interpoly_Xsize];
	int GrayOrder[test_vec_num];
	int uncom_elem_interpoint[eta][3];


	//Initialize
	for (i = 0; i<test_vec_num; i++)
		for (u = 0; u<interpoly_Ysize; u++)	//rs
			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
				Q_interpoly_BF[i][u][v] = 0;

	//Interpolation for first path

	//set common element interpoint (xi,ri)
	for (i = 0; i<n; ++i)
	{
		com_elem_interpoint[i][0] = x_ordered[i];
		com_elem_interpoint[i][1] = large_vec[0][(int)test_set_ordered[1][i]];
	}

	//Initialization
	for (j = 0; j<init_polyNum; ++j)
		for (u = 0; u<interpoly_Ysize; ++u)
			for (v = 0; v<interpoly_Xsize; ++v)
				poly[j][u][v] = 0;

	for (i = 0; i<(lm + 1); ++i)
			poly[i][i][0] = 1;

#ifdef _PrintLCCLod_
	printf("\nBF-LCC result:");
#endif

	for (i = 0; i < n; i++)
	{
		//Calculate each polynomial's leading order
		for (j = 0; j<init_polyNum; j++)	//num_poly
		{
			lod_temp = 0;
			lod[j] = 0;

			for (u = 0; u<interpoly_Ysize; u++)	//max(deg_y)+1
			{
				for (v = 0; v<interpoly_Xsize; v++)	//w+1
				{
					if (poly[j][u][v] != 0)
					{
						lod_temp = mono_order[u][v];
						if (lod_temp>lod[j])
							lod[j] = lod_temp;
					}
				}
			}
		}

		//Calculate the hasse derivative mapping of each polynomials
		j_min = -1;
		for (flag_lod_min=1 ,j = 0; j<init_polyNum; j++)	//wrt each polynomial
		{
			J[j] = 0;
			delta[j] = 0;

			//Hasse derivative mapping
			for (u = 0; u<interpoly_Ysize; u++)	//deg_z
			{
				delta_temp1 = 0;

				for (v = 0; v<interpoly_Xsize; v++)	//deg_y
				{
					if (poly[j][u][v] != 0)	//coputation num concerning coefficients
					{
						delta_temp = power(com_elem_interpoint[i][0], v);
						delta_temp = mul(delta_temp, poly[j][u][v]);
						//Hasse derivative mapping
						delta_temp1 = add(delta_temp1, delta_temp);
					}
				}

				if (u == 0)	//deg_z==0
				{
					delta[j] = delta_temp1;
				}
				else if (u>0)	//deg_z>0
				{
					delta_temp1 = mul(delta_temp1, power(com_elem_interpoint[i][1], u));
					delta[j] = add(delta[j], delta_temp1);
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
			if (J[j] != 0 && lod[j]<lod_min)
			{
				lod_min = lod[j];
				j_min = j;
			}
		//printf("\nj_min=%d\n", j_min);

		//record the j_min seq for every poly
		if (i - n + eta >= 0)
			poly_min_order[i - n + eta] = j_min;

		poly_min_order_all[i] = j_min;


#ifdef _PrintLCCLod_
		printf("\nthe %d iter lod = {%d\t%d\t%d}", i, lod[0], lod[1], lod[2]);
		printf("\nj_min = %d", j_min);
		printf("\nJ =\t");
		for (u = 0; u < init_polyNum; u++)
			if (J[u] != 0)
				printf("%d\t", u);
		printf("\n");
#endif

		//Polynomial modification
		if (j_min != -1)
		{
			//f=g[j_min]
			for (u = 0; u < interpoly_Ysize; u++)	//rs
				for (v = 0; v < interpoly_Xsize; v++)	//max(deg_y)+1
					f[u][v] = poly[j_min][u][v];

			//Modify nonzero polynomials
			for (j = 0; j < init_polyNum; j++)	//num of polys
			{
				int temp1, temp2;

				if (J[j] == 1)
				{
					if (j != j_min)
					{
						//delta*poly_k+delta_k*f
						for (u = 0; u < interpoly_Ysize; u++)	//rs
							for (v = 0; v < interpoly_Xsize; v++)	//max(deg_y)+1
							{
								if (poly[j][u][v] != 0)
								{
									temp1 = mul(delta[j_min], poly[j][u][v]);
								}
								else
									temp1 = 0;

								if (f[u][v] != 0)
								{
									temp2 = mul(delta[j], f[u][v]);
								}
								else
									temp2 = 0;

								if (temp1 != 0 || temp2 != 0)
								{
									poly[j][u][v] = add(temp1, temp2);
								}
								else
									poly[j][u][v] = 0;
							}
					}
					else if (j == j_min)
					{
						for (u = 0; u < interpoly_Ysize; u++)	//max(deg_y)+1
						{
							for (v = 0; v < (interpoly_Xsize + 1); v++)	//w+2
								g1[u][v] = 0;
							for (v = 0; v < interpoly_Xsize; v++)	//w+1
								g2[u][v] = 0;
						}

						//g1=x*f
						for (u = 0; u < interpoly_Ysize; u++)	//rs
							for (v = 0; v < interpoly_Xsize; v++)	//max(deg_y)+1
								if (poly[j][u][v] != 0)
								{
									g1[u][v + 1] = poly[j][u][v];
								}

						//g2=xi*f
						for (u = 0; u < interpoly_Ysize; u++)	//rs
							for (v = 0; v < interpoly_Xsize; v++)	//max(deg_y)+1
								if (poly[j][u][v] != 0)
								{
									g2[u][v] = mul(com_elem_interpoint[i][0], poly[j][u][v]);
								}
						//poly=g1+g2
						for (u = 0; u < interpoly_Ysize; u++)	//rs
							for (v = 0; v < interpoly_Xsize; v++)	//max(deg_y)+1
								if (g1[u][v] != 0 || g2[u][v] != 0)
								{
									poly[j][u][v] = add(g1[u][v], g2[u][v]);
								}
								else
									poly[j][u][v] = 0;
					}
				}
			}
		}      
	}

	j_min = ChooseTheMinPoly(poly, init_polyNum);
	for (u = 0; u < interpoly_Ysize; ++u)
		for (v = 0; v < interpoly_Xsize; ++v)
			Q_interpoly_BF[0][u][v] = poly[j_min][u][v];

	//Produce GrayCode used by the ordering of BFinterpolation
	for (i = 0; i < test_vec_num; ++i)
	{
		GrayOrder[i] = i ^ (i >> 1);
	}

	//produce interpoint group for Backward interpolation
	for (i = 0; i<eta; ++i)
	{
		uncom_elem_interpoint[i][0] = x_ordered[n - eta + i];
		uncom_elem_interpoint[i][1] = large_vec[0][(int)test_set_ordered[1][n - eta + i]];
		uncom_elem_interpoint[i][2] = large_vec[1][(int)test_set_ordered[1][n - eta + i]];
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
		int Bpoint[2] = { uncom_elem_interpoint[location][0], uncom_elem_interpoint[location][1 + Blocation] };
		int Flocation = (result & GrayOrder[i]) ? 1 : 0;
		int Fpoint[2] = { uncom_elem_interpoint[location][0], uncom_elem_interpoint[location][1 + Flocation] };

		//Calculate each poly's leading order
		for (j = 0; j<init_polyNum; ++j)
		{
			lod_temp = 0;
			lod[j] = 0;
			for (u = 0; u<interpoly_Ysize; ++u)
				for (v = 0; v<interpoly_Xsize; ++v)
				{
					if (poly[j][u][v] != 0)
					{
						lod_temp = mono_order[u][v];
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
			for (u = 0; u<interpoly_Ysize; ++u)
				for (v = 0; v<interpoly_Xsize; ++v)
				{
					if (poly[j][u][v] != 0)
					{
						lod_temp = mono_order[u][v];
						if (lod_temp > lod[j])
							lod[j] = lod_temp;
					}
				}
		}


		//debug: print the lod after BF
#ifdef _PrintGroupLod_
		printf("\nAfter %dst BI, the lod =\t{%d,\t%d,\t%d}-->", i, lod[0], lod[1], lod[2]);
#endif

		//Forward Interpolation
		InterOnce_2(poly, Fpoint, poly_min_order, location, poly_min_order_all, location_all);

		//Calculate each poly's leading order
		for (j = 0; j<init_polyNum; ++j)
		{
			lod_temp = 0;
			lod[j] = 0;
			for (u = 0; u<interpoly_Ysize; ++u)
				for (v = 0; v<interpoly_Xsize; ++v)
				{
					if (poly[j][u][v] != 0)
					{
						lod_temp = mono_order[u][v];
						if (lod_temp > lod[j])
							lod[j] = lod_temp;
					}
				}
		}

		//debug: print the lod after BF
#ifdef _PrintGroupLod_
		printf("After %dst FI, the lod =\t{%d,\t%d,\t%d}", i, lod[0], lod[1], lod[2]);
#endif

		//find the min polynomial
		j_min = -1;
		j_min = ChooseTheMinPoly(poly, init_polyNum);
		for (u = 0; u<interpoly_Ysize; ++u)
			for (v = 0; v<interpoly_Xsize; ++v)
				Q_interpoly_BF[GrayOrder[i]][u][v] = poly[j_min][u][v];
	}
#ifdef _PrintGroupLod_
	printf("\n\n");
#endif


}

void InterOnce_2(int poly[][interpoly_Ysize][interpoly_Xsize], int interpoint[2], int poly_min_order[eta], int location, int poly_min_order_all[n], int location_all=0)
{
	int j, j_min, u, v, z;
	int lod[init_polyNum], lod_temp, lod_min, delta[init_polyNum], delta_temp, flag_lod_min;
	int J[init_polyNum];
	int f[interpoly_Ysize][interpoly_Xsize];
	int g1[interpoly_Ysize][interpoly_Xsize + 1];
	int g2[interpoly_Ysize][interpoly_Xsize];

	//Calculate each poly's leading order
	for (j = 0; j<init_polyNum; ++j)
	{
		lod_temp = 0;
		lod[j] = 0;
		for (u = 0; u<interpoly_Ysize; ++u)
			for (v = 0; v<interpoly_Xsize; ++v)
			{
				if (poly[j][u][v] != 0)
				{
					lod_temp = mono_order[u][v];
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
		for (u = 0; u<interpoly_Ysize; ++u)
		{
			delta_temp = 0;

			for (v = 0; v<interpoly_Xsize; ++v)
					if (poly[j][u][v] != 0)
					{
						int temp = power(interpoint[0], v);
						temp = mul(temp, poly[j][u][v]);
						//Hasse derivative mapping
						delta_temp = add(delta_temp, temp);
					}

			if (u == 0)
				delta[j] = delta_temp;
			else if (u>0)
			{
				delta_temp = mul(delta_temp, power(interpoint[1], u));
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
		for (u = 0; u<interpoly_Ysize; ++u)
			for (v = 0; v<interpoly_Xsize; ++v)
				f[u][v] = poly[j_min][u][v];

		//Modify nonzero polynomials
		for (j = 0; j<init_polyNum; ++j)
		{
			int temp1, temp2;

			if (J[j] != 0)
			{
				if (j != j_min)
				{
					//delta*poly_k+delta_k*f
					for (u = 0; u < interpoly_Ysize; ++u)
						for (v = 0; v < interpoly_Xsize; ++v)
						{
							if (poly[j][u][v] != 0)
								temp1 = mul(delta[j_min], poly[j][u][v]);
							else
								temp1 = 0;

							if (f[u][v] != 0)
								temp2 = mul(delta[j], f[u][v]);
							else
								temp2 = 0;

							if (temp1 != 0 || temp2 != 0)
								poly[j][u][v] = add(temp1, temp2);
							else
								poly[j][u][v] = 0;
						}
				}
				else if (j == j_min)
				{
					for (u = 0; u < interpoly_Ysize; u++)	//max(deg_y)+1
					{
						for (v = 0; v < (interpoly_Xsize + 1); v++)	//w+2
							g1[u][v] = 0;
						for (v = 0; v < interpoly_Xsize; v++)	//w+1
							g2[u][v] = 0;
					}


					//g1 = x * f
					for (u = 0; u < interpoly_Ysize; ++u)
						for (v = 0; v < interpoly_Xsize; ++v)
							if (poly[j][u][v] != 0)
								g1[u][v + 1] = poly[j][u][v];

					//g2 = xi * f
					for (u = 0; u < interpoly_Ysize; ++u)
						for (v = 0; v < interpoly_Xsize; ++v)
								if (poly[j][u][v] != 0)
									g2[u][v] = mul(interpoint[0], poly[j][u][v]);

					//poly = g1 + g2
					for (u = 0; u < interpoly_Ysize; ++u)
						for (v = 0; v < interpoly_Xsize; ++v)
							if (g1[u][v] != 0 || g2[u][v] != 0)
								poly[j][u][v] = add(g1[u][v], g2[u][v]);
							else
								poly[j][u][v] = 0;
				}
			}
		}
	}

}

int ChooseTheMinPoly(int polySour[][interpoly_Ysize][interpoly_Xsize], int PolyNum)
{
	int temp, index_temp;
	int *degree_temp = (int *)malloc(sizeof(int)*PolyNum);


	for (int j = 0; j < PolyNum; ++j)
	{
		temp = -1;
		degree_temp[j] = 0;
		for (int u = 0; u < interpoly_Ysize; ++u)
			for (int v = 0; v < interpoly_Xsize; ++v)
				if (polySour[j][u][v] != 0 && temp < mono_order[u][v])
					temp = mono_order[u][v];

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