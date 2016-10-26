#include <stdio.h>
#include "FiniteFieldBasisGF(4).h"
#include "main.h"
#include "NewInter.h"
#include "BackInter.h"




void test_LM_validity(int poly_input[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int point_BF[3], int point_FI[3], int poly_min_order[eta], int location)
{
	int i, j, u, v, z;
	int lod_temp, lod[init_polyNum];
	int empty1[n], empty2;
	int poly_result[init_polyNum][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize];

	for (i = 0; i < init_polyNum; ++i)
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
					poly_result[i][u][v][z] = poly_input[i][u][v][z];

	//before BF, calculate the lod
	for (j = 0; j<init_polyNum; ++j)
	{
		lod_temp = 0;
		lod[j] = 0;
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
				{
					if (poly_result[j][u][v][z] != 0)
					{
						lod_temp = MonoOrderConvert(u, v, z);
						if (lod_temp > lod[j])
							lod[j] = lod_temp;
					}
				}
	}	


	//BF point_BF
	/*poly_min_order[0] = 2;*/
	BackInterpolation_2(poly_result, point_BF, poly_min_order, location, lod, empty1, &empty2);	//Todo 修改形参表

	//calculate the order
	for (j = 0; j<init_polyNum; ++j)
	{
		lod_temp = 0;
		lod[j] = 0;
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
				{
					if (poly_result[j][u][v][z] != 0)
					{
						lod_temp = MonoOrderConvert(u, v, z);
						if (lod_temp > lod[j])
							lod[j] = lod_temp;
					}
				}
	}	



	//FI point_FI
	InterOnce_2(poly_result, point_FI, poly_min_order, location, empty1, empty2);

	//printf the LM of every poly in Group
	for (j = 0; j<init_polyNum; ++j)
	{
		lod_temp = 0;
		lod[j] = 0;
		for (u = 0; u<New_interpoly_Zsize; ++u)
			for (v = 0; v<New_interpoly_Ysize; ++v)
				for (z = 0; z<New_interpoly_Xsize; ++z)
				{
					if (poly_result[j][u][v][z] != 0)
					{
						lod_temp = MonoOrderConvert(u, v, z);
						if (lod_temp > lod[j])
							lod[j] = lod_temp;
					}
				}
	}

	printf("\ndebug_lod = [%d, %d, %d, %d]\n", lod[0], lod[1], lod[2], lod[3]);


}