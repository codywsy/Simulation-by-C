#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "FiniteFieldBasisGF(8).h"
#include "FiniteFieldBasisCal.h"
#include "BackInter.h"

//from Herm(8,4))_BackwardInter_v1.0.cpp
extern int x_ordered[2][n];
extern int large_vec[choose_num][n];
extern float test_set_ordered[2][n];
extern int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
extern unsigned long int seq_num_Now;

//function: realize the Backinterpolation for Hermition code
void BackInterpolation(int Q_com_poly[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int backPoint[3])
{
	int j;
//	int backPoint[3];
	int max = w-1;
	int A[init_polyNum];
	int inverA[init_polyNum];
	int g[init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	memset(inverA, 0, sizeof(inverA));
	/*memset(A, 1, sizeof(int)*init_polyNum);*/
	for (int u = 0; u < init_polyNum; ++u)
		A[u] = 1;

	//Equivalent Grobner basis transformation
	for(int a=0; a<(lm+1); ++a)
		for(int b=0; b<=max; ++b)
			if(ItsSize(A, sizeof(A)/sizeof(int))>0)
			{
				int *qResult = ComputeResult(Q_com_poly, A, backPoint[0], lm-a, max-b, init_polyNum);	//reverse order
				//int *qResult = ComputeResult(Q_com_poly, A, backPoint[0], a, b, init_polyNum);	//order
				int temp_debug = ItsSize(qResult, init_polyNum);
				if (temp_debug>0 && ItsSize(A, sizeof(A) / sizeof(int))>0)
				{
					int minIndex;
					int LeadingOrder[init_polyNum];
					////caculate the LeadingOrder of each poly
					//for(j=0; j<init_polyNum; ++j)
					//{
					//	LeadingOrder[j] = -1;
					//	if(qResult[j]>0 && A[j]>0)
					//	{	
					//		LeadingOrder[j] = 0;
					//		for(int u=0; u<New_interpoly_Zsize; ++u)
					//			for(int v=0; v<New_interpoly_Ysize; ++v)
					//				for (int z = 0; z<New_interpoly_Xsize; ++z)
					//					if (Q_com_poly[j][u][v][z] != 0 && MonoOrderConvert(u,v,z)>LeadingOrder[j])
					//					{
					//						LeadingOrder[j] = MonoOrderConvert(u, v, z);
					//					}
					//	}
					//}

					//Calculate each poly's leading order
					//1st, poly converts to g
					for (j = 0; j<init_polyNum; ++j)
						for (int u = 0; u<interpoly_Zsize; ++u)
							for (int v = 0; v<interpoly_Ysize; ++v)
								for (int z = 0; z<interpoly_Xsize; ++z)
									g[j][u][v][z] = 0;

					ConvertX2Y(Q_com_poly, g, init_polyNum, w);

					//2nd, use g to calculate the lod[]
					for (j = 0; j<init_polyNum; ++j)
					{
						int lod_temp = 0;
						LeadingOrder[j] = 0;
						for (int u = 0; u<New_interpoly_Zsize; ++u)
							for (int v = 0; v<New_interpoly_Ysize; ++v)
								for (int z = 0; z<New_interpoly_Xsize; ++z)
								{
									if (g[j][u][v][z] != 0)
									{
										lod_temp = mono_order[u][v][z];
										if (lod_temp > LeadingOrder[j])
											LeadingOrder[j] = lod_temp;
									}
								}
					}


					//find the index with the min leading order
					minIndex = -1;
					for(j=0; j<init_polyNum; ++j)
					{
						if(qResult[j]>0 && A[j]>0)
						{
							if(minIndex==-1)
							{
								minIndex = j;
								continue;
							}
							else if(LeadingOrder[j] < LeadingOrder[minIndex])
							{
								minIndex = j;
							}
						}
					}
					if(minIndex==-1)
						printf("\n\nBackInterpolation has error\n\n");
					
					//update A and inverA
					A[minIndex] = 0;
					inverA[minIndex] = 1;
					
					//update the poly group
					for (j = 0; j < init_polyNum; ++j)
					{
						if(qResult[j]>0 && A[j]>0)
						{
							//update poly
							BF_Update(qResult[j], Q_com_poly[j], qResult[minIndex], Q_com_poly[minIndex]);
						}
					}
				}
			}

				
	//Generalized backward interpolation for m=1
	if (ItsSize(A, sizeof(A) / sizeof(int))>0)
	{
		for(j=0; j<init_polyNum; ++j)
			if(A[j]>0)
			{
				UpdatePolyWithDivision(Q_com_poly[j], backPoint[0]);
			}
	}
}

/******************
Function(todo):
BackInterpolation with choosing only one polynomial
********************/
void BackInterpolation_2(int Q_com_poly[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int backPoint[3], int poly_min_order[eta], int location, int lod[init_polyNum], int poly_min_order_all[n], int *location_all)
{
	int j;
	//	int backPoint[3];
	int max = w - 1;
	int A[init_polyNum];
	int inverA[init_polyNum];
	int g[init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	memset(inverA, 0, sizeof(inverA));
	/*memset(A, 1, sizeof(int)*init_polyNum);*/
	for (int u = 0; u < init_polyNum; ++u)
		A[u] = 1;

	//Equivalent Grobner basis transformation
	for (int a = 0; a<(lm + 1); ++a)
		for (int b = 0; b <= max; ++b)
			if (ItsSize(A, sizeof(A) / sizeof(int))>0)
			{
				int *qResult = ComputeResult(Q_com_poly, A, backPoint[0], lm - a, max - b, init_polyNum);	//reverse order
				//int *qResult = ComputeResult(Q_com_poly, A, backPoint[0], a, b, init_polyNum);	//order
				int temp_debug = ItsSize(qResult, init_polyNum);
				if (temp_debug>0 && ItsSize(A, sizeof(A) / sizeof(int))>0)
				{
					int minIndex;
					int LeadingOrder[init_polyNum];
					////caculate the LeadingOrder of each poly
					//for(j=0; j<init_polyNum; ++j)
					//{
					//	LeadingOrder[j] = -1;
					//	if(qResult[j]>0 && A[j]>0)
					//	{	
					//		LeadingOrder[j] = 0;
					//		for(int u=0; u<New_interpoly_Zsize; ++u)
					//			for(int v=0; v<New_interpoly_Ysize; ++v)
					//				for (int z = 0; z<New_interpoly_Xsize; ++z)
					//					if (Q_com_poly[j][u][v][z] != 0 && MonoOrderConvert(u,v,z)>LeadingOrder[j])
					//					{
					//						LeadingOrder[j] = MonoOrderConvert(u, v, z);
					//					}
					//	}
					//}

					//Calculate each poly's leading order
					//1st, poly converts to g
					for (j = 0; j<init_polyNum; ++j)
						for (int u = 0; u<interpoly_Zsize; ++u)
							for (int v = 0; v<interpoly_Ysize; ++v)
								for (int z = 0; z<interpoly_Xsize; ++z)
									g[j][u][v][z] = 0;

					ConvertX2Y(Q_com_poly, g, init_polyNum, w);

					//2nd, use g to calculate the lod[]
					for (j = 0; j<init_polyNum; ++j)
					{
						int lod_temp = 0;
						LeadingOrder[j] = 0;
						for (int u = 0; u<New_interpoly_Zsize; ++u)
							for (int v = 0; v<New_interpoly_Ysize; ++v)
								for (int z = 0; z<New_interpoly_Xsize; ++z)
								{
									if (g[j][u][v][z] != 0)
									{
										lod_temp = mono_order[u][v][z];
										if (lod_temp > LeadingOrder[j])
											LeadingOrder[j] = lod_temp;
									}
								}
					}


					//find the index with the min leading order
					minIndex = -1;
					for (j = 0; j<init_polyNum; ++j)
					{
						if (qResult[j]>0 && A[j]>0)
						{
							if (minIndex == -1)
							{
								minIndex = j;
								continue;
							}
							else if (LeadingOrder[j] < LeadingOrder[minIndex])
							{
								minIndex = j;
							}
						}
					}
					if (minIndex == -1)
						printf("\n\nBackInterpolation has error\n\n");

					//update A and inverA
					A[minIndex] = 0;
					inverA[minIndex] = 1;

					//update the poly group
					for (j = 0; j < init_polyNum; ++j)
					{
						//if (seq_num_Now != 13)
						//{
							if (qResult[j] > 0 && A[j]>0)		//here, we choose the poly in A with qResult!=0 to update
							{
								//update poly
								BF_Update(qResult[j], Q_com_poly[j], qResult[minIndex], Q_com_poly[minIndex]);
							}
						//}
						//else if (seq_num_Now == 13)
						//{
						//	if (A[j] > 0)
						//	{
						//		//update poly
						//		BF_Update(qResult[j], Q_com_poly[j], qResult[minIndex], Q_com_poly[minIndex]);
						//	}
						//}
					}
				}
			}


	//Generalized backward interpolation for m=1
	int count_temp = ItsSize(A, sizeof(A) / sizeof(int));
	if (count_temp > 0)
		if (count_temp == 1)	//Normally, there is only one poly in A
		{
			*location_all = 0;
			for (j = 0; j < init_polyNum; ++j)
				if (A[j]>0)
				{
					UpdatePolyWithDivision(Q_com_poly[j], backPoint[0]);
					break;
				}
		}
		else if (count_temp > 1)	//When BF applying in Herm code, it's possible to have more than one poly in A
		{
			int flag = 0;
			*location_all = 0;
			//plan1, use poly_min_order[location]
			for (j = 0; j < init_polyNum; ++j)
				if (A[j]>0 && j == poly_min_order[location])
				{
					UpdatePolyWithDivision(Q_com_poly[j], backPoint[0]);
					flag = 1;
					break;
				}
			/*	else if (j==init_polyNum-1)
					printf("\n[errorInBI]\tseq_num=%d, There is no j == poly_min_order[location]\n", seq_num_Now);*/
			
			//plan2, choose the poly in A with max lod to UpdatePolyWithDivision
			if (flag != 1)
			{
				int lod_max = -1;
				int index_max = -1;
				for (j = 0; j < init_polyNum; ++j)
					if (A[j]>0 && lod[j] > lod_max)
					{
						lod_max = lod[j];
						index_max = j;
					}

				UpdatePolyWithDivision(Q_com_poly[index_max], backPoint[0]);
			}

			//plan3, choose the poly in A with current j_min record
			//int flag_temp = 0;
			//for (int i = n-1; i >= 0; --i)
			//{
			//	for (j = 0; j < init_polyNum; ++j)
			//		if (A[j] != 0 && j == poly_min_order_all[i])
			//		{
			//			*location_all = i;
			//			UpdatePolyWithDivision(Q_com_poly[j], backPoint[0]);
			//			flag_temp = 1;
			//			break;
			//		}

			//	if (flag_temp == 1)
			//		break;
			//}
		}
}



/******************
Function:
	Judge if array is empty
********************/
int ItsSize(int *A, int sizeA)
{
	int count = 0;
	for(int i=0; i<sizeA; ++i)
		if(A[i]>0)
			++count;

	return count;
}

int *ComputeResult(const int Q[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int *A, int alpha, int Zindex, int Yindex, int polyNum)
{
	int *result = (int *)malloc(polyNum * sizeof(int));
	for(int i=0; i<polyNum; ++i)
	{
		result[i] = -1;
		if(A[i]==1)
		{
			result[i] = 0;
			for(int z=0; z < New_interpoly_Xsize ; ++z)
				if (Q[i][Zindex][Yindex][z]!=0)
				{
					int temp = power(alpha,z);
					if(temp!=0)
					{
						temp = mul(Q[i][Zindex][Yindex][z],temp);
						result[i] = add(result[i],temp);
					}
				}
		}
	}
	return result;
}

void BF_Update(const int alpha, int Q1[][New_interpoly_Ysize][New_interpoly_Xsize], const int beta, const int Q2[][New_interpoly_Ysize][New_interpoly_Xsize])
{
	//calculate alpha * Q2 + beta * Q1
	int temp1, temp2;
//	int g[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	//start
	for(int u=0; u<New_interpoly_Zsize; ++u)
		for(int v=0; v<New_interpoly_Ysize; ++v)
			for(int z=0; z<New_interpoly_Xsize; ++z)
			{
				if(Q1[u][v][z]!=0)
					temp1 = mul(beta,Q1[u][v][z]);
				else
					temp1 = 0;

				if(Q2[u][v][z]!=0)
					temp2 = mul(alpha,Q2[u][v][z]);
				else
					temp2 = 0;

				if(temp1!=0 || temp2!=0)
					Q1[u][v][z] = add(temp1,temp2);
				else
					Q1[u][v][z] = 0;
			}	

}

void UpdatePolyWithDivision(int Q[][New_interpoly_Ysize][New_interpoly_Xsize], int alpha)
{
	int X_poly[New_interpoly_Xsize], div_poly[New_interpoly_Xsize], poly[2];
	int divisor[2], dividend[2], quo[2];
	int valid_flag;
	int result[New_interpoly_Xsize];

	poly[0] = alpha;
	poly[1] = 1;

	for(int j=0; j<New_interpoly_Zsize; ++j)
		for(int i=0; i<New_interpoly_Ysize; ++i)
		{
			int v;
			for(v=0; v<New_interpoly_Xsize; ++v)
			{
				X_poly[v] = Q[j][i][v];
				result[v] = 0;
			}
			
			divisor[0] = 1;
			divisor[1] = 1;

			valid_flag = 0;
			quo[0] = quo[1] = 0;
			for(int u=0; u<New_interpoly_Xsize-1; ++u)
			{
				for(v=0; v<New_interpoly_Xsize; ++v)
					div_poly[v] = 0;

				//find the dividend[2]
				dividend[0] = 0;
				dividend[1] = 0;
				for(v=New_interpoly_Xsize-1; v>=0; --v)
					if(X_poly[v]!=0)
					{
						dividend[0] = X_poly[v];
						dividend[1] = v;
						break;
					}
				
				//effcetive solution, no eemainder
				if(v==-1 && X_poly[0]==0)
				{
					valid_flag = 1;
					break;
				}
				
				//uneffective solution, remainder
				if(dividend[1] < divisor[1])
				{
					valid_flag = 0;
					break;
				}

				//find the quotient
				quo[0] = mul(dividend[0], inv(divisor[0]));
				quo[1] = dividend[1] - divisor[1];

				if((0<=quo[1]) && (quo[1]<=New_interpoly_Xsize-1))
				{
					result[quo[1]] = quo[0];	
					
					//div_poly = divisor*quotient
					for(v=0; v<2; ++v)
						if(poly[v]!=0)
							div_poly[v+quo[1]] = mul(quo[0],poly[v]);

					//dividend-div_poly
					for(v=0; v<New_interpoly_Xsize; ++v)
						X_poly[v] = add(X_poly[v],div_poly[v]);
				}
				else if(quo[1]<0 || quo[1]>(New_interpoly_Xsize-1))
				{
					valid_flag = 0;
					break;
				}
			}
			if(valid_flag==0)
				printf("\n\ndivision (x+ai) failed\n\n");
			
			//if division is effective, then update the Q
			for(v=0; v<New_interpoly_Xsize; v++)
				Q[j][i][v] = result[v];
		}

}

/***********************************************************************
Function:	
	used to dectect the polynomial of group who has factor (x+a_i) 
Argument:	
	Q_poly	---	polynomial group
	polyNum	---	the number of polynomial
	ai		--- the x value of interpoint 
************************************************************************/
void DectectIfFactorInPoly(int Q_poly[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int polyNum, int xi)
{
	for (int i = 0; i < polyNum; i++)
	{
		for (int u = 0; u < interpoly_Zsize; ++u)
			for (int v = 0; v < interpoly_Ysize; ++v)
			{
				int temp = Q_poly[i][u][v][0];
				for (int z = 1; z < interpoly_Xsize; ++z)
					if (Q_poly[i][u][v][z] != 0)
					{
						int temp1 = power(xi, z);
						temp1 = mul(temp1, Q_poly[i][u][v][z]);
						temp = add(temp, temp1);
					}

				if (temp != 0)
				{
					printf("\nQ[%d] does not have factor", i);
					goto out;
				}
			}

		printf("\nQ[%d] have factor", i);

		out: ;
	}

	printf("\n");
}