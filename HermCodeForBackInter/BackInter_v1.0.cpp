#include "main.h"

void BackInterpolation(int &Q_com_poly[][intepoly_Ysize][interpoly_Xsize], int position)
{
	int i, j;
	int backPoint[3];
	int max = interpoly_Ysize-1;
	int A[init_polyNum];
	int inverA[init_polyNum];

	
	backPoint[0] = x_order[0][position];
	backPoint[1] = x_order[1][position];
	backPoint[2] = large_vec[1][(int)test_set_ordered[1][position]];
	memset(A, 1, sizeof(A));
	memset(inverA, 0, sizeof(A));
	//Equivalent Grobner basis transformation
	for(int a=0; a<=lm; ++a)
		for(int b=0; b<=max; ++b)
			if(ItsSize(A)>0)
			{
				int *qResult = ComputeResult(A, backPoint[0], lm-a, max-b, init_polyNum); 
				if(ItsSize(qResult)>0 && ItsSize(A)>0)
				{
					int minIndex;
					int LeadingOrder[init_polyNum];
					//caculate the LeadingOrder of each poly
					for(j=0; j<init_polyNum; ++j)
					{
						LeadingOrder[j] = -1;
						if(qResult[j]>0 && A[j]>0)
						{	
							LeadingOrder[j] = 0;
							for(int u=0; u<interpoly_Zsize; ++u)
								for(int v=0; v<interpoly_Ysize; ++v)
									for(int z=0; z<interpoly_Xsize; ++z)
										if(Q_com_poly[j][u][v][z]!=0 && mono_order[u][v][z]>LeadingOrder[j])
										{
											LeadingOrder[j] = mono_order[u][v][z];
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
					for(j=0; j<init_polyNum; ++j)
						if(qResult[j]>0 && A[j]>0)
						{
							//update poly
							Update2(qResult[j], Q_com_poly[j], qResult[minIndex], Q_com_poly[minIndex]);
						}
				}
			}

				
	//Generalized backward interpolation for m=1
	if(ItsSize(A)>0)
	{
		for(j=0; j<init_polyNum; ++j)
			if(A[j]>0)
			{
				UpdatePolyWithDivision(Q_com_poly[j], backPoint[0]);
			}
	}
}


int ItsSize(int *A)	//judge array is empty?
{
	int sizeA = sizeof(A)/sizeof(int);
	int count = 0;
	for(int i=0; i<sizeA; ++i)
		if(A[i]>0)
			++count;

	return count;
}

int *ComputeResult(int *A, int alpha, int Zindex, int Yindex, int polyNum)
{
	int *result = (int *)malloc(polyNum * sizeof(int));
	for(int i=0; i<polyNum; ++i)
	{
		result[i] = -1;
		if(A[i]!=1)
		{
			result[i] = 0;
			for(int z=0; z<=w ; ++z)
			{
				int temp = pow(alpha,z);
				if(temp!=0)
				{
					temp = mul(Q_com_poly[i][Zindex][Yindex][z],temp);
					if(temp!=0)
						result[i] = add(result[i],temp);
				}
			}
		}
	}
	return result;
}

void Update2(const int alpha, int &Q1[][interpoly_Ysize][interpoly_Xsize], const int beta, const int &Q2[][interpoly_Ysize][interpoly_Xsize])
{
	//calculate alpha * Q2 + beta * Q1
	int temp1, temp2;
	int g[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	//start
	for(int u=0; u<interpoly_Zsize; ++u)
		for(int v=0; v<interpoly_Ysize; ++v)
			for(int z=0; z<interpoly_Xsize; ++z)
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

void UpdatePolyWithDivision(int &Q[][interpoly_Ysize][interpoly_Xsize], int alpha);
{
	int X_poly[interpoly_Xsize], div_poly[interpoly_Xsize], poly[2];
	int divisor[2], dividend[2], quo[2];
	int valid_flag;
	int result[interpoly_Xsize];

	poly[0] = alpha;
	poly[1] = 1;

	for(int j=0; j<interpoly_Zsize; ++j)
		for(int i=0; i<interpoly_Ysize; ++i)
		{
			int v;
			for(v=0; v<interpoly_Xsize; ++v)
			{
				X_poly[v] = Q[j][i][v];
				result[v] = 0;
			}
			
			divisor[0] = 1;
			divisor[1] = 1;

			quo[0]=quo[1]=0;
			valid_flag=0;
			for(int u=0; u<interpoly_Xsize-1; ++u)
			{
				for(v=0; v<interpoly_Xsize; ++v)
					div_poly[v] = 0;

				//find the dividend[2]
				dividend[0] = 0;
				dividend[1] = 0;
				for(v=interpoly_Xsize-1; v>=0; --v)
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
				quo[0] = mul(dividend[0], inv(divisor[0]);
				quo[1] = dividend[1] - divisor[1];

				if((0<=quo[1]) && (quo[1]<=interpoly_Xsize-1))
				{
					result[quo[1]] = quo[0];	
					
					//div_poly = divisor*quotient
					for(v=0; v<2; ++v)
						if(poly[v]!=0)
							div_poly[v+quo[1]] = mul(quo[0],poly[v]);

					//dividend-div_poly
					for(v=0; v<interpoly_Xsize; ++v)
						X_poly[v] = add(X_poly[v],div_poly[v]);
				}
				else if(quo[1]<0 || quo[1]>(interpoly_Xsize-1))
				{
					valid_flag = 0;
					break;
				}
			}
			if(valid_flag==0)
				printf("\n\ndivision (x+ai) failed\n\n");
			
			//if division is effective, then update the Q
			for(v=0; v<interpoly_Xsize; v++)
				Q[j][i][v] = result[v];
		}

}

