#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "main.h"
#include "FiniteFieldBasisGF(8).h"
#include "FiniteFieldBasisCal.h"
#include "BackInter.h"
#include "BF_LCC_Inter.h"

#define _NoReductionCom_
#define _NoReductionUncom_
#define OpenFile fp=fopen("BFLCC_RS(7,2).txt","a")
#define FrameError 309

//**************basic var*****************
float start, finish;
unsigned long int seq_num;	//number of input binary sequences
float SNR;
double BER, FER;
double BER_BF, FER_BF;
float N0, sgm;
const float pi = 3.141593;	// pai

int bi_message[k*p], message[k];	//transmitted messge[k]
int codeword[n], bi_codeword[n*p]; //codewords
//demodulation()
float tx_symbol[p*n][2], rx_symbol[p*n][2];
float RM[q][n];

//mono_table()
int mono_order[monoTable_Ysize][monoTable_Xsize];

//generator()
int gmatrix[k][n];	//generator matrix[k][n]
int dec_message[k], dec_codeword[n], dec_bicodeword[n*p]; //decoding result

//test_vec_construction()
int large_vec[choose_num][n];
float test_set_ordered[2][n];
int x_ordered[n];

//interpolation()
int Q_interpoly[test_vec_num][interpoly_Ysize][interpoly_Xsize];	//chosen interpolated polynomial [z-deg+1][max(deg_y)+1][w+1]



//factorization() and rcs()
//int Q[k][facpoly_Ysize][facpoly_Xsize];	//sequtial deduction polynomial
int uu;	//factorisation step index
int l, listNum[test_vec_num];	//candidate output index
int output[lm + 1][k], outputList[test_vec_num][lm + 1][k];	//the list of candidate message [expeced number of candidate message, >rs][length of message, k]

//choose()
int epcount1, testCount1[test_vec_num], testCount2, testCount1_com;

//**************debug var****************
unsigned long int seq_num_Now;	//record the seq_num at the moment
//some var about computation complexity
unsigned long int addNum, mulNum, addNum_temp, mulNum_temp;
int flag_addNum = 0, flag_mulNum = 0;

int degree_test_temp[test_vec_num][init_polyNum];  //store the degree of leading monomial of polynomial
int degree_test[test_vec_num];		//store the degree of the choosen polynomial

//**************BF var********************
//use for BackInter
int Q_interpoly_BF[test_vec_num][interpoly_Ysize][interpoly_Xsize];
int dec_codeword_BF[n], dec_bicodeword_BF[n*p];
int outputList_BF[test_vec_num][lm + 1][k];
int listNum_BF[test_vec_num];



//***************my function***************** 
void mono_table(void);
void test_vec_construction(void);
void interpolation(void);
void com_elem_interpolation(int g[][interpoly_Ysize][interpoly_Xsize], int interpoint[][n - eta]);
void uncom_elem_interpolation(int interpoint[3], int inGroup[][interpoly_Ysize][interpoly_Xsize], int outGroup1[][interpoly_Ysize][interpoly_Xsize], int outGroup2[][interpoly_Ysize][interpoly_Xsize]);
//void factorisation(int Q_input[test_vec_num][interpoly_Ysize][interpoly_Xsize], int output_list[test_vec_num][lm + 1][k], int list_num[test_vec_num]);
void rcs(int);
void choose(int output_codeword[n], int output_bicodeword[n*p], int output_list[test_vec_num][lm + 1][k], int list_num[test_vec_num], int flag_decoding_alg);
int result_compare(int output_list[test_vec_num][lm + 1][k], int output_codeword[n], int output_list_BF[test_vec_num][lm + 1][k], int output_codeword_BF[n]);
//*********************************

//JIONGYUE's factorization
#define poly_X_size interpoly_Xsize//larger than C
#define poly_Y_size interpoly_Ysize//Ysize=lm+1
#define fac_size (lm+1)*k//the probable numbers of factorisation polynomials, bigger than lm*k
#define coeff_size fac_size//the probable numbers of coefficients of factorisation polynomials, bigger than lm*k

//factorisation
int Q[fac_size][poly_Y_size][poly_X_size + lm + 1];//because Q_(s+1)(x,y)=Q_s(x,xy+p_s), the X_size should increase to C+lm
int times;//factorisation times
int pai[fac_size];//pai[u]:the parent vertex of u
int deg[fac_size];//deg[u]:degree of u=distance from root-1
int coeff[coeff_size];//the coefficients of factorisation polynomials
int fx[lm + 1][k];//the factorisation polynomials, the number of it is not sure, maybe overflow
int f_num;//the number of factorisation polynomials, it is not sure
int fac_test_vec_mark;

//factorisation step
void factorisation(int Q_input[poly_Y_size][poly_X_size]);
void DFS(int u);
void out_fac_poly();
int find_root(int u,int i);
void update_fac_poly(int v,int u,int alpha);
void output_fac_poly(int u);
int combination(int a, int b);

//****basic function statement****
void generator(void);
void encoder(int message_temp[], int codeword_temp[]);
void modulation(void);
void channel(void);
void demodulation(void);
float PDF(int, int);
//*********************************

void main()
{
	int i, u, m, num, value;
	unsigned int mask = 1;
	long int error, ferror;
	long int error_BF, ferror_BF;
	double progress, b;
	unsigned long int j;
	double addNum_count, mulNum_count, totalNum_count;

	FILE *fp;
	if ((OpenFile) == NULL)
	{
		printf("Can't open the data.txt\n");
		exit(0);
	}
	fclose(fp);

	srand(1977);

	mono_table();

	generator();

	//*****input data from basic_input.txt*********
	FILE * fp_input = fopen("basic_input.txt", "r");
	if (fp_input != NULL) {
		fscanf(fp_input, "start SNR:%f\n", &start);
		fscanf(fp_input, "finish SNR:%f\n", &finish);
		fscanf(fp_input, "seq_num:%ld", &seq_num);

		fclose(fp_input);
	}
	else{
		printf("\n\ncan't read the basic_input file\n\n");
	}

	printf("start SNR: %0.2f\n\n", start);
	printf("finish SNR: %0.2f\n\n", finish);
	printf("seq_num: %d\n", seq_num);
	//*******************************

	for (SNR = start; SNR <= finish; SNR += interval)
	{
		N0 = (1.0 / (float(k) / float(n))) / pow(10.0, SNR / 10.0);
		sgm = sqrt(N0 / 2);
		b = 1.0;

		error = 0;
		ferror = 0;
		error_BF = 0;
		ferror_BF = 0;

		addNum_count = 0.0;
		mulNum_count = 0.0;
		totalNum_count = 0.0;

		for (j = 1; j <= seq_num; j++)
		{
			//*****debug*****
			seq_num_Now = j;
			//***************

			//generate binary input sequence
			for (u = 0; u < k*p; u++)	//k*4
				bi_message[u] = rand() % 2;

			//convert to decimal input sequence
			for (i = 0; i < k; i++)
			{
				num = 1;
				message[i] = 0;
				for (u = 0; u < p; u++)
				{
					message[i] = message[i] + (bi_message[p*i + u] * num);
					num = num * 2;
				}
			}

			encoder(message, codeword);

			//convert to binary 
			for (u = 0; u < n*p; u++)	//n*4
				bi_codeword[u] = 0;

			for (u = 0; u < n; u++)	//n
			{
				value = codeword[u];
				mask = 1;
				for (m = 0; m < p; m++)
				{
					if ((value & mask) > 0)
						bi_codeword[p*u + m] = 1;
					else
						bi_codeword[p*u + m] = 0;
					mask = mask << 1;
				}
			}


			modulation();

			channel();

			demodulation(); // output soft information

			test_vec_construction();


			//BF_LCC interpolation
			BF_LCC_inter();
			//BF_LCC factorisation
			for (u = 0; u < test_vec_num; u++)
			{
				factorisation(Q_interpoly_BF[u]);
				for (int h = 0; h<(lm + 1); h++)
					for (int z = 0; z<k; z++)
						outputList_BF[u][h][z] = fx[h][z];

				listNum_BF[u] = f_num;
			}
			//choose
			choose(dec_codeword_BF, dec_bicodeword_BF, outputList_BF, listNum_BF, 1);

			//bit error rate calculation, BF-LCC
			int temp_BF = error_BF;
			for (u = 0; u<n*p; u++)
			{
				if (bi_codeword[u] != dec_bicodeword_BF[u])
					error_BF++;
			}

			//frame error rate calculation - BF-LCC
			if (error_BF > temp_BF)	//BF-LCC decoding error
				ferror_BF++;

			BER_BF = (double)(error_BF) / (double)(n*p*j);
			FER_BF = (double)(ferror_BF) / (double)(j);

			//********LCC interpolation*****************
			//interpolation
			interpolation();
			//factorisation
			for (u = 0; u < test_vec_num; u++)
			{
				factorisation(Q_interpoly[u]);
				for (int h = 0; h<(lm + 1); h++)
					for (int z = 0; z<k; z++)
						outputList[u][h][z] = fx[h][z];

				listNum[u] = f_num;

			}
			//choose
			choose(dec_codeword, dec_bicodeword, outputList, listNum, 0);

			//bit error rate calculation
			int temp = error;
			for (i = 0; i < p*n; i++)
			{
				if (bi_codeword[i] != dec_bicodeword[i])
				{
					error++;
				}
			}

			//frame error rate calculation
			if (error > temp)
				ferror++;

			progress = (double)(j * 100) / (double)seq_num;
			BER = (double)(error) / (double)(n*p*j);
			FER = (double)(ferror) / (double)(j);




			//compare result
			if (int compare_result = result_compare(outputList, dec_codeword, outputList_BF, dec_codeword_BF))
				printf("\nseq_num=%d:\tBF_interpolation has error code %d\n", seq_num_Now, compare_result);


			printf("Progress=%0.1f, seq_num_Now=%d, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, BER_BF=%E, FER_BF=%E\r", progress, seq_num_Now, SNR, error, BER, ferror, FER, BER_BF, FER_BF);

			if (ferror > FrameError)
				break;
		}

		if (ferror > FrameError)
		{
			BER = (double)error / (double)(n*p*j);
			FER = (double)(ferror) / (double)(j);
		}
		else
		{
			BER = (double)error / (double)(n*p*seq_num);
			FER = (double)(ferror) / (double)(seq_num);
		}

		printf("Progress=%0.1f, seq_num_Now=%d, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, BER_BF=%E, FER_BF=%E\n", progress, seq_num_Now, SNR, error, BER, ferror, FER, BER_BF, FER_BF);

		OpenFile;
		fprintf(fp, "Progress=%0.1f, seq_num_Now=%d, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, BER_BF=%E, FER_BF=%E\n", progress, seq_num_Now, SNR, error, BER, ferror, FER, BER_BF, FER_BF);
		fclose(fp);
	}

}


void interpolation()
{
	int i, j, u, v, z, num, temp, temp_index;
	int com_elem_interpoint[2][(n - eta)];
	int Q_com_elem[init_polyNum][interpoly_Ysize][interpoly_Xsize];
	int degree_temp[test_vec_num][init_polyNum];
	int Q_uncom_elem[test_vec_num][init_polyNum][interpoly_Ysize][interpoly_Xsize];

	//set common element interpoint (xi,yi)
	for (i = 0; i<n - eta; i++)
	{
		com_elem_interpoint[0][i] = x_ordered[i];
		com_elem_interpoint[1][i] = large_vec[0][(int)test_set_ordered[1][i]];
	}


	for (i = 0; i<init_polyNum; i++)	//num_polys
			for (u = 0; u<interpoly_Ysize; u++)	//max(deg_y)+1
				for (v = 0; v<interpoly_Xsize; v++)	//w+1
					Q_com_elem[i][u][v] = 0;

	com_elem_interpolation(Q_com_elem, com_elem_interpoint);

#ifndef _GS_Normal_
	//uncommon element interpolation
	int uncom_elem_interpoint[eta][3];
	//set uncomon element interpolation
	for (i = n - eta; i<n; i++)
	{
		uncom_elem_interpoint[i - n + eta][0] = x_ordered[i];
		uncom_elem_interpoint[i - n + eta][1] = large_vec[0][(int)test_set_ordered[1][i]];
		uncom_elem_interpoint[i - n + eta][2] = large_vec[1][(int)test_set_ordered[1][i]];
	}

	//initialize
	for (z = 0; z<test_vec_num; z++)
		for (i = 0; i<init_polyNum; i++)
				for (u = 0; u<interpoly_Ysize; u++)
					for (v = 0; v<interpoly_Xsize; v++)
						Q_uncom_elem[z][i][u][v] = 0;

	for (i = 0; i<init_polyNum; i++)
			for (u = 0; u<interpoly_Ysize; u++)
				for (v = 0; v<interpoly_Xsize; v++)
					Q_uncom_elem[0][i][u][v] = Q_com_elem[i][u][v];

	//uncommon elem interpolation
	for (i = 0; i<eta; i++)
	{
		num = (int)pow(2.0, i);
		for (u = num - 1; u >= 0; u--)
			uncom_elem_interpolation(uncom_elem_interpoint[i], Q_uncom_elem[u], Q_uncom_elem[2 * u + 0], Q_uncom_elem[2 * u + 1]);
	}

	
	//choose the poly to factorization

	//initialize
	for (i = 0; i<test_vec_num; i++)
		for (u = 0; u<interpoly_Ysize; u++)	//max(deg_y)+1
			for (v = 0; v<interpoly_Xsize; v++)	//w+1
				Q_interpoly[i][u][v] = 0;

	//eta>0
#ifdef _PrintGroupLod_
	printf("\n");
#endif
	for (i = 0; i < test_vec_num; i++)
	{
		//calculate the degree of poly
		for (j = 0; j < init_polyNum; j++)
		{
			temp = -1;
			degree_temp[i][j] = 0;
			//*******debug*******
			degree_test_temp[i][j] = 0;
			//******************
			for (u = 0; u < interpoly_Ysize; u++)	//rs
				for (v = 0; v < interpoly_Xsize; v++)	//max(deg_y)+1
					if (Q_uncom_elem[i][j][u][v] != 0)
						if (temp < mono_order[u][v])
						{
							temp = mono_order[u][v];

							//***debug**************
							degree_test_temp[i][j] = v + u*(k - 1);
							//***********************
						}
			degree_temp[i][j] = temp;
		}

		//choose the min degree
		temp_index = 0;
		temp = degree_temp[i][0];
		for (j = 1; j<init_polyNum; j++)
			if (temp>degree_temp[i][j] && degree_temp[i][j] <= iterNum)
			{
				temp = degree_temp[i][j];
				temp_index = j;

			}

		//****debug*******
		degree_test[i] = degree_test_temp[i][temp_index];

#ifdef _PrintGroupLod_
		printf("\ntest_vector[%d]'s lod = [%d,\t%d,\t%d]", i, degree_temp[i][0], degree_temp[i][1], degree_temp[i][2]);
#endif
		//******************

		//assignment
		for (u = 0; u < interpoly_Ysize; u++)	//rs
			for (v = 0; v < interpoly_Xsize; v++)	//max(deg_y)+1
				Q_interpoly[i][u][v] = Q_uncom_elem[i][temp_index][u][v];
	}
#ifdef _PrintGroupLod_
	printf("\n");
#endif
	//******************************
#else
	//eta=0
	for (i = 0; i<test_vec_num; i++)
	{
		//calculate the degree of poly
		for (j = 0; j<init_polyNum; j++)
		{
			temp = -1;
			degree_temp[i][j] = 0;
			for (u = 0; u<interpoly_Ysize; u++)	//rs
				for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
					if (Q_com_elem[j][u][v] != 0)
						if (temp < mono_order[u][v])
							temp = mono_order[u][v];

			degree_temp[i][j] = temp;
		}

		//choose the min degree
		temp_index = 0;
		temp = degree_temp[i][0];
		for (j = 1; j<init_polyNum; j++)
			if (temp>degree_temp[i][j])
			{
				temp = degree_temp[i][j];
				temp_index = j;
			}

		degree_test[i] = degree_temp[i][temp_index];

		//assignment
		for (u = 0; u<interpoly_Ysize; u++)	//rs
			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
				Q_interpoly[i][u][v] = Q_com_elem[temp_index][u][v];
	}
	//********************
#endif

	//********debug*************
	//used to prove the polynomials choosen for fac equals to 0 over all the interpoint
	//prove the efficiencies of the interpolation
	//printf("\ninterpolation check procesing result:\n");
	int temp_x, temp_y;
	for (i = 0; i<test_vec_num; i++)
	{
		for (j = 0; j<n - eta; j++)
		{
			temp = 0;
			for (u = 0; u<interpoly_Ysize; u++)	//rs
				for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
					if (Q_interpoly[i][u][v] != 0)
					{
						temp_x = power(com_elem_interpoint[0][j], v);
						temp_y = power(com_elem_interpoint[1][j], u);
						temp = add(temp, mul(Q_interpoly[i][u][v], mul(temp_y, temp_x)));
					}
			if (temp != 0)
				printf("\n[interpolation_error]Q_interpoly[%d] in point[%d](%d,%d) = %d", i, j, com_elem_interpoint[0][j], com_elem_interpoint[1][j], temp);
		}

		unsigned int mask = (int)pow(2.0, eta-1);
		int value = 0;

		for (j = 0; j<eta; j++)
		{
			value = (i&mask)>0 ? 1 : 0;
			mask = mask >> 1;

			temp = 0;
			for (u = 0; u<interpoly_Ysize; u++)	//max(deg_y)+1
				for (v = 0; v<interpoly_Xsize; v++)	//w+1
					if (Q_interpoly[i][u][v] != 0)
					{
						temp_x = power(uncom_elem_interpoint[j][0], v);
						temp_y = power(uncom_elem_interpoint[j][value+1], u);
						temp = add(temp, mul(Q_interpoly[i][u][v], mul(temp_y, temp_x)));
					}
			if (temp != 0)
				printf("\n[interpolation_error]Q_interpoly[%d] in point[%d](%d,%d) = %d\n", i, j+n-eta, uncom_elem_interpoint[j][0], uncom_elem_interpoint[j][value+1], temp);
		}

	}
//	printf("\n");
	//**************************

}

void com_elem_interpolation(int g[][interpoly_Ysize][interpoly_Xsize], int interpoint[][n - eta])
{
	int i, j, u, v, z, delta_temp, delta_temp1, lod_temp, temp1, temp2;
	int delta[init_polyNum], J[init_polyNum], act[init_polyNum], lod[init_polyNum], lod_min, j_min;	//(delta, J, act, lod)[num of polys]
	int f[interpoly_Ysize][interpoly_Xsize];  //g[num_polys][z-deg+1][max(deg_y)+1][w+1], g1[z-deg+1][max(deg_y)+1][w+2], (g2, f)[z-deg+1][max(deg_y)+1][w+1]	
	int g1[interpoly_Ysize][interpoly_Xsize + 1], g2[interpoly_Ysize][interpoly_Xsize];

	//Initialisation
	for (i = 0; i<init_polyNum; i++)	//num_polys
		for (u = 0; u<interpoly_Ysize; u++)	//max(deg_y)+1
			for (v = 0; v<interpoly_Xsize; v++)	//w+1
				g[i][u][v] = 0;

	for (i = 0; i<(lm + 1); i++)	//rs
		g[i][i][0] = 1;	//j+w*i


#ifdef _PrintLCCLod_
	printf("\nLCC result:");
#endif

	//Interpolation
	for (i = 0; i<n - eta; i++)	//wrt each point
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
						if (g[j][u][v] != 0)
						{
							lod_temp = mono_order[u][v];
							if (lod_temp>lod[j])
								lod[j] = lod_temp;
						}
					}
				}
		}

#ifndef _NoReductionCom_
		//Initialise the eliminator array act[num_poly]
		for (j = 0; j<init_polyNum; j++)	//num_poly
		{
			if (lod[j] <= iterNum)	//C=n when multiplicity = 1
				act[j] = 1;
			else
				act[j] = 0;
		}
#else
		for (j = 0; j<init_polyNum; j++)	//num_poly
			act[j] = 1;
#endif


		//Calculate the hasse derivative mapping of each polynomials
		j_min = -1;
		for (j = 0; j<init_polyNum; j++)	//wrt each polynomial
		{
			J[j] = 0;
			delta[j] = 0;
			if (act[j] == 1)	//for polynomials with leading order less of equal to C
			{
				//Hasse derivative mapping
				for (u = 0; u<interpoly_Ysize; u++)	//deg_z
				{
					delta_temp1 = 0;

					for (v = 0; v<interpoly_Xsize; v++)	//deg_y
					{
						if (g[j][u][v] != 0)	//coputation num concerning coefficients
						{
							delta_temp = power(interpoint[0][i], v);
							delta_temp = mul(delta_temp, g[j][u][v]);
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
						delta_temp1 = mul(delta_temp1, power(interpoint[1][i], u));
						delta[j] = add(delta[j], delta_temp1);
					}

				}

				if (delta[j] != 0)
				{
					J[j] = 1;	//record those polynomial with a nonzero hasse derivative mapping
					lod_min = lod[j];
					j_min = j;
				}
			}
		}

		//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
		for (j = 0; j<init_polyNum; j++)	//num_polys
		{
			if (J[j] == 1 && lod[j]<lod_min)
			{
				lod_min = lod[j];
				j_min = j;
			}
		}
		//printf("\nj_min=%d\n", j_min);

#ifdef _PrintLCCLod_
		printf("\nthe %d iter lod = {%d\t%d\t%d}", i, lod[0], lod[1], lod[2]);
		printf("\nj_min = %d", j_min);
		printf("\nJ =\t");
		for (u = 0; u < init_polyNum; u++)
			if (J[u] != 0)
				printf("%d\t", u);
		printf("\n");
#endif

		if (j_min != -1)
		{
			//f=g[j_min]
			for (u = 0; u<interpoly_Ysize; u++)	//rs
				for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
					f[u][v] = g[j_min][u][v];

			//Modify nonzero polynomials
			for (j = 0; j<init_polyNum; j++)	//num of polys
			{
				if (J[j] == 1)
				{
					if (j != j_min)
					{
						//delta*g_k+delta_k*f
						for (u = 0; u<interpoly_Ysize; u++)	//rs
							for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
							{
								if (g[j][u][v] != 0)
								{
									temp1 = mul(delta[j_min], g[j][u][v]);
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
									g[j][u][v] = add(temp1, temp2);
								}
								else
									g[j][u][v] = 0;
							}
					}
					else if (j == j_min)
					{
							for (u = 0; u<interpoly_Ysize; u++)	//max(deg_y)+1
							{
								for (v = 0; v<(interpoly_Xsize + 1); v++)	//w+2
									g1[u][v] = 0;
								for (v = 0; v<interpoly_Xsize; v++)	//w+1
									g2[u][v] = 0;
							}

						//g1=x*f
						for (u = 0; u<interpoly_Ysize; u++)	//rs
							for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
								if (g[j][u][v] != 0)
								{
									g1[u][v + 1] = g[j][u][v];
								}

						//g2=xi*f
						for (u = 0; u<interpoly_Ysize; u++)	//rs
							for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
									if (g[j][u][v] != 0)
									{
										g2[u][v] = mul(interpoint[0][i], g[j][u][v]);
									}
						//g=g1+g2
						for (u = 0; u<interpoly_Ysize; u++)	//rs
							for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
								if (g1[u][v] != 0 || g2[u][v] != 0)
								{
									g[j][u][v] = add(g1[u][v], g2[u][v]);
								}
								else
									g[j][u][v] = 0;
					}
				}
			}
		}
		//debug: dectect the Q_com_elem[i] who has factor (x+a_i)
		//		DectectIfFactorInPoly(g, init_polyNum, interpoint[0][i]);


		//****debug**********
		//Calculate each polynomial's leading order
		//for (j = 0; j<init_polyNum; j++)	//num_poly
		//{
		//	lod_temp = 0;
		//	lod[j] = 0;
		//	for (u = 0; u<interpoly_Ysize; u++)	//rs
		//	{
		//		for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
		//		{
		//			if (g[j][u][v] != 0)
		//			{
		//				lod_temp = mono_order[u][v];
		//				if (lod_temp>lod[j])
		//					lod[j] = lod_temp;
		//			}
		//		}
		//	}
		//}

#ifdef _Check_Com_Inter_
		//int temp_x, temp_y, temp;

		//for (j = 0; j<init_polyNum; j++)
		//{
		//	temp = 0;
		//	for (u = 0; u<interpoly_Ysize; u++)	//rs
		//		for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
		//			if (g[j][u][v] != 0)
		//			{
		//				temp_x = power(interpoint[0][i], v);
		//				temp_y = power(interpoint[1][i], u);
		//				temp = add(temp, mul(g[j][u][v], mul(temp_y, temp_x)));
		//			}

		//	if ((lod[j] <= iterNum) && temp != 0)
		//		printf("\n[common_interpolation_error]g[%d] with lod[%d]=%d in point[%d](%d,%d) = %d", j, j, lod[j], i, interpoint[0][i], interpoint[1][i], temp);

		//}
		//	printf("\n");
#endif
	}


}

void uncom_elem_interpolation(int interpoint[3], int inGroup[][interpoly_Ysize][interpoly_Xsize], int outGroup1[][interpoly_Ysize][interpoly_Xsize], int outGroup2[][interpoly_Ysize][interpoly_Xsize])
{
	int i, j, u, v, z, J[2][init_polyNum], act[init_polyNum], lod_temp, lod[init_polyNum], lod_min[2], j_min[2];	//(delta, J, act, lod)[num of polys]
	int g1[interpoly_Ysize][interpoly_Xsize + 1], g2[interpoly_Ysize][interpoly_Xsize];
	int inGroup_temp[init_polyNum][interpoly_Ysize][interpoly_Xsize];
	int c[init_polyNum][lm + 1], result_temp[init_polyNum][2];
	int index_min, poly_temp1, poly_temp2;

	//initialization
	for (j = 0; j<init_polyNum; j++)
		for (i = 0; i<2; i++)
		{
			result_temp[j][i] = 0;
		}

	for (j = 0; j<init_polyNum; j++)
		for (u = 0; u<interpoly_Ysize; u++)	//rs
			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
			{
				inGroup_temp[j][u][v] = inGroup[j][u][v];
			}

	//Interpolation

	//Calculate each polynomial's leading order
	for (j = 0; j<init_polyNum; j++)	//num_poly
	{
		lod_temp = 0;
		lod[j] = 0;
		for (u = 0; u<interpoly_Ysize; u++)	//rs
		{
			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
			{
				if (inGroup_temp[j][u][v] != 0)
				{
					lod_temp = mono_order[u][v];
					if (lod_temp>lod[j])
						lod[j] = lod_temp;
				}
			}
		}
	}


	//Initialise the eliminator array act[num_poly]
#ifndef _NoReductionUncom_
	for (j = 0; j<init_polyNum; j++)	//num_poly
	{
		if (lod[j] <= iterNum)	//C=n when multiplicity = 1
			act[j] = 1;
		else
			act[j] = 0;
	}
#else
	for (j = 0; j<init_polyNum; j++)
		act[j] = 1;
#endif

	//Calculate the hasse derivative mapping of each polynomials
	for (j = 0; j<init_polyNum; j++)
		if (act[j] == 1)
			for (u = 0; u<interpoly_Ysize; u++)
			{
				c[j][u] = 0;

				for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
					if (inGroup_temp[j][u][v] != 0)
					{
						poly_temp1 = power(interpoint[0], v);
						poly_temp1 = mul(inGroup_temp[j][u][v], poly_temp1);
						c[j][u] = add(c[j][u], poly_temp1);
						//	c[j][u] = add( c[j][u],mul( inGroup_temp[j][u][v][z],mul( power(interpoint[0],z),power(interpoint[1],v) ) ) );
					}
			}

	for (j = 0; j<init_polyNum; j++)
		if (act[j] == 1)
			for (i = 0; i<2; i++)
			{
				result_temp[j][i] = 0;

				for (u = 0; u<interpoly_Ysize; u++)
				{
					poly_temp1 = power(interpoint[i + 1], u);
					poly_temp1 = mul(c[j][u], poly_temp1);
					result_temp[j][i] = add(result_temp[j][i], poly_temp1);
					//	result_temp[j][i] = add( result_temp[j][i],mul( c[j][u],power(interpoint[i+2],u) ) );
				}
			}

	for (i = 0; i<2; i++)
	{
		j_min[i] = -1;
		for (j = 0; j<init_polyNum; j++)
		{
			J[i][j] = 0;
			if (act[j] == 1 && result_temp[j][i] != 0)
			{
				J[i][j] = 1;	//record those polynomial with a nonzero hasse derivative mapping
				lod_min[i] = lod[j];
				j_min[i] = j;
			}
		}
	}
	//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
	for (i = 0; i<2; i++)
		for (j = 0; j<init_polyNum; j++)	//num_polys
		{
			if (J[i][j] == 1 && lod[j]<lod_min[i])
			{
				lod_min[i] = lod[j];
				j_min[i] = j;
			}
		}
	//printf("\nj_min=%d\n", j_min);

	//initialize outGroup1 and outGroup2
	for (j = 0; j<init_polyNum; j++)
		for (u = 0; u<interpoly_Ysize; u++)	//rs
			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
			{
				outGroup1[j][u][v] = inGroup_temp[j][u][v];
				outGroup2[j][u][v] = inGroup_temp[j][u][v];
			}

	//update the poly of outGroup1 
	if (j_min[0] != -1)
	{
		index_min = j_min[0];	//index_min = j'

		//Modify nonzero polynomials
		for (j = 0; j<init_polyNum; j++)	//num of polys
		{
			if (J[0][j] == 1)
			{
				if (j != j_min[0])
				{
					//delta*g_k+delta_k*f
					for (u = 0; u<interpoly_Ysize; u++)	//rs
						for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
						{
							if (inGroup_temp[j][u][v] != 0)
							{
								poly_temp1 = mul(result_temp[index_min][0], inGroup_temp[j][u][v]);
							}
							else
								poly_temp1 = 0;

							if (inGroup_temp[index_min][u][v] != 0)
							{
								poly_temp2 = mul(result_temp[j][0], inGroup_temp[index_min][u][v]);
							}
							else
								poly_temp2 = 0;

							if (poly_temp1 != 0 || poly_temp2 != 0)
							{
								outGroup1[j][u][v] = add(poly_temp1, poly_temp2);
							}
							else
								outGroup1[j][u][v] = 0;

							//	poly_temp1 = mul( result_temp[index_min][0],inGroup_temp[j][u][v][z] );
							//	poly_temp2 = mul( result_temp[j][0],inGroup_temp[index_min][u][v][z] ) ;
							//	outGroup1[j][u][v][z]=add(poly_temp1,poly_temp2);
							//	outGroup1[j][u][v][z]=add(2,2);
						}
				}
				else if (j == j_min[0])
				{
					for (u = 0; u<interpoly_Ysize; u++)	//rs
					{
						for (v = 0; v<(interpoly_Xsize + 1); v++)	//w+2
							g1[u][v] = 0;
						for (v = 0; v<interpoly_Xsize; v++)	//w+1
							g2[u][v] = 0;
					}

					//g1=x*f
					for (u = 0; u<interpoly_Ysize; u++)	//rs
						for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
								g1[u][v + 1] = inGroup_temp[index_min][u][v];

					//g2=xi*f
					for (u = 0; u<interpoly_Ysize; u++)	//rs
						for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
							if (inGroup_temp[index_min][u][v]!= 0)
							{
								g2[u][v] = mul(interpoint[0], inGroup_temp[index_min][u][v]);
							}
					//g=g1+g2
					for (u = 0; u<interpoly_Ysize; u++)	//rs
						for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
								if (g1[u][v] != 0 || g2[u][v] != 0)
								{
									outGroup1[index_min][u][v] = add(g1[u][v], g2[u][v]);
								}
								else
									outGroup1[index_min][u][v] = 0;
				}
			}
		}
	}

	//update the poly of outGroup2 
	if (j_min[1] != -1)
	{
		index_min = j_min[1];

		//Modify nonzero polynomials
		for (j = 0; j<init_polyNum; j++)	//num of polys
		{
			if (J[1][j] == 1)
			{
				if (j != j_min[1])
				{
					//delta*g_k+delta_k*f
					for (u = 0; u<interpoly_Ysize; u++)	//rs
						for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
						{
							if (inGroup_temp[j][u][v] != 0)
							{
								poly_temp1 = mul(result_temp[index_min][1], inGroup_temp[j][u][v]);
							}
							else
								poly_temp1 = 0;

							if (inGroup_temp[index_min][u][v] != 0)
							{
								poly_temp2 = mul(result_temp[j][1], inGroup_temp[index_min][u][v]);
							}
							else
								poly_temp2 = 0;

							if (poly_temp1 != 0 || poly_temp2 != 0)
							{
								outGroup2[j][u][v] = add(poly_temp1, poly_temp2);
							}
							else
								outGroup2[j][u][v] = 0;

							//	poly_temp1 = mul( result_temp[index_min][1],inGroup_temp[j][u][v][z] );
							//	poly_temp2 = mul( result_temp[j][1],inGroup_temp[index_min][u][v][z] );
							//	outGroup2[j][u][v][z]=add( poly_temp1,poly_temp2 );	
						}
				}
				else if (j == j_min[1])
				{
					for (u = 0; u<interpoly_Ysize; u++)	//rs
					{
						for (v = 0; v<(interpoly_Xsize + 1); v++)	//w+2
							g1[u][v] = 0;
						for (v = 0; v<interpoly_Xsize; v++)	//w+1
							g2[u][v] = 0;
					}

					//g1=x*f
					for (u = 0; u<interpoly_Ysize; u++)	//rs
						for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
							g1[u][v + 1] = inGroup_temp[index_min][u][v];

					//g2=xi*f
					for (u = 0; u<interpoly_Ysize; u++)	//rs
						for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
							if (inGroup_temp[index_min][u][v] != 0)
							{
								g2[u][v] = mul(interpoint[0], inGroup_temp[index_min][u][v]);
							}

					//g=g1+g2
					for (u = 0; u<interpoly_Ysize; u++)	//rs
						for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
							if (g1[u][v] != 0 || g2[u][v] != 0)
							{
								outGroup2[index_min][u][v] = add(g1[u][v] , g2[u][v]);
							}
							else
								outGroup2[index_min][u][v] = 0;
				}
			}
		}
	}

#ifdef _Check_Uncom_Inter_
	//*******debug*********
	//outGroup1
	for (j = 0; j<init_polyNum; j++)
	{
		lod_temp = 0;
		lod[j] = 0;
		for (u = 0; u<interpoly_Ysize; u++)	//rs
		{
			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
			{
				if (outGroup1[j][u][v] != 0)
				{
					lod_temp = mono_order[u][v];
					if (lod_temp>lod[j])
						lod[j] = lod_temp;
				}
			}
		}
	}

	int temp_x, temp_y, temp;

	for (j = 0; j<init_polyNum; j++)
	{
		temp = 0;
		for (u = 0; u<interpoly_Ysize; u++)	//rs
			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
				if (outGroup1[j][u][v] != 0)
				{
					temp_x = power(interpoint[0], v);
					temp_y = power(interpoint[1], u);
					temp = add(temp, mul(outGroup1[j][u][v], mul(temp_y, temp_x)));
				}
		if (temp!=0)
			printf("\noutGroup1[%d] with lod[%d]=%d in point[up](%d,%d) = %d\n", j, j, lod[j], interpoint[0], interpoint[1], temp);
	}
//	printf("\n");

	//outGroup2
	for (j = 0; j<init_polyNum; j++)
	{
		lod_temp = 0;
		lod[j] = 0;
		for (u = 0; u<interpoly_Ysize; u++)	//rs
		{
			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
			{
				if (outGroup2[j][u][v] != 0)
				{
					lod_temp = mono_order[u][v];
					if (lod_temp>lod[j])
						lod[j] = lod_temp;
				}
			}
		}
	}

	for (j = 0; j<init_polyNum; j++)
	{
		temp = 0;
		for (u = 0; u<interpoly_Ysize; u++)	//rs
			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
				if (outGroup2[j][u][v] != 0)
				{
					temp_x = power(interpoint[0], v);
					temp_y = power(interpoint[2], u);
					temp = add(temp, mul(outGroup2[j][u][v], mul(temp_y, temp_x)));
				}
		if (temp != 0)
			printf("\noutGroup2[%d] with lod[%d]=%d in point[down](%d,%d) = %d\n", j, j, lod[j], interpoint[0], interpoint[2], temp);
	}
	//printf("\n");
#endif
	//******************

}

//void factorisation(int Q_input[test_vec_num][interpoly_Ysize][interpoly_Xsize], int output_list[test_vec_num][lm + 1][k], int list_num[test_vec_num])
//{
//	int i, j, u, v, z;
//
//	for (i = 0; i<test_vec_num; i++)
//	{
//		//Initialisation
//		for (j = 0; j<k; j++)	//number of fac steps=k
//			for (u = 0; u<facpoly_Ysize; u++)	//rs
//				for (v = 0; v<facpoly_Xsize; v++)	//y_size
//						Q[j][u][v] = 0;
//
//		//		for(u=0;u<k;u++)	//number of fac steps=k
//		//			for(v=0;v<lm+1;v++)	//5>rs=expected number of roots
//		//				rootlist[u][v]=-1;	
//
//		for (u = 0; u<lm + 1; u++)	//5>(rs-1)=expected number of output lists
//			for (v = 0; v<k; v++)	//k
//				output[u][v] = -1;
//
//		//Initialisation of factorisation
//		uu = 0;	//recursive deduction index
//		l = 0;	//candidate output index
//		//q_0(z)=itp(z)
//		for (u = 0; u<interpoly_Ysize; u++)	//rs
//			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
//					Q[uu][u][v] = Q_input[i][u][v];
//
//		//recursive coefficient search
//		rcs(uu);
//
//		//store the necessary data
//		for (u = 0; u<lm + 1; u++)	//5>(rs-1)=expected number of output lists
//			for (v = 0; v<k; v++)	//k
//				output_list[i][u][v] = output[u][v];
//
//		list_num[i] = l;
//
//	}
//
//}

void choose(int output_codeword[n], int output_bicodeword[n*p], int output_list[test_vec_num][lm + 1][k], int list_num[test_vec_num], int flag_decoding_alg)
{
	int i, j, u, v, z, value, flag, min_index_1, min_index_2;
	unsigned int mask = 1;
	//	float proba[lm*test_vec_num], proba_temp, temp;
	double proba[lm*test_vec_num], proba_temp, temp;
	int codeword_temp[n];

	//normal mode
	//Initialise hamming distance counter
	for (u = 0; u<lm*test_vec_num; u++)
		proba[u] = -1.0;

	proba_temp = 0.0;
	flag = 0;
	min_index_1 = -1;
	min_index_2 = -1;

	for (i = 0; i<test_vec_num; i++)
		if (list_num[i] != 0)
		{
			for (j = 0; j<list_num[i]; j++)
			{
				//reencoding
				encoder(output_list[i][j], codeword_temp);
				//calculate the posteriori probablity
				temp = 1.0;
				for (u = 0; u<n; u++)
				{
					for (v = 0; v<q; v++)
						if (codeword_temp[u] == root[v])
							temp = temp*RM[v][u];
				}

				if (proba_temp < temp)	//< or <=
				{
					proba_temp = temp;
					min_index_1 = i;
					min_index_2 = j;

					for (v = 0; v<n; v++)
						output_codeword[v] = codeword_temp[v];

					flag = 1;	//exist at less one valid solution
				}
			}
		}

	//output the decoding result 
	if (flag == 0)	// not exist a valid solution
	{
		for (u = 0; u<n; u++)
			output_codeword[u] = large_vec[0][u];

		//nonbinary --> binary
		for (u = 0; u<n; u++)
		{
			value = output_codeword[u];
			mask = 1;
			for (v = 0; v<p; v++)
			{
				if ((value & mask)>0)
					output_bicodeword[p*u + v] = 1;
				else
					output_bicodeword[p*u + v] = 0;
				mask = mask << 1;
			}
		}
	}
	else if (flag == 1)	//exist a valid solution
	{

		//encode message[min_index]
		//		encoder(output_list[min_index_1][min_index_2],output_codeword);
		//nonbinary --> binary
		for (u = 0; u<n; u++)
		{
			value = output_codeword[u];
			mask = 1;
			for (v = 0; v<p; v++)
			{
				if ((value & mask)>0)
					output_bicodeword[p*u + v] = 1;
				else
					output_bicodeword[p*u + v] = 0;
				mask = mask << 1;
			}
		}
	}
	//*********************


	//Check Herm(64, 39)'s working perperty
	int epcount2 = 0;
	for (u = 0; u<n; u++)
		if (codeword[u] != output_codeword[u])
			epcount2++;

	//******debug*******
	//caculate the degree_test[i]
	for (i = 0; i < test_vec_num; ++i)
	{
		int mono_temp = -1;
		degree_test[i] = 0;
		for (u = 0; u<interpoly_Ysize; u++)	//rs
			for (v = 0; v<interpoly_Xsize; v++)	//max(deg_y)+1
				if (Q_interpoly[i][u][v] != 0 && mono_temp < mono_order[u][v])
				{
					mono_temp = mono_order[u][v];
					degree_test[i] = v + u*(k-1);
				}
	}

	//calculate the testCount1
	testCount1[0] = 0;
	for (i = 0; i < n - eta; ++i)
	{
		int index_temp = (int)test_set_ordered[1][i];
		if (codeword[index_temp] != large_vec[0][index_temp])
			testCount1[0]++;
	}

	for (i = 1; i < test_vec_num; ++i)
		testCount1[i] = testCount1[0];

	for (i = 0; i < test_vec_num; i++)
	{
		int index_temp;
		mask = 1;
		for (int j = n - 1; j >= n - eta - 1; j--)
		{
			//bit caculate
			if ((i & mask) > 0)
				index_temp = 1;
			else
				index_temp = 0;

			if (codeword[(int)test_set_ordered[1][j]] != large_vec[index_temp][(int)test_set_ordered[1][j]])
				testCount1[i]++;

			mask = mask << 1;
		}
	}

	//calculate the testCount2
	//for (i = 0; i<test_vec_num; i++)
	//{
	//	testCount2[i] = 0;
	//	for (j = 0; j<k; j++)
	//		if (message[j] != output_list[i][0][j])
	//			testCount2[i]++;
	//}

	if (flag_decoding_alg == 0)
	{
		for (i = 0; i < test_vec_num; i++)
		{
			if ((n - testCount1[i])>degree_test[i])
			{
				flag = 0;
				for (j = 0; j < lm + 1; j++)
				{
					testCount2 = 0;
					for (u = 0; u < k; u++)
						if (output_list[i][j][u] != -1 && message[u] != output_list[i][j][u])
							testCount2++;

					if (testCount2 == 0)
					{
						flag = 1;
						break;
					}
				}

				if (flag==0)
				{
					printf("\n\nseq_num_Now=%d, No.%d test vector factorization has failed!!, testCount1=%d, testCount2=%d, degreeOfPolynomail=%d, list_num=%d\n", seq_num_Now, i, testCount1[i], testCount2, degree_test[i], list_num[i]);
					//				printf("\nNo.%d test vector has error!!, testCount1=%d, testCount2=%d, testCount1_uncom=%d, list_num=%d\n\n", i, testCount1[i], testCount2[i], (testCount1[i]-testCount1_com), list_num[i]); 
					//				printf("\nseq_num_Now=%d, This sequence is decoding failed£¬ and error num is %d!\n\n", seq_num_Now, epcount1);
					//********debug*************
					int temp_x, temp_y, temp_z;
					int temp, h, z;

					//**************
					printf("\n");
					for (u = 0; u < interpoly_Ysize; u++)	//rs
						for (v = 0; v < interpoly_Xsize; v++)	//max(deg_y)+1
							for (z = 0; z < interpoly_Xsize; z++)	//w+1
								if (Q_interpoly[i][u][v] != 0)
								{
									printf("Q_interpoly[%d][%d][%d]=%d\n", i, u, v, Q_interpoly[i][u][v]);
								}
					printf("\n");

					//************
					printf("\nTx message is:\n");
					for (u = 0; u < k; u++)	//k
						printf("%d\t", message[u]);

					for (u = 0; u < list_num[i]; u++)
					{
						printf("\nlist_%d:\n", u);
						for (v = 0; v < k; v++)
							printf("%d\t", output_list[i][u][v]);
					}
					printf("\n\n");
					//***************
					printf("\nTx message is:\n");
					for (u = 0; u < k; u++)	//k
						printf("%d ", message[u]);
					printf("\nTx code word is:\n");
					for (u = 0; u < n; u++)	//n
						printf("%d ", codeword[u]);

					printf("\n\n");
					for (v = 0; v < n*p / 4; v++)
						printf("\t(%f,%f,%f,%f)", rx_symbol[v * 2 + 0][0], rx_symbol[v * 2 + 0][1], rx_symbol[v * 2 + 1][0], rx_symbol[v * 2 + 1][0]);

					printf("\n\n");
					for (u = 0; u < q; u++)
					{
						printf("\n");
						printf("\t%d", root[u]);
						for (v = 0; v < n; v++)
							printf("\t%f", RM[u][v]);
					}
					printf("\n\n");

					printf("\n\n");
					printf("large_vec[]\t");
					for (j = 0; j < n; j++)
						printf("\t%d", large_vec[0][j]);

					printf("\n\n");
					printf("large_vec2[]\t");
					for (j = 0; j < n; j++)
						printf("\t%d", large_vec[1][j]);


					printf("\n\n");
					printf("x_order[]\t");
					for (i = 0; i < n; i++)
						printf("\t%d", x_ordered[i]);


					printf("\n\n");
					printf("large_vec_order[]");
					for (j = 0; j < n; j++)
						printf("\t%d", large_vec[0][(int)test_set_ordered[1][j]]);

					printf("\n\n");
					printf("large_vec2_order[]");
					for (j = 0; j < n; j++)
						printf("\t%d", large_vec[1][(int)test_set_ordered[1][j]]);

					printf("\n\n");
					printf("codeword_order[]");
					for (j = 0; j < n; j++)
						printf("\t%d", codeword[(int)test_set_ordered[1][j]]);

					printf("\n%d errors in received!", epcount1);
					printf("\noutputlist:");

					for (i = 0; i < test_vec_num; i++)
					{
						printf("\n\nlist_num[%d]=%d", i, list_num[i]);
						for (u = 0; u < list_num[i]; u++)
						{
							printf("\n[%d][%d]", i, u);
							for (v = 0; v < k; v++)
								printf("\t%d", output_list[i][u][v]);
							printf("\tproba[%d][%d]=%f", i, u, proba[i*(lm + 1) + u]);

						}
					}
					printf("\nDecoded code word is:\n");
					for (u = 0; u < n; u++)	//n
						printf("%d ", output_codeword[u]);
					printf("\n%d errors after decoding!\n\n", epcount2);
					//***********
				}
			}
		}
		//*******************
	}


}

int result_compare(int output_list[test_vec_num][lm + 1][k], int output_codeword[n], int output_list_BF[test_vec_num][lm + 1][k], int output_codeword_BF[n])
{
	int result = 0;
	for (int i = 0; i < n; ++i)
		if (output_codeword[i] != output_codeword_BF[i])
		{
			result = 1;
			break;
		}

	for (int i = 0; i < test_vec_num; ++i)
		for (int j = 0; j < lm + 1; ++j)
			for (int u = 0; u < k; ++u)
				if (output_list[i][j][u] != output_list_BF[i][j][u])
				{
					result += 2;
					return result;
				}

	return result;
}

void test_vec_construction(void)
{	//construct the test_vec[test_vec_num]
	int i, j, temp_index;
	float temp;
	float test_set[2][n]; // [0] stores value, [1] stores index
	float large_vec_proba[choose_num][n]; // largest index is 0, second largest index is 1

	//cal test_set
	for (i = 0; i<n; i++)
	{
		//finde the largest vec
		large_vec_proba[0][i] = 0;
		for (j = 0; j<q; j++)
			if (large_vec_proba[0][i]<RM[j][i])
			{
				large_vec_proba[0][i] = RM[j][i];
				large_vec[0][i] = root[j];
			}
		//find the second largest vec
		large_vec_proba[1][i] = 0;
		for (j = 0; j<q; j++)
			if ((large_vec_proba[1][i]<RM[j][i]) && (root[j] != large_vec[0][i]))
			{
				large_vec_proba[1][i] = RM[j][i];
				large_vec[1][i] = root[j];
			}

		test_set[0][i] = large_vec_proba[1][i] / large_vec_proba[0][i];
		test_set[1][i] = i;
	}


	//order test_vec
	//key part!!
	for (i = 0; i<n; i++)
	{
		test_set_ordered[0][i] = test_set[0][i];
		test_set_ordered[1][i] = test_set[1][i];
	}

	for (i = 0; i<(n - 1); i++)
		for (j = 0; j<(n - (i + 1)); j++)
			if (test_set_ordered[0][j]>test_set_ordered[0][j + 1])
			{
				temp = test_set_ordered[0][j];
				test_set_ordered[0][j] = test_set_ordered[0][j + 1];
				test_set_ordered[0][j + 1] = temp;

				temp_index = test_set_ordered[1][j];
				test_set_ordered[1][j] = test_set_ordered[1][j + 1];
				test_set_ordered[1][j + 1] = temp_index;
			}

	//contruction test_vec
	//for (i = 0; i<n - eta; i++)
	//	test_vec_com[i] = large_vec[0][(int)test_set_ordered[1][i]];

	//construct x_ordered
	for (i = 0; i<n; i++)
		x_ordered[i] = mularray[(int)test_set_ordered[1][i]];

#ifdef checkPrint	
	//************debug***************
	printf("\n\n");
	printf("xi[]\t");
	for (j = 0; j<n; j++)
		printf("\t%d", mularray[j]);

	printf("\n\n");
	printf("large_vec[]\t");
	for (j = 0; j<n; j++)
		printf("\t%d", large_vec[0][j]);

	printf("\n\n");
	printf("large_vec2[]\t");
	for (j = 0; j<n; j++)
		printf("\t%d", large_vec[1][j]);

	printf("\n\n");
	printf("x_order[]\t");
	for (j = 0; j<n; j++)
		printf("\t%d", x_ordered[j]);


	printf("\n\n");
	printf("large_vec_order[]");
	for (j = 0; j<n; j++)
		printf("\t%d", large_vec[0][(int)test_set_ordered[1][j]]);

	printf("\n\n");
	printf("large_vec2_order[]");
	for (j = 0; j<n; j++)
		printf("\t%d", large_vec[1][(int)test_set_ordered[1][j]]);

	printf("\n\n");
	printf("codeword[]\t");
	for (j = 0; j<n; j++)
		printf("\t%d", codeword[(int)test_set_ordered[1][j]]);
	//********************************
#endif
}

void mono_table(void)
{
	int i, j, v, cal_flag;
	int weight_degree[weight_Ysize][weight_Xsize], mono_order_temp[weight_Ysize][weight_Xsize];
	
	//Initialization
	for (i = 0; i < weight_Ysize; i++)
		for (j = 0; j < weight_Xsize; j++)
			mono_order_temp[i][j] = -1;

	for (i = 0; i < weight_Ysize; i++)
		for (j = 0; j < weight_Xsize; j++)
			weight_degree[i][j] = j + i*(k - 1);

	//calculate RS's monoTable
	cal_flag = 0;
	for (v = 0; v <= target_weight_degree; v++)
		for (i = 0; i < weight_Ysize; i++)
			for (j = 0; j < weight_Xsize; j++)
			{
				if (weight_degree[i][j] == v)
				{
					mono_order_temp[i][j] = cal_flag;
					cal_flag++;
					break;
				}
			}

	//set mono_table
	for (i = 0; i < monoTable_Ysize; i++)
		for (j = 0; j < monoTable_Xsize; j++)
			if (mono_order_temp[i][j] != -1)
				mono_order[i][j] = mono_order_temp[i][j];

}

void generator(void)
{
	int i, j;
	for (i = 0; i<k; i++)
		for (j = 0; j<n; j++)	//n
			gmatrix[i][j] = 0;

	//generator matrix
	for (i = 0; i<k; i++)	//k
		for (j = 0; j<n; j++)	//n
			gmatrix[i][j] = power(mularray[j], i);

}

void encoder(int message_temp[], int codeword_temp[])
{
	int i, j, temp;

	//Encoding
	for (i = 0; i<n; i++)	//n
	{
		codeword_temp[i] = 0;
		for (j = 0; j<k; j++)	//k
		{
			temp = mul(message_temp[j], gmatrix[j][i]);
			codeword_temp[i] = add(codeword_temp[i], temp);
		}
	}

}

void modulation(void)
{
	int i;

	//BPSK
	for (i = 0; i<n*p; i++)
	{
		tx_symbol[i][0] = -2 * bi_codeword[i] + 1;	//0-->(1,0)
		tx_symbol[i][1] = 0;  					//1-->(-1,0)
	}
}

void channel(void)
{
	int i, j;
	float u, r, g;

	//Add AWGN 
	for (i = 0; i<n*p; i++)
	{
		for (j = 0; j<2; j++)	//Inphase + Quadrature
		{
			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			r = sgm*sqrt(2.0*log(1.0 / (1.0 - u)));

			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			g = (float)r*cos(2 * pi*u);

			rx_symbol[i][j] = tx_symbol[i][j] + g;
		}
	}
}

void demodulation(void)  //reconed constellation diagram
{	//soft decision
	int i, j, u, v;
	float sum, proba[pointNum];
	float Pr[p][2];

	//initialize
	for (u = 0; u<q; u++)
		for (v = 0; v<n; v++)
			RM[u][v] = 0;

	//cal PDF
	for (i = 0; i<n; i++)
	{
		for (j = 0; j<p; j++)
		{
			sum = 0;
			for (u = 0; u<pointNum; u++)
			{
				proba[u] = PDF(i*p + j, u);
				sum += proba[u];
			}

			Pr[j][0] = proba[0] / sum;	//proba-->00	0
			Pr[j][1] = proba[1] / sum;	//proba-->01	1

		}

		j = 0;

		for (int x = 0; x<pointNum; x++)
			for (int y = 0; y<pointNum; y++)
				for (int h = 0; h<pointNum; h++)
				{
					RM[j][i] = (Pr[0][h] * Pr[1][y] * Pr[2][x]);
					j++;
				}

	}

}

float PDF(int i, int u)
{
	float temp1, temp2, temp3;
	float constell_point[2][2] = { { 1, 0 }, { -1, 0 } }; // S1{1,0}-->0,S2{-1,0}-->1

	temp3 = 1.0 / (pi*N0);
	temp1 = pow((rx_symbol[i][0] - constell_point[u][0]), 2) + pow((rx_symbol[i][1] - constell_point[u][1]), 2);
	temp2 = temp3*exp(-temp1 / N0);  // 0.5 is the normalization coefficient

	return temp2;
}

void factorisation(int Q_input[poly_Y_size][poly_X_size])
{
	/**********************/
	//factorisation step
	/**********************/
	int i, j, h, u;
	int min_dist_pos;//the position of minimun Euclidean distance
	//initialization Q(x,y)
	for (h = 0; h<fac_size; h++)
		for (j = 0; j<poly_Y_size; j++)
			for (i = 0; i<poly_X_size; i++)
				Q[h][j][i] = 0;

	//initialization pai[],deg[]
	for (i = 0; i<fac_size; i++)
	{
		pai[i] = -1;
		deg[i] = -1;
	}

	//initialization coeff[]
	for (i = 0; i<coeff_size; i++)
		coeff[i] = -1;

	//initialization the factorisation polynomials
	for (i = 0; i<(lm + 1); i++)
		for (j = 0; j<k; j++)
			fx[i][j] = -1;

	//Q_0(x,y)=<<Q(x,y)>>=<<inter_poly_minlod(x,y)>>, <<f(x)>> means f(x)/x^m, where m is the biggest number such that x^m|f(x)
	h = -1;
	for (i = 0; i<poly_X_size; i++)
	{
		for (j = 0; j<poly_Y_size; j++)
		{
			if (Q_input[j][i] != 0)
			{
				h = i;
				break;
			}
		}
		if (h == i)
			break;
	}
	for (j = 0; j<poly_Y_size; j++)
		for (i = 0; i<(poly_X_size - h); i++)
			Q[0][j][i] = Q_input[j][i + h];

	//begin factorization
	f_num = 0;
	pai[0] = -1;
	deg[0] = -1;
	times = 1;
	u = 0;
	//deep first search in vertex u
	DFS(u);
}

void DFS(int u)
{
	/**********************/
	//deep first search in vertex u
	/**********************/
	int temp, list_size;
	int i, v;
	int rootlist[lm];//the root of polynomial
	//initialization
	for (i = 0; i<lm; i++)
		rootlist[i] = -1;

	temp = 0;
	for (i = 0; i<(poly_X_size + lm + 1); i++)
	{
		if (Q[u][0][i] != 0)
			temp++;
	}
	//if Q_u(x,0)=0, output the factorisation polynomial
	//causion: deg[u]%(k-1) is to ensure to DFS the last degree
	if (temp == 0 && deg[u] % (k - 1) == 0 && deg[u] != 0)
	{
		output_fac_poly(u);
	}
	//explore edges from vertex u
	else if (deg[u]<(k - 1))//deg f(x)<D=k-1
	{
		list_size = 0;
		//find the root of Q_u(0,y)
		for (i = 0; i<q; i++)
		{
			temp = find_root(u, i);
			if (temp == 0)
			{
				rootlist[list_size] = root[i];
				list_size++;
			}
		}
		//for each root, deep search the factorisation coefficient
		for (i = 0; i<list_size; i++)
		{
			v = times;
			times = times + 1;
			pai[v] = u;//the parent of v is u
			deg[v] = deg[u] + 1;
			coeff[v] = rootlist[i];;
			update_fac_poly(v, u, rootlist[i]);//update the factorisation polynomial, Q_v(x,y)=<<Q_u(x,xy+alpha)>>
			DFS(v);
		}
	}
}
void update_fac_poly(int v, int u, int alpha)
{
	/**********************/
	//update the factorisation polynomial, Q_v(x,y)=<<Q_u(x,xy+alpha)>>
	/**********************/
	int temp_fac_poly[poly_Y_size][poly_X_size + lm + 1];//the x size of polynomial should plus (lm+1) since the second part of Q_u(x,xy+alpha)
	int i, j, h;
	int beta, sign, temp;
	int poly_deg = 0;

	//initialize temp_fac_poly
	for (j = 0; j<poly_Y_size; j++)
		for (i = 0; i<(poly_X_size + lm + 1); i++)
			temp_fac_poly[j][i] = 0;
	//update Q_u(x,xy+alpha), its coeffecients Q_u[j][i] can be expressed as D_j*g_(i-j)(alpha)=(beta_up,j_down)*Q_(beta,(i-j))*alpha^(beta-j)
	for (j = 0; j<poly_Y_size; j++)
	{
		for (i = j; i<(poly_X_size + lm + 1); i++)
		{
			temp = 0;
			for (beta = j; beta<poly_Y_size; beta++)
			{
				sign = combination(beta, j);
				if ((sign % 2) != 0)
				{
					temp = add(temp, mul(Q[u][beta][i - j], power(alpha, (beta - j))));
					//if(Q[u][beta][i-j]!=0)
					//	fac_mul_com=fac_mul_com+2;
					//if(mul(Q[u][beta][i-j],power(alpha,(beta-j)))!=0)
					//	fac_add_com++;
				}
			}
			temp_fac_poly[j][i] = temp;
		}
	}
	//Q_v(x,y)=<<Q_u(x,y)>>
	h = -1;
	for (i = 0; i<(poly_X_size + lm + 1); i++)
	{
		for (j = 0; j<poly_Y_size; j++)
		{
			if (temp_fac_poly[j][i] != 0)
			{
				h = i;
				break;
			}
		}
		if (h == i)
			break;
	}
	//refresh Q_v(x,y)
	for (j = 0; j<poly_Y_size; j++)
		for (i = 0; i<(poly_X_size + lm + 1); i++)
			Q[v][j][i] = 0;
	//update Q_v(x,y)
	for (j = 0; j<poly_Y_size; j++)
		for (i = 0; i<(poly_X_size + lm + 1 - h); i++)
			Q[v][j][i] = temp_fac_poly[j][i + h];
}
void output_fac_poly(int u)
{
	/**********************/
	//output the factorisation polynomial
	/**********************/
	int i, temp, deg_u;

	deg_u = deg[u];
	//f_u(x)=coeff[u]*x^deg[u]+coeff[pai[u]]*x^deg[pai[u]]+...
	while (coeff[u] != -1 && deg[u] != -1 && pai[u] != -1)
	{
		//outputList[fac_test_vec_mark][f_num][deg_u] = coeff[u];
		fx[f_num][deg_u] = coeff[u];
		u = pai[u];
		deg_u = deg[u];
	}
	//validate the f(x) 
	temp = 0;
	for (i = 0; i<k; i++)
		if (fx[f_num][i] != -1)
			temp++;
	//if f(x) has k bits, f(x) is a valid message polynomial
	if (temp == k)
		f_num++;
}
int find_root(int u, int i)
{
	/**********************/
	//find the root of Q_u(0,y)
	/**********************/
	int j, temp;
	temp = 0;
	for (j = 0; j<poly_Y_size; j++)
	{
		temp = add(temp, mul(Q[u][j][0], power(root[i], j)));
		//if(Q[u][j][0]!=0)
		//	fac_mul_com=fac_mul_com+2;
		//if(mul(Q[u][j][0],power(root[i],j))!=0)
		//	fac_add_com++;
	}
	return temp;
}

int combination(int a, int b)
{
	/**********************/
	//judge whether the combination number of two numbers (a_up,b_dowm) is even or odd
	/**********************/
	int value, sign;
	if (a >= b)
	{
		value = a&b;

		if (value == b)
			sign = 1;
		else
			sign = 2;
	}
	else
	{
		sign = -1;
	}
	return sign;
}