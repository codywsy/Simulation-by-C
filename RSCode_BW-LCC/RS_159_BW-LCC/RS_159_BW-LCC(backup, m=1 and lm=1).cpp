#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

//mode_EuclideanDisdance	mode = 1
//mode_PostProba			mode = 2
//mode_MLrule				mode = 3
//mode_genius				mode = 4
#define mode 2

#define OpenFile fp=fopen("RS159_BW-LCC_η=4.txt","a")

#define finiteField 16
#define FrameError 309
//#define cheatingDec
#define cheatingEncoding
//#define ChaseGameDebug
//#define checkInterpolation
//#define checkChoose	
#define eta_bit 4
#define tv_num_bit 16 


//************variable define***************
#define k 9  //length of message
#define able_correct_num ((n-k)/2)
#define m 1 // the multiplicity of LCC algorithm is 1
#define lm 1 // LCC's lm is equal to 1, solid!!
#define V 1 //projected basic function BPSK=1, QPSK=2
#define pointNum 2 // the number of constell point == 2^V
#define choose_num 2 //the num of codeword choice, choose 2 in BASIC LCC algorithm
#define interval 1
#define fac_size (lm*k+4)//the probable numbers of factorisation polynomials, bigger than lm*k
#define coeff_size (lm*k+4)//the probable numbers of coefficients of factorisation polynomials, bigger than lm*k

//****modify variable**********
#define WeightDegree_Ysize 2
#define WeightDegree_Xsize 200
#define mono_ordinSize WeightDegree_Xsize	// larger than (monoTable_Ysize * monoTable_Xsize)
#define monoTable_Ysize 2	//designed to cover fully, larger than monoOrder_Ysize
#define monoTable_Xsize 200	// designed to cover fully, larger than monoOrder_Xsize and (n+k+1)
#define interpoly_Ysize (lm+1)	//value+1, fix value is lm 
#define interpoly_Xsize (n+1)	//value+1, max value is iterNum in this algorithm
//*********************************************

//************************
#if finiteField == 16
//(15,X)
#define n 15	//length of codeword
#define q 16  //GF(q) q=n+1
#define p 4 //GF(2^p) = GF(q)
int mularray[] = { 1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9 }; 	//this array is used to the infinite field element mul
int root[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };			//n+1		be used to the factorization
int logNum[] = { -1, 0, 1, 4, 2, 8, 5, 10, 3, 14, 9, 7, 6, 13, 11, 12 };	//used to locate the degree of finitefield element through the value of finitefield element

#elif finiteField == 32
//(31,X)
#define n 31	//length of codeword
#define q 32  //GF(q) q=n+1
#define p 5	//GF(2^p) = GF(q)
int mularray[] = { 1, 2, 4, 8, 16, 5, 10, 20, 13, 26, 17, 7, 14, 28, 29, 31, 27, 19, 3, 6, 12, 24, 21, 15, 30, 25, 23, 11, 22, 9, 18 };
//this array is used to the infinite field element mul
int root[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 };
//n+1		be used to the factorization
int logNum[] = { -1, 0, 1, 18, 2, 5, 19, 11, 3, 29, 6, 27, 20, 8, 12, 23,
4, 10, 30, 17, 7, 7, 22, 28, 26, 21, 25, 9, 16, 13, 14,
15 };	//used to locate the degree of finitefield element through the value of finitefield element

#elif finiteField == 64
//(63,X)
#define n 63	//length of codeword
#define q 64  //GF(q) q=n+1
#define p 6	//GF(2^p) = GF(q)
int mularray[] = { 1, 2, 4, 8, 16, 32, 3, 6, 12, 24, 48, 35, 5, 10, 20, 40,
19, 38, 15, 30, 60, 59, 53, 41, 17, 34, 7, 14, 28, 56, 51, 37,
9, 18, 36, 11, 22, 44, 27, 54, 47, 29, 58, 55, 45, 25, 50, 39,
13, 26, 52, 43, 21, 42, 23, 46, 31, 62, 63, 61, 57, 49, 33 };  //this array is used to the infinite field element mul

int root[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 };	//n+1		be used to the factorization

int logNum[] = { -1, 0, 1, 6, 2, 12, 7, 26, 3, 32, 13, 35, 8, 48, 27, 18,
4, 24, 33, 16, 14, 52, 36, 54, 9, 45, 49, 38, 28, 41, 19,
56, 5, 62, 25, 11, 34, 31, 17, 47, 15, 23, 53, 51, 37, 44,
55, 40, 10, 61, 46, 30, 50, 22, 39, 43, 29, 60, 42, 21, 20,
59, 57, 58 };	//used to locate the degree of finitefield element through the value of finitefield element

#elif finiteField == 256
//(63,X)
#define n 255	//length of codeword
#define q 256  //GF(q) q=n+1
#define p 8	//GF(2^p) = GF(q)
int mularray[] = { 1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 152, 45, 90, 180, 117, 234,
201, 143, 3, 6, 12, 24, 48, 96, 192, 157, 39, 78, 156, 37, 74, 148, 53, 106, 212, 181, 119, 238, 193, 159, 35, 70,
140, 5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111, 222, 161, 95, 190, 97, 194, 153, 47, 94, 188, 101, 202,
137, 15, 30, 60, 120, 240, 253, 231, 211, 187, 107, 214, 177, 127, 254, 225, 223, 163, 91, 182, 113, 226, 217, 175,
67, 134, 17, 34, 68, 136, 13, 26, 52, 104, 208, 189, 103, 206, 129, 31, 62, 124, 248, 237, 199, 147, 59, 118, 236, 197,
151, 51, 102, 204, 133, 23, 46, 92, 184, 109, 218, 169, 79, 158, 33, 66, 132, 21, 42, 84, 168, 77, 154, 41, 82, 164, 85,
170, 73, 146, 57, 114, 228, 213, 183, 115, 230, 209, 191, 99, 198, 145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 255,
227, 219, 171, 75, 150, 49, 98, 196, 149, 55, 110, 220, 165, 87, 174, 65, 130, 25, 50, 100, 200, 141, 7, 14, 28, 56, 112,
224, 221, 167, 83, 166, 81, 162, 89, 178, 121, 242, 249, 239, 195, 155, 43, 86, 172, 69, 138, 9, 18, 36, 72, 144, 61, 122,
244, 245, 247, 243, 251, 235, 203, 139, 11, 22, 44, 88, 176, 125, 250, 233, 207, 131, 27, 54, 108, 216, 173, 71, 142 };  //this array is used to the infinite field element mul

int root[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88,
89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114,
115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138,
139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162,
163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186,
187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210,
211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234,
235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255 };	//n+1		be used to the factorization

int logNum[] = { -1, 0, 1, 25, 2, 50, 26, 198, 3, 223, 51, 238, 27, 104, 199, 75, 4, 100, 224, 14, 52, 141, 239, 129, 28,
193, 105, 248, 200, 8, 76, 113, 5, 138, 101, 47, 225, 36, 15, 33, 53, 147, 142, 218, 240, 18, 130, 69, 29, 181, 194, 125,
106, 39, 249, 185, 201, 154, 9, 120, 77, 228, 114, 166, 6, 191, 139, 98, 102, 221, 48, 253, 226, 152, 37, 179, 16, 145, 34,
136, 54, 208, 148, 206, 143, 150, 219, 189, 241, 210, 19, 92, 131, 56, 70, 64, 30, 66, 182, 163, 195, 72, 126, 110, 107, 58,
40, 84, 250, 133, 186, 61, 202, 94, 155, 159, 10, 21, 121, 43, 78, 212, 229, 172, 115, 243, 167, 87, 7, 112, 192, 247, 140,
128, 99, 13, 103, 74, 222, 237, 49, 197, 254, 24, 227, 165, 153, 119, 38, 184, 180, 124, 17, 68, 146, 217, 35, 32, 137, 46,
55, 63, 209, 91, 149, 188, 207, 205, 144, 135, 151, 178, 220, 252, 190, 97, 242, 86, 211, 171, 20, 42, 93, 158, 132, 60, 57,
83, 71, 109, 65, 162, 31, 45, 67, 216, 183, 123, 164, 118, 196, 23, 73, 236, 127, 12, 111, 246, 108, 161, 59, 82, 41, 157, 85,
170, 251, 96, 134, 177, 187, 204, 62, 90, 203, 89, 95, 176, 156, 169, 160, 81, 11, 245, 22, 235, 122, 117, 44, 215, 79, 174, 213,
233, 230, 231, 173, 232, 116, 214, 244, 234, 168, 80, 88, 175 };	//used to locate the degree of finitefield element through the value of finitefield element

#endif

//**************basic var*****************
unsigned long int seq_num;	//number of input binary sequences
float SNR;
double BER, FER;

int bi_message[k*p];	//binary input message
int message[k];			//decimal input message
int codeword[n];		//encoded code word
int bi_codeword[n*p];   //encoded binary codeword
int dec_bicodeword[n*p];  //decoded binary codeword
int dec_codeword[n];    //decoded codeword
float pi = 3.141593;
float tx_symbol[n*p][2];	//modulated tranmitted signal
float rx_symbol[n*p][2];	//received destorted signal
unsigned long int seq_num_Now;
int mono_order[monoTable_Ysize][monoTable_Xsize];	//multiplexity matrix
float N0, sgm;
int gmatrix[k][n];
//***************

//************demodulation var*************
float RM[q][n];   // reliability matrix
int T_poly[k]; // store T(x)
int hard_decision[n];	//store hard decision codeword
int reen_codeword[k][3];	//store k symbols' to re-encoding
int tv_common[n - k - 1][3];	//store the symbols to common interpolation
int flip_sym_count;	//store numbers of flipping symbols
int tv_uncommon_y[tv_num_bit][eta_bit];	//store the test-vector with uncommon elements
int tv_uncommon_x[eta_bit];
int test_vec_bit[tv_num_bit][n];	
float rx_relb[2][n*p];	//store the probility of every recieve bit
int testCount[tv_num_bit];

//*******factorization by JiongYue**********
int Q[fac_size][interpoly_Ysize][interpoly_Xsize+lm+1];	//because Q_(s+1)(x,y)=Q_s(x,xy+p_s), the X_size should increase to C+lm
int times;//factorisation times
int pai[fac_size];//pai[u]:the parent vertex of u
int deg[fac_size];//deg[u]:degree of u=distance from root-1
int coeff[coeff_size];//the coefficients of factorisation polynomials
int fx[lm + 1][k];//the factorisation polynomials, the number of it is not sure, maybe overflow
int f_num;//the number of factorisation polynomials, it is not sure
int output[tv_num_bit][lm + 1][k];
int output_list_num[tv_num_bit];

int psi_codeword_temp[n]; 	// store psi 
int v_poly[k + 1]; //store v(x)
int Q_interpoly[tv_num_bit][interpoly_Ysize][interpoly_Xsize];
//int message_temp[tv_num_bit][k];
int valid_flag[tv_num_bit]; //store if the message_temp[j] of fac is vaild
int Q_uncom_elem[tv_num_bit][lm + 1][interpoly_Ysize][interpoly_Xsize];

#ifdef checkChoose
int ChosenWrong_SeqNum;
int DecSucc_flag;
double ChosenWrong_Rate;
#endif

//*******debug*************
unsigned long int addNum, mulNum, addNum_temp, mulNum_temp;
int flag_addNum = 0, flag_mulNum = 0;
//***************************

//***************basic function statement**************
int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
void mono_table(void);
void encoder(int message_temp[], int codeword_temp[]);
void modulation(void);
void channel(void);
void generator(void);
//**********************************************
//***************my function***************** 
void demodulation(void);
float PDF(int, int);
void test_set_construction(void);
void reencoding(void);
void cal_T_poly(void);
void cal_v_poly(void);
void interpolation(void);
void com_elem_interpolation(int Q_temp[][lm + 1][n + 1], int interpoint[][2]);
void uncom_tv_interpolation(int Q_temp[][interpoly_Ysize][interpoly_Xsize], int interpoint_x[], int interpoint_y[]);
int cal_delta(int Q_temp[interpoly_Ysize][interpoly_Xsize], int, int, int, int);
float comb(float, float);
void factorisation(int);
//void factorization(void);
void choose(void);

//*********factorisation by JiongYue**************
void fac(void);
void DFS(int);
void update_fac_poly(int, int, int);
void output_fac_poly(int);
int find_root(int, int);

//***************testing test_vec construction
void chase_formulation(void);



#ifdef cheatingDec
void cheatDecoding(void);
#endif


void main()
{
	int i, u, x, value, num;
	float start, finish;
	unsigned long int j;
	unsigned int mask = 1;
	long int error, ferror;
	double progress, b;
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

#ifdef cheatingEncoding
	//以读入文件方式进行预编码，而不进行真正的编码
	int file_index;
	FILE *fin1;
	if ((fin1 = fopen("bi_message.txt", "r")) == NULL)
	{
		printf("bi_message.txt open filed");
	}
	else
	{
		file_index = 0;
		while (fscanf(fin1, "%d", &bi_message[file_index]) != -1)
		{
			++file_index;
		}
	}
	fclose(fin1);

	if ((fin1 = fopen("message.txt", "r")) == NULL)
	{
		printf("message.txt open filed");
	}
	else
	{
		file_index = 0;
		while (fscanf(fin1, "%d", &message[file_index]) != -1)
		{
			++file_index;
		}
	}
	fclose(fin1);

	if ((fin1 = fopen("codeword.txt", "r")) == NULL)
	{
		printf("codeword.txt open filed");
	}
	else
	{
		file_index = 0;
		while (fscanf(fin1, "%d", &codeword[file_index]) != -1)
		{
			++file_index;
		}
	}
	fclose(fin1);

	if ((fin1 = fopen("bi_codeword.txt", "r")) == NULL)
	{
		printf("bi_codeword.txt open filed");
	}
	else
	{
		file_index = 0;
		while (fscanf(fin1, "%d", &bi_codeword[file_index]) != -1)
		{
			++file_index;
		}
	}
	fclose(fin1);

#endif

	for (SNR = start; SNR <= finish; SNR += interval)
	{
		N0 = (1.0 / (float(k) / float(n))) / pow(10.0, SNR / 10.0);
		sgm = sqrt(N0 / 2);
		b = 1.0;
		error = 0;
		ferror = 0;

		addNum_count = 0.0;
		mulNum_count = 0.0;
		totalNum_count = 0.0;

#ifdef checkChoose
		ChosenWrong_SeqNum = 0;
		ChosenWrong_Rate = 0.0;
#endif

		for (j = 1; j <= seq_num; j++)
		{
			flag_mulNum = 0;
			flag_addNum = 0;

			addNum = 0;
			mulNum = 0;

			seq_num_Now = j;

#ifndef cheatingEncoding
			//printf("\n\n*************For the %dth frame*************:\n", j);
			//Generate binary message
			for (i = 0; i<k*p; i++)
				bi_message[i] = rand() % 2;

			//Convert to nonbinary

			for (i = 0; i<k; i++)
			{
				num = 1;
				message[i] = 0;
				for (x = 0; x<p; x++)
				{
					message[i] = message[i] + (bi_message[p*i + x] * num);
					num = num * 2;
				}
			}


			encoder(message, codeword);

			//Convert the codeword into binary
			for (i = 0; i<n; i++)
			{
				value = codeword[i];
				mask = 1;
				for (x = 0; x<p; x++) //for(m=p-1;m>=0;m--)
				{
					if ((value & mask)>0)
						bi_codeword[p*i + x] = 1;
					else
						bi_codeword[p*i + x] = 0;
					mask = mask << 1;
				}
			}

#endif

			modulation();
			channel();
			demodulation(); // output soft information
			chase_formulation();

#ifdef cheatingDec
			cheatDecoding();

			//frame error rate calculation
			int temp = ferror;
			for (u = 0; u<n; u++)
				if (dec_codeword[u] != codeword[u])
				{
					ferror++;
					break;
				}

#else

			reencoding();
			interpolation();
			fac();
			//factorization();
			choose();

			//frame error rate calculation
			int temp = ferror;
			for (u = 0; u<n; u++)
				if (dec_codeword[u] != codeword[u])
				{
					ferror++;
					break;
				}
#endif

#ifdef checkChoose
			//calculate the chooseError num 
			if ((DecSucc_flag == 1) && (error>temp))
			{
				ChosenWrong_SeqNum++;
			}

			ChosenWrong_Rate = (double)ChosenWrong_SeqNum / (double)seq_num;
#endif

			addNum_count = addNum_count + (addNum - addNum_count) / (double)(j);
			mulNum_count = mulNum_count + (mulNum - mulNum_count) / (double)(j);
			totalNum_count = addNum_count + mulNum_count;

			progress = (double)(j * 100) / (double)seq_num;

			//BER = (double)(error) / (double)(n*p*j);
			FER = (double)(ferror) / (double)(j);
			//printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\r", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);
			printf("Progress=%0.1f, SNR=%2.2f, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\r", progress, SNR, ferror, FER, addNum_count, mulNum_count, totalNum_count);


			if (ferror>FrameError)
				break;
		}


		if (ferror>FrameError)
		{
			BER = (double)error / (double)(n*p*j);
			FER = (double)(ferror) / (double)(j);
		}
		else
		{
			BER = (double)error / (double)(n*p*seq_num);
			FER = (double)(ferror) / (double)(seq_num);
		}

		//printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);
		printf("Progress = %0.1f, SNR = %2.2f, Frame Errors = %2.1d, FER = %E, addNum = %0.2f, mulNum = %0.2f, total_num = %0.2f\n", progress, SNR, ferror, FER, addNum_count, mulNum_count, totalNum_count);

		OpenFile;
		//fprintf(fp, "Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);
		fprintf(fp, "Progress=%0.1f, SNR=%2.2f, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n", progress, SNR, ferror, FER, addNum_count, mulNum_count, totalNum_count);
		fclose(fp);

	}

}

void mono_table()
{
	int weight_degree[WeightDegree_Ysize][WeightDegree_Xsize];
	int i, j, count, count_temp, v;

	// initialise the weight degree table
	for (j = 0; j < WeightDegree_Ysize; j++)
		for (i = 0; i < WeightDegree_Xsize; i++)
			weight_degree[j][i] = i - j;

	//initiallise the weight degree table
	for (j = 0; j < monoTable_Ysize; j++)
		for (i = 0; i < monoTable_Xsize; i++)
			mono_order[j][i] = -1;

	count = 0;
	for (v = -(WeightDegree_Ysize - 1); v < mono_ordinSize; v++)
	{
		for (j = 0; j < WeightDegree_Ysize; j++)
		{
			for (i = 0; i < WeightDegree_Xsize; i++)
				if (weight_degree[j][i] == v)
				{
					mono_order[j][i] = count;
					count++;
					break;
				}
		}
	}

}

int power(int a, int b)
{
	int temp, pow_result = -1;

	//mulNum increasing one
	if (flag_mulNum == 1)
	{
		if (b>1)
			mulNum += (b - 1);
	}

	if (b == 0)
	{
		pow_result = 1;
		return pow_result;
	}
	else if (a == 0 && b != 0)
	{
		pow_result = 0;
		return pow_result;
	}
	else if (a>0 && b != 0)
	{
		if (flag_mulNum == 1)
		{
			mulNum += b;
		}
		temp = (logNum[a] * b) % (q - 1);
		pow_result = mularray[temp];
		return pow_result;
	}

	if (a<0)
	{
		printf("\n\n power has error!!");
	}

}

int mul(int fac1, int fac2)
{
	int mulresult = 0, temp;

	//mulNum increasing one
	if (flag_mulNum == 1)
	{
		mulNum++;
	}

	if (fac1 == 0 || fac2 == 0)
	{
		mulresult = 0;
		return mulresult;
	}
	else
	{
		temp = (logNum[fac1] + logNum[fac2]) % (q - 1);
		mulresult = mularray[temp];
		return mulresult;
	}

}

int add(int fac1, int fac2)
{
	//addNum increasing one
	if (flag_addNum == 1)
	{
		addNum++;
	}
	return (fac1 ^ fac2);
}

int inv(int fac)
{
	int i, invresult;

	if (fac == 0)
		invresult = 0;
	else
	{
		i = logNum[fac] % n;
		invresult = mularray[(n - i) % n];
	}
	return invresult;
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

void fastfadingchannel(void)
{

}

void demodulation(void)  //reconed constellation diagram
{	//soft decision
	int i, j, u, v;
	float sum, proba[pointNum];
	float Pr[p][2];

	//initialize
	for (u = 0; u < q; u++)
		for (v = 0; v < n; v++)
			RM[u][v] = 0;

	//cal PDF
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < p; j++)
		{
			sum = 0;
			for (u = 0; u < pointNum; u++)
			{
				proba[u] = PDF(i*p + j, u);
				sum += proba[u];
			}

			Pr[j][0] = proba[0] / sum;	//proba-->00	0
			Pr[j][1] = proba[1] / sum;	//proba-->01	1

			rx_relb[0][i*p + j] = Pr[j][0];
			rx_relb[1][i*p + j] = Pr[j][1];

		}

		j = 0;
#if finiteField == 16
		for (int z3 = 0; z3 < pointNum; z3++)
			for (int z2 = 0; z2 < pointNum; z2++)
				for (int z1 = 0; z1 < pointNum; z1++)
					for (int z0 = 0; z0 < pointNum; z0++)
					{
						RM[j][i] = (Pr[0][z0] * Pr[1][z1] * Pr[2][z2] * Pr[3][z3]);
						j++;
					}
#elif finiteField == 32
		for (int z4 = 0; z4<pointNum; z4++)
			for (int z3 = 0; z3 < pointNum; z3++)
				for (int z2 = 0; z2 < pointNum; z2++)
					for (int z1 = 0; z1 < pointNum; z1++)
						for (int z0 = 0; z0 < pointNum; z0++)
						{
							RM[j][i] = (Pr[0][z0] * Pr[1][z1] * Pr[2][z2] * Pr[3][z3] * Pr[4][z4]);
							j++;
						}
#elif finiteField == 64
		for (int z5 = 0; z5<pointNum; z5++)
			for (int z4 = 0; z4<pointNum; z4++)
				for (int z3 = 0; z3 < pointNum; z3++)
					for (int z2 = 0; z2 < pointNum; z2++)
						for (int z1 = 0; z1 < pointNum; z1++)
							for (int z0 = 0; z0 < pointNum; z0++)
							{
								RM[j][i] = (Pr[0][z0] * Pr[1][z1] * Pr[2][z2] * Pr[3][z3] * Pr[4][z4] * Pr[5][z5]);
								j++;
							}

#elif finiteField == 256
		for (int z7 = 0; z7<pointNum; z7++)
			for (int z6 = 0; z6<pointNum; z6++)
				for (int z5 = 0; z5<pointNum; z5++)
					for (int z4 = 0; z4<pointNum; z4++)
						for (int z3 = 0; z3<pointNum; z3++)
							for (int z2 = 0; z2<pointNum; z2++)
								for (int z1 = 0; z1<pointNum; z1++)
									for (int z0 = 0; z0 < pointNum; z0++)
									{
										RM[j][i] = (Pr[0][z0] * Pr[1][z1] * Pr[2][z2] * Pr[3][z3] * Pr[4][z4] * Pr[5][z5] * Pr[6][z6] * Pr[7][z7]);
										j++;
									}

#endif

	}


	////print codeword
	//printf("\ncodeword:\n\t");
	//for (i = 0; i < n; i++)
	//	printf("\t%d", codeword[i]);

	////print RM
	//printf("\nRM:");
	//for (i = 0; i < q; i++)
	//{
	//	printf("\n");
	//	printf("\t%d", i);
	//	for (j = 0; j < n; j++)
	//		printf("\t%0.2f", RM[i][j]);
	//}

	//printf("\n");


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

void chase_formulation(void)  //The LCC formulation approach (symbol wise sorting)
{
	int i, j, v, t, pos, count, s[p], theta[k], theta_comp[n-k], sort_sym_index[n], sort_bit_index[p*n], flip_bit_sym[2][eta_bit], temp_table[2][p], value, flip_sym[eta_bit], sym_reals[q][eta_bit], sym_reals_1[q][eta_bit], reals_num[eta_bit], location, location_count, num, bingo;
	float sym_relb[n], sym_relb_1[n], bit_llr[p*n], bit_llr_1[p*n], temp;
	unsigned int mask=1;

	//Determine the bit llrs
	for(i=0;i<n*p;i++)
	{
		bit_llr[i]=log(rx_relb[0][i]/rx_relb[1][i]);
		bit_llr_1[i]=bit_llr[i];
	}
	//sort bit |LLR|s in an ascending order
	for(i=0;i<n*p;i++)
	{
		temp=0;
		for(j=0;j<n*p;j++)
		{
			if(abs(bit_llr_1[j])>temp)
			{
				temp=abs(bit_llr_1[j]);
				v=j;
			}
		}
		sort_bit_index[n*p-1-i]=v;
		bit_llr_1[v]=0;
	}

#ifdef ChaseGameDebug
	printf("\nThe bit wise LLRs are:\n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
			printf("%0.3f ", bit_llr[p*i+j]);
		printf("\t");
	}

	printf("\nThe sorted bit indices (in ascending order) are:\n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
			printf("%d ", sort_bit_index[p*i+j]);
		printf("\t");
	}
#endif

#ifdef ChaseGameDebug
	printf("\nThe frozen sorted bit indices are:\n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
			printf("%d ", sort_bit_index[p*i+j]);
		printf("\t");
	}
#endif

	//Mark the eta bits and symbols that will be flipped
	v=0;
	for(i=0;i<n*p;i++)
	{
		if(sort_bit_index[i]!=-1)
		{
			flip_bit_sym[0][v]=sort_bit_index[i];
			flip_bit_sym[1][v]=(int)floor(sort_bit_index[i]/(float)p);
			v++;
		}
		if(v==eta_bit)
			break;
	}

#ifdef ChaseGameDebug
	printf("\nThe %d bits and symbols to be flipped are:\n", eta);
	for(j=0;j<2;j++)
	{
		for(i=0;i<eta;i++)
			printf("%d ", flip_bit_sym[j][i]);
		printf("\n");
	}
#endif

	//Determine the flipped symbols' realisations
	for(i=0;i<q;i++)
	{
		for (j = 0; j<eta_bit; j++)
		{
			sym_reals[i][j]=-1;
			sym_reals_1[i][j]=-1;
		}
	}
	for (i = 0; i<eta_bit; i++)
	{
		flip_sym[i]=-1;
		reals_num[i]=-1;
	}
	v=0; //index for the flipped symbols
	for (i = 0; i<eta_bit; i++)
	{
		bingo=0;
		for (j = 0; j<eta_bit; j++)
		{
			if (flip_bit_sym[1][j] != -1)
			{
				value=flip_bit_sym[1][j];
				flip_sym[v]=value; //Mark the flipped symbol
				bingo=1;
				break;
			}
		}

		//Start to fulfill the temp_table that leads to at most q realisations of the flipped symbol
		if(bingo==1)
		{
			//Table that remembers the realisations of the flipped bits
			for(j=0;j<2;j++)
				for(t=0;t<p;t++)
					temp_table[j][t]=-1;

			// //For the flipped bits
			for (j = 0; j<eta_bit; j++)
			{
				if(flip_bit_sym[1][j]==value) 
				{
					pos=flip_bit_sym[0][j]%p;
					temp_table[0][pos]=0;
					temp_table[1][pos]=1;
					flip_bit_sym[1][j]=-1; //Mark the bit (symbol) has been considered
				}
			}

			// //For the unflipped bits
			for(j=0;j<p;j++)
			{
				if(temp_table[0][j]==-1)
				{
					if(rx_relb[0][p*value+j]>rx_relb[1][p*value+j])
						temp_table[0][j]=0;
					else
						temp_table[1][j]=1;
				}
			}

			//Start to translate the temp_table to sym_reals
			for(t=0;t<q;t++)
			{
				//convert t into binary sequence
				mask=1;
				for(j=0;j<p;j++)
				{
					if((t & mask)>0)
						s[j]=1;
					else
						s[j]=0;
					mask=mask<<1;
				}
				//Check if realisation t will be realised
				count=0;
				for(j=0;j<p;j++)
					if(temp_table[s[j]][j]!=-1)
						count++;
				if(count==p)
					sym_reals[t][v]=t;
			}

			v++; //update the flipped symbol index
		}
	}

#ifdef ChaseGameDebug
	printf("\nThe flipped symbols are:\n");
	for(i=0;i<eta;i++)
		printf("%d\t", flip_sym[i]);
	printf("\nTheir realisations are:\n");
	for (j = 0; j < q; j++)
	{
		for (i = 0; i < eta; i++)
			printf("%d\t", sym_reals[j][i]);
		printf("\n");
	}
#endif

	//Generate the tvnum (2^eta) test-vectors
	//Initialisation the tv to be the hard-decision 
	for (i = 0; i < tv_num_bit; i++)
		for (j = 0; j < eta_bit; j++)
			tv_uncommon_y[i][j] = -1;

	for (j = 0; j < n; j++)
	{
		float relb_temp = 0.0;
		hard_decision[j] = -1;
		for (i = 0; i < q; i++)
			if (RM[i][j]>relb_temp)
			{
				relb_temp = RM[i][j];
				hard_decision[j] = i;
			}
	}

	//Generate sym_reals_1[q][eta] and reals_num[eta]
	for(j=0;j<eta_bit;j++)
	{
		v=0;
		for(i=0;i<q;i++)
		{
			if(sym_reals[i][j]!=-1)
			{
				sym_reals_1[v][j]=sym_reals[i][j];
				v++;
			}
		}
		reals_num[j]=v; //remember the number of realisations for each flipped symbol
	}

#ifdef ChaseGameDebug
	printf("\nThe refreshed realisations are:\n");
	for(j=0;j<q;j++)
	{
		for(i=0;i<eta;i++)
			printf("%d\t", sym_reals_1[j][i]);
		printf("\n");
	}
	printf("\nRealisation number for each of them is:\n");
	for(i=0;i<eta;i++)
		printf("%d\t", reals_num[i]);
#endif	

	//Count the number of flipped symbols
	flip_sym_count=0;
	for(j=0;j<eta_bit;j++)
		if(reals_num[j]>0)
			flip_sym_count++;
	//Construct the tvnum (2^eta) test-vectors
	for(i=0;i<tv_num_bit;i++)
	{
		value=i;
		mask = tv_num_bit >> 1;
		for (j = 0; j<flip_sym_count; j++)
		{
			if(flip_sym[j]!=-1)
			{
				location_count=log((float)reals_num[j])/log(2.0);
				location=0;
				num=reals_num[j]>>1;	//reals_num[]，记录每个symbol位置的翻转符号数

				for(v=0;v<location_count;v++)
				{
					if((value&mask)>0)
						location+=num;

					mask=mask>>1;
					num=num>>1;
				}

				tv_uncommon_y[i][j] = sym_reals_1[location][j];
			}
		}
	}

	//initialize tv_uncommon_x
	for (i = 0; i < eta_bit; i++)
	{
		if (i < flip_sym_count)
			tv_uncommon_x[i] = mularray[flip_sym[i]];
		else
			tv_uncommon_x[i] = -1;
	}

#ifdef cheatingDec
	//calculate the num of diffs between codeword and test_vec_bit[i]
	for (i = 0; i < tv_num_bit; i++)
	{
		testCount[i] = 0;
		for (j = 0; j < n; j++)
			if (codeword[j] != test_vec_bit[i][j])
			{
				testCount[i]++;
			}
	}
#endif

	//choose k symbols' to re-encoding
	int temp_codeword[n];

	//initialization
	j = 0;
	for (i = 0; i < n; i++)
		temp_codeword[i] = 1;
	for (i = 0; i < flip_sym_count; i++)
		temp_codeword[flip_sym[i]] = -1;

	//calculate k symbols
	for (i = 0; i < n;i++)
		if (temp_codeword[i] == 1)
		{
			temp_codeword[i] = -1;
			reen_codeword[j][0] = i;	//location
			reen_codeword[j][1] = mularray[i];	//xi
			reen_codeword[j][2] = hard_decision[i];	//yi
			j++;

			if (j >= k)
				break;
		}

	//choose the symbols to common element interpolation
	j = 0;
	for (i = 0; i < n;i++)
		if (temp_codeword[i] == 1)
		{
			temp_codeword[i] = -1;
			tv_common[j][0] = i;
			tv_common[j][1] = mularray[i];
			tv_common[j][2] = hard_decision[i];
			j++;
		}

	if (j >= n - k - 1)
		printf("\n tv_common construction has error。\n");

}



#ifdef cheatingDec
void cheatDecoding(void)
{
	int i, u, v, value, flag_judge;
	unsigned int mask = 1;
	int genius, errorCorrectionNum;

	//calculate the error correction ability
	errorCorrectionNum = (n - k) / 2;

	//start judgment
	flag_judge = 0;	//used to judge if there is correct codeword, 0-->No, 1-->Yes
	for (i = 0; i<tv_num_bit; i++)
		if (testCount[i] <= errorCorrectionNum)
		{
			flag_judge = 1;
			for (u = 0; u<n; u++)
			{
				dec_codeword[u] = codeword[u];
			}

			break;
		}


	if (!flag_judge)	//flag_judge==0 means that there is not correct codeword
	{
		for (u = 0; u<n; u++)
		{
			dec_codeword[u] = large_vec[0][u];
		}
	}


}
#endif



void reencoding()
{
	int i, u, temp1;

	flag_mulNum = 1;
	flag_addNum = 1;

	cal_T_poly();

	//1st, cal psi, the codeword_temp for re-encodeword part
	for (i = 0; i < k; i++)
	{
		//temp1 = 1;
		psi_codeword_temp[i] = 0;
		for (u = 0; u < k; u++)
		{
			//psi_codeword_temp[i] = add(psi_codeword_temp[i], mul(T_poly[u], temp1));
			//temp1 = mul(temp1, reen_codeword[i][1]);

			temp1 = mul(T_poly[u], power(reen_codeword[i][1], u));
			psi_codeword_temp[i] = add(psi_codeword_temp[i], temp1);
		}
	}

	//cal y'=y+psi for re-encodeword part
	for (i = 0; i<k; i++)
		reen_codeword[i][2] = add(reen_codeword[i][2], psi_codeword_temp[i]);

	//2nd, cal psi, the codeword_temp for common element part
	for (i = 0; i < n - k - flip_sym_count; i++)
	{
		//temp1 = 1;
		psi_codeword_temp[i + k] = 0;
		for (u = 0; u < k; u++)
		{
			//psi_codeword_temp[i] = add(psi_codeword_temp[i], mul(T_poly[u], temp1));
			//temp1 = mul(temp1, x_ordered[i]);

			temp1 = mul(T_poly[u], power(tv_common[i][1], u));
			psi_codeword_temp[i + k] = add(psi_codeword_temp[i + k], temp1);
		}
	}

	//cal y'=y+psi for common element part
	for (i = 0; i < n - k - flip_sym_count; i++)
		tv_common[i][2] = add(tv_common[i][2], psi_codeword_temp[i + k]);	//caution: the start location of psi_codeword_temp[i]

	//3rd, cal psi for uncommon element part
	for (i = 0; i < flip_sym_count; i++)
	{
		//temp1 = 1;
		psi_codeword_temp[i + n - flip_sym_count] = 0;
		for (u = 0; u < k; u++)
		{
			//psi_codeword_temp[i] = add(psi_codeword_temp[i], mul(T_poly[u], temp1));
			//temp1 = mul(temp1, x_ordered[i]);

			temp1 = mul(T_poly[u], power(tv_uncommon_x[i], u));
			psi_codeword_temp[i + n - flip_sym_count] = add(psi_codeword_temp[i + n - flip_sym_count], temp1);
		}
	}

	flag_mulNum = 0;
	flag_addNum = 0;

}

void cal_T_poly(void)
{
	/************************
	Editor: CodyWu
	Description: calculate the T_poly[] of reencoding.
				 calculate T_poly[] = T_poly[0] * 1 + T_poly[1] * x + T_poly[2] * x^2 + ··· + T_poly[k-1] * x^(k-1).
	************************/

	int i, u, v;
	int temp1, temp2;
	int poly[k], poly_temp1[k + 1], poly_temp2[k];

	//initialize
	for (i = 0; i<k; i++)
		T_poly[i] = 0;

	for (u = 0; u<k; u++)
	{
		//initialize
		for (i = 0; i<k; i++)
			poly[i] = poly_temp2[i] = 0;
		for (i = 0; i<k + 1; i++)
			poly_temp1[i] = 0;

		poly[0] = 1;
		temp1 = 1;
		for (v = 0; v<k; v++)
			if (v != u)
			{
				for (i = 0; i<k + 1; i++)
					poly_temp1[i] = 0;
				for (i = 0; i<k; i++)
					poly_temp2[i] = 0;

				//calculate poly
				for (i = 0; i<k; i++)
					if (poly[i] != 0)
					{
						poly_temp1[i + 1] = poly[i];	//x*poly[]
						poly_temp2[i] = mul(reen_codeword[v][1], poly[i]);	//x_order[v] * poly[]
					}

				for (i = 0; i<k; i++)
					poly[i] = add(poly_temp1[i], poly_temp2[i]);	//(x+x_order[v]) * poly[]

				//calculate the factor of denominator
				temp1 = mul(temp1, add(reen_codeword[u][1], reen_codeword[v][1]));

			}

		temp2 = mul(reen_codeword[u][2], inv(temp1));
		for (i = 0; i<k; i++)
			poly[i] = mul(temp2, poly[i]);

		for (i = 0; i<k; i++)
			T_poly[i] = add(T_poly[i], poly[i]);

	}

	/*	//***************debug***************
	for(i=0;i<k;i++)
	{
	temp1=0;
	for(v=0;v<k;v++)
	if(T_poly[v]!=0)
	{
	temp2=mul( power( x_ordered[i],v ),T_poly[v] );
	temp1=add( temp1,temp2 );
	}

	temp1=add( temp1,test_vec[0][i] );

	if( temp1!=0 )
	printf("\n\n\tT_poly is error");
	}
	*/	//***********************************

}

void cal_T_poly_2(void)
{
	/************************
	Editor: Jiongyue Xing
	Description: calculate the T_poly[] of reencoding
	************************/
	int i, j, h, temp, temp1, temp2;
	int reliable_i, reliable_j;
	int unreliable_i;
	int re_encoding_poly[n], Lagrange_basis_poly[k][n], temp1_poly[n + 1], temp2_poly[n];

	//initialisation
	for (i = 0; i < n + 1; i++)
		temp1_poly[i] = 0;

	for (i = 0; i < n; i++)
		temp2_poly[i] = 0;

	for (i = 0; i < n; i++)
		re_encoding_poly[i] = 0;

	for (i = 0; i < k; i++)
		T_poly[i] = 0;

	for (i = 0; i<k; i++)
	{
		Lagrange_basis_poly[i][0] = 1;
		for (j = 1; j<n; j++)
			Lagrange_basis_poly[i][j] = 0;
	}

	//construct re-encoding polynomial
	for (j = 0; j<k; j++)
	{
		reliable_j = j;
		for (i = 0; i<k; i++)
		{
			reliable_i = i;
			if (reliable_i != reliable_j)
			{
				//poly * x
				for (h = 0; h<n; h++)
					temp1_poly[h + 1] = Lagrange_basis_poly[j][h];
				//poly * alpha_i
				for (h = 0; h<n; h++)
					temp2_poly[h] = mul(reen_codeword[reliable_i][1], Lagrange_basis_poly[j][h]);
				//poly * (x+alpha_i)
				for (h = 0; h<n; h++)
					Lagrange_basis_poly[j][h] = add(temp1_poly[h], temp2_poly[h]);
				for (h = 0; h<n; h++)
					Lagrange_basis_poly[j][h] = mul(Lagrange_basis_poly[j][h], inv(add(reen_codeword[reliable_j][1], reen_codeword[reliable_i][1])));
				//initialization
				for (h = 0; h<n; h++)
				{
					temp1_poly[h] = 0;
					temp2_poly[h] = 0;
				}
				temp1_poly[n] = 0;
			}
		}
		for (i = 0; i<n; i++)
			re_encoding_poly[i] = add(re_encoding_poly[i], mul(reen_codeword[reliable_j][2], Lagrange_basis_poly[j][i]));
	}

	for (i = 0; i < k; i++)
		T_poly[i] = re_encoding_poly[i];
}

void cal_v_poly(void)
{
	int i, u;
	int poly_temp1[k + 2], poly_temp2[k + 1];

	//initialize
	for (u = 0; u<k + 1; u++)
		v_poly[u] = 0;


	v_poly[0] = 1;
	for (i = 0; i<k; i++)
	{
		//initialize
		for (u = 0; u<k + 1; u++)
		{
			poly_temp1[u] = 0;
			poly_temp2[u] = 0;
		}
		poly_temp1[k + 1] = 0;

		for (u = 0; u<k + 1; u++)
			if (v_poly[u] != 0)
			{
				poly_temp1[u + 1] = v_poly[u];
				poly_temp2[u] = mul(reen_codeword[i][1], v_poly[u]);
			}

		for (u = 0; u<k + 1; u++)
			v_poly[u] = add(poly_temp1[u], poly_temp2[u]);
	}


	/*	//************debug****************
	int temp,temp1;
	for(i=0;i<n-eta-k;i++)
	{
	temp1=0;
	//cal v(ai)
	for(u=0;u<k+1;u++)
	{
	temp=mul( v_poly[u], power( x_ordered[i],u ) );
	temp1=add( temp1,temp );
	}

	if(temp1!=0)
	printf("\n\n\tv(x) is error ");
	}

	*/	//***************************************

}



void interpolation(void)
{
	int i, j, u, v, num, temp, temp1, index_temp;
	int com_elem_interpoint[n - k - 1][2];
	int Q_com_elem[lm + 1][lm + 1][n + 1];
	int uncom_tv_temp[tv_num_bit][eta_bit];
	int degree_temp[tv_num_bit][lm+1];

	flag_mulNum = 1;
	flag_addNum = 1;

	//common element interpolation

	//initialize
	for (j = 0; j<lm + 1; j++)
		for (u = 0; u<interpoly_Ysize; u++)
			for (v = 0; v<interpoly_Xsize; v++)
				Q_com_elem[j][u][v] = 0;

	//set common element interpoint (xi,ri)
	cal_v_poly();

	//initialization
	for (i = 0; i < n - k - flip_sym_count; i++)
	{
		temp1 = 1;
		temp = 0;
		//cal v(ai)
		for (u = 0; u < k + 1; u++)
		{
			temp1 = power(tv_common[i][1], u);
			temp1 = mul(v_poly[u], temp1);
			temp = add(temp, temp1);
		}

		//set common interpolation point 
		com_elem_interpoint[i][0] = tv_common[i][1];
		com_elem_interpoint[i][1] = mul(tv_common[i][2], inv(temp));
	}
	//start
	com_elem_interpolation(Q_com_elem, com_elem_interpoint);

#ifdef checkInterpolation

	flag_mulNum = 0;
	flag_addNum = 0;
	//*********debug**********************
	for (i = 0; i < n - k - flip_sym_count; i++)
	{
		for (j = 0; j < lm + 1; j++)
		{
			temp = 0;
			for (u = 0; u < interpoly_Ysize; u++)
				for (v = 0; v < interpoly_Xsize; v++)
				{
					temp1 = mul(power(com_elem_interpoint[i][0], v), power(com_elem_interpoint[i][1], u));
					temp1 = mul(temp1, Q_com_elem[j][u][v]);
					temp = add(temp, temp1);
				}

			if (temp != 0)
				printf("\ncommon_interpoint[%d] has the common interpolation error in poly[%d]!!\n", i, j);
		}
	}
	//**************************************
	flag_mulNum = 1;
	flag_addNum = 1;

#endif

	//uncommon element interpolation
	//initialise the uncommon element interpolation
	for (i = 0; i < tv_num_bit; i++)
		for (j = 0; j < eta_bit; j++)
			uncom_tv_temp[i][j] = -1;

	for (i = 0; i < tv_num_bit; i++)
		for (j = 0; j < flip_sym_count; j++)
			uncom_tv_temp[i][j] = tv_uncommon_y[i][j];

	//set uncomon element interpolation
	for (j = 0; j < flip_sym_count;j++)
		for (i = 0; i<tv_num_bit; i++)
		{
			temp1 = 1;
			temp = 0;
			//cal v(ai)
			for (u = 0; u<k + 1; u++)
			{
				temp1 = power(tv_uncommon_x[j], u);
				temp1 = mul(v_poly[u], temp1);
				temp = add(temp, temp1);
			}
			uncom_tv_temp[i][j] = add(uncom_tv_temp[i][j], psi_codeword_temp[j + n - flip_sym_count]);
			uncom_tv_temp[i][j] = mul(uncom_tv_temp[i][j], inv(temp));
		}

	//initialize
	for (i = 0; i < tv_num_bit; i++)
		for (j = 0; j<lm + 1; j++)
			for (u = 0; u<interpoly_Ysize; u++)
				for (v = 0; v<interpoly_Xsize; v++)
					Q_uncom_elem[i][j][u][v] = Q_com_elem[j][u][v];

	//start
	for (i = 0; i < tv_num_bit; i++)
		uncom_tv_interpolation(Q_uncom_elem[i], tv_uncommon_x, uncom_tv_temp[i]);

#ifdef checkInterpolation

	flag_mulNum = 0;
	flag_addNum = 0;
	//************debug*****************
	for (i = 0; i < tv_num_bit; i++)
	{
		for (int h = 0; h < flip_sym_count; h++)
		{
			for (j = 0; j < lm + 1; j++)
			{
				temp = 0;
				for (u = 0; u < interpoly_Ysize; u++)
					for (v = 0; v < interpoly_Xsize; v++)
					{
						temp1 = mul(power(tv_uncommon_x[h], v), power(uncom_tv_temp[i][h], u));
						temp1 = mul(temp1, Q_uncom_elem[i][j][u][v]);
						temp = add(temp, temp1);
					}
				if (temp != 0)
					printf("\nTV[%d]'s uncommon_interpoint[%d] has the common interpolation error in poly[%d]!!\n", i, h, j);
			}
		}
	}
	flag_mulNum = 1;
	flag_addNum = 1;
#endif

	//initialise
	for (i = 0; i<tv_num_bit; i++)
		for (u = 0; u<interpoly_Ysize; u++)
			for (v = 0; v<interpoly_Xsize; v++)
				Q_interpoly[i][u][v] = 0;

	//choose the poly to factorization
	for (i = 0; i < tv_num_bit; i++)
	{
		//calculate the degree of poly
		for (j = 0; j < lm+1; j++)
		{
			temp = -1;
			degree_temp[i][j] = 0;
			for (u = 0; u < interpoly_Ysize; u++)
				for (v = 0; v < interpoly_Xsize; v++)
					if (Q_uncom_elem[i][j][u][v] != 0)
						if (temp < mono_order[u][v])
						{
							temp = mono_order[u][v];
						}
			degree_temp[i][j] = temp;	//deg[1,-1](x,y)
		}

		//choose the min degree
		index_temp = -1;
		temp = degree_temp[i][0];
		for (j = 1; j < lm + 1; j++)
			if (temp>degree_temp[i][j])
			{
				temp = degree_temp[i][j];
				index_temp = j;
			}
		
		//assigment
		for(u=0;u<interpoly_Ysize;u++)
			for(v=0;v<interpoly_Xsize;v++)
				Q_interpoly[i][u][v] = Q_uncom_elem[i][index_temp][u][v];

		// Q_interpoly[] update processing: Q_interpoly(x,y) = q_0(x) + y * q_1(x) --> q_0(x) * V(x) + y * q_1(x)
		int X_temp[k + 1 + interpoly_Xsize];
		int X_poly[interpoly_Xsize];
		
		for (u = 0; u < interpoly_Xsize; u++)
			X_poly[u] = 0;

		for (u = 0; u < k + 1;u++)
			if (v_poly[u] != 0)
			{
				for (v = 0; v < k + 1 + interpoly_Xsize; v++)
					X_temp[v] = 0;

				for (v = 0; v < interpoly_Xsize; v++)
					if (Q_interpoly[i][0][v] != 0)
						X_temp[v+u] = mul(v_poly[u], Q_interpoly[i][0][v]);

				for(v=0;v<interpoly_Xsize;v++)
					X_poly[v] = add(X_poly[v], X_temp[v]);
			}

		for (v = 0; v < interpoly_Xsize; v++)
			Q_interpoly[i][0][v] = X_poly[v];	//update the q_0(x)

	}
	//****************

	flag_mulNum = 0;
	flag_addNum = 0;
}


void com_elem_interpolation(int Q_temp[][interpoly_Ysize][interpoly_Xsize], int interpoint[][2])
{
	int i, j, u, v, alpha, beta, lod_temp, temp1, temp2;
	int delta[lm + 1], lod[lm + 1], J[lm + 1], lod_min, j_min;
	int f[interpoly_Ysize][interpoly_Xsize];
	int g1[interpoly_Ysize][interpoly_Xsize + 1], g2[interpoly_Ysize][interpoly_Xsize];

	//initialize the G0
	for (i = 0; i < lm + 1; i++)
		Q_temp[i][i][0] = 1;
	
	//interpolation
	for (i = 0; i<n - k - flip_sym_count; i++)
	{
		for (alpha = 0; alpha < m; alpha++)
			for (beta = 0; beta < (m - alpha); beta++)
			{
				//calculate lod[]
				for (j = 0; j < lm + 1; j++)
				{
					lod_temp = 0;
					lod[j] = -1;
					for (u = 0; u < interpoly_Ysize; u++)
						for (v = 0; v < interpoly_Xsize; v++)
							if (Q_temp[j][u][v] != 0)
							{
								lod_temp = mono_order[u][v];
								if (lod_temp > lod[j])
									lod[j] = lod_temp;
							}
				}

				//compute delta[]
				for (j = 0; j < lm + 1; j++)
				{
					delta[j] = cal_delta(Q_temp[j], interpoint[i][0], interpoint[i][1], alpha, beta);
				}

				//choose
				j_min = -1;
				lod_min = INT_MAX;
				for (j = 0; j < lm + 1; j++)
				{
					J[j] = 0;
					if (delta[j] != 0)
					{
						J[j] = 1;
						if (lod_min > lod[j])
						{
							lod_min = lod[j];
							j_min = j;
						}
					}
				}

				//update
				if (j_min != -1)
				{
					// f = Q_temp[j_min]
					for (u = 0; u < interpoly_Ysize; u++)
						for (v = 0; v < interpoly_Xsize; v++)
							f[u][v] = Q_temp[j_min][u][v];

					//g(x,y)
					for (j = 0; j < lm + 1; j++)
						if (J[j] == 1)
							if (j != j_min)	//delta==0 don't need to be update
							{
								//delta[min]*g_j - delta_j*f
								for (u = 0; u < interpoly_Ysize; u++)
									for (v = 0; v < interpoly_Xsize; v++)
									{
										if (Q_temp[j][u][v] != 0)
										{
											temp1 = mul(delta[j_min], Q_temp[j][u][v]);
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
											Q_temp[j][u][v] = add(temp1, temp2);
										}
										else
											Q_temp[j][u][v] = 0;
									}
							}
							else if (j == j_min)
							{
								//initialise the poly temp
								for (u = 0; u < interpoly_Ysize; u++)
								{
									for (v = 0; v < interpoly_Xsize + 1; v++)
										g1[u][v] = 0;
									for (v = 0; v < interpoly_Xsize; v++)
										g2[u][v] = 0;
								}

								//g1 = x * f
								for (u = 0; u < interpoly_Ysize; u++)
									for (v = 0; v < interpoly_Xsize; v++)
										if (Q_temp[j][u][v] != 0)
										{
											g1[u][v + 1] = Q_temp[j][u][v];
										}

								//g2 = xi * f
								for (u = 0; u < interpoly_Ysize; u++)
									for (v = 0; v < interpoly_Xsize; v++)
										if (Q_temp[j][u][v] != 0)
										{
											g2[u][v] = mul(interpoint[i][0], Q_temp[j][u][v]);
										}

								//Q_temp = g1 + g2
								for (u = 0; u < interpoly_Ysize; u++)
									for (v = 0; v < interpoly_Xsize; v++)
										if(g1[u][v]!=0 || g2[u][v] != 0)
										{
											Q_temp[j][u][v] = add(g1[u][v], g2[u][v]);
										}
										else
											Q_temp[j][u][v] = 0;
							}
				}
			}
	}

}

void uncom_tv_interpolation(int Q_temp[][interpoly_Ysize][interpoly_Xsize], int interpoint_x[], int interpoint_y[])
{
	int i, j, u, v, alpha, beta, lod_temp, temp1, temp2;
	int delta[lm + 1], lod[lm + 1], J[lm + 1], lod_min, j_min;
	int f[interpoly_Ysize][interpoly_Xsize];
	int g1[interpoly_Ysize][interpoly_Xsize + 1], g2[interpoly_Ysize][interpoly_Xsize];

	//interpolation
	for (i = 0; i<flip_sym_count; i++)
	{
		for (alpha = 0; alpha < m; alpha++)
			for (beta = 0; beta < (m - alpha); beta++)
			{
				//calculate lod[]
				for (j = 0; j < lm + 1; j++)
				{
					lod_temp = 0;
					lod[j] = -1;
					for (u = 0; u < interpoly_Ysize; u++)
						for (v = 0; v < interpoly_Xsize; v++)
							if (Q_temp[j][u][v] != 0)
							{
								lod_temp = mono_order[u][v];
								if (lod_temp > lod[j])
									lod[j] = lod_temp;
							}
				}

				//compute delta[]
				for (j = 0; j < lm + 1; j++)
				{
					delta[j] = cal_delta(Q_temp[j], interpoint_x[i], interpoint_y[i], alpha, beta);
				}

				//choose
				j_min = -1;
				lod_min = INT_MAX;
				for (j = 0; j < lm + 1; j++)
				{
					J[j] = 0;
					if (delta[j] != 0)
					{
						J[j] = 1;
						if (lod_min > lod[j])
						{
							lod_min = lod[j];
							j_min = j;
						}
					}
				}

				//update
				if (j_min != -1)
				{
					// f = Q_temp[j_min]
					for (u = 0; u < interpoly_Ysize; u++)
						for (v = 0; v < interpoly_Xsize; v++)
							f[u][v] = Q_temp[j_min][u][v];

					//g(x,y)
					for (j = 0; j < lm + 1; j++)
						if (J[j] == 1)
							if (j != j_min)	//delta==0 don't need to be update
							{
								//delta[min]*g_j - delta_j*f
								for (u = 0; u < interpoly_Ysize; u++)
									for (v = 0; v < interpoly_Xsize; v++)
									{
										if (Q_temp[j][u][v] != 0)
										{
											temp1 = mul(delta[j_min], Q_temp[j][u][v]);
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
											Q_temp[j][u][v] = add(temp1, temp2);
										}
										else
											Q_temp[j][u][v] = 0;
									}
							}
							else if (j == j_min)
							{
								//initialise the poly temp
								for (u = 0; u < interpoly_Ysize; u++)
								{
									for (v = 0; v < interpoly_Xsize + 1; v++)
										g1[u][v] = 0;
									for (v = 0; v < interpoly_Xsize; v++)
										g2[u][v] = 0;
								}

								//g1 = x * f
								for (u = 0; u < interpoly_Ysize; u++)
									for (v = 0; v < interpoly_Xsize; v++)
										if (Q_temp[j][u][v] != 0)
										{
											g1[u][v + 1] = Q_temp[j][u][v];
										}

								//g2 = xi * f
								for (u = 0; u < interpoly_Ysize; u++)
									for (v = 0; v < interpoly_Xsize; v++)
										if (Q_temp[j][u][v] != 0)
										{
											g2[u][v] = mul(interpoint_x[i], Q_temp[j][u][v]);
										}

								//Q_temp = g1 + g2
								for (u = 0; u < interpoly_Ysize; u++)
									for (v = 0; v < interpoly_Xsize; v++)
										if (g1[u][v] != 0 || g2[u][v] != 0)
										{
											Q_temp[j][u][v] = add(g1[u][v], g2[u][v]);
										}
										else
											Q_temp[j][u][v] = 0;
							}
				}
			}
	}

}

//int cal_delta(int Q_temp[interpoly_Ysize][interpoly_Xsize], int x, int y, int alpha, int beta)
//{
//	int u, v;
//	int temp, temp1;
//
//	temp = 0;
//	for (u = 0; u < interpoly_Ysize; u++)
//		for (v = 0; v < interpoly_Xsize; v++)
//		{
//			temp1 = mul(power(x, v), power(y, u));
//			temp1 = mul(temp1, Q_temp[u][v]);
//			temp = add(temp, temp1);
//		}
//
//	return temp;
//}

int cal_delta(int Q_temp[interpoly_Ysize][interpoly_Xsize], int x, int y, int alpha, int beta)
{
	/**************
	calculate the delta of interpolation polynomial Q[h] of (x,y)
	****************/
	int i, j, flag, delta, delta_temp;
	int a, b;

	delta = 0;

	for (j = 0; j < interpoly_Ysize; j++)
		for (i = 0; i < interpoly_Xsize; i++)
		{
			if ((Q_temp[j][i] != 0) && (j >= beta) && (i >= alpha))
			{
				a = i;
				b = j;
				flag = ((int)comb(a, alpha)) * ((int)comb(b, beta));
				// if flag is even number , delta_temp is equal to 0
				// if flag is odd number, delta_temp can be calculated to the delta
				if ((flag % 2) != 0)
				{
					delta_temp = mul(Q_temp[j][i], power(x, (a - alpha)));
					delta_temp = mul(delta_temp, power(y, (b - beta)));
					delta = add(delta, delta_temp);
				}
			}
		}

	return delta;
}

float comb(float a, float b)  
{
	/***********
	calculate combination(a_up,b_down)
	************/
	int i;
	float temp1, temp2 = 1.0;

	for (i = 0; i<(int)b; i++)
	{
		temp1 = (a - (float)i) / (b - (float)i);
		temp2 *= temp1;
	}

	return temp2;
}

void choose(void)
{
	int i, j, v, u, value;
	unsigned int mask = 1;
	int temp, min_i, min_j;
	int codeword_temp[tv_num_bit][lm + 1][n], bicodeword_temp[tv_num_bit][lm + 1][n*p];
	//int bicodeword_coor[n*p];

	flag_mulNum = 1;
	flag_addNum = 1;

	//Encoding
	//calculate the real output[][][]
	for (i = 0; i<tv_num_bit; i++)
		if (valid_flag[i] == 1)
		{
			for (j = 0; j < output_list_num[i];j++)
				for (u = 0; u<k; u++)
					output[i][j][u] = add(output[i][j][u], T_poly[u]);
		}
	
	//initialise the codeword_temp[][][]
	for (i = 0; i < tv_num_bit; i++)
		for (j = 0; j < lm + 1; j++)
			for (u = 0; u < n; u++)
				codeword_temp[i][j][u] = 0;

	for (i = 0; i<tv_num_bit; i++)
		if (valid_flag[i] == 1)
		{
			for (j = 0; j < output_list_num[i];j++)
				encoder(output[i][j], codeword_temp[i][j]);
		}

	////Convert the codeword into binary
	////initialise the bicodeword_temp[][][]
	//for (i = 0; i < tv_num_bit; i++)
	//	for (j = 0; j < lm + 1; j++)
	//		for (u = 0; u < n*p; u++)
	//			bicodeword_temp[i][j][u] = 0;

	//for (i = 0; i < tv_num_bit; i++)
	//	if (valid_flag[i] == 1)
	//		for (j = 0; j < output_list_num[i]; j++)
	//		{
	//			for (u = 0; u < n; u++)
	//			{
	//				value = codeword_temp[i][j][u];
	//				mask = 1;
	//				for (v = 0; v < p; v++)
	//				{
	//					if ((value & mask) > 0)
	//						bicodeword_temp[i][j][p*i + v] = 1;
	//					else
	//						bicodeword_temp[i][j][p*i + v] = 0;
	//					mask = mask << 1;
	//				}
	//			}
	//		}
	////rx_bicodeword, bicodeword_coor

	//for (i = 0; i<n; i++)
	//{
	//	value = large_vec[0][i];
	//	mask = 1;
	//	for (v = 0; v<p; v++)
	//	{
	//		if ((value & mask)>0)
	//			bicodeword_coor[p*i + v] = 1;
	//		else
	//			bicodeword_coor[p*i + v] = 0;
	//		mask = mask << 1;
	//	}
	//}


#if mode == 1
	//mode1: calculate Euclidean distance
	//cal symbol_temp
	for (j = 0; j<test_vec_num; j++)
		if (valid_flag[j] == 1)
			for (i = 0; i<(n*p) / 2; i++)
			{
				tx_symbol_temp[j][i][0] = -2 * bicodeword_temp[j][i * 2 + 0] + 1;	//0-->1
				tx_symbol_temp[j][i][1] = -2 * bicodeword_temp[j][i * 2 + 1] + 1;  	//1-->-1
			}


	//calculate the distance
	float temp_x, temp_y;
	float distance_temp[test_vec_num];
	min_num = -1;
	float dis_temp = 1000.0;
	for (j = 0; j<test_vec_num; j++)
		if (valid_flag[j] == 1)
		{
			distance_temp[j] = 0.0;
			for (i = 0; i<(n*p) / 2; i++)
			{
				temp_x = tx_symbol_temp[j][i][0] - rx_symbol[i][0];
				temp_y = tx_symbol_temp[j][i][1] - rx_symbol[i][1];
				distance_temp[j] += sqrt((temp_x*temp_x) + (temp_y*temp_y));
			}

			if (dis_temp>distance_temp[j])
			{
				dis_temp = distance_temp[j];
				min_num = j;
			}

		}
	//***********************
#elif mode == 2		
	//mode2: posteriori probablility
	double proba[tv_num_bit][lm+1];
	double proba_temp = 0.0;

	for (i = 0; i < tv_num_bit; i++)
		for (j = 0; j < lm + 1; j++)
			proba[i][j] = 0;

	min_i = -1;
	min_j = -1;
	for (i = 0; i<tv_num_bit; i++)
		if (valid_flag[i] == 1)
		{
			for (j = 0; j < output_list_num[i]; j++)
			{
				proba[i][j] = 1.0;
				for (v = 0; v < n; v++)
				{
					for (u = 0; u < q; u++)
						if (codeword_temp[i][j][v] == root[u])
						{
							proba[i][j] = proba[i][j] * RM[u][v];
							break;
						}
				}

				if (proba_temp < proba[i][j])
				{
					proba_temp = proba[i][j];
					min_i = i;
					min_j = j;
				}
			}

		}
	//***********************************

#elif mode == 3	
	//mode3: ML rule
	int index_c = 0;
	int index_r = 0;
	float total_temp[test_vec_num];
	for (j = 0; j<test_vec_num; j++)
		total_temp[j] = 10000.0;

	for (j = 0; j<test_vec_num; j++)
		if (valid_flag[j] == 1)
		{
			float temp_c = 0.0;
			float temp_r = 0.0;
			for (i = 0; i<n; i++)
				if (codeword_temp[j][i] != large_vec[0][i])
				{
					for (u = 0; u<q; u++)
						if (codeword_temp[j][i] == root[u])
							index_c = u;

					temp_c += log(RM[index_c][i]);

					for (u = 0; u<q; u++)
						if (large_vec[0][i] == root[u])
							index_r = u;

					temp_r += log(RM[index_r][i]);
				}

			total_temp[j] = temp_r - temp_c;
		}

	min_num = -1;
	float temp1 = 1000.0;
	for (j = 0; j<test_vec_num; j++)
		if (valid_flag[j] == 1)
			if (temp1>total_temp[j])
			{
				temp1 = total_temp[j];
				min_num = j;
			}
	//*****************************************

	/*		//*****debug*************
	printf("\n\nmin_num=%d\n\n",min_num);
	for(j=0;j<test_vec_num;j++)
	if(valid_flag[j]==1)
	printf("\ntotal_temp[%d]=%f\n",j,total_temp[j]);
	*/		//***********************
#endif

#ifdef checkChoose	
	//**********genius: optimal choose***********
	DecSucc_flag = -1;	//1--DecSucc  0--DecError but have choose  -1--DecError and no choose
	for (j = 0; j<test_vec_num; j++)
		if (valid_flag[j] == 1)
		{
			int temp = 0;
			for (i = 0; i<n; i++)
				if (message[i] == message_temp[j][i])
					temp++;

			if (temp == k)
			{
				DecSucc_flag = 1;
				break;
			}
			else if (temp != k)
			{
				DecSucc_flag = 0;
			}
		}
	//***************************
#endif

	//judge
	if (min_i == -1)	//there is no available codeword
	{
		//for (i = 0; i<n*p; i++)
		//	dec_bicodeword[i] = bicodeword[i];

		for (i = 0; i<n; i++)
			dec_codeword[i] = hard_decision[i];
	}
	else if (min_i >= 0 && min_j>=0)
	{
		//for (u = 0; u<n*p; u++)
		//	dec_bicodeword[u] = bicodeword_temp[min_i][min_j][u];

		for (u = 0; u<n; u++)
			dec_codeword[u] = codeword_temp[min_i][min_j][u];

	}
	else printf("\n\nchoose error, min_i=%d, min_j=%d\n", min_i, min_j);

	flag_mulNum = 1;
	flag_addNum = 1;
}

void factorization(void)
{
	/**********************
	Editor: CodyWu
	Description: the root factorization realised by polynomial division. It only can be applied in multiplicity of one.
	Parameters:
		1. Q_interpoly[][][]	the minimum interpolation polynomial to factorization
		2. message_temp[][]		the factorization output result
		3. valid_flag[]			the flag to show whether the tv[i] is able to factorization
	***********************/

	int i, j, u, v, index_temp;
	int X_poly_size, X_temp_size;
	int message_temp[tv_num_bit][k];
	int X_poly[tv_num_bit][n + k + 1], X_temp[n + k + 1 + k];	//X_poly[][n+k+1] = q0 * v(x)
	int divisor[2], dividend[2], quo[2];
	int div_poly[n + k + 1];

	X_poly_size = n + k + 1;
	X_temp_size = n + k + 1 + k;

	//initialize
	for (j = 0; j<tv_num_bit; j++)
		for (v = 0; v<k; v++)
			message_temp[j][v] = 0;

	for (j = 0; j<tv_num_bit; j++)
	{
		valid_flag[j] = -1;
	}

	//calculate message_temp poly
	for (j = 0; j<tv_num_bit; j++)
	{
		for (v = 0; v<X_poly_size; v++)
			X_poly[j][v] = 0;


		//calculate the real q0(x,y) * v(x)
		for (i = 0; i<k + 1; i++)
			if (v_poly[i] != 0)
			{
				for (v = 0; v<X_temp_size; v++)
					X_temp[v] = 0;

				for (v = 0; v<n + 1; v++)
					if (Q_interpoly[j][0][v] != 0)
						X_temp[v + i] = mul(v_poly[i], Q_interpoly[j][0][v]);

				for (v = 0; v<X_poly_size; v++)
					X_poly[j][v] = add(X_poly[j][v], X_temp[v]);
			}

		/**********debug***********
		printf("\n\nX_poly[%d]",j);
		for(v=0;v<X_poly_size;v++)
		printf("\t%d",X_poly[j][v]);
		*****************************/

		//calculate the div of poly
		//find the divisor[2]
		divisor[0] = divisor[1] = 0;
		for (v = n; v >= 0; v--)
			if (Q_interpoly[j][1][v] != 0)
			{
				divisor[0] = Q_interpoly[j][1][v];	//coefficient 
				divisor[1] = v;						//index
				break;
			}

		//processing
		quo[0] = 0;
		quo[1] = 0;
		dividend[0] = 0;
		dividend[1] = 0;
		valid_flag[j] = -1;
		for (i = 0; i<k + 1; i++)
		{
			//initialize
			for (v = 0; v<X_poly_size; v++)
				div_poly[v] = 0;


			//			index_temp = dividend[1];

			//find the dividend[2]
			for (v = X_poly_size - 1; v >= 0; v--)
				if (X_poly[j][v] != 0)
				{
					dividend[0] = X_poly[j][v];	//coefficient
					dividend[1] = v;				//v
					break;
				}
			/*
			if( i!=0 && index_temp != dividend[1]+1 )
			{
			while( index_temp != dividend[1]+1 )
			{
			index_temp--;
			message_temp[j][index_temp-divisor[1]] = 0;
			}
			}
			*/

			//effective solution, no remainder
			if (v == -1 && X_poly[j][0] == 0)
			{
				valid_flag[j] = 1;
				break;
			}
			//uneffective solution, remainder
			if (dividend[1]<divisor[1])
			{
				valid_flag[j] = 0;
				break;
			}
			//find the quotient
			quo[0] = mul(dividend[0], inv(divisor[0]));	//fac message's coefficient
			quo[1] = dividend[1] - divisor[1];				//fac message's index

			if ((0 <= quo[1]) && (quo[1] <= (k - 1)))
			{
				message_temp[j][quo[1]] = quo[0];

				// div_poly=divisor*quotient
				for (v = 0; v<n + 1; v++)
					if (Q_interpoly[j][1][v] != 0)
						div_poly[v + quo[1]] = mul(quo[0], Q_interpoly[j][1][v]);

				//dividend-div_poly
				for (v = 0; v<X_poly_size; v++)
					X_poly[j][v] = add(X_poly[j][v], div_poly[v]);
			}
			else if ((quo[1]<0) || (quo[1]>(k - 1)))
			{
				valid_flag[j] = 0;
				//				printf("\ninvalid message_temp[%d],quo[0]=%d,quo[1]=%d\n",j,quo[0],quo[1]);
				//				printf("\nmessage_temp[%d] has reminder",j);
				break;
			}
		}

	}

	for (j = 0; j<tv_num_bit; j++)
		if (valid_flag[j] == 1)
		{
			for (v = 0; v<k; v++)
				message_temp[j][v] = add(message_temp[j][v], T_poly[v]);
		}


	printf("\ntesting\n");
}

void fac()
{
	int i, j, u;

	flag_mulNum = 1;
	flag_addNum = 1;

	//initialise the output[][][]
	for (i = 0; i < tv_num_bit; i++)
		for (j = 0; j < lm + 1; j++)
			for (u = 0; u < k; u++)
				output[i][j][u] = -1;

	//start to factorisation
	for (i = 0; i < tv_num_bit; i++)
	{
		f_num = 0;
		factorisation(i);

		output_list_num[i] = f_num;	//TODO: Here f_num is not larger than lm+1!!

		valid_flag[i] = 0;
		if (f_num != 0)
			valid_flag[i] = 1;

		for (j = 0; j < f_num; j++)
			for (u = 0; u < k; u++)
				output[i][j][u] = fx[j][u];
	}

	flag_mulNum = 0;
	flag_addNum = 0;
}
void factorisation(int index_temp)
{
	/**********************
	Editor: Jiongyue Xing
	Description: the root factorization realised by polynomial division. It can be applied in any multiplicity value. 
	***********************/
	int i, j, h, u;
	int min_dist_pos;//the position of minimun Euclidean distance
	//initialization Q(x,y)
	for (h = 0; h<fac_size; h++)
		for (j = 0; j<interpoly_Ysize; j++)
			for (i = 0; i<interpoly_Xsize+lm+1; i++)
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
	for (i = 0; i<interpoly_Xsize; i++)
	{
		for (j = 0; j<interpoly_Ysize; j++)
		{
			if (Q_interpoly[index_temp][j][i] != 0)
			{
				h = i;
				break;
			}
		}
		if (h == i)
			break;
	}
	for (j = 0; j<interpoly_Ysize; j++)
		for (i = 0; i<(interpoly_Xsize - h); i++)
			Q[0][j][i] = Q_interpoly[index_temp][j][i + h];

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
	for (i = 0; i<(interpoly_Xsize + lm + 1); i++)
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
			coeff[v] = rootlist[i];
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
	int temp_fac_poly[interpoly_Ysize][interpoly_Xsize + lm + 1];//the x size of polynomial should plus (lm+1) since the second part of Q_u(x,xy+alpha)
	int i, j, h;
	int beta, sign, temp, temp1;

	//initialize temp_fac_poly
	for (j = 0; j<interpoly_Ysize; j++)
		for (i = 0; i<(interpoly_Xsize + lm + 1); i++)
			temp_fac_poly[j][i] = 0;
	//update Q_u(x,xy+alpha), its coeffecients Q_u[j][i] can be expressed as D_j*g_(i-j)(alpha)=(beta_up,j_down)*Q_(beta,(i-j))*alpha^(beta-j)
	for (j = 0; j<interpoly_Ysize; j++)
	{
		for (i = j; i<(interpoly_Xsize + lm + 1); i++)
		{
			temp = 0;
			for (beta = j; beta<interpoly_Ysize; beta++)
			{
				sign = (int)comb(beta, j);
				if ((sign % 2) != 0)
				{
					//temp=add(temp,mul(Q[u][beta][i-j],power(alpha,(beta-j))));
					if (Q[u][beta][i - j] != 0)
					{
						temp1 = power(alpha, (beta - j));
						temp1 = mul(Q[u][beta][i - j], temp1);
						if (temp1 != 0)
						{
							temp = add(temp, temp1);
						}
					}
				}
			}
			temp_fac_poly[j][i] = temp;
		}
	}
	//Q_v(x,y)=<<Q_u(x,y)>>
	h = -1;
	for (i = 0; i<(interpoly_Xsize + lm + 1); i++)
	{
		for (j = 0; j<interpoly_Ysize; j++)
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
	for (j = 0; j<interpoly_Ysize; j++)
		for (i = 0; i<(interpoly_Xsize + lm + 1); i++)
			Q[v][j][i] = 0;
	//update Q_v(x,y)
	for (j = 0; j<interpoly_Ysize; j++)
		for (i = 0; i<(interpoly_Xsize + lm + 1 - h); i++)
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
		fx[f_num][deg_u] = coeff[u];
		u = pai[u];
		deg_u = deg[u];
	}
	//validate the f(x) 
	temp = 0;
	for (i = 0; i<k; i++)
	{
		if (fx[f_num][i] != -1)
			temp++;
	}
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
	for (j = 0; j<interpoly_Ysize; j++)
	{
		temp = add(temp, mul(Q[u][j][0], power(root[i], j)));
	}
	return temp;
}


