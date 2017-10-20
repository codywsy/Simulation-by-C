#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/********change element with different code************
1.n,k,q,p			m-->iter_num,tm,lm,mono_size,fa_size
2.a[n],mul_array[n]
3.demodulation()--> RM[][] contruction
*****************************************************/

//mode_EuclideanDisdance	mode = 1
//mode_PostProba			mode = 2
//mode_MLrule				mode = 3
//mode_genius				mode = 4
#define mode 2

//RamdonMode
//LCC_Mode		0
//Method2-1		21
//Method2-2		22	//fixedly flip the j_{n-1} unreliable symbols
//Method2-3		23	//pick up the first Ω test_vec with the largest proba 
#define RamdonMode 23

#if RamdonMode == 0
#define OpenFile fp=fopen("[16QAM]RS255239_LCC_η=4.txt","a")
#elif RamdonMode == 21
#define OpenFile fp=fopen("[16QAM]RS255239_RPLCC-21_Ω=16.txt","a")
#elif RamdonMode == 22
#define OpenFile fp=fopen("[16QAM]RS255239_RPLCC-22_Ω=16.txt","a")
#elif RamdonMode == 23
#define OpenFile fp=fopen("[16QAM]RS255239_RPLCC-23_Ω=16.txt","a")
#endif

#define finiteField 256
#define FrameError 209
#define cheatingDec
#define cheatingEncoding
#define _16QAM_

#define eta 4 // chosen num
#define test_vec_num 64	//the num of test_vec is 2^(n-k)
#define eta_max 11	//is equals to log_{2}(test_vec_num)
#define choosen_tv_num 16


//#define checkFac
//#define checkDecoding
//#define checkPrint
//#define checkInterpolation
//#define checkChoose	

//************variable define***************
#define k 239  //length of message
//#define eta 7 // chosen num
//#define test_vec_num 128 //the num of test_vec is 2^eta
#define able_correct_num ((n-k)/2)
#define m 1 // the multiplicity of LCC algorithm is 1
#define lm 1 // LCC's lm is equal to 1, solid!!
#define V 1 //projected basic function BPSK=1, QPSK=2
#define pointNum 2 // the number of constell point == 2^V
#define choose_num 2 //the num of codeword choice, choose 2 in BASIC LCC algorithm

#define interval 0.2

#ifdef _16QAM_
//16QAM modulation
#define SMALLPROB 1.0e-12
#define len_input 4
#define len_output 2
#define num_symbol 16

double output_symbol[16][2] = {
	{ 0.9486832981, 0.9486832981 }, { 0.3162277660, 0.9486832981 }, { -0.9486832981, 0.9486832981 }, { -0.3162277660, 0.9486832981 },
	{ 0.9486832981, 0.3162277660 }, { 0.3162277660, 0.3162277660 }, { -0.9486832981, 0.3162277660 }, { -0.3162277660, 0.3162277660 },
	{ 0.9486832981, -0.9486832981 }, { 0.3162277660, -0.9486832981 }, { -0.9486832981, -0.9486832981 }, { -0.3162277660, -0.9486832981 },
	{ 0.9486832981, -0.3162277660 }, { 0.3162277660, -0.3162277660 }, { -0.9486832981, -0.3162277660 }, { -0.3162277660, -0.3162277660 }
};
int input_symbol[16][4] = {
	{ 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 0, 0, 1, 0 }, { 0, 0, 1, 1 }, { 0, 1, 0, 0 }, { 0, 1, 0, 1 }, { 0, 1, 1, 0 }, { 0, 1, 1, 1 },
	{ 1, 0, 0, 0 }, { 1, 0, 0, 1 }, { 1, 0, 1, 0 }, { 1, 0, 1, 1 }, { 1, 1, 0, 0 }, { 1, 1, 0, 1 }, { 1, 1, 1, 0 }, { 1, 1, 1, 1 }
};

double m_symbol_prob[num_symbol];
#endif


//****modify variable**********
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
//******************************
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
double RM[q][n];   // reliability matrix
int test_vec[test_vec_num][n]; //store test vectors
int test_vec_com[n - eta];	//store first n-eta symbols
float test_set_ordered[2][n]; // [0] stores value, [1] stores index
int large_vec[choose_num][n]; // largest index is 0, second largest index is 1
int x_ordered[n];
int T_poly[k]; // store T(x)
#ifdef checkDecoding
int test_vec_check[test_vec_num][eta];
#endif


int psi_codeword_temp[n]; 	// store psi 
int reen_codeword_ordered[n - eta];
int v_poly[k + 1]; //store v(x)
int Q_interpoly[test_vec_num][lm + 1][n + 1];
int message_temp[test_vec_num][k];
int valid_flag[test_vec_num]; //store if the message_temp[j] of fac is vaild
int Q_uncom_elem[test_vec_num][lm + 1][lm + 1][n + 1];

#ifdef checkFac
int checkFac_flag[test_vec_num];
#endif

#ifdef checkChoose
int ChosenWrong_SeqNum;
int DecSucc_flag;
double ChosenWrong_Rate;
#endif

#ifdef cheatingDec
int testCount1[test_vec_num];
#endif

//*******debug*************
unsigned long int addNum, mulNum, addNum_temp, mulNum_temp;
int flag_addNum = 0, flag_mulNum = 0;

#if RamdonMode == 23
int serial_with_most_proba;
#endif
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
void uncom_elem_interpolation(int interpoint[], int inGroup[][lm + 1][n + 1], int outGroup1[][lm + 1][n + 1], int outGroup2[][lm + 1][n + 1]);
void factorization(void);
void choose(void);
//***************testing test_vec construction
void test_set_construction_method2_1(void);
void test_set_construction_method2_2(void);
void test_set_construction_method2_3(void);

#ifdef _16QAM_
//16QAM
void Mapping(int *bin_cc, double *xx, int num_bit_in_blk);
void Demapping(double *symRin, int num_bit_in_blk);
void ProbClip(double *xx, int len_xx);
void channel_16QAM(double *xx, int xx_len, double *symRin, int symRin_len);
#endif


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

	flag_addNum = 1;
	flag_mulNum = 1;

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
#ifndef _16QAM_
		//N0=(1.0/(float(k)/float(n)))/pow(10.0, SNR/10.0);	//BPSK
#endif
#ifdef _16QAM_
		N0 = ((1.0 / (float(k) / float(n))) / pow(10.0, SNR / 10.0)) / ((float)(len_input));	//16QAM
#endif
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

#ifndef _16QAM_
			modulation();
			channel();
			demodulation(); // output soft information
#endif

#ifdef _16QAM_
			double xx[n*p / len_input * 2];	//输出比特序列
			Mapping(bi_codeword, xx, n*p);
			double symRin[num_symbol*(n*p / len_input)];	//输出符号概率
			channel_16QAM(xx, n*p / len_input * 2, symRin, num_symbol*(n*p / len_input));
			Demapping(symRin, n*p);


#endif


#if RamdonMode == 0
			test_set_construction();
#elif RamdonMode == 21
			test_set_construction_method2_1();
#elif RamdonMode == 22
			test_set_construction_method2_2();
#elif RamdonMode == 23
			test_set_construction_method2_3();
#endif

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

			factorization();

			choose();


			//bit error rate calculation
			int temp = error;
			for (i = 0; i<p*n; i++)
			{
				if (bi_codeword[i] != dec_bicodeword[i])
				{
					error++;
				}
			}

			//frame error rate calculation
			if (error>temp)
				ferror++;

#endif
			////////////////////////////////


#ifdef checkChoose
			//calculate the chooseError num 
			if ((DecSucc_flag == 1) && (error>temp))
			{
				ChosenWrong_SeqNum++;
			}

			ChosenWrong_Rate = (double)ChosenWrong_SeqNum / (double)seq_num;
#endif


			/*
			if(f_error_temp==1)
			printf("\nerror=%d\n",error);
			f_error_temp=0;
			*/			/*
			if(error!=0)
			ferror++;
			*/
			addNum_count = addNum_count + (addNum - addNum_count) / (double)(j);
			mulNum_count = mulNum_count + (mulNum - mulNum_count) / (double)(j);
			totalNum_count = addNum_count + mulNum_count;

			progress = (double)(j * 100) / (double)seq_num;

			BER = (double)(error) / (double)(n*p*j);
			FER = (double)(ferror) / (double)(j);
			printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\r", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);

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

		printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);

		OpenFile;
		fprintf(fp, "Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);
		fclose(fp);

	}


	getchar();


}

void mono_table()
{
	int i;

	// 1st row
	for (i = 0; i<monoTable_Xsize; i++)
		mono_order[0][i] = 2 * i + 1;

	// 2nd row
	mono_order[1][0] = 0;
	for (i = 1; i<monoTable_Xsize; i++)
		mono_order[1][i] = 2 * i;

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
		/*		for(i=0;i<n;i++)
		{
		if(fac==mularray[i])
		invresult=mularray[(n-i)%n];
		}
		*/
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

void demodulation(void)  //reconed constellation diagram
{	//soft decision
	int i, j, u, v;
	double sum, proba[pointNum];
	double Pr[p][2];

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
				proba[u] = (double)PDF(i*p + j, u);
				sum += proba[u];
			}

			Pr[j][0] = proba[0] / sum;	//proba-->00	0
			Pr[j][1] = proba[1] / sum;	//proba-->01	1

		}

		j = 0;
#if finiteField == 16
		for (v = 0; v<pointNum; v++)
			for (int x = 0; x<pointNum; x++)
				for (int y = 0; y<pointNum; y++)
					for (int h = 0; h<pointNum; h++)
					{
						RM[j][i] = (Pr[0][h] * Pr[1][y] * Pr[2][x] * Pr[3][v]);
						j++;
					}
#elif finiteField == 32
		for (u = 0; u<pointNum; u++)
			for (v = 0; v<pointNum; v++)
				for (int x = 0; x<pointNum; x++)
					for (int y = 0; y<pointNum; y++)
						for (int h = 0; h<pointNum; h++)
						{
							RM[j][i] = (Pr[0][h] * Pr[1][y] * Pr[2][x] * Pr[3][v] * Pr[4][u]);
							j++;
						}
#elif finiteField == 64
		for (int z = 0; z<pointNum; z++)
			for (u = 0; u<pointNum; u++)
				for (v = 0; v<pointNum; v++)
					for (int x = 0; x<pointNum; x++)
						for (int y = 0; y<pointNum; y++)
							for (int h = 0; h<pointNum; h++)
							{
								RM[j][i] = (Pr[0][h] * Pr[1][y] * Pr[2][x] * Pr[3][v] * Pr[4][u] * Pr[5][z]);
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

#if RamdonMode == 0
void test_set_construction(void)
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
	for (i = 0; i<n - eta; i++)
		test_vec_com[i] = large_vec[0][(int)test_set_ordered[1][i]];

	//construct x_ordered
	for (i = 0; i<n; i++)
		x_ordered[i] = mularray[(int)test_set_ordered[1][i]];

#ifdef cheatingDec
	//calculate the num of diffs between codeword and test_vec[i]
	testCount1[0] = 0;
	for (i = 0; i<n - eta; i++)
		if (codeword[(int)test_set_ordered[1][i]] != test_vec_com[i])
		{
			testCount1[0]++;
		}

	for (j = 1; j<test_vec_num; j++)
	{
		testCount1[j] = testCount1[0];
	}

	//the num of loop equals to eta
	unsigned long int mask;
	for (i = 0; i<test_vec_num; i++)
	{
		mask = 1;
		for (j = n - 1; j>n - eta - 1; j--)
		{
			//bit calculate
			if ((i & mask)>0)
			{
				temp_index = 1;
			}
			else
			{
				temp_index = 0;
			}

			if (codeword[(int)test_set_ordered[1][j]] != large_vec[temp_index][(int)test_set_ordered[1][j]])
			{
				testCount1[i]++;
			}

			mask = mask << 1;
		}
	}
#endif

}
#endif

#if RamdonMode == 21
void test_set_construction_method2_1(void)
{	//construct the test_vec[test_vec_num]
	int i, j, temp_index;
	float temp;
	float test_set[2][n]; // [0] stores value, [1] stores index
	float large_vec_proba[choose_num][n]; // largest index is 0, second largest index is 1
	unsigned long int mask;

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
	//construct 2^{n-k} = test_vec_num test_vector
	for (i = 0; i < test_vec_num; i++)
		for (j = 0; j < n; j++)
			test_vec[i][j] = large_vec[0][(int)test_set_ordered[1][j]];

	for (i = 0; i < test_vec_num; i++)
	{
		mask = 1;
		for (j = n - 1; j > n - eta_max - 1; j--)
		{
			//bit calculate
			int bit_index;
			if ((i & mask) > 0)
			{
				bit_index = 1;
			}
			else
			{
				bit_index = 0;
			}

			test_vec[i][j] = large_vec[bit_index][(int)test_set_ordered[1][j]];

			mask = mask << 1;
		}
	}

	//ramdonly choose \Omega test_vec
	int ramdon_index[choosen_tv_num];
	for (i = 0; i < choosen_tv_num; i++)
	{
		ramdon_index[i] = -1;
		int flag_repeated = 0;
		if (i == 0)
		{
			ramdon_index[i] = rand() % ((test_vec_num - 1) + 1);
		}
		else {
			while (flag_repeated == 0)
			{
				ramdon_index[i] = rand() % ((test_vec_num - 1) + 1);	//index range in [0, choose_num-1]
				for (j = 0; j < i; j++)
					if (ramdon_index[j] == ramdon_index[i])
					{
						flag_repeated = 1;
						break;
					}

				if (flag_repeated == 0)	//if flag_repeated == 1, indicates there are repeated ramdon number
					break;
				else if (flag_repeated == 1)
					flag_repeated = 0;
			}
		}

	}
	//ordering operation for ramdon_index
	for (i = 1; i < choosen_tv_num; i++)
	{
		int temp = ramdon_index[i];
		j = i - 1;
		while (ramdon_index[j] >= temp && j >= 0)
		{
			//check if some repeated elements are in ramdon_index
			if (ramdon_index[j] == temp)
				printf("\nsome repeated elements in ramdon_index\n");

			ramdon_index[j + 1] = ramdon_index[j];
			j--;
		}
		ramdon_index[j + 1] = temp;
	}



#ifdef cheatingDec
	//calculate the num of diffs between codeword and test_vec[i]
	for (i = 0; i < test_vec_num; i++)
		testCount1[i] = -1;
	
	for (i = 0; i < choosen_tv_num; i++)
		if (ramdon_index[i] != -1)
		{
			temp_index = ramdon_index[i];
			testCount1[temp_index] = 0;
			for (j = 0; j < n; j++)
			if (codeword[(int)test_set_ordered[1][j]] != test_vec[temp_index][j])
				{
					testCount1[temp_index]++;
				}
		} 
		else {
			printf("\nramdon_index has error!!\n");
		}

#endif

}
#endif

#if RamdonMode == 22
void test_set_construction_method2_2(void)	
{	//construct the test_vec[test_vec_num]
	int i, j, temp_index;
	float temp;
	float test_set[2][n]; // [0] stores value, [1] stores index
	float large_vec_proba[choose_num][n]; // largest index is 0, second largest index is 1
	unsigned long int mask;

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
	//construct 2^{n-k} = test_vec_num test_vector
	for (i = 0; i < test_vec_num; i++)
		for (j = 0; j < n; j++)
			test_vec[i][j] = large_vec[0][(int)test_set_ordered[1][j]];

	for (i = 0; i < test_vec_num; i++)
	{
		mask = 1;
		for (j = n - 1; j > n - eta_max - 1; j--)
		{
			//bit calculate
			int bit_index;
			if ((i & mask) > 0)
			{
				bit_index = 1;
			}
			else
			{
				bit_index = 0;
			}

			test_vec[i][j] = large_vec[bit_index][(int)test_set_ordered[1][j]];

			mask = mask << 1;
		}
	}

	//ramdonly choose \Omega test_vec
	int ramdon_index[choosen_tv_num];
	for (i = 0; i < choosen_tv_num; i++)
	{
		ramdon_index[i] = -1;
		int flag_repeated = 0;
		if (i == 0)
		{
			ramdon_index[i] = rand() % ((test_vec_num/2 - 1) + 1);
		}
		else {
			while (flag_repeated == 0)
			{
				ramdon_index[i] = rand() % ((test_vec_num/2 - 1) + 1);	//index range in [0, test_vec_num/2-1]
				for (j = 0; j < i; j++)
					if (ramdon_index[j] == ramdon_index[i])
					{
						flag_repeated = 1;
						break;
					}

				if (flag_repeated == 0)	//if flag_repeated == 1, indicates there are repeated ramdon number
					break;
				else if (flag_repeated == 1)
					flag_repeated = 0;
			}
		}

	}

	//KeyPart: correct the index for the test_vecs had flipped j_{n-1}
	for (i = 0; i < choosen_tv_num; i++)
		if (i<choosen_tv_num/2)
		{
			ramdon_index[i] = ramdon_index[i] * 2 + 1;	//half of indexes for the test_vecs had flipped j_{n-1}
		}
		else
		{
			ramdon_index[i] = ramdon_index[i] * 2;	//another half for the test_vecs had not flipped j_{n-1}
		}

	//ordering operation for ramdon_index
	for (i = 1; i < choosen_tv_num; i++)
	{
		int temp = ramdon_index[i];
		j = i - 1;
		while (ramdon_index[j] >= temp && j >= 0)
		{
			//check if some repeated elements are in ramdon_index
			if (ramdon_index[j] == temp)
				printf("\nsome repeated elements in ramdon_index\n");

			ramdon_index[j + 1] = ramdon_index[j];
			j--;
		}
		ramdon_index[j + 1] = temp;
	}


#ifdef cheatingDec
	//calculate the num of diffs between codeword and test_vec[i]
	for (i = 0; i < test_vec_num; i++)
		testCount1[i] = -1;

	for (i = 0; i < choosen_tv_num; i++)
		if (ramdon_index[i] != -1)
		{
			temp_index = ramdon_index[i];
			testCount1[temp_index] = 0;
			for (j = 0; j < n; j++)
				if (codeword[(int)test_set_ordered[1][j]] != test_vec[temp_index][j])
				{
					testCount1[temp_index]++;
				}
		}
		else {
			printf("\nramdon_index has error!!\n");
		}

#endif

}
#endif

#if RamdonMode == 23
void test_set_construction_method2_3(void)
{	//construct the test_vec[test_vec_num]
	int i, j, temp_index;
	float temp;
	float test_set[2][n]; // [0] stores value, [1] stores index
	float large_vec_proba[choose_num][n]; // largest index is 0, second largest index is 1
	unsigned long int mask;

	//cal test_set
	for (i = 0; i < n; i++)
	{
		//finde the largest vec
		large_vec_proba[0][i] = 0;
		for (j = 0; j < q; j++)
			if (large_vec_proba[0][i] < RM[j][i])
			{
				large_vec_proba[0][i] = RM[j][i];
				large_vec[0][i] = root[j];
			}
		//find the second largest vec
		large_vec_proba[1][i] = 0;
		for (j = 0; j < q; j++)
			if ((large_vec_proba[1][i] < RM[j][i]) && (root[j] != large_vec[0][i]))
			{
				large_vec_proba[1][i] = RM[j][i];
				large_vec[1][i] = root[j];
			}

		test_set[0][i] = large_vec_proba[1][i] / large_vec_proba[0][i];
		test_set[1][i] = i;
	}


	//order test_vec
	//key part!!
	for (i = 0; i < n; i++)
	{
		test_set_ordered[0][i] = test_set[0][i];
		test_set_ordered[1][i] = test_set[1][i];
	}

	for (i = 0; i < (n - 1); i++)
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
	//construct 2^{n-k} = test_vec_num test_vector
	for (i = 0; i < test_vec_num; i++)
		for (j = 0; j < n; j++)
			test_vec[i][j] = large_vec[0][(int)test_set_ordered[1][j]];

	for (i = 0; i < test_vec_num; i++)
	{
		mask = 1;
		for (j = n - 1; j > n - eta_max - 1; j--)
		{
			//bit calculate
			int bit_index;
			if ((i & mask) > 0)
			{
				bit_index = 1;
			}
			else
			{
				bit_index = 0;
			}

			test_vec[i][j] = large_vec[bit_index][(int)test_set_ordered[1][j]];

			mask = mask << 1;
		}
	}
	//ordering test_vec by the value of the posteriori probablility
	double log_proba[2][test_vec_num];
	for (i = 0; i < test_vec_num; i++)
	{
		log_proba[1][i] = i;
		log_proba[0][i] = 0.0;
		for (j = 0; j < n; j++)
		{
			for (int u = 0; u < q; u++)
				if (test_vec[i][j] == root[u])
				{
					//the porba in test_vec[i] is a ordered element's proba, so the column of RM shall be a ordered serial num
					log_proba[0][i] = log_proba[0][i] + log(RM[u][(int)test_set_ordered[1][j]]);	
					break;
				}
		}
	}
	//ordering proba
	for (i = 1; i < test_vec_num; i++)
	{
		double temp = log_proba[0][i];
		j = i - 1;
		while (log_proba[0][j] <= temp && j >= 0)
		{
			log_proba[0][j + 1] = log_proba[0][j];
			log_proba[1][j + 1] = log_proba[1][j];
			j--;
		}
		log_proba[0][j + 1] = temp;
		log_proba[1][j + 1] = i;
	}
	serial_with_most_proba = (int)log_proba[1][0];

#ifdef cheatingDec
	//calculate the num of diffs between codeword and test_vec[i]
	for (i = 0; i < test_vec_num; i++)
		testCount1[i] = -1;

	for (i = 0; i < choosen_tv_num; i++)
	{
		temp_index = (int)log_proba[1][i];
		testCount1[temp_index] = 0;
		for (j = 0; j < n; j++)
			if (codeword[(int)test_set_ordered[1][j]] != test_vec[temp_index][j])
			{
				testCount1[temp_index]++;
			}
	}

#endif

	////************debug***************
	//printf("\n\n");
	//printf("normal_sequence\t");
	//for (j = 0; j<n; j++)
	//	printf("\t%d", j);

	//printf("\n\n");
	//printf("large_vec[]\t");
	//for (j = 0; j<n; j++)
	//	printf("\t%d", large_vec[0][j]);

	//printf("\n\n");
	//printf("large_vec2[]\t");
	//for (j = 0; j<n; j++)
	//	printf("\t%d", large_vec[1][j]);

	//printf("\n\n");
	//printf("ordered_seq\t");
	//for (j = 0; j<n; j++)
	//	printf("\t%d", (int)test_set_ordered[1][j]);

	//printf("\n\n");
	//printf("large_vec_order[]");
	//for(j=0;j<n;j++)
	//printf("\t%d",large_vec[0][(int)test_set_ordered[1][j]]);

	//printf("\n\n");
	//printf("large_vec2_order[]");
	//for(j=0;j<n;j++)
	//printf("\t%d",large_vec[1][(int)test_set_ordered[1][j]]);

	////printf("\n\n");
	////printf("x_order[]\t");
	////for(i=0;i<n;i++)
	////printf("\t(%d,%d)",x_ordered[0][i],x_ordered[1][i]);

	//printf("\n\n");
	//for(int u=0;u<q;u++)
	//{
	//	printf("\n");
	//	printf("\t%d",root[u]);
	//	for(int v=0;v<n;v++)
	//		printf("\t%f",RM[u][v]);
	//}

	////********************************	


}
#endif

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
	for (i = 0; i<test_vec_num; i++)
		if (testCount1[i] != -1 && testCount1[i] <= errorCorrectionNum)
		{
			flag_judge = 1;
			for (u = 0; u<n; u++)
			{
				dec_codeword[u] = codeword[u];
			}

			//nonbinary --> binary
			//for(u=0;u<n;u++)
			//{	
			//	value=dec_codeword[u];
			//	mask=1;
			//	for(v=0;v<p;v++)
			//	{
			//		if((value & mask)>0)
			//			dec_bicodeword[p*u+v]=1;
			//		else
			//			dec_bicodeword[p*u+v]=0;
			//		mask=mask<<1;
			//	}
			//}	

			break;
		}


#if RamdonMode < 23
	if (!flag_judge)	//flag_judge==0 means that there is not correct codeword
	{
		//select the dec_codeword by the posteriori probablility
		double log_proba[test_vec_num];
		double log_proba_temp = 0.0;
		int flag_proba = 0;
		int max_num = -1;
		for (int j = 0; j<test_vec_num; j++)
			if (testCount1[j] > -1)
			{
				log_proba[j] = 0.0;
				for (i = 0; i<n; i++)
				{
					for (u = 0; u<q; u++)
						if (test_vec[j][i] == root[u])
						{
							//the porba in test_vec[j] is a ordered element's proba, so the column of RM shall be a ordered serial num
							log_proba[j] = log_proba[j] + log(RM[u][(int)test_set_ordered[1][i]]);
							break;
						}
				}

				if (flag_proba == 0)
				{
					log_proba_temp = log_proba[j];
					max_num = j;
					flag_proba = 1;
				}
				else {
					if (log_proba_temp < log_proba[j])
					{
						log_proba_temp = log_proba[j];
						max_num = j;
					}
				}

			}

		if (max_num != -1)
		{
			for (u = 0; u < n; u++)
			{
				dec_codeword[(int)test_set_ordered[1][u]] = test_vec[max_num][u];
			}
		}
		else {
			printf("\nThe selecting process in checkDecoding has error!\n");
		}



		//nonbinary --> binary
		//for(u=0;u<n;u++)
		//{	
		//	value=dec_codeword[u];
		//	mask=1;
		//	for(v=0;v<p;v++)
		//	{
		//		if((value & mask)>0)
		//			dec_bicodeword[p*u+v]=1;
		//		else
		//			dec_bicodeword[p*u+v]=0;
		//		mask=mask<<1;
		//	}
		//}	
	}
#elif RamdonMode >= 23
	if (!flag_judge)	//flag_judge==0 means that there is not correct codeword
	{
		for (u = 0; u < n; u++)
		{
			dec_codeword[(int)test_set_ordered[1][u]] = test_vec[serial_with_most_proba][u];
		}
	}
#endif

}
#endif



void reencoding(void)
{
	int i, u, temp1;

	cal_T_poly();

	//cal psi, the codeword_temp
	for (i = 0; i<n; i++)
	{
		temp1 = 1;
		psi_codeword_temp[i] = 0;
		for (u = 0; u<k; u++)
		{
			psi_codeword_temp[i] = add(psi_codeword_temp[i], mul(T_poly[u], temp1));
			temp1 = mul(temp1, x_ordered[i]);

			//			temp1=mul( T_poly[u],power( x_ordered[i],u ) );
			//			psi_codeword_temp[i]=add( psi_codeword_temp[i],temp1 );

		}
	}


	//cal y'=y+psi
	for (i = 0; i<n - eta; i++)
		reen_codeword_ordered[i] = add(test_vec_com[i], psi_codeword_temp[i]);

	/*	//*******debug*****************
	printf("\n\npsi_codeword_temp\t");
	for(i=0;i<n;i++)
	printf("\t%d",psi_codeword_temp[i]);

	printf("\n\nreen_codeword_ordered");
	for(i=0;i<n-eta;i++)
	printf("\t%d",reen_codeword_ordered[i]);

	printf("\n\n");
	*/	//*********************************

}


void cal_T_poly(void)
{
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
						poly_temp1[i + 1] = poly[i];
						poly_temp2[i] = mul(x_ordered[v], poly[i]);
					}

				for (i = 0; i<k; i++)
					poly[i] = add(poly_temp1[i], poly_temp2[i]);

				//calculate the factor of denominator
				temp1 = mul(temp1, add(x_ordered[u], x_ordered[v]));

			}

		temp2 = mul(test_vec_com[u], inv(temp1));
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
				poly_temp2[u] = mul(x_ordered[i], v_poly[u]);
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
	int i, j, u, v, num, temp, temp1;
	int com_elem_interpoint[n - k - eta][2];
	int Q_com_elem[lm + 1][lm + 1][n + 1];
	//	int uncom_elem_interpoint[eta][3];
	//	int Q_uncom_elem[test_vec_num][lm+1][lm+1][n+1];
	int degree_temp[test_vec_num][2];

	//common element interpolation

	//initialize
	for (j = 0; j<lm + 1; j++)
		for (u = 0; u<interpoly_Ysize; u++)
			for (v = 0; v<interpoly_Xsize; v++)
				Q_com_elem[j][u][v] = 0;

	//set common element interpoint (xi,ri)
	cal_v_poly();

	for (i = k; i<n - eta; i++)
	{
		temp1 = 1;
		temp = 0;
		//cal v(ai)
		for (u = 0; u<k + 1; u++)
		{
			//			temp =add( temp,mul(v_poly[u],temp1) );
			//			temp1=mul( temp1,x_ordered[i] );

			temp1 = power(x_ordered[i], u);
			temp1 = mul(v_poly[u], temp1);
			temp = add(temp, temp1);

		}

		//set common interpolation point 
		com_elem_interpoint[i - k][0] = x_ordered[i];
		com_elem_interpoint[i - k][1] = mul(reen_codeword_ordered[i], inv(temp));
	}

	/*	//*********debug***********
	printf("\n\ncom_elem_interpoint\n");
	for(i=0;i<n-k-eta;i++)
	printf("\t%d",com_elem_interpoint[i][0]);
	printf("\n");
	for(i=0;i<n-k-eta;i++)
	printf("\t%d",com_elem_interpoint[i][1]);

	*/	//***************************


	com_elem_interpolation(Q_com_elem, com_elem_interpoint);
	//com_elem interpolation finish

	/*	//*********debug**********************

	for(i=0;i<lm+1;i++)
	{
	printf("\n\nQ_com_elem[%d]:",i);
	for(u=0;u<lm+1;u++)
	{
	printf("\n");
	for(v=0;v<n+1;v++)
	printf("\t%d",Q_com_elem[i][u][v]);
	}
	}
	*/	//**************************************

	//uncommon element interpolation
	int uncom_elem_interpoint[eta][3];
	//set uncomon element interpolation
	for (i = n - eta; i<n; i++)
	{
		temp1 = 1;
		temp = 0;
		//cal v(ai)
		for (u = 0; u<k + 1; u++)
		{
			//			temp= add( temp,mul(v_poly[u],temp1) );
			//			temp1=mul( temp1,x_ordered[i] ); 

			temp1 = power(x_ordered[i], u);
			temp1 = mul(v_poly[u], temp1);
			temp = add(temp, temp1);

			//			temp=mul( v_poly[u], power( x_ordered[i],u ) );
			//			temp1=add( temp1,temp );

		}

		//****debug***
		//		printf("\nv(a%d)=%d",i,temp1);
		//***********

		uncom_elem_interpoint[i - n + eta][0] = x_ordered[i];
		uncom_elem_interpoint[i - n + eta][1] = mul(add(large_vec[0][(int)test_set_ordered[1][i]], psi_codeword_temp[i]), inv(temp));
		uncom_elem_interpoint[i - n + eta][2] = mul(add(large_vec[1][(int)test_set_ordered[1][i]], psi_codeword_temp[i]), inv(temp));
	}

	/*	//******debug********
	printf("\n\nuncom_elem_interpoint\n");
	for(i=0;i<eta;i++)
	printf("\t%d",uncom_elem_interpoint[i][0]);
	printf("\n");
	for(i=0;i<eta;i++)
	printf("\t%d",uncom_elem_interpoint[i][1]);
	printf("\n");
	for(i=0;i<eta;i++)
	printf("\t%d",uncom_elem_interpoint[i][2]);

	*/	//*********************

	//initialize
	for (j = 0; j<test_vec_num; j++)
		for (i = 0; i<lm + 1; i++)
			for (u = 0; u<interpoly_Ysize; u++)
				for (v = 0; v<interpoly_Xsize; v++)
					Q_uncom_elem[j][i][u][v] = 0;

	for (i = 0; i<lm + 1; i++)
		for (u = 0; u<interpoly_Ysize; u++)
			for (v = 0; v<interpoly_Xsize; v++)
				Q_uncom_elem[0][i][u][v] = Q_com_elem[i][u][v];



	//interpolation
	for (i = 0; i<eta; i++)
	{
		num = (int)pow(2.0, i);
		for (u = num - 1; u >= 0; u--)
			uncom_elem_interpolation(uncom_elem_interpoint[i], Q_uncom_elem[u], Q_uncom_elem[2 * u + 0], Q_uncom_elem[2 * u + 1]);
	}


	/*	//********debug***********
	//calculate Q_uncom_elem[0]
	for(i=0;i<eta;i++)
	{

	for(j=0;j<lm+1;j++)
	{
	temp1=0;
	for(u=0;u<lm+1;u++)
	for(v=0;v<n+1;v++)
	if(Q_uncom_elem[0][j][u][v]!=0)
	{
	temp2=mul( power(uncom_elem_interpoint[i][0],v),power(uncom_elem_interpoint[i][1],u) );
	temp2=mul( temp2,Q_uncom_elem[0][j][u][v] );
	temp1=add( temp1,temp2 );
	}

	if(temp1!=0)
	printf("\ninterpoint[%d],poly[%d] is error",i,j);
	//				else if(temp1==0)
	//					printf("\nsucceed");
	}
	}

	//calculate Q_uncom_elem[1]
	for(j=0;j<lm+1;j++)
	{
	temp1=0;
	for(u=0;u<lm+1;u++)
	for(v=0;v<n+1;v++)
	if(Q_uncom_elem[1][j][u][v]!=0)
	{
	temp2=mul( power(uncom_elem_interpoint[0][0],v),power(uncom_elem_interpoint[0][1],u) );
	temp2=mul( temp2,Q_uncom_elem[1][j][u][v] );
	temp1=add( temp1,temp2 );
	}

	if(temp1!=0)
	printf("\ninterpoint[%d],poly[%d] is error",i,j);
	//				else if(temp1==0)
	//					printf("\nsucceed");
	}

	for(j=0;j<lm+1;j++)
	{
	temp1=0;
	for(u=0;u<lm+1;u++)
	for(v=0;v<n+1;v++)
	if(Q_uncom_elem[1][j][u][v]!=0)
	{
	temp2=mul( power(uncom_elem_interpoint[1][0],v),power(uncom_elem_interpoint[1][2],u) );
	temp2=mul( temp2,Q_uncom_elem[1][j][u][v] );
	temp1=add( temp1,temp2 );
	}

	if(temp1!=0)
	printf("\ninterpoint[%d],poly[%d] is error",i,j);
	//				else if(temp1==0)
	//					printf("\nsucceed");
	}

	//calculate Q_uncom_elem[2]
	for(j=0;j<lm+1;j++)
	{
	temp1=0;
	for(u=0;u<lm+1;u++)
	for(v=0;v<n+1;v++)
	if(Q_uncom_elem[2][j][u][v]!=0)
	{
	temp2=mul( power(uncom_elem_interpoint[0][0],v),power(uncom_elem_interpoint[0][2],u) );
	temp2=mul( temp2,Q_uncom_elem[2][j][u][v] );
	temp1=add( temp1,temp2 );
	}

	if(temp1!=0)
	printf("\ninterpoint[%d],poly[%d] is error",i,j);
	//				else if(temp1==0)
	//					printf("\nsucceed");
	}

	for(j=0;j<lm+1;j++)
	{
	temp1=0;
	for(u=0;u<lm+1;u++)
	for(v=0;v<n+1;v++)
	if(Q_uncom_elem[2][j][u][v]!=0)
	{
	temp2=mul( power(uncom_elem_interpoint[1][0],v),power(uncom_elem_interpoint[1][1],u) );
	temp2=mul( temp2,Q_uncom_elem[2][j][u][v] );
	temp1=add( temp1,temp2 );
	}

	if(temp1!=0)
	printf("\ninterpoint[%d],poly[%d] is error",i,j);
	//				else if(temp1==0)
	//					printf("\nsucceed");
	}

	//calculate Q_uncom_elem[0]
	for(i=0;i<eta;i++)
	{

	for(j=0;j<lm+1;j++)
	{
	temp1=0;
	for(u=0;u<lm+1;u++)
	for(v=0;v<n+1;v++)
	if(Q_uncom_elem[3][j][u][v]!=0)
	{
	temp2=mul( power(uncom_elem_interpoint[i][0],v),power(uncom_elem_interpoint[i][2],u) );
	temp2=mul( temp2,Q_uncom_elem[3][j][u][v] );
	temp1=add( temp1,temp2 );
	}

	if(temp1!=0)
	printf("\ninterpoint[%d],poly[%d] is error",i,j);
	//				else if(temp1==0)
	//					printf("\nsucceed");
	}
	}

	*/	//************************

	//initialize
	for (j = 0; j<test_vec_num; j++)
		for (u = 0; u<interpoly_Ysize; u++)
			for (v = 0; v<interpoly_Xsize; v++)
				Q_interpoly[j][u][v] = 0;

	//choose the poly to factorization
	for (j = 0; j<test_vec_num; j++)
	{
		//calculate the degree of poly
		for (i = 0; i<2; i++)
		{
			temp = -1;
			degree_temp[j][i] = 0;
			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					if (Q_uncom_elem[j][i][u][v] != 0)
						if (temp < mono_order[u][v])
						{
							temp = mono_order[u][v];
						}
			degree_temp[j][i] = temp;	//deg[1,-1](x,y)
		}
		//choose the poly with the min degree
		if (degree_temp[j][0]<degree_temp[j][1])
		{
			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					if (Q_uncom_elem[j][0][u][v] != 0)
						Q_interpoly[j][u][v] = Q_uncom_elem[j][0][u][v];
		}
		else if (degree_temp[j][0]>degree_temp[j][1])
		{
			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					if (Q_uncom_elem[j][1][u][v] != 0)
						Q_interpoly[j][u][v] = Q_uncom_elem[j][1][u][v];
		}
		else if (degree_temp[j][0] == degree_temp[j][1])
		{
			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					if (Q_uncom_elem[j][1][u][v] != 0)
						Q_interpoly[j][u][v] = Q_uncom_elem[j][1][u][v];

			//			printf("in choose, deg[%d][0]==deg[%d][1]",j,j);
		}

	}
	//****************

}


void com_elem_interpolation(int Q_temp[][interpoly_Ysize][interpoly_Xsize], int interpoint[][2])
{
	int i, j, u, v;
	int temp1, temp2, temp3, result_temp[2], deg_temp[2], update_poly_index, index_temp[2];
	int poly_temp1[lm + 1][n - k - eta + 2], poly_temp2[lm + 1][n - k - eta + 2];


	//initialize the G0
	Q_temp[0][0][0] = 1;
	Q_temp[1][1][0] = 1;

	//interpolation
	for (i = 0; i<n - k - eta; i++)
	{
		//compute
		for (j = 0; j<2; j++)
		{
			temp2 = 0;
			for (u = 0; u<interpoly_Ysize; u++)
				for (v = 0; v<n - k - eta + 1; v++)
					if (Q_temp[j][u][v] != 0)
					{
						temp1 = mul(power(interpoint[i][0], v), power(interpoint[i][1], u));
						temp1 = mul(Q_temp[j][u][v], temp1);
						temp2 = add(temp2, temp1);
					}
			result_temp[j] = temp2;
		}

		//choose
		update_poly_index = -1;

		for (j = 0; j<2; j++)
		{
			index_temp[0] = -1;
			index_temp[1] = -1;
			temp1 = -1;
			for (u = 0; u<interpoly_Ysize; u++)
				for (v = 0; v<n - k - eta + 1; v++)
					if (Q_temp[j][u][v] != 0)
						if (temp1 < mono_order[u][v])
						{
							temp1 = mono_order[u][v];
							index_temp[1] = u;
							index_temp[0] = v;
						}

			deg_temp[j] = temp1;	//deg[1,-1](x,y)

		}

		if (result_temp[0] != 0 && result_temp[1] != 0)
		{
			if (deg_temp[0]<deg_temp[1])
				update_poly_index = 0;
			else if (deg_temp[0]>deg_temp[1])
				update_poly_index = 1;
			else if (deg_temp[0] == deg_temp[1])
				update_poly_index = 0;  // arbitrarily choose

		}
		else if (result_temp[0] != 0 && result_temp[1] == 0)
			update_poly_index = 0;
		else if (result_temp[0] == 0 && result_temp[1] != 0)
			update_poly_index = 1;
		else if (result_temp[0] == 0 && result_temp[1] == 0)
			update_poly_index = -1;

		//update
		if (update_poly_index != -1)
		{
			//g(x,y)
			for (j = 0; j<2; j++)
				if (j != update_poly_index && result_temp[j] != 0)	//delta==0 don't need to be update
				{

					temp1 = mul(result_temp[j], inv(result_temp[update_poly_index]));

					if (temp1 != 0)
					{
						for (u = 0; u<interpoly_Ysize; u++)
							for (v = 0; v<(n - k - eta + 1); v++)
								if (Q_temp[update_poly_index][u][v] != 0 || Q_temp[j][u][v] != 0)
								{
									temp2 = mul(Q_temp[update_poly_index][u][v], temp1);
									Q_temp[j][u][v] = add(Q_temp[j][u][v], temp2);
								}
					}
				}

			//f(x,y)
			if (result_temp[update_poly_index] != 0)
			{
				for (u = 0; u<interpoly_Ysize; u++)
					for (v = 0; v<n - k - eta + 2; v++)
					{
						poly_temp1[u][v] = 0;
						poly_temp2[u][v] = 0;
					}

				//xi*f(x,y)
				for (u = 0; u<interpoly_Ysize; u++)
					for (v = 0; v<(n - k - eta + 1); v++)
						if (Q_temp[update_poly_index][u][v] != 0)
							poly_temp1[u][v] = mul(interpoint[i][0], Q_temp[update_poly_index][u][v]);

				poly_temp2[0][0] = 0;
				poly_temp2[1][0] = 0;

				for (u = 0; u<interpoly_Ysize; u++)
					for (v = 0; v<(n - k - eta + 1); v++)
						poly_temp2[u][v + 1] = Q_temp[update_poly_index][u][v];

				for (u = 0; u<interpoly_Ysize; u++)
					for (v = 0; v<(n - k - eta + 1); v++)
						Q_temp[update_poly_index][u][v] = add(poly_temp1[u][v], poly_temp2[u][v]);
			}

		}
#ifdef checkInterpolation
		//********debug***************
		for (int h = 0; h<lm + 1; h++)
		{
			temp1 = 0;
			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					if (Q_temp[h][u][v] != 0)
					{
						temp2 = mul(power(interpoint[i][0], v), power(interpoint[i][1], u));
						temp2 = mul(temp2, Q_temp[h][u][v]);
						temp1 = add(temp1, temp2);
					}

			if (temp1 != 0)
				printf("\ninterpoint[%d], poly[%d] is error", i, h);
			//			else if(temp1==0)
			//				printf("\nsucceed");

		}
		//***************************
#endif
	}

}

void uncom_elem_interpolation(int interpoint[3], int inGroup[][lm + 1][n + 1], int outGroup1[][lm + 1][n + 1], int outGroup2[][lm + 1][n + 1])
{
	int i, j, u, v, a, y1_HD, y2_HD;
	int c0_0, c0_1, c1_0, c1_1;
	int result_temp[2][2], index_temp[2], deg_temp[2], update_poly_index[2];
	int temp1;
	int poly_temp1[lm + 1][n + 2], poly_temp2[lm + 1][n + 2];
	int inGroup_temp[lm + 1][lm + 1][n + 1];



	//initialize
	a = interpoint[0];
	y1_HD = interpoint[1];
	y2_HD = interpoint[2];


	for (j = 0; j<lm + 1; j++)
		for (u = 0; u<lm + 1; u++)
			for (v = 0; v<n + 1; v++)
				inGroup_temp[j][u][v] = inGroup[j][u][v];

	//compute
	c0_0 = 0;
	c0_1 = 0;
	c1_0 = 0;
	c1_1 = 0;
	for (v = 0; v<n + 1; v++)
	{
		temp1 = power(a, v);

		if (inGroup_temp[0][0][v] != 0)
		{
			c0_0 = add(c0_0, mul(inGroup_temp[0][0][v], temp1));
		}

		if (inGroup_temp[0][1][v] != 0)
		{
			c0_1 = add(c0_1, mul(inGroup_temp[0][1][v], temp1));
		}

		if (inGroup_temp[1][0][v] != 0)
		{
			c1_0 = add(c1_0, mul(inGroup_temp[1][0][v], temp1));
		}

		if (inGroup_temp[1][1][v] != 0)
		{
			c1_1 = add(c1_1, mul(inGroup_temp[1][1][v], temp1));
		}

	}

	result_temp[0][0] = add(c0_0, mul(c0_1, y1_HD));
	result_temp[0][1] = add(c1_0, mul(c1_1, y1_HD));
	result_temp[1][0] = add(c0_0, mul(c0_1, y2_HD));
	result_temp[1][1] = add(c1_0, mul(c1_1, y2_HD));

	//choose 
	update_poly_index[0] = -1;
	update_poly_index[1] = -1;

	for (j = 0; j<2; j++)
	{
		index_temp[0] = -1;
		index_temp[1] = -1;
		temp1 = -1;
		for (u = 0; u<lm + 1; u++)
			for (v = 0; v<n - k - eta + 1; v++)
				if (inGroup_temp[j][u][v] != 0)
					if (temp1 < mono_order[u][v])
					{
						temp1 = mono_order[u][v];
						index_temp[1] = u;
						index_temp[0] = v;
					}

		deg_temp[j] = temp1;

	}

	for (j = 0; j<2; j++)
	{
		if (result_temp[j][0] != 0 && result_temp[j][1] != 0)
		{
			if (deg_temp[0]<deg_temp[1])
				update_poly_index[j] = 0;
			else if (deg_temp[0]>deg_temp[1])
				update_poly_index[j] = 1;
			else if (deg_temp[0] == deg_temp[1])
				update_poly_index[j] = 0;
		}
		else if (result_temp[j][0] != 0 && result_temp[j][1] == 0)
			update_poly_index[j] = 0;
		else if (result_temp[j][0] == 0 && result_temp[j][1] != 0)
			update_poly_index[j] = 1;
		else if (result_temp[j][0] == 0 && result_temp[j][1] == 0)
			update_poly_index[j] = -1;
	}

	//update y1_HD
	if (update_poly_index[0] != -1)
	{
		//update g_1(x,y)
		for (i = 0; i<2; i++)
			if (i != update_poly_index[0])
				if (result_temp[0][i] != 0)
				{
					temp1 = mul(result_temp[0][i], inv(result_temp[0][update_poly_index[0]]));
					if (temp1 != 0)
					{
						for (u = 0; u<lm + 1; u++)
							for (v = 0; v<n + 1; v++)
								outGroup1[i][u][v] = add(inGroup_temp[i][u][v], mul(temp1, inGroup_temp[update_poly_index[0]][u][v]));
					}
					else if (temp1 == 0)
					{
						for (u = 0; u<lm + 1; u++)
							for (v = 0; v<n + 1; v++)
								outGroup1[i][u][v] = inGroup_temp[i][u][v];
					}
				}
				else if (result_temp[0][i] == 0)		//when result_temp=0, outGroup1[]=inGroup[]
				{
					for (u = 0; u<lm + 1; u++)
						for (v = 0; v<n + 1; v++)
							outGroup1[i][u][v] = inGroup_temp[i][u][v];
				}



		//update f_1(x,y),update_poly_index
		if (result_temp[0][update_poly_index[0]] != 0)
		{
			//initialization
			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					poly_temp1[u][v] = poly_temp2[u][v] = 0;


			for (u = 0; u<(lm + 1); u++)
				for (v = 0; v<(n + 1); v++)
					if (inGroup_temp[update_poly_index[0]][u][v] != 0)
						poly_temp1[u][v] = mul(a, inGroup_temp[update_poly_index[0]][u][v]);

			poly_temp2[0][0] = 0;
			poly_temp2[1][0] = 0;

			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<(n + 1); v++)
					poly_temp2[u][v + 1] = inGroup_temp[update_poly_index[0]][u][v];

			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					outGroup1[update_poly_index[0]][u][v] = add(poly_temp1[u][v], poly_temp2[u][v]);
		}
		else if (result_temp[0][update_poly_index[0]] == 0)	//when result_temp=0, outGroup1[]=inGroup[]
		{
			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					outGroup1[update_poly_index[0]][u][v] = inGroup_temp[update_poly_index[0]][u][v];
		}

	}

	if (update_poly_index[1] != -1)
	{
		//update g_2(x,y)
		for (i = 0; i<2; i++)
			if (i != update_poly_index[1])
				if (result_temp[1][i] != 0)
				{
					temp1 = mul(result_temp[1][i], inv(result_temp[1][update_poly_index[1]]));
					if (temp1 != 0)
					{
						for (u = 0; u<lm + 1; u++)
							for (v = 0; v<n + 1; v++)
								outGroup2[i][u][v] = add(inGroup_temp[i][u][v], mul(temp1, inGroup_temp[update_poly_index[1]][u][v]));
					}
					else if (temp1 == 0)
					{
						for (u = 0; u<lm + 1; u++)
							for (v = 0; v<n + 1; v++)
								outGroup2[i][u][v] = inGroup_temp[i][u][v];
					}
				}
				else if (result_temp[1][i] == 0)
				{
					for (u = 0; u<lm + 1; u++)
						for (v = 0; v<n + 1; v++)
							outGroup2[i][u][v] = inGroup_temp[i][u][v];
				}


		//update f_2(x,y)
		if (result_temp[1][update_poly_index[1]] != 0)
		{
			//initialization
			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					poly_temp1[u][v] = poly_temp2[u][v] = 0;


			for (u = 0; u<(lm + 1); u++)
				for (v = 0; v<(n + 1); v++)
					if (inGroup_temp[update_poly_index[1]][u][v] != 0)
						poly_temp1[u][v] = mul(a, inGroup_temp[update_poly_index[1]][u][v]);

			poly_temp2[0][0] = 0;
			poly_temp2[1][0] = 0;

			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<(n + 1); v++)
					poly_temp2[u][v + 1] = inGroup_temp[update_poly_index[1]][u][v];

			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					outGroup2[update_poly_index[1]][u][v] = add(poly_temp1[u][v], poly_temp2[u][v]);
		}
		else if (result_temp[1][update_poly_index[1]] == 0)
		{
			for (u = 0; u<lm + 1; u++)
				for (v = 0; v<n + 1; v++)
					outGroup2[update_poly_index[1]][u][v] = inGroup_temp[update_poly_index[1]][u][v];
		}




	}

#ifdef checkInterpolation
	//*******debug************
	for (j = 0; j<lm + 1; j++)
	{
		temp1 = 0;
		for (u = 0; u<lm + 1; u++)
			for (v = 0; v<n + 1; v++)
				if (outGroup1[j][u][v] != 0)
				{
					temp2 = mul(power(a, v), power(y1_HD, u));
					temp2 = mul(outGroup1[j][u][v], temp2);
					temp1 = add(temp1, temp2);
				}

		if (temp1 != 0)
			printf("\ninterpoint(%d,%d),Group1 is error");
		//		else if(temp1==0)
		//			printf("\nsucceed");
	}


	for (j = 0; j<lm + 1; j++)
	{
		temp1 = 0;
		for (u = 0; u<lm + 1; u++)
			for (v = 0; v<n + 1; v++)
				if (outGroup2[j][u][v] != 0)
				{
					temp2 = mul(power(a, v), power(y2_HD, u));
					temp2 = mul(outGroup2[j][u][v], temp2);
					temp1 = add(temp1, temp2);
				}

		if (temp1 != 0)
			printf("\ninterpoint(%d,%d), Group2 is error");
		//		else if(temp1==0)
		//			printf("\nsucceed");
	}
	//*********************
#endif

}


void factorization(void)
{
	int i, j, u, v, index_temp;
	int X_poly_size, X_temp_size;
	//	int Vpoly[k+1],Vpoly_temp1[k+2],Vpoly_temp2[k+1];
	int X_poly[test_vec_num][n + k + 1], X_temp[n + k + 1 + k];	//X_poly[][n+k+1] = q0 * v(x)
	int divisor[2], dividend[2], quo[2];
	int div_poly[n + k + 1];

	X_poly_size = n + k + 1;
	X_temp_size = n + k + 1 + k;

	//initialize
	for (j = 0; j<test_vec_num; j++)
		for (v = 0; v<k; v++)
			message_temp[j][v] = 0;

	for (j = 0; j<test_vec_num; j++)
	{
		valid_flag[j] = -1;
	}

	//calculate message_temp poly
	for (j = 0; j<test_vec_num; j++)
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

		/*		//**********debug***********
		printf("\n\nX_poly[%d]",j);
		for(v=0;v<X_poly_size;v++)
		printf("\t%d",X_poly[j][v]);
		*/		//*****************************
#ifdef checkFac
		//*********debug*************
		int test_poly[test_vec_num][lm + 1][n + k + 1];
		int errorCheck_temp;

		for (u = 0; u<lm + 1; u++)
			for (v = 0; v<X_poly_size; v++)
				test_poly[j][u][v] = 0;

		for (v = 0; v<X_poly_size; v++)
			test_poly[j][0][v] = X_poly[j][v];
		for (v = 0; v<n + 1; v++)
			test_poly[j][1][v] = Q_interpoly[j][1][v];

		//cal A(p,R)
		errorCheck_temp = 0;
		for (v = 0; v<n - eta; v++)
			if (large_vec[0][(int)test_set_ordered[1][v]] == codeword[(int)test_set_ordered[1][v]])
			{
				errorCheck_temp++;
			}

		if (large_vec[j][(int)test_set_ordered[1][n - 1]] == codeword[(int)test_set_ordered[1][n - 1]])
		{
			errorCheck_temp++;
		}

		//deg[1,k-1](Q)
		int temp = -1;
		int degreeFac_temp = 0;
		for (u = 0; u<test_vec_num; u++)
			for (v = 0; v<X_poly_size; v++)
				if (test_poly[j][u][v] != 0)
				{
					temp = u*(k - 1) + v;

					if (degreeFac_temp<temp)
					{
						degreeFac_temp = temp;
					}
				}

		checkFac_flag[j] = 0;
		if (errorCheck_temp > degreeFac_temp)
		{
			checkFac_flag[j] = 1;
		}
		//***************************
#endif


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

	/*	//************debug***********
	printf("\n\nv(x)\t");
	for(i=0;i<k+1;i++)
	printf("\t%d",v_poly[i]);


	for(j=0;j<test_vec_num;j++)
	{
	printf("\n\nmessage[%d]",j);
	for(v=0;v<k;v++)
	printf("\t%d",message_temp[j][v]);
	}

	*/	//******************************

	/*	//*******debug****************
	int test_codeword[n];
	int temp,temp1,degree_temp;



	for(i=0;i<n;i++)
	test_codeword[i]=add(codeword[i],psi_codeword_temp[i]);

	for(j=0;j<test_vec_num;j++)
	{
	temp=0;
	for(i=0;i<n;i++)
	if(test_codeword[i]==reen_codeword_ordered[j][i])
	temp++;

	temp1=-1;
	degree_temp=0;
	for(u=0;u<lm+1;u++)
	for(v=0;v<X_poly_size;v++)
	if(test_poly[j][u][v]!=0 )
	if(temp1 < mono_order[u][v] )
	{
	temp1=mono_order[u][v];
	degree_temp=v-u;
	}

	//calculate the codeword score
	if( (temp>degree_temp) && ( valid_flag[j]!=1))
	{
	printf("\n\n%d factorization error, valid_flag[%d]=%d, temp=%d, degree_temp=%d\n",j,j,valid_flag[j],temp,degree_temp);

	printf("\nmessage=\t");
	for(i=0;i<k;i++)
	printf("%d\t",message[i]);
	printf("\n");
	}

	if(valid_flag[j]==-1)
	{
	printf("\n\nvalid_flag[%d] error,temp=%d, degree_temp=%d\n",j,temp,degree_temp);
	printf("\n\nQ_interpoly[%d]\n",j);
	for(u=0;u<lm+1;u++)
	{
	printf("\n");
	for(v=0;v<X_poly_size;v++)
	printf("\t%d",test_poly[j][u][v]);
	}
	printf("\n\n");

	}

	}


	*/	//****************************





}

void choose(void)
{
	int i, j, v, u, value;
	unsigned int mask = 1;
	int temp, min_num = -1;
	int codeword_temp[test_vec_num][n], bicodeword_temp[test_vec_num][n*p], bicodeword_coor[n*p];

#ifdef checkFac
	for (i = 0; i<test_vec_num; i++)
		if (checkFac_flag[i] == 1 && valid_flag[i] != 1)
			printf("\n\ntest_vec[%d] factorization failed, flagForCondition=%d, flagForFac=%d\n\n", i, checkFac_flag[i], valid_flag[i]);
#endif



	min_num = -1;
	//Encoding
	//calculate the real message_temp
	for (j = 0; j<test_vec_num; j++)
		if (valid_flag[j] == 1)
		{
			for (v = 0; v<k; v++)
				message_temp[j][v] = add(message_temp[j][v], T_poly[v]);
		}
	/*	//************debug***********
	printf("\n\nT(x)\t");
	for(v=0;v<k;v++)
	printf("\t%d",T_poly[v]);

	for(j=0;j<test_vec_num;j++)
	{
	printf("\n\nreal_message[%d]",j);
	for(v=0;v<k;v++)
	printf("\t%d",message_temp[j][v]);
	}

	printf("\n\nmessage\t");
	for(v=0;v<k;v++)
	printf("\t%d",message[v]);

	*/	//******************************


	//encoding		
	for (j = 0; j<test_vec_num; j++)
		if (valid_flag[j] == 1)
		{
			encoder(message_temp[j], codeword_temp[j]);
		}


	//Convert the codeword into binary
	for (j = 0; j<test_vec_num; j++)
		if (valid_flag[j] == 1)
			for (i = 0; i<n; i++)
			{
				value = codeword_temp[j][i];
				mask = 1;
				for (v = 0; v<p; v++)
				{
					if ((value & mask)>0)
						bicodeword_temp[j][p*i + v] = 1;
					else
						bicodeword_temp[j][p*i + v] = 0;
					mask = mask << 1;
				}
			}

	//rx_bicodeword, bicodeword_coor

	for (i = 0; i<n; i++)
	{
		value = large_vec[0][i];
		mask = 1;
		for (v = 0; v<p; v++)
		{
			if ((value & mask)>0)
				bicodeword_coor[p*i + v] = 1;
			else
				bicodeword_coor[p*i + v] = 0;
			mask = mask << 1;
		}
	}


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
	float proba[test_vec_num];
	float proba_temp = 0.0;
	min_num = -1;
	for (j = 0; j<test_vec_num; j++)
		if (valid_flag[j] == 1)
		{
			proba[j] = 1.0;
			for (i = 0; i<n; i++)
			{
				for (u = 0; u<q; u++)
					if (codeword_temp[j][i] == root[u])
					{
						proba[j] = proba[j] * RM[u][i];
						break;
					}
			}

			if (proba_temp < proba[j])
			{
				proba_temp = proba[j];
				min_num = j;
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
	if (min_num == -1)	//there is no available codeword
	{
		//		printf("\ndecoding failed\n");

		for (i = 0; i<n*p; i++)
			dec_bicodeword[i] = bicodeword_coor[i];

		for (i = 0; i<n; i++)
			dec_codeword[i] = large_vec[0][i];

	}
	else if ((min_num >= 0) && (min_num<test_vec_num))
	{
		for (i = 0; i<n*p; i++)
			dec_bicodeword[i] = bicodeword_temp[min_num][i];

		for (i = 0; i<n; i++)
			dec_codeword[i] = codeword_temp[min_num][i];

	}
	else printf("\n\nchoose error, min_num=%d\n", min_num);

#ifdef checkDecoding
	//**********debug******************
	flag_addNum = 0;
	flag_mulNum = 0;
	//	unsigned int mask = 1;
	int flag = 0, index_temp, choosen_temp;
	int epcount1[test_vec_num], epcount2[test_vec_num];	// the num of uncorrect codeword

	//initialize
	for (i = 0; i<test_vec_num; i++)
	{
		for (j = 0; j<eta; j++)
			test_vec_check[i][j] = -1;

		epcount1[i] = 0;
		epcount2[i] = 0;
	}

	//start
	for (i = 0; i<test_vec_num; i++)
	{
		mask = 1;
		for (j = 0; j<eta; j++)
		{
			if ((i&mask) > 0)
				choosen_temp = 1;
			else
				choosen_temp = 0;

			index_temp = n - j - 1;
			test_vec_check[i][eta - 1 - j] = large_vec[choosen_temp][(int)test_set_ordered[1][index_temp]];
			mask = mask << 1;
		}
	}

	//error count
	for (j = 0; j<n - eta; j++)
		if (codeword[(int)test_set_ordered[1][j]] != test_vec_com[j])
			epcount1[0]++;

	for (i = 1; i<test_vec_num; i++)
		epcount1[i] = epcount1[0];

	for (i = 0; i<test_vec_num; i++)
	{
		for (j = 0; j<eta; j++)
		{
			if (codeword[(int)test_set_ordered[1][n - eta + j]] != test_vec_check[i][j])
				epcount1[i]++;
		}
	}

	for (j = 0; j<test_vec_num; j++)
	{
		for (i = 0; i<k; i++)
			if (message[i] != message_temp[j][i])
			{
				epcount2[j]++;	//epcount!=0 means unable to decoding
			}
	}

	//************strange part****************
	for (j = 0; j<test_vec_num; j++)
	{
		flag = 0;
		for (i = 0; i<k; i++)
			if (message[i] == message_temp[j][i])
				flag++;

		if (flag == k)
			break;
	}
	//****************************

	for (i = 0; i<test_vec_num; i++)
	{
		if ((epcount1[i] <= able_correct_num) && (epcount2[j] != 0))
		{
			printf("\n\ntest_vec[%d] decoding failed, errorNum[%d]=%d\n", i, i, epcount1[i]);

			printf("codeword_order:\n");
			for (j = 0; j<n; j++)
				printf("\t%d", codeword[(int)test_set_ordered[1][j]]);
			printf("\n\n");

			printf("test_vec[%d]:\n", i);
			for (j = 0; j<n - eta; j++)
				printf("\t%d", large_vec[0][(int)test_set_ordered[1][j]]);
			for (j = 0; j<eta; j++)
				printf("\t%d", test_vec_check[i][j]);
			printf("\n\n");
		}

	}
	flag_addNum = 1;
	flag_mulNum = 1;
	//*************************************
#endif
}

#ifdef _16QAM_
void Mapping(int *bin_cc, double *xx, int num_bit_in_blk)
{
	int i, j, u, v;
	int num_point_in_blk;
	int symbol;

	num_point_in_blk = num_bit_in_blk / len_input;//序列表示星座点个数

	u = 0;
	v = 0;
	for (i = 0; i < num_point_in_blk; i++)
	{
		symbol = 0;
		for (j = 0; j < len_input; j++)
		{
			//symbol = (symbol << 1) | bin_cc[u];
			symbol = (symbol << 1) + bin_cc[u];
			u++;
		}
		for (j = 0; j < len_output; j++)
		{
			xx[v] = output_symbol[symbol][j];
			v++;
		}
	}

	//printf("\ncodeword:\n");
	//for (i = 0; i < n; i++)
	//	printf("\t%d", codeword[i]);
	//printf("\n");

	//printf("\nbitcc:\n");
	//for (i = 0; i < n; i++)
	//{
	//	printf("\t");
	//	for (j = 0; j < p; j++)
	//		printf("%d", bin_cc[i*p + j]);
	//	//printf("\t");
	//}
	//printf("\n");

	//printf("\nMappingSymbol:\n");
	//for (i = 0; i < v / 2; i++)
	//	printf("  (%.3f, %.3f)", xx[i * 2 + 0], xx[i * 2 + 1]);
	//printf("\n\n");

	return;
}

void Demapping(double *symRin, int num_bit_in_blk)
{
	int i, j, u, v;
	double bitLout[n*p], bitLout_temp[n*p][2];
	int num_point_in_blk;
	double P_0, P_1;
	double sum;

	num_point_in_blk = num_bit_in_blk / len_input;

	ProbClip(symRin, num_symbol*num_point_in_blk);

	for (i = 0; i < num_point_in_blk; i++) //计算 prob(x=0 | y);
	{
		for (u = 0; u < num_symbol; u++)
		{
			m_symbol_prob[u] = 1.0;
		}

		//compute symbol full message
		sum = 0.0;
		v = i * num_symbol;
		for (u = 0; u < num_symbol; u++)
		{
			m_symbol_prob[u] *= symRin[v];
			sum += m_symbol_prob[u];
			v++;
		}

		//normalized
		for (u = 0; u < num_symbol; u++)
			m_symbol_prob[u] /= sum;

		//compute bit extrinsic message for output
		v = i * len_input;
		for (j = 0; j < len_input; j++)
		{
			P_0 = 0;
			P_1 = 0;
			for (u = 0; u < num_symbol; u++)
			{
				if (input_symbol[u][j] == 0)
					P_0 += m_symbol_prob[u];
				else
					P_1 += m_symbol_prob[u];
			}
			//P_0 /= bitLin[v];
			//P_1 /= (1.0-bitLin[v]);
			bitLout[v] = P_0 / (P_0 + P_1);
			v++;
		}
	} //end of for (i = 0; i < num_point_in_blk; i++)

	ProbClip(bitLout, num_bit_in_blk);

	//bitLout[] transforms to RM[][]
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < p; j++)
		{
			bitLout_temp[i*p + j][0] = bitLout[i*p + j];
			bitLout_temp[i*p + j][1] = 1 - bitLout[i*p + j];
		}

		j = 0;
#if finiteField == 16
		for (int z3 = 0; z3 < pointNum; z3++)
			for (int z2 = 0; z2 < pointNum; z2++)
				for (int z1 = 0; z1 < pointNum; z1++)
					for (int z0 = 0; z0 < pointNum; z0++)
					{
						RM[j][i] = bitLout_temp[i*p + 0][z0] * bitLout_temp[i*p + 1][z1] * bitLout_temp[i*p + 2][z2] * bitLout_temp[i*p + 3][z3];
						j++;
					}

#elif finiteField == 64
		for (int z5 = 0; z5 < pointNum; z5++)
			for (int z4 = 0; z4 < pointNum; z4++)
				for (int z3 = 0; z3 < pointNum; z3++)
					for (int z2 = 0; z2 < pointNum; z2++)
						for (int z1 = 0; z1 < pointNum; z1++)
							for (int z0 = 0; z0 < pointNum; z0++)
							{
								RM[j][i] = bitLout_temp[i*p + 0][z0] * bitLout_temp[i*p + 1][z1] * bitLout_temp[i*p + 2][z2] * bitLout_temp[i*p + 3][z3] * bitLout_temp[i*p + 4][z4] * bitLout_temp[i*p + 5][z5];
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
										RM[j][i] = bitLout_temp[i*p + 0][z0] * bitLout_temp[i*p + 1][z1] * bitLout_temp[i*p + 2][z2] * bitLout_temp[i*p + 3][z3]
											* bitLout_temp[i*p + 4][z4] * bitLout_temp[i*p + 5][z5] * bitLout_temp[i*p + 6][z6] * bitLout_temp[i*p + 7][z7];
										j++;
									}

#endif
	}

	//printf("\nRM:\n");
	//for (j = 0; j < q; j++)
	//{
	//	printf("\t%d", j);
	//	for (i = 0; i < n; i++)
	//		printf("\t%0.3f", RM[j][i]);
	//	printf("\n");
	//}
	//printf("\n");

	return;
}

void ProbClip(double *xx, int len_xx)
{
	int i;

	for (i = 0; i < len_xx; i++){
		if (xx[i] < SMALLPROB)
			xx[i] = SMALLPROB;
		else if (xx[i] > 1.0 - SMALLPROB)
			xx[i] = 1.0 - SMALLPROB;
	}

	return;
}

void channel_16QAM(double *xx, int xx_len, double *symRin, int symRin_len)
{
	int i, j, v;
	float u, r, g;
	double temp1;

	//Add AWGN 
	for (i = 0; i < xx_len; i++)
	{
		u = (float)rand() / (float)RAND_MAX;
		if (u == 1.0)
			u = 0.999999;
		r = sgm*sqrt(2.0*log(1.0 / (1.0 - u)));

		u = (float)rand() / (float)RAND_MAX;
		if (u == 1.0)
			u = 0.999999;
		g = (float)r*cos(2 * pi*u);

		xx[i] = xx[i] + g;
	}


	//printf("\nMappingSymbol + Noise:\n");
	//for (i = 0; i < xx_len / 2; i++)
	//	printf("  (%.3f, %.3f)", xx[i * 2 + 0], xx[i * 2 + 1]);
	//printf("\n\n");


	//calculate symRin
	v = 0;
	for (i = 0; i < xx_len / 2; i++)
	{
		for (j = 0; j < num_symbol; j++)
		{
			temp1 = pow(xx[i * 2 + 0] - output_symbol[j][0], 2) + pow(xx[i * 2 + 1] - output_symbol[j][1], 2);
			symRin[v] = exp(-temp1 / N0);
			v++;
		}
	}

	//printf("\nsymRin:\n");
	//for (i = 0; i < num_symbol; i++)
	//	printf("\t%d", i);
	//printf("\n");

	//for (i = 0; i < xx_len/2; i++)
	//{
	//	printf("%d", i);
	//	for (j = 0; j < num_symbol; j++)
	//		printf("\t%0.3f",symRin[i*num_symbol + j]);
	//	printf("\n");
	//}
}
#endif