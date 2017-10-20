#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "FiniteFieldBasisGF(16).h"
#include "FiniteFieldBasisCal.h"

#define _NoReduction_
//#define _PolyCoeffNumUncom_
//#define _PolyCoeffNumFac_
#define myWay
#define OpenFile fp=fopen("Herm(64,49)_GS_optimal.txt","a")
#define FrameError 229

#define cheatingEncoding

//#define printCoeffTable
//#define printDemodulation
//#define printfMonotable

//global define
#define w 4
#define weiz 54
#define multiplicity 1 //multiplicty 
#define lm 1	//design length
#define able_correct 4
#define pointNum 2
#define interval 0.5


/*
//tgorder()
#define tg_size 275	//more than k, tg_size represent the probably used pole basis num, (w+1)*w/2 + (w+1) * ( (w+interpoly_Ysize)-w+1 )
//mono_table()
#define weight_Zsize 4
#define weight_XYsize 275	//equals to the tg_size
#define mono_ordinSize 700	//choose the value with debugging
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 67	//large than interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize 350 //>index of term[max(deg_y)][w] in the pole basis + 1, nearly equals to tg_size
//polebasis()
*/

//tgorder()
#define tg_size 815	//more than k, tg_size represent the probably used pole basis num, (w+1)*w/2 + (w+1) * ( (w+interpoly_Ysize)-w+1 ) 
//mono_table()
#define weight_Zsize 35
#define weight_XYsize 815	//equals to the tg_size
#define mono_ordinSize 855	//choose the value with debugging
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 166	//large than interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize 815 //>index of term[max(deg_y)][w] in the pole basis + 1, nearly equals to tg_size
//*******************

//******varible*********

//main()
unsigned long int seq_num;	//number of input binary sequences
float SNR, N0;
double BER, FER;
float pi = 3.141593;	// pai
int iterNum = 0; //iteration num
int bi_message[k*p], message[k];	//transmitted messge[k]
int codeword[n], bi_codeword[n*p]; //codewords
float tx_symbol[p*n][2], rx_symbol[p*n][2], sgm;
int rxcodeword[n], bi_rxcodeword[n*p];	//receive codeword

float RM[q][n];	//Reliability Matrix
int MM[q][n];	//Multiplicity Matrix
//polebasis()

//findpoint()
int point[2][n];	//rational points[2][w^3], [0]->x, [1]->y
//tgorder()
int tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
//mono_table()
int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
int gmatrix[k][n];	//generator matrix[k][n]
int dec_message[k], dec_codeword[n], dec_bicodeword[n*p]; //decoding result
//MatrixConvert()
int s = 0;	//large enough to make the circle of MatrixConvert() keep moving
//interpolation
int epcount1, epcount2, testCount1, testCount2, testCount1_com;


//****debug**********
#ifdef checkFac
int Deg_iterNum, Deg_poly, codewordScore;
#endif

//some var about choosing incorrectly
unsigned long int seq_num_Now;	//record the seq_num at the moment
unsigned long int DecSucc_SeqNum, ChosenWrong_SeqNum;	//record the seq_num should be decoding success in genius mode
double CWR;	//record the rate of chosen uncorrectly

//some var about computation complexity
unsigned long int addNum, mulNum, addNum_temp, mulNum_temp;
int flag_addNum = 0, flag_mulNum = 0;

//*******************


//**********function**********
void findpoints(void);
void mono_table(void);
//void polyexp(int, int, int);
//void zerobasis(void);
//void polebasis(void);
void tgorder(void);
//void coefficient(void);
//void test_vec_contruction(void);
void interpolation(void);
void factorisation(void);
void rcs(int);
void choose(void);
void generator(void);
void encoder(int message_temp[], int codeword_temp[]);
void modulation(void);
void channel(void);
void demodulation(void);
void demodulation_GS(void);
float PDF(int, int);
//*****************
//void MatrixConvert(void);
//int cal_max_a();
//int gama(int);
//void zerobasis(int, int LM_zb[][zb_alpha]);
//void coefficientSearch();
//****************
int cal_delta(int interpoint[], int, int, int);
float comb(float, float);

void decoding(void);





void main()
{
	int i, u, m, num, value;
	float start, finish;
	unsigned long int j, v;
	unsigned int mask = 1;
	long int error, ferror;
	double progress;
	double channelError_count, successError_count;
	double addNum_count, mulNum_count, totalNum_count;

	FILE *fp;
	if ((OpenFile) == NULL)
	{
		printf("Can't open the data.txt\n");
		exit(0);
	}
	fclose(fp);


	findpoints();
	/*	//***debug******
	printf("findpoint function:\n");
	for(i=0;i<n;i++)
	printf("affine point[%d]=(%d,%d)\n",i,point[0][i],point[1][i]);
	*/	//***************
	tgorder();
	/*	//****debug*****
	printf("\n\ntgorder of polebasis:\n");
	for(i=0;i<n;i++)
	printf("fai_%d=x^%d*y^%d\n",i,tg_order[i][0],tg_order[i][1]);
	*/	//**************
	mono_table();

	generator();
	//****debug******
	/*	printf("\n\ngenerator:\n");
	for(i=0;i<n;i++)
	printf("\tp%d",i);
	for(u=0;u<k;u++)
	{
	printf("\nfai%d",u);
	for(i=0;i<n;i++)
	printf("\t%d",gmatrix[u][i]);
	}
	printf("\n\n");
	*/	//***************


	srand(1977);

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

	for (SNR = start; SNR <= finish; SNR = SNR + interval)
	{
		N0 = (1.0 / ((float)k / (float)n)) / pow(10.0, SNR / 10.0);
		sgm = sqrt(N0 / 2);

		error = 0;
		ferror = 0;

		channelError_count = 0.0;
		successError_count = 0.0;
		v = 0;

		DecSucc_SeqNum = 0;
		ChosenWrong_SeqNum = 0;
		CWR = 0.0;

		addNum_count = 0.0;
		mulNum_count = 0.0;
		totalNum_count = 0.0;

		flag_addNum = 0;
		flag_mulNum = 0;

		for (j = 1; j <= seq_num; j++)
		{
			addNum = 0;
			mulNum = 0;
			//*****debug*****
			seq_num_Now = j;
			//***************

#ifndef cheatingEncoding
			//generate binary input sequence
			for (u = 0; u<k*p; u++)	//k*4
				bi_message[u] = rand() % 2;

			//convert to decimal input sequence
			for (i = 0; i<k; i++)
			{
				num = 1;
				message[i] = 0;
				for (u = 0; u<p; u++)
				{
					message[i] = message[i] + (bi_message[p*i + u] * num);
					num = num * 2;
				}
			}

			//****debug*****
			//			message[0]=0;	message[1]=0;	message[2]=0;	message[3]=0;
			//**************

			encoder(message, codeword);

			//convert to binary 
			for (u = 0; u<n*p; u++)	//
				bi_codeword[u] = 0;

			for (u = 0; u<n; u++)	//n
			{
				value = codeword[u];
				mask = 1;
				for (m = 0; m<p; m++)
				{
					if ((value & mask)>0)
						bi_codeword[p*u + m] = 1;
					else
						bi_codeword[p*u + m] = 0;
					mask = mask << 1;
				}
			}
#endif
			//modulation
			modulation();

			//channel
			channel();

			//demodulation
			//demodulation();	//这是KV版的demodulation
			demodulation_GS();	//这是GS版的

			//decoding by optimal bound
			decoding();

			//frame error rate calculation
			int temp = ferror;
			for (u = 0; u<n; u++)
				if (dec_codeword[u] != codeword[u])
				{
					ferror++;
					break;
				}

			if (temp == ferror)
				v++;



			//channelError_count = channelError_count + (epcount1-channelError_count)/(double)(j);
			addNum_count = addNum_count + (addNum - addNum_count) / (double)(j);
			mulNum_count = mulNum_count + (mulNum - mulNum_count) / (double)(j);
			totalNum_count = addNum_count + mulNum_count;

			progress = (double)(j * 100) / (double)seq_num;
			BER = (double)(error) / (double)(n*p*j);
			FER = (double)(ferror) / (double)(j);

			CWR = (double)(ChosenWrong_SeqNum) / (double)(DecSucc_SeqNum);

			printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f, Choice Errors=%2.1d, CWR=%E\r", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count, ChosenWrong_SeqNum, CWR);

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

		printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f, Choice Errors=%2.1d, CWR=%E\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count, ChosenWrong_SeqNum, CWR);

		OpenFile;
		fprintf(fp, "Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f, Choice Errors=%2.1d, CWR=%E\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count, ChosenWrong_SeqNum, CWR);
		fclose(fp);


		flag_addNum = 0;
		flag_mulNum = 0;

	}

}


void modulation(void)
{
	int i;

	//BPSK
	for (i = 0; i<n*p; i++)
	{
		tx_symbol[i][0] = -2 * bi_codeword[i] + 1;
		tx_symbol[i][1] = 0;
	}
}

void channel(void)
{
	int i, j;
	float u, r, g;

	//add noise

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

void demodulation_GS(void)
{
	int i;
	float d1, d2;
	for (i = 0; i<n*p; i++)
	{
		d1 = (rx_symbol[i][0] - 1)*(rx_symbol[i][0] - 1) + rx_symbol[i][1] * rx_symbol[i][1];   // a*a= the square of a 
		d2 = (rx_symbol[i][0] + 1)*(rx_symbol[i][0] + 1) + rx_symbol[i][1] * rx_symbol[i][1];
		if (d1<d2)
			bi_rxcodeword[i] = 0;
		else
			bi_rxcodeword[i] = 1;

	}

	//Convert rx_bicodeword to nonbinary's rx_codeword
	for (i = 0; i < n; i++)
		rxcodeword[i] = bi_rxcodeword[p*i] + 2 * bi_rxcodeword[p*i + 1] + 4 * bi_rxcodeword[p*i + 2] + 8 * bi_rxcodeword[p*i + 3];
}

void findpoints()
{
	int i, j, u, a1, a2, a3, x, y;

	//Initialisation
	for (i = 0; i<2; i++)
		for (j = 0; j<n; j++)	//w^3
			point[i][j] = 0;

	//find points over x^3-y^2-y=0
	u = 0;
	for (i = 0; i<q; i++)	//number of elements in GF(q)
	{
		x = root[i];
		for (j = 0; j<q; j++)	//number of elements in GF(q)
		{
			y = root[j];
			a1 = power(x, w + 1);
			a2 = power(y, w);
			a3 = y;
			if (add(a3, add(a1, a2)) == 0)
			{
				point[0][u] = x;
				point[1][u] = y;
				u++;
			}
		}
	}

}

void tgorder()
{
	int i, j, u, index_temp;

	for (j = 0; j<tg_size; j++)
		for (i = 0; i<2; i++)
			tg_order[j][i] = -1;

	//judge the index's scale of coresponding tgsize
	/*pesudo code
	if(  0 <= index < w ) tgsize = (index+2) * (index+1) / 2
	else if( w <= index ) tgsize = w*(w+1)/2 + (w+1)*(index-w+1)
	*/
	index_temp = 164;	//according to the last formualtion to calculate out

	j = 0;
	for (i = 0; i <= index_temp; i++)
		for (u = i; u >= 0; u--)
			if (u <= w)
			{
				tg_order[j][0] = u;	// 0<deg_x<=w
				tg_order[j][1] = i - u;	// 0<deg_y
				j++;
			}

}

void generator()
{
	int i, j;
	for (i = 0; i<k; i++)
		for (j = 0; j<n; j++)	//n
			gmatrix[i][j] = 0;

	//generator matrix
	for (i = 0; i<k; i++)	//k
		for (j = 0; j<n; j++)	//n
			gmatrix[i][j] = mul(power(point[0][j], tg_order[i][0]), power(point[1][j], tg_order[i][1]));

}

void mono_table()
{
	int i, j, u, v, z, weight[weight_Zsize][weight_XYsize], mono_order_1[monoTable_Zsize][weight_XYsize];	//26=(100/weight(z))+1, 100=pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][max(deg_y)+1][w+1]

	//Initialisation
	for (i = 0; i<monoTable_Zsize; i++)	//(230/weight(z))+1
	{
		for (j = 0; j<weight_XYsize; j++)	//pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][y_size][w+1].
		{
			weight[i][j] = -1;
			mono_order_1[i][j] = -1;
		}
	}

	for (i = 0; i<monoTable_Zsize; i++)	// rs
		for (j = 0; j<monoTable_Ysize; j++)	//max(deg_y)+1
			for (u = 0; u<monoTable_Xsize; u++)	//w+1
				mono_order[i][j][u] = -1;




	//weight of monomial over function x^3-y^2-y=0
	for (i = 0; i<weight_Zsize; i++)
	{
		if (i == 0)
			for (j = 0; j<weight_XYsize; j++)
				weight[i][j] = tg_order[j][0] * w + tg_order[j][1] * (w + 1);
		else if (i>0)
			for (j = 0; j<weight_XYsize; j++)
				weight[i][j] = weiz*i + weight[0][j];
	}


	/*
	for(i=0;i<weight_Zsize;i++)
	{
	weight[i][0]=weiz*i;
	for(j=1;j<weight_XYsize;j++)
	weight[i][j]=weiz*i+j+1;
	}
	*/
	//construction 2-dimensional mono table
	z = 0;
	for (v = 0; v<mono_ordinSize; v++)	//for each possible weight until weight(term[rs-1][max(deg_y)+1][w-1])+1, note term[rs-1][max(deg_y)+1][w-1] is the term next to term[rs-1][max(deg_y)][w]
	{
		for (i = 0; i<monoTable_Zsize; i++)	//230/weight(z)
		{
			for (j = 0; j<weight_XYsize; j++)	//230, >the index of the largest term with deg(z)=0, and weight=weight(term[rs-1][max(deg_y)+1][w-1])
			{
				if (weight[i][j] == v)
				{
					mono_order_1[i][j] = z;
					z++;
					break;
				}
			}
		}
	}

	// 2-dimensional table transform into 3-dimensional table
	for (u = 0; u<monoTable_Zsize; u++)	//rs
		for (z = 0; z<monoTable_totalSize; z++)	//>index of term[max(deg_y)][w] in the pole basis + 1
			mono_order[u][tg_order[z][1]][tg_order[z][0]] = mono_order_1[u][z];

#ifdef printfMonotable
	//**********debug********
	for (j = 0; j<weight_XYsize; j++)
		printf("\t%d", j);
	printf("\n");
	for (i = 0; i<weight_Zsize; i++)
	{
		printf("\n");
		for (j = 0; j<weight_XYsize; j++)
			printf("\t%d", weight[i][j]);
	}

	printf("\n\n");
	for (i = 0; i<monoTable_Zsize; i++)
	{
		printf("\n");
		for (j = 0; j<weight_XYsize; j++)
			printf("\t%d", mono_order_1[i][j]);
	}

	printf("\nMonmial basis is:\n");
	for (i = 0; i<monoTable_Zsize; i++)
	{
		printf("\n\nZ=%d", i);

		for (z = 0; z<monoTable_Xsize; z++)
			printf("\t%d ", z);
		printf("\n\n");
		for (j = 0; j<monoTable_Ysize; j++)
		{
			printf("\n%d", j);
			for (z = 0; z<monoTable_Xsize; z++)
				printf("\t%d ", mono_order[i][j][z]);
		}
	}

	//******************
#endif
}

void encoder(int message_temp[], int codeword_temp[])
{
	int i, j;

	for (i = 0; i<n; i++)	//n
	{
		codeword_temp[i] = 0;
		for (j = 0; j<k; j++)	//k
			codeword_temp[i] = add(codeword_temp[i], mul(message_temp[j], gmatrix[j][i]));
	}

}


void decoding(void)
{
	int i, j, value;
	int d_min, genius, GS_bound;
	int count;

	//caculate GS optimal bound
	genius = w * (w - 1) / 2;
	d_min = n - k - genius + 1;
	GS_bound = n - (int)sqrt((float)(n*(n - d_min))) - 1;


	//decoding
	count = 0;
	for (i = 0; i < n; i++)
		if (codeword[i] != rxcodeword[i])
			count++;

	if (count <= GS_bound)	//decdoing correctly
	{
		//codeword
		for (i = 0; i<n; i++)
			dec_codeword[i] = codeword[i];

		//bi_codeword
		for (i = 0; i<n*p; i++)
			dec_bicodeword[i] = bi_codeword[i];
	}
	else if (count > GS_bound)	//decoding incorrectly
	{
		for (i = 0; i < n; i++)
			dec_codeword[i] = rxcodeword[i];

	}

}








