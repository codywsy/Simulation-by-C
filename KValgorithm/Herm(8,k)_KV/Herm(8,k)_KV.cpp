﻿#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "FiniteFieldBasisGF(4).h"
#include "FiniteFieldBasisCal.h"

//#define _NoReduction_
//#define _PolyCoeffNumUncom_
//#define _PolyCoeffNumFac_
#define myWay
#define OpenFile fp=fopen("Herm(8,4)_KV_l=2.txt","a")
#define FrameError 209

#define checkInter
#define checkFac
#define cheatingDecoding

//#define printCoeffTable
//#define printDemodulation
//#define printfMonotable

//global define
#define w 2
#define weiz 4
#define lm 2	//design length
#define able_correct 14
#define pointNum 2
#define interval 1

//*****need to be modify**************
//coefficientSearch()
#define max_m 3	//design large enough to cover, dependent on the value of lm 
#define max_alpha 8	//alpha < max_m, designing alpha large enough to cover all the m_ij of Mulplicity Matrix
#define tableSize_a 20	//more than max_a, can be achieved by the hint of coefficientSearch()
#define tableSize_alpha (max_alpha+30) 	//more than max_alpha, can be achieved by the hint of coefficientSearch()
//zerobasie()
#define zb_alpha 20	//size+1, larger than tableSize_alpha, according to the maxDeg_x and maxDeg_y of PbMax, we can judge the max num Zb we need
#define zb_Ysize (zb_alpha/(w+1)+1)	//this term must be changed when alpha greater than (w+1)
#define zb_Xsize (w+1)
//cal_max_a()
#define N 500
//*****************************************
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
#define tg_size 297	//more than k, tg_size represent the probably used pole basis num, (w+1)*w/2 + (w+1) * ( (w+interpoly_Ysize)-w+1 ) 
//mono_table()
#define weight_Zsize 75
#define weight_XYsize 297	//equals to the tg_size
#define mono_ordinSize 700	//choose the value with debugging
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 100	//large than interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize 297 //>index of term[max(deg_y)][w] in the pole basis + 1, nearly equals to tg_size
//*******************

//*****need to be modify_third**************
//interpolation()
#define init_polyNum (w*(lm+1))	//change with diff lm, the poly num of the init polyGroup
#define interpoly_Zsize (lm+1)	//maxValue of z is lm=1, so the Zsize should add 1 more.
#define interpoly_Ysize 52	//maxdeg of y is (w-1) + w*(iterNum_max/(w+1)), so the Ysize should add 1 more.
#define interpoly_Xsize	(w+1)	//maxdeg of x is w, so the Xsize should add 1 more.
//factorization()
#define facpoly_Zsize (lm+1)// same as the interpoly_Zsize
#define facpoly_Ysize 52	// same as the interpoly_Ysize
#define facpoly_Xsize (w+1)	// same as the interpoly_Xsize
//rcs()
#define rcspoly_Ysize 60	//degY+1, degY = interpoly_Ysize + max(j_1) + w*( (max(i_1)+w)/(w+1) ), rcspoly_Ysize = (expoly_Ysize-1) + (facpoly_Ysize-1) + 1
#define rcspoly_Xsize 8	    //degX+1, degX = max(i_1) + w, rcspoly_Xsize = (expoly_Xsize-1) + (facpoly_Xsize-1) + 1
//#define rcspoly_Ysize 67	//degY+1, degY = interpoly_Ysize + max(j_1) + w*( (max(i_1)+w)/(w+1) )
//#define rcspoly_Xsize 9	//degX+1, degX = max(i_1) + w
//expoly()-->expanded polynomial
#define expoly_Ysize (lm*faiMax_Ysize+1)	
#define expoly_Xsize (lm*faiMax_Xsize+1)
//polebasis
#define faiMax_Ysize 7	//7change with diff w and k, the max degY of probably used polebasis
#define faiMax_Xsize 4	//4change with diff w and k, the max degX of probably used polebasis


//******varible*********

//main()
unsigned long int seq_num;	//number of input binary sequences
float SNR,N0;
double BER, FER;
float pi=3.141593;	// pai
int iterNum=0; //iteration num
int bi_message[k*p], message[k];	//transmitted messge[k]
int codeword[n], bi_codeword[n*p]; //codewords
float tx_symbol[p*n][2], rx_symbol[p*n][2], sgm;
float RM[q][n];	//Reliability Matrix
int MM[q][n];	//Multiplicity Matrix
//polebasis()

//zerobasis
int zb[zb_alpha][zb_Ysize][zb_Xsize];
//coefficientSearch
int coeff_table[n][tableSize_alpha][tableSize_a];	//n points
//findpoint()
int point[2][n];	//rational points[2][w^3], [0]->x, [1]->y
//tgorder()
int tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
//mono_table()
int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
int gmatrix[k][n];	//generator matrix[k][n]
int dec_message[k],dec_codeword[n],dec_bicodeword[n*p]; //decoding result
//MatrixConvert()
int s=0;	//large enough to make the circle of MatrixConvert() keep moving
//interpolation
int bi_rxword[p*n], rxword[n]; //received wordss
int epcount1, epcount2, testCount1, testCount2, testCount1_com;
int Q_interpoly[init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];	//chosen interpolated polynomial [z-deg+1][max(deg_y)+1][w+1]
int uu;	//factorisation step index
int l, listNum;	//candidate output index
//factorization() and rcs()
int Q[k][facpoly_Zsize][facpoly_Ysize][facpoly_Xsize];	//sequtial deduction polynomial [number of fac steps=k][rs][y_size][w+1], y_size> maxdeg_y]+rs*(deg_¦µ(k-1))
int rootlist[k][lm+1];	//the list of roots [number of fac steps=k][expected number of roots in each solution, 5>rs]
int output[lm+1][k];	//the list of candidate message [expeced number of candidate message, >rs][length of message, k]
int expoly[2][expoly_Ysize][expoly_Xsize];	//expanded polynomial in [z+f_k-1-u*pb_k-1-u]^rs, expoly[rs][3>(max(deg_y) in encoding functions)*(rs-1)][3>(max(deg_x) in encoding functions)*(rs-1)]


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
	int flag_addNum=0, flag_mulNum=0;

	//*******************


//**********function**********
void findpoints(void);
void mono_table(void);
//void polyexp(int, int, int);
void polyexp1(int, int, int, int, int poly[][expoly_Ysize][expoly_Xsize]);
//void zerobasis(void);
void polebasis(void);
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
float PDF(int,int);
//*****************
void MatrixConvert(void);
int cal_max_a();
int gama(int);
void zerobasis(int, int LM_zb[][zb_alpha]);
void coefficientSearch();
//****************
int cal_delta(int interpoint[], int, int, int);
float comb(float, float);

#ifdef cheatingDecoding
void decoding(void);
#endif



void main()
{
	int i, u, m, num, value;	
	float start, finish;
	unsigned long int j,v;
	unsigned int mask=1;
	long int error,ferror;
	double progress;
	double channelError_count, successError_count;
	double addNum_count, mulNum_count, totalNum_count;

	FILE *fp;
	if((OpenFile)==NULL)
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

//**the part of coefficient search
	//initilaization
	for(i=0;i<n;i++)
		for(u=0;u<tableSize_alpha;u++)
			for(m=0;m<tableSize_a;m++)
			{
				coeff_table[i][u][m] = -1;
			}

	coefficientSearch();

//********************************


	srand(1977);

	//*****input data from basic_input.txt*********
	FILE * fp_input = fopen("basic_input.txt","r");
	if(fp_input!=NULL) {
		fscanf(fp_input,"start SNR:%f\n", &start);
		fscanf(fp_input,"finish SNR:%f\n", &finish);
		fscanf(fp_input,"seq_num:%ld", &seq_num);

		fclose(fp_input);
	}else{
		printf("\n\ncan't read the basic_input file\n\n");
	}
	
	printf("start SNR: %0.2f\n\n", start);
	printf("finish SNR: %0.2f\n\n", finish);
	printf("seq_num: %d\n", seq_num);
	//*******************************

	for(SNR=start; SNR<=finish; SNR=SNR+interval)
	{
		N0=(1.0/((float)k/(float)n))/pow(10.0, SNR/10.0);
		sgm=sqrt(N0/2);
		
		error=0;
		ferror=0;

		channelError_count=0.0;
		successError_count=0.0;
		v=0;

		DecSucc_SeqNum=0;
		ChosenWrong_SeqNum=0;
		CWR=0.0;

		addNum_count=0.0;
		mulNum_count=0.0;
		totalNum_count=0.0;

		flag_addNum=1; 
		flag_mulNum=1;

		for(j=1;j<=seq_num;j++)
		{
			addNum=0;
			mulNum=0;
			//*****debug*****
			seq_num_Now=j;
			//***************

			//generate binary input sequence
			for(u=0;u<k*p;u++)	//k*4
				bi_message[u]=rand()%2;

			//convert to decimal input sequence
			for(i=0;i<k;i++)
			{
				num=1;
				message[i]=0;
				for(u=0;u<p;u++)
				{
					message[i]=message[i]+(bi_message[p*i+u]*num);
					num=num*2;
				}
			}

			//****debug*****
//			message[0]=0;	message[1]=0;	message[2]=0;	message[3]=0;
			//**************

			encoder(message,codeword);

			//convert to binary 
			for(u=0;u<n*p;u++)	//n*4
				bi_codeword[u]=0;

			for(u=0;u<n;u++)	//n
			{
				value=codeword[u];
				mask=1;
				for(m=0;m<p;m++)
				{
					if((value & mask)>0)
						bi_codeword[p*u+m]=1;
					else
						bi_codeword[p*u+m]=0;
					mask=mask<<1;
				}
			}

			//modulation
			modulation();

			//channel
			channel();

			//demodulation
			demodulation();

			//MatrixConvert
			iterNum = 0;
			MatrixConvert();

			//LIST DECODER
			//test vector construction
//			test_vec_contruction();

			//interpolation
			interpolation();

#ifndef cheatingDecoding
			//factorisation
			factorisation();

			//choose
			choose();
#endif
			
#ifdef cheatingDecoding
			//decoding by cheat
			decoding();
#endif
			//bit error rate calculation
			int temp=error;
			for(u=0;u<n*p;u++)	//n*4
			{
				if(dec_bicodeword[u]!=bi_codeword[u])
					error++;
			}

			//frame error rate calculation
			if( error>temp )
				ferror++;
			else if(error==temp) //decoding success
			{
				//calculate the aver error of decoding success situation
				v++;
			//	successError_count = successError_count + (epcount1-successError_count)/(double)(v);
			}

			//channelError_count = channelError_count + (epcount1-channelError_count)/(double)(j);
			addNum_count = addNum_count + (addNum-addNum_count)/(double)(j);
			mulNum_count = mulNum_count + (mulNum-mulNum_count)/(double)(j);
			totalNum_count = addNum_count + mulNum_count;

			progress=(double)(j*100)/(double)seq_num;
			BER=(double)(error)/(double)(n*p*j);
			FER=(double)(ferror)/(double)(j);
			
			CWR=(double)(ChosenWrong_SeqNum)/(double)(DecSucc_SeqNum);
			
			printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f, Choice Errors=%2.1d, CWR=%E\r", progress, SNR, error, BER, ferror, FER,	addNum_count, mulNum_count, totalNum_count, ChosenWrong_SeqNum, CWR);

			if(ferror>FrameError)
				break;
		}

		if(ferror>FrameError)
		{
			BER=(double)error/(double)(n*p*j);
			FER=(double)(ferror)/(double)(j);
		}
		else
		{
			BER=(double)error/(double)(n*p*seq_num);
			FER=(double)(ferror)/(double)(seq_num);
		}

		printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f, Choice Errors=%2.1d, CWR=%E\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count, ChosenWrong_SeqNum, CWR);

		OpenFile;
		fprintf(fp,"Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f, Choice Errors=%2.1d, CWR=%E\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count, ChosenWrong_SeqNum, CWR);
		fclose(fp);


		flag_addNum=0; 
		flag_mulNum=0;

	}

}


void modulation(void)
{
	int i;

	//BPSK
	for(i=0;i<n*p;i++)	
	{
		tx_symbol[i][0]=-2*bi_codeword[i]+1;
		tx_symbol[i][1]=0;
	}
}

void channel(void)
{
	int i, j;
	float u, r, g;

	//add noise

	for(i=0;i<n*p;i++)	
	{
		for(j=0;j<2;j++)	//Inphase + Quadrature
		{
			u=(float)rand()/(float)RAND_MAX;
			if(u==1.0)
				u=0.999999;
			r=sgm*sqrt(2.0*log(1.0/(1.0-u)));

			u=(float)rand()/(float)RAND_MAX;
			if(u==1.0)
				u=0.999999;
			g=(float)r*cos(2*pi*u);

			rx_symbol[i][j]=tx_symbol[i][j]+g;
		}
	}
}

void demodulation(void)
{
	int i,j,u,v;
	float sum,proba[pointNum];
	float Pr[p][2];

	//initializ
	for(u=0;u<q;u++)
		for(v=0;v<n;v++)
			RM[u][v]=0;

	//cal PDF
	for(i=0;i<n;i++)
	{	
		for(j=0;j<p;j++)
		{
			sum=0;
			for(u=0;u<pointNum;u++)
			{	
				proba[u]=PDF(i*p+j,u);
				sum+=proba[u];
			}

			Pr[j][0]=proba[0]/sum;	//proba-->00	0
			Pr[j][1]=proba[1]/sum;	//proba-->01	1

		}

		j=0;
//		for(int z=0;z<pointNum;z++)
		for(u=0;u<pointNum;u++)
		for(v=0;v<pointNum;v++)
//		for(int x=0;x<pointNum;x++)
//		for(int y=0;y<pointNum;y++)
//		for(int h=0;h<pointNum;h++)
		{
			RM[j][i]=(Pr[0][v]*Pr[1][u]);  
			j++;
		}
	
	}

#ifdef printDemodulation
	//*****debug**********
	printf("\n\ncodeword");
	for(v=0;v<n;v++)
		printf("\t%d",codeword[v]);

/*	printf("\n\nrxword\t");
	for(v=0;v<n;v++)
		printf("\t%d",rxword[v]);
	
	printf("\n\n");
	for(v=0;v<n*p/4;v++)
		printf("\t(%f,%f,%f,%f)",rx_symbol[v*2+0][0],rx_symbol[v*2+0][1],rx_symbol[v*2+1][0],rx_symbol[v*2+1][0]);
*/
	printf("\n\n");
	for(u=0;u<q;u++)
	{
		printf("\n");
		printf("\t%d",root[u]);
		for(v=0;v<n;v++)
			printf("\t%f",RM[u][v]);
	}
	printf("\n\n");
	//*******************
#endif

}

float PDF(int i,int u)
{
	float temp1,temp2,temp3;
	float constell_point[2][2]={ {1,0},{-1,0} }; // S1{1,0}-->0,S2{-1,0}-->1

	temp3= 1.0/(pi*N0);
	temp1= pow( (rx_symbol[i][0]-constell_point[u][0]) ,2) + pow( (rx_symbol[i][1]-constell_point[u][1]) ,2);
	temp2= temp3*exp(-temp1/N0);  // 0.5 is the normalization coefficient
	
	return temp2;
}

void MatrixConvert( )
{
	int i, j, u, v, index_i, index_j, lM, flag_iterNum;
	unsigned long int degree_iterNum, iterNum_temp;
	float RM_temp[q][n], tempRM, max_temp;
	int MulpMatrix[q][n], tempMM;

	//initialization
	for(i=0;i<n;i++)
		for(j=0;j<q;j++)
		{
			RM_temp[j][i] = RM[j][i];
			MulpMatrix[j][i] = 0; 
			MM[j][i] = 0;
		}

	s = 3000;	//important part, must be check
	lM = 0;

	//start
	while( s>0 && lM<=lm )
	{
		max_temp = 0.0;
		index_j = -1;
		index_i = -1;


		//find the maximal entry factor in RM
		for(i=0;i<n;i++)
			for(j=0;j<q;j++)
			{
				if( max_temp < RM_temp[j][i] )
				{
//					sumNum1 += 4;
					max_temp = RM_temp[j][i];
					index_j = j;
					index_i = i;

#ifdef myWay
//					sumNum2 += 5;
					//because the sum of a cloumn of RM is equal to 1. so 
					//for a cloumn, if there is a factor larger than 0.5, that must
					//be the largest one of the cloumn
					if( max_temp > 0.5 )
					{
//						sumNum2 += 4;
						max_temp = RM_temp[j][i];
						index_j = j;
						index_i = i;
						break; //jump out the cloumn circle
					}
#endif
				}
			}

		//update
		RM_temp[index_j][index_i] = RM_temp[index_j][index_i] / (MulpMatrix[index_j][index_i]+2);
		MulpMatrix[index_j][index_i] += 1;
		s -= 1; 

/*	
	fp=fopen("TestResult.txt","a");
	fprintf(fp, "\nindex_j = %d\tindex_i = %d\ts=%d", index_j, index_i, s);

	fprintf(fp, "\n\nMM:\n");
	for(j=0;j<q;j++)
	{
		for(i=0;i<n;i++)
			fprintf(fp, "%d\t",MulpMatrix[j][i]);
		fprintf(fp,"\n");
	}

	fprintf(fp,"\n");
//	printf("\n\nsumNum1=%d\tsumNum2=%d", sumNum1, sumNum2);
	fclose(fp);
*/
		//cal CM
		iterNum_temp = 0;
		for(j=0;j<q;j++)
			for(i=0;i<n;i++)
			{
				iterNum_temp +=( MulpMatrix[j][i] * (MulpMatrix[j][i]+1));
			}
		iterNum_temp = iterNum_temp/2;

		//cal the deg(1,wz) corresponding to CM
		degree_iterNum = 0;
		flag_iterNum = 1;	//make the search part more effciency
		for(u=0; u<(lm+1) && flag_iterNum ;u++)
			for(j=0; j<monoTable_Ysize && flag_iterNum ;j++)
				for(i=0; i<monoTable_Xsize && flag_iterNum ;i++)
					if ( mono_order[u][j][i] == iterNum_temp )
					{
						degree_iterNum = i*w + j*(w+1) + u*weiz;
						flag_iterNum = 0;
					}

		//cal lM
		lM = degree_iterNum/weiz;
		
	}

	//transform to the MM[q][n]
	for(j=0;j<q;j++)
		for(i=0;i<n;i++)
		{
			MM[j][i] = MulpMatrix[j][i];
		}

	iterNum = iterNum_temp;	//CM

#ifdef checkFac
	//cal degree_iterNum
	Deg_iterNum = 10000;
	Deg_iterNum = degree_iterNum;

	//cal codewordScore
	codewordScore = 0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<q;j++)
			if( codeword[i]==root[j] )
				{
					codewordScore += MulpMatrix[j][i];
				}
	}
#endif



	
}

void interpolation()
{
	int i, j, u, v, z, temp, temp1, temp2, lod_temp, index_temp, degree_temp[init_polyNum];
	int lod[init_polyNum], delta[init_polyNum], J[init_polyNum], act[init_polyNum];
	int lod_min, j_min, mulpi, alpha, beta, xi;
	int interpoint[q*n][3], interpoint_num;	//(pointIndex,ri,mij),the maximum num of interpoints is n*q
	int f[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];	//temp for poly[j_min]
	int g1[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize+1], g2[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	
	//initiliazation
	for(j=0;j<q*n;j++)
		for(i=0;i<3;i++)
		{
			interpoint[j][i] = 0;
		}

	for(j=0;j<init_polyNum;j++)
		for(i=0;i<interpoly_Zsize;i++)
			for(u=0;u<interpoly_Ysize;u++)
				for(v=0;v<interpoly_Xsize;v++)
				{
					Q_interpoly[j][i][u][v] = 0;
				}

			
	//start to interpolation
	//set the interpoint
	interpoint_num = 0;
	for(i=0;i<n;i++)
		for(j=0;j<q;j++)
			if( MM[j][i] != 0 )
			{
				interpoint[interpoint_num][0] = i;	//pointIndex
				interpoint[interpoint_num][1] = root[j];	//ri
				interpoint[interpoint_num][2] = MM[j][i];	//mij
				interpoint_num ++;
			}

/*	//****debug***********
	interpoint_num = 0;
	
	for(i=0;i<n;i++)
		for(j=0;j<q;j++)
			if( MM[j][i] != 0 && i!=1 )
			{
				interpoint[interpoint_num][0] = i;	//pointIndex
				interpoint[interpoint_num][1] = root[j];	//ri
				interpoint[interpoint_num][2] = MM[j][i];	//mij
				interpoint_num ++;
			}
	
	interpoint[interpoint_num][0] = 1;
	interpoint[interpoint_num][1] = root[3];
	interpoint[interpoint_num][2] = 2;
	interpoint_num++;
*/	//************************

	//set Group initialization
	for(i=0;i<(lm+1);i++)	//rs
		for(j=0;j<w;j++)	//w
			Q_interpoly[w*i+j][i][j][0]=1;	//j+w*i

	//interpolation main progress
	for(i=0;i<interpoint_num;i++)
	{
		mulpi = interpoint[i][2];	//multiplicity
		for(alpha=0;alpha<mulpi;alpha++)
			for(beta=0;beta<(mulpi-alpha);beta++)
			{
				//Calculate each polynomial's leading order
				for(j=0;j<init_polyNum;j++)
				{
					lod_temp=0;
					lod[j]=0;
					for(u=0;u<interpoly_Zsize;u++)
						for(v=0;v<interpoly_Ysize;v++)
							for(z=0;z<interpoly_Xsize;z++)
								if( Q_interpoly[j][u][v][z] != 0 )
								{
									lod_temp = mono_order[u][v][z];
									if(lod_temp>lod[j])
									{
										lod[j]=lod_temp;
									}
		
								}
				}

#ifndef _NoReduction_
				//Initialise the eliminator array act[num_poly]
				for(j=0;j<init_polyNum;j++)	//num_poly
				{
					if(lod[j]<=iterNum)	//C=n when multiplicity = 1
						act[j]=1;
					else
						act[j]=0;
				}
#else		
				for(j=0;j<init_polyNum;j++)	//num_poly
					act[j] = 1;		
#endif		

				//Calculate the hasse derivative mapping of each polynomials
				j_min = -1;
				for(j=0;j<init_polyNum;j++)	//wrt each polynomial
				{
					J[j] = 0;
					delta[j] = 0;
					if( act[j]==1 )	//for polynomials with leading order less of equal to C
					{
						delta[j] = cal_delta(interpoint[i],j,alpha,beta);	//(interpoint, polyIndex, alpha, beta)
					}

					if(delta[j]!=0)
					{
						J[j]=1;	//record those polynomial with a nonzero hasse derivative mapping
						lod_min = lod[j];
						j_min = j;
					}
				}

				//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
				for(j=0;j<init_polyNum;j++)	//num_polys
				{
					if(J[j]==1 && lod[j]<lod_min)
					{
						lod_min=lod[j];
						j_min=j;
					}
				}
					//printf("\nj_min=%d\n", j_min);
				
				//update poly group
				if(j_min!=-1)
				{
					//f = Q_interpoly[j_min]
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
							{
								f[u][v][z] = Q_interpoly[j_min][u][v][z];
							}

					//Modify nonzero polynomials
					for(j=0;j<init_polyNum;j++)
					{
						if(J[j]==1)
						{
							if( j!=j_min )
							{
								//delta*g_k+delta_k*f
								for(u=0;u<interpoly_Zsize;u++)	//rs
									for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
										for(z=0;z<interpoly_Xsize;z++)	//w+1
										{
											if(	Q_interpoly[j][u][v][z]!=0 )
											{
												temp1 = mul(delta[j_min],Q_interpoly[j][u][v][z]);
											}
											else
												temp1 = 0;
		
											if(	f[u][v][z]!=0 )
											{
												temp2 = mul(delta[j],f[u][v][z]);
											}
											else
												temp2 = 0;
											
											if( temp1!=0 || temp2!=0 )
											{
												Q_interpoly[j][u][v][z]=add(temp1,temp2);	
											}
											else
												Q_interpoly[j][u][v][z] = 0;
										}
							}
							else if(j==j_min)
							{
								for(u=0;u<interpoly_Zsize;u++)	//rs
								{
									for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
									{
										for(z=0;z<(interpoly_Xsize+1);z++)	//w+2
										{
											g1[u][v][z]=0;
										}

										for(z=0;z<interpoly_Xsize;z++)	//w+1
										{
											g2[u][v][z]=0;
										}
									}
								}

								//g1=x*f
								for(u=0;u<interpoly_Zsize;u++)	//rs
									for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
										for(z=0;z<interpoly_Xsize;z++)	//w+1
											if(Q_interpoly[j][u][v][z]!=0)
											{
												g1[u][v][z+1]=Q_interpoly[j][u][v][z];
											}

								//convert x^w+1=y^w+y, difference with diff code
								for(u=0;u<interpoly_Zsize;u++)	//rs
								{
									for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
									{
										if(g1[u][v][w+1]!=0)	//refering to two finitefield add 
										{
											g1[u][v+1][0]=add(g1[u][v+1][0],g1[u][v][w+1]);
											g1[u][v+w][0]=add(g1[u][v+w][0],g1[u][v][w+1]);
											g1[u][v][w+1]=0;
										}
									}
								}

								//g2=xi*f
								index_temp = interpoint[i][0];
								xi = point[0][index_temp];
								for(u=0;u<interpoly_Zsize;u++)	//rs
									for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
										for(z=0;z<interpoly_Xsize;z++)	//w+1
											if(Q_interpoly[j][u][v][z]!=0)
											{
												g2[u][v][z] = mul( xi,Q_interpoly[j][u][v][z] );
											}
								//Q_interpoly=g1+g2
								for(u=0;u<interpoly_Zsize;u++)	//rs
									for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
										for(z=0;z<interpoly_Xsize;z++)	//w+1	
											if( g1[u][v][z]!=0 || g2[u][v][z]!=0 )
											{
												Q_interpoly[j][u][v][z] = add( g1[u][v][z],g2[u][v][z] );
											}
											else
												Q_interpoly[j][u][v][z] = 0;								

							}
						}
					}
				}

			}

	}


	//find out the poly for factorization
	//calculate the lod of poly
	for(j=0;j<init_polyNum;j++)
	{
		temp = -1;
		degree_temp[j] = 0;
	
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					if( Q_interpoly[j][u][v][z]!=0 )
						if( temp < mono_order[u][v][z] )
						{
							temp = mono_order[u][v][z];
						}

		degree_temp[j] = temp ;
	}

	//choose the min dgree
	index_temp = 0;
	temp = degree_temp[0];
	for(j=1;j<init_polyNum;j++)
		if( (temp>degree_temp[j]) && (degree_temp[j]<=iterNum) )
		{
			temp = degree_temp[j];
			index_temp = j;
		}

	//initialization poly for factorization
	for(j=0;j<k;j++)	//number of fac steps=k
		for(u=0;u<facpoly_Zsize;u++)	//rs
			for(v=0;v<facpoly_Ysize;v++)	//y_size
				for(z=0;z<facpoly_Xsize;z++)	//w+1
					Q[j][u][v][z]=0;

	for(u=0;u<interpoly_Zsize;u++)	//rs
		for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
			for(z=0;z<interpoly_Xsize;z++)	//w+1
				Q[0][u][v][z]=Q_interpoly[index_temp][u][v][z];

#ifdef checkFac
	//cal the deg(1,wz) corresponding to Q_min
	Deg_poly = -1;
	int	flag_iterNum = 1;	//make the search part more effciency
		for(u=0; u<(lm+1) && flag_iterNum ;u++)
			for(j=0; j<monoTable_Ysize && flag_iterNum ;j++)
				for(i=0; i<monoTable_Xsize && flag_iterNum ;i++)
					if ( mono_order[u][j][i] == degree_temp[index_temp] )
					{
						Deg_poly = i*w + j*(w+1) + u*weiz;
						flag_iterNum = 0;
					}

	if(Deg_poly < 0)
	{
		printf("\n\nDeg_poly is error\n\n");
	}
#endif

#ifdef checkInter
	int delta_temp;

	flag_addNum=0; 
	flag_mulNum=0;

	j = index_temp;
	for(i=0;i<interpoint_num;i++)
	{
		mulpi = interpoint[i][2];	//multiplicity
		for(alpha=0;alpha<mulpi;alpha++)
			for(beta=0;beta<(mulpi-alpha);beta++)
			{	
				delta_temp = -1;
				delta_temp = cal_delta(interpoint[i],j,alpha,beta);
				
				if(delta_temp!=0)
				{
					printf("\n\nQ_interpoly in point[%d](alpha=%d,beta=%d) has interpolation error\n\n", i, alpha, beta);
				}
			}
	}

	flag_addNum=1; 
	flag_mulNum=1;

#endif




}

void factorisation()	//output: output[lm+1][k], listNum
{
	int i, j, u, v, z;


	//initialization
	for(u=0;u<k;u++)	//number of fac steps=k
		for(v=0;v<lm+1;v++)	//5>rs=expected number of roots
			rootlist[u][v]=-1;	

	for(u=0;u<lm+1;u++)	//5>(rs-1)=expected number of output lists
		for(v=0;v<k;v++)	//k
			output[u][v]=-1;

	//Initialisation of factorisation
	uu=0;	//recursive deduction index
	l=0;	//candidate output index
		
		
	#ifdef _PolyCoeffNumFac_
		int NumTemp0=0, NumTemp1=0;
		for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
			for(z=0;z<interpoly_Xsize;z++)	//w+1
				if(Q_interpoly[i][0][v][z]!=0)
						NumTemp0++;

		for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
			for(z=0;z<interpoly_Xsize;z++)	//w+1
				if(Q_interpoly[i][1][v][z]!=0)
						NumTemp1++;

		NumTemp0 += 0;
		NumTemp1 += 0;
		//*********************
	#endif

	//recursive coefficient search
	rcs(uu);

	listNum = l;

}

void rcs(int uu)
{
	int i, j, u, v, m, z, t, r, i_1, j_1, i_2, j_2, a, b, leadMono, leadMono_temp, alpha, act, temp;
	int lc[lm+1], q_temp[lm+1][rcspoly_Ysize][rcspoly_Xsize];	//q_temp[z-deg+1][y_size>max(deg_y)+1+(rs-1)*(max(deg_y) in encoding functions)][14>w+(rs-1)*w], lc[rs]--leading coefficient polynomial
	int index_flag;
	int rcspoly_Ysize_1,rcspoly_Xsize_1;	//the q_temp size of the first step of factorization

	//array size initialization
	rcspoly_Ysize_1= faiMax_Ysize + facpoly_Ysize + 1;
	rcspoly_Xsize_1= faiMax_Xsize + facpoly_Xsize + 1;	

	leadMono=0; leadMono_temp=0;	//leading monomial index
	act=0;	//judge value for recursive search of each f_k-1-u

	//find pb_k-1-uu
	j_1 = tg_order[k-1-uu][1];
	i_1 = tg_order[k-1-uu][0];

	//initialise q_temp
	for(i=0;i<lm+1;i++)	//rs
		for(j=0;j<rcspoly_Ysize;j++)	
			for(u=0;u<rcspoly_Xsize;u++)	
				q_temp[i][j][u]=0;

	//Calculate q_temp[pb_k-1-u]=q[u][pb_k-1-u]
	for(i=0;i<facpoly_Zsize;i++)	//rs
	{
		for(j=0;j<facpoly_Ysize;j++)	//y_size=max(deg_y) in encoding function*(rs-1), and 31>=max(deg_y)+1
		{
			for(u=0;u<facpoly_Xsize;u++)	//
			{
				q_temp[i][j+i*j_1][u+i*i_1] = Q[uu][i][j][u];	//this time, q_temp size is [rs+1][maxY(j_1)+maxY(Q)][maxX(i_1)+maxX(Q)]
			}
		}
	}

	//*****debug*****
	if( rcspoly_Xsize_1 < (w+1) )
		printf("\n\nrcs() has error\n\n");
	//**************

	//x^(w+1)=y^w+y
	for(i=0;i<lm+1;i++)	//rs
	{
		index_flag = rcspoly_Xsize_1-1;	//the max deg_x for q_temp
		while( index_flag>w )
		{
			temp = index_flag-(w+1);	// deg_x - (w+1)
			for(j=0;j<rcspoly_Ysize_1;j++)
				if( q_temp[i][j][index_flag]!=0 )
				{
					q_temp[i][j+1][temp] = add( q_temp[i][j+1][temp],q_temp[i][j][index_flag] );	//y^1
					q_temp[i][j+w][temp] = add( q_temp[i][j+w][temp],q_temp[i][j][index_flag] );	//y^(w+1)
					q_temp[i][j][index_flag] = 0;
				}
			index_flag = index_flag-1;
		}
	}
	//*****************
	
	//find the leading monomial in q_temp
	leadMono=-1;
	for(i=0;i<lm+1;i++)	//rs
	{
		for(j=0;j<rcspoly_Ysize;j++)	//y_size
		{
			for(u=0;u<w+1;u++)	//w+1
			{
				if(q_temp[i][j][u]!=0)
				{
					if( leadMono <= mono_order[0][j][u] )
					{
						leadMono = mono_order[0][j][u];
						j_2=j;	//record the leading monomial's y degree
						i_2=u;	//record the leading monomial's x degree
					}
				}
			}
		}
	}

	//find the leading coefficient from each polynomial lc[rs]
	for(i=0;i<lm+1;i++)	//rs
	{
		lc[i]=0;
		for(j=0;j<rcspoly_Ysize;j++)	//y_size
		{
			for(u=0;u<w+1;u++)	//w+1
			{
				if(q_temp[i][j][u]!=0 && j==j_2 && u==i_2)
					lc[i]=q_temp[i][j][u];
			}
		}
	}

	u = 0;	//root index
	act = 0;
	//find the roots of leading coefficient polynomial lc[rs]
	for(i=0;i<q;i++)	//number of elements in GF(4)
	{
		b=0;
		for(j=0;j<lm+1;j++)	//rs
		{
			a=mul(lc[j], power(root[i], j));
			b=add(a, b);
		}
		if(b==0)	//the root is found out, and set in rootlist[uu]
		{
			act=1;
			rootlist[uu][u]=root[i];
			u++;
		}
	}		

	//For each distinct root of rootlist[u]
	if(act==1)	//act==1 means there is at least one root in rootlist[u];
	{
		for(i=0;i<lm+1;i++)	//2>rs
		{
			if(rootlist[uu][i]!=-1)
			{	
				alpha=rootlist[uu][i];
				
				output[l][k-1-uu]=alpha;	//output[l][k-1-uu]
					
				if(uu==(k-1))	//when u=k-1, output candidate polynomial
				{
					for(j=0;j<k;j++)	//k
					{
						if(output[l][j]==-1)	//this judgement is used to make sure the outputList is ever not be used 
							output[l][j]=output[l-1][j];
					}

					l++;	//locate next candidate output
				}
				else  //update the q[uu+1]
				{	
					r=uu+1;

					//calculate Q[uu+1]
					for(j=0;j<lm+1;j++)	//rs
						if(j==0)
						{
							for(m=0;m<facpoly_Ysize;m++)	//y_size
								for(z=0;z<facpoly_Xsize;z++)	//w+1
									Q[r][j][m][z] = Q[uu][j][m][z];

							polyexp1(alpha, i_1, j_1, j, expoly);

						}
						else if(j>0)
						{
							//Initialise q_temp
							for(u=0;u<lm+1;u++)	//rs
								for(m=0;m<facpoly_Ysize;m++)	//y_size
									for(z=0;z<facpoly_Xsize;z++)	//w+1
										q_temp[u][m][z]=0;	

							//calculate (z+f_k-1-u*pb_k-1-u)^j
							polyexp1(alpha, i_1, j_1, j, expoly);
							
							//calculate q_temp
							for(u=0;u<lm+1;u++)	//rs
							{
								for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
								{
									for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
									{
										if(expoly[u][m][z]!=0)
										{
											for(v=0;v<facpoly_Ysize;v++)	//y_size
											{
												for(t=0;t<facpoly_Xsize;t++)	//w+1
												{
													if(Q[uu][j][v][t]!=0)
													{
														temp = mul(expoly[u][m][z], Q[uu][j][v][t]);

														//if(q_temp[u][v+m][t+z]!=0)
														//	printf("\nq_temp[u][v+m][t+z]!=0\n");

														q_temp[u][v+m][t+z] = temp; 													}
													// caution overflow problem of q_temp
													// here q_temp, ysize=exploy_Ysize+facpoly_Ysize-2+1=9, xsize=expoly_Xsize+facpoly_Xsize-2+1=6
												}
											}
										}
									}
								}

								//convert x^w+1=y^w+y, difference with diff code
								index_flag = rcspoly_Xsize-1;	//the max deg_x for q_temp, rcspoly_Xsize = (expoly_Xsize-1) + (facpoly_Xsize-1) + 1
								while( index_flag>w )
								{
									temp = index_flag-(w+1);	// deg_x - (w+1)
									for(m=0;m<rcspoly_Ysize;m++)
										if( q_temp[u][m][index_flag]!=0 )
										{
											q_temp[u][m+1][temp] = add( q_temp[u][m+1][temp],q_temp[u][m][index_flag] );	//y^1
											q_temp[u][m+w][temp] = add( q_temp[u][m+w][temp],q_temp[u][m][index_flag] );	//y^(w+1)
											q_temp[u][m][index_flag] = 0;
										}
									index_flag = index_flag-1;
								}

							}

							//Q[r] = add( Q[r],q_temp )
							for(u=0;u<lm+1;u++)	//rs
								for(m=0;m<facpoly_Ysize;m++)	//y_size
									for(z=0;z<facpoly_Xsize;z++)	//w+1
										Q[r][u][m][z]=add( Q[r][u][m][z],q_temp[u][m][z] );
						}
				
					//next coefficient searching
					rcs(r);
				}
			}
		}
	}

}

void polyexp1(int c, int i, int j, int deg_z, int poly[][expoly_Ysize][expoly_Xsize])
{
	int u, m, z;	//p1[rs+1][18>(max(deg_y) in encoding functions)*(rs-1)][10>(max(deg_x) in encoding functions)*(rs-1)], p2[rs][26=18+max(deg_y) in encoding functions][14=10+max(deg_x) in encoding functions]
	int poly_temp[lm+1+1][expoly_Ysize+faiMax_Ysize][expoly_Xsize+faiMax_Xsize];
	int temp, temp_Ysize, temp_Xsize;

	temp_Ysize = expoly_Ysize+faiMax_Ysize;
	temp_Xsize = expoly_Xsize+faiMax_Xsize;

	//start
	if(deg_z==0)
	{
		for(u=0;u<lm+1;u++)	//rs
			for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
				for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
					poly[u][m][z]=0;

		poly[0][0][0]=1;
	}
	else if(deg_z==1)	//because deg_z is equal to or less than 1, so deg_z>0 <-> deg_z==1
	{
		poly[0][0][0]=0;

		poly[1][0][0]=1;
		poly[0][j][i]=c;
	}
	else if(deg_z>1)
	{
		for(u=0;u<lm+1+1;u++)	//rs
			for(m=0;m<temp_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
				for(z=0;z<temp_Ysize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
				{
					poly_temp[u][m][z]=0;
				}

		//calculate poly_temp * poly
		for(u=0;u<lm+1;u++)	//rs
			for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
				for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
					if( poly[u][m][z]!=0 )
					{
						poly_temp[u+1][m][z] = poly[u][m][z];
					}
	
		//calculate y^j*x^i*poly
		for(u=0;u<lm+1;u++)	//rs
			for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
				for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
					if( poly[u][m][z]!=0 )
					{
						temp = mul( c,poly[u][m][z] );
						poly_temp[u][j+m][i+z] = add( poly_temp[u][j+m][i+z],temp );
					}

		for(u=0;u<lm+1;u++)	//rs
			for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
				for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
				{
					poly[u][m][z] = poly_temp[u][m][z];
				}
	}
}

void choose()
{
	int i, j, u, v, value, flag, min_index;
	unsigned int mask=1;
	float proba[lm+1], proba_temp, temp;
	int codeword_temp[n];


	//normal mode
	//Initialise hamming distance counter
	for(u=0;u<lm+1;u++)
		proba[u]=-1.0;

	proba_temp=0.0;
	flag = 0;
	min_index = -1;

	if( listNum!=0 )
	{
		for(j=0;j<listNum;j++)
		{
			//reencoding
			encoder(output[j],codeword_temp);
			//calculate the posteriori probablity
			temp = 1.0;
			for(u=0;u<n;u++)
			{
				for(v=0;v<q;v++)
					if(codeword_temp[u]==root[v])
						temp = temp*RM[v][u];
			}

			if( proba_temp <= temp )
			{
				proba_temp = temp;
				min_index = j;

				for(v=0;v<n;v++)
					dec_codeword[v] = codeword_temp[v];

				flag = 1;	//exist at less one valid solution
			}

		}

	}

	//output the decoding result
	if(flag==0)	//listNum == 0
	{
		//search dec_codeword[n] from RM
		for(v=0;v<n;v++)
		{
			temp = 0.0;
			for(u=0;u<q;u++)
				if( temp<RM[u][v] )
				{
					temp = RM[u][v];
					dec_codeword[v] = root[u];
				}
		}

		//nonbinary --> binary
		for(u=0;u<n;u++)
		{	
			value=dec_codeword[u];
			mask=1;
			for(v=0;v<p;v++)
			{
				if((value & mask)>0)
					dec_bicodeword[p*u+v]=1;
				else
					dec_bicodeword[p*u+v]=0;
				mask=mask<<1;
			}
		}	
	}
	else if(flag==1)	//exist a valid solution
	{
		//nonbinary --> binary
		for(u=0;u<n;u++)
		{	
			value=dec_codeword[u];
			mask=1;
			for(v=0;v<p;v++)
			{
				if((value & mask)>0)
					dec_bicodeword[p*u+v]=1;
				else
					dec_bicodeword[p*u+v]=0;
				mask=mask<<1;
			}
		}		

	}

#ifdef checkFac
	if( flag==0 && (codewordScore>Deg_iterNum) )	//Sm(c) > delta(Cm)
		printf("\n\n this frame fac error by iterNum!!\\n\n");

	if( flag==0 && (codewordScore>Deg_poly) )	//Sm(C) > deg(Q)
		printf("\n\n this frame fac error by poly!!\\n\n");
#endif


	//******************
	flag = 0;
	int count_temp;

	for(j=0;j<lm+1;j++)
	{
		count_temp=0;
		for(u=0;u<k;u++)
			if(output[j][u]==message[u])
			{
				count_temp++;
			}	

		if(count_temp==k)
		{
			flag=1;
			DecSucc_SeqNum++;

			if( min_index>=0 )	//min_index==-1 represent there is no valid solution
				if( min_index!=j )
				{
					ChosenWrong_SeqNum++;
				}

			break;
		}

	}
			
	//*********************


	//Check Herm(64, 39)'s working perperty
	epcount2=0;
	for(u=0;u<n;u++)
		if(codeword[u]!=dec_codeword[u])
			epcount2++;
	//printf("\nepcount2=%d\n", epcount2);

	//******debug*******
	//caculate the testCount2

	testCount2=0;
	for(j=0;j<k;j++)
		if( message[j]!=output[0][j] )
			testCount2++;

/*	//****debug**********
	if( epcount1<=able_correct && epcount2!=0)	//this seq_num has chosen the wrong one
	{
		printf("\n\nDecoding is error\n\n");
	}
*/
}


void findpoints()
{
	int i, j, u, a1, a2, a3, x, y;

	//Initialisation
	for(i=0;i<2;i++)
		for(j=0;j<n;j++)	//w^3
			point[i][j]=0;

	//find points over x^3-y^2-y=0
	u=0;
	for(i=0;i<q;i++)	//number of elements in GF(q)
	{
		x=root[i];
		for(j=0;j<q;j++)	//number of elements in GF(q)
		{
			y=root[j];
			a1=power(x, w+1);
			a2=power(y, w);
			a3=y;
			if(add(a3, add(a1, a2))==0)
			{
				point[0][u]=x;
				point[1][u]=y;
				u++;
			}
		}
	}

}

void tgorder()
{
	int i, j, u, index_temp;

	for(j=0;j<tg_size;j++)
		for(i=0;i<2;i++)
			tg_order[j][i]=-1;

	//judge the index's scale of coresponding tgsize
	/*pesudo code
	if(  0 <= index < w ) tgsize = (index+2) * (index+1) / 2
	else if( w <= index ) tgsize = w*(w+1)/2 + (w+1)*(index-w+1)
	*/
	index_temp = 99;	//according to the last formualtion to calculate out

	j=0;
	for(i=0;i<=index_temp;i++)
		for(u=i;u>=0;u--)
			if(u<=w)
			{
				tg_order[j][0]=u;	// 0<deg_x<=w
				tg_order[j][1]=i-u;	// 0<deg_y
				j++;
			}

}

void generator()
{
	int i, j;
	for(i=0;i<k;i++)
		for(j=0;j<n;j++)	//n
			gmatrix[i][j]=0;

	//generator matrix
	for(i=0;i<k;i++)	//k
		for(j=0;j<n;j++)	//n
			gmatrix[i][j]=mul(power(point[0][j], tg_order[i][0]), power(point[1][j], tg_order[i][1]));

}

void mono_table()
{
	int i, j, u, v, z, weight[weight_Zsize][weight_XYsize], mono_order_1[monoTable_Zsize][weight_XYsize];	//26=(100/weight(z))+1, 100=pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][max(deg_y)+1][w+1]

	//Initialisation
	for(i=0;i<monoTable_Zsize;i++)	//(230/weight(z))+1
	{
		for(j=0;j<weight_XYsize;j++)	//pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][y_size][w+1].
		{
			weight[i][j]=-1;
			mono_order_1[i][j]=-1;
		}
	}

	for(i=0;i<monoTable_Zsize;i++)	// rs
		for(j=0;j<monoTable_Ysize;j++)	//max(deg_y)+1
			for(u=0;u<monoTable_Xsize;u++)	//w+1
				mono_order[i][j][u]=-1;




	//weight of monomial over function x^3-y^2-y=0
	for(i=0;i<weight_Zsize;i++)	
	{
		if(i==0)
			for(j=0;j<weight_XYsize;j++)
				weight[i][j]=tg_order[j][0]*w + tg_order[j][1]*(w+1);
		else if(i>0)
			for(j=0;j<weight_XYsize;j++)
				weight[i][j]= weiz*i + weight[0][j];
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
	z=0;
	for(v=0;v<mono_ordinSize;v++)	//for each possible weight until weight(term[rs-1][max(deg_y)+1][w-1])+1, note term[rs-1][max(deg_y)+1][w-1] is the term next to term[rs-1][max(deg_y)][w]
	{
		for(i=0;i<monoTable_Zsize;i++)	//230/weight(z)
		{
			for(j=0;j<weight_XYsize;j++)	//230, >the index of the largest term with deg(z)=0, and weight=weight(term[rs-1][max(deg_y)+1][w-1])
			{
				if(weight[i][j]==v)
				{
					mono_order_1[i][j]=z;
					z++;
					break;
				}
			}
		}
	}

	// 2-dimensional table transform into 3-dimensional table
	for(u=0;u<monoTable_Zsize;u++)	//rs
		for(z=0;z<monoTable_totalSize;z++)	//>index of term[max(deg_y)][w] in the pole basis + 1
			mono_order[u][tg_order[z][1]][tg_order[z][0]]=mono_order_1[u][z];

#ifdef printfMonotable
	//**********debug********
	for(j=0;j<weight_XYsize;j++)
		printf("\t%d",j);
	printf("\n");
	for(i=0;i<weight_Zsize;i++)
	{
		printf("\n");
		for(j=0;j<weight_XYsize;j++)
			printf("\t%d",weight[i][j]);
	}

	printf("\n\n");
	for(i=0;i<monoTable_Zsize;i++)
	{
		printf("\n");
		for(j=0;j<weight_XYsize;j++)
			printf("\t%d",mono_order_1[i][j]);
	}	
	
	printf("\nMonmial basis is:\n");
	for(i=0;i<monoTable_Zsize;i++)
	{
		printf("\n\nZ=%d",i);

		for(z=0;z<monoTable_Xsize;z++)
			printf("\t%d ",z);
		printf("\n\n");
		for(j=0;j<monoTable_Ysize;j++)
		{
			printf("\n%d",j);
			for(z=0;z<monoTable_Xsize;z++)
				printf("\t%d ", mono_order[i][j][z]);
		}
	}

	//******************
#endif
}

void encoder(int message_temp[], int codeword_temp[])
{
	int i, j;

	for(i=0;i<n;i++)	//n
	{
		codeword_temp[i]=0;
		for(j=0;j<k;j++)	//k
			codeword_temp[i]=add(codeword_temp[i], mul(message_temp[j], gmatrix[j][i]));
	}

}

int cal_delta(int point_temp[3], int j,int alpha,int beta)
{
	int i, u, v, indexSum, flag;
	int a, b, ri;
	int delta, delta_temp1, delta_temp2;

	ri = point_temp[1];

	delta = 0;
	for(i=0;i<interpoly_Zsize;i++)
		for(u=0;u<interpoly_Ysize;u++)
			for(v=0;v<interpoly_Xsize;v++)
				if( (Q_interpoly[j][i][u][v]!=0) && (i>=beta) )
				{
					//search a from x^v and y^u
					indexSum = u+v;
					if( indexSum < w )
					{
						a = (indexSum+1) * indexSum / 2 + u + 1;
						a = a - 1;	//index should begin from 0
					}
					else if( indexSum >= w)
					{
						a = (w+1)*w/2 + (indexSum-w)*(w+1) + (w-v+1);
						a = a - 1;
					}

					//search the coefficient 
					if( coeff_table[point_temp[0]][alpha][a] != 0 )
					{
						b = i;
						flag = (int)comb(b,beta);
						// if flag is even number , delta_temp is equal to 0
						// if flag is odd number, delta_temp can be calculated to the delta
						if( (flag%2)!=0 )
						{
							delta_temp1 = power( ri,(b-beta) );
							delta_temp2 = mul( Q_interpoly[j][i][u][v],coeff_table[point_temp[0]][alpha][a] );
							delta_temp2 = mul( delta_temp2,delta_temp1 );
							delta 		= add( delta,delta_temp2 );
						}
					}
				}

	return delta;
}

float comb(float a,float b)	//calculate combination(a_up,b_down)
{
	int i;
	float temp1,temp2=1.0;

	if( (int)a >= (int)b )
	{
		if( (int)a==(int)b && (int)a==0 )	//C(0,0)
		{
			temp2=1.0;
		}

		for(i=0;i<(int)b;i++)
		{
			temp1=(a-(float)i)/(b-(float)i);
			temp2*=temp1;
		}
	}
	else if( (int)a < (int)b )
	{
		printf("\n\ncomb() error");
	}

	return temp2;	

}


//****coefficientSearch*********
int cal_max_a()
{
	unsigned long int u, gene, iter_num_temp, u_temp, deg_Q, a_temp;
	float temp;
	int lm_temp,tm;

	//calculate gene, C=iter_num, 
	gene=w*(w-1)/2;

	iter_num_temp=n*(max_m+1)*max_m/2;

	//calculate lm
	u_temp=-1;
	for(u=0;u<N;u++)
	{
		temp=((float)(weiz*u)/2.0-gene)*(u-1);
		if(temp<=iter_num_temp)
		{
			u_temp=u;
		}
		else if(temp>iter_num_temp)
		{
			break;
		}
	}

	lm_temp=u_temp-1;

	//calculate tm
	u_temp=-1;
	for(u=0;u<N;u++)
	{
		temp=(lm_temp+1)*u + weiz*(lm_temp+1)*lm_temp/2.0 - lm_temp*gene - gama(u);
		if(temp<=iter_num_temp)
		{	
			u_temp=u;
		}
		if(temp>iter_num_temp)
		{
			break;
		}
	}

	tm=u_temp;

	//calculate deg_Q
	deg_Q = lm_temp*weiz + tm;

	//search a
	a_temp=-1;
	for(u=0;u<tg_size;u++)
	{
		u_temp= tg_order[u][0]*w + tg_order[u][1]*(w+1);
		if( u_temp>=deg_Q )	//key part, determine can be "<" or "<="
		{
			a_temp=u;
			break;
		}
	}

	if(a_temp==-1)
	{
		printf("\n\ncalculate a is erro\n");
	}

	return a_temp;

}

int gama(int u)
{
	int i,j,num,temp;

	j=0;
	num=0;
	for(i=0;i<=u;i++)
	{
		temp=tg_order[j][0]*w + tg_order[j][1]*(w+1);
		if(temp!=i)
			num++;
		else if(temp==i)
			j++;
	}

	return num;
}

//if(alpha<(w+1))
void zerobasis(int pointIndex, int LM_zb[][zb_alpha])
{
	int i,j,u,v;
	int temp,flag,xi,yi,index_flag,temp2;
	int zb_temp1[zb_Xsize], zb_temp2[zb_Ysize][zb_Xsize+1];
	int poly_temp1[zb_Ysize+w][zb_Xsize+w],poly_temp2[zb_Ysize+1][zb_Xsize];

	//Initialisation
	for(i=0;i<zb_alpha;i++)
		for(u=0;u<zb_Ysize;u++)		
			for(v=0;v<zb_Xsize;v++)
			{
				zb[i][u][v] = 0;
			}

	xi = point[0][pointIndex];
	yi = point[1][pointIndex];
	flag = 0;

	//calculate Zb
	for(i=0;i<zb_alpha;i++)
	{

		//set poly_temp1 to be 0
		for(u=0;u<zb_Ysize+w;u++)
			for(v=0;v<zb_Xsize+w;v++)
			{
				poly_temp1[u][v]=0;				
			}
		for(u=0;u<zb_Ysize+1;u++)
			for(v=0;v<zb_Xsize;v++)
			{
				poly_temp2[u][v]=0;
			}

		//zb_a
		temp=i%(w+1);
		if(temp==0)
		{
			//initialization
			for(v=0;v<zb_Xsize;v++)
			{
				zb_temp1[v]=0;
			}

			zb_temp1[0]=1;
			flag=1;	//change the zb_temp2
		}
		else if(temp==1)
		{
			zb_temp1[0]=xi;
			zb_temp1[1]=1;
		}
		else if(temp>1 && temp<=w)
		{
			//x*zb
			poly_temp1[0][0]=0;	//only the first line of poly_temp1 will be uesd
			for(v=0;v<zb_Xsize;v++)
			{
				poly_temp1[0][v+1]=zb_temp1[v];
			}

			//xi*zb+x*zb
			for(v=0;v<zb_Xsize;v++)
			{
				zb_temp1[v]=mul( xi,zb_temp1[v] );
				zb_temp1[v]=add( zb_temp1[v],poly_temp1[0][v] );				
			}

		}


		//zb_b;
		if(flag==1)
		{
			temp=i/(w+1);
			if(temp==0)
			{
				//initialization
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize+1;v++)
					{
						zb_temp2[u][v]=0;
					}

				//calculate zb
				zb_temp2[0][0] = 1;
			}
			else if(temp==1)
			{
				zb_temp2[1][0] = 1;
				zb_temp2[0][1] = power( xi,w );
				zb_temp2[0][0] = mul( zb_temp2[0][1],xi );
				zb_temp2[0][0] = add( zb_temp2[0][0],yi );
			}
			else if(temp>1)
			{
				//cal xi^w*x*zb_temp2
				temp2 = power( xi,w );
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb_temp2[u][v]!=0 )
						{
							poly_temp1[u][v+1]=mul( temp2,zb_temp2[u][v] );
						}
				//x^(w+1)=y^w+y
				for(u=0;u<zb_Ysize;u++)
					if( poly_temp1[u][w+1]!=0 )
					{
						poly_temp1[u+1][0] = add( poly_temp1[u+1][0],poly_temp1[u][w+1] );	//y^1
						poly_temp1[u+w][0] = add( poly_temp1[u+w][0],poly_temp1[u][w+1] );	//y^2
						poly_temp1[u][w+1] = 0;
					}

					
				//cal y+zb_temp2
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb_temp2[u][v]!=0 )
						{
							poly_temp2[u+1][v]=zb_temp2[u][v];
						}

				//real zb_temp2
				temp = power( xi,(w+1) );
				temp = add( temp,yi );
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						{
							zb_temp2[u][v] = mul( zb_temp2[u][v],temp );
							zb_temp2[u][v] = add( zb_temp2[u][v],poly_temp1[u][v] );
							zb_temp2[u][v] = add( zb_temp2[u][v],poly_temp2[u][v] );
						}
			}

			//close the changing of zb_temp2
			flag=0;
		}

		//zb = zb_temp1 * zb_temp2
		for(j=0;j<zb_Xsize;j++)
			if(	zb_temp1[j]!=0 )
			{

				//set poly_temp1 to be 0
				for(u=0;u<zb_Ysize+w;u++)
					for(v=0;v<zb_Xsize+w;v++)
					{
						poly_temp1[u][v]=0;
					}
				
				//cal zb_temp1[j]*zb_temp2
				temp = zb_temp1[j];
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb_temp2[u][v]!=0 )
						{
							poly_temp1[u][v+j] = mul( temp,zb_temp2[u][v] );
						}

				//x^(w+1)=y^w+y
				index_flag = zb_Xsize+j-1;
				while( index_flag>w )
				{
					temp = index_flag-(w+1);
					for(u=0;u<zb_Ysize;u++)
						if( poly_temp1[u][index_flag]!=0 )
						{
							poly_temp1[u+1][temp] = add( poly_temp1[u+1][temp],poly_temp1[u][index_flag] );	//y^1
							poly_temp1[u+w][temp] = add( poly_temp1[u+w][temp],poly_temp1[u][index_flag] );	//y^2
							poly_temp1[u][index_flag] = 0;
						}
					index_flag = index_flag-1;
				}

				
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb[i][u][v]!=0 || poly_temp1[u][v]!=0 )
						{
							zb[i][u][v] = add( zb[i][u][v],poly_temp1[u][v] );
						}				
			}

		//calculate LM_zb and LC_zb
		LM_zb[0][i] = i%(w+1);
		LM_zb[1][i] = i/(w+1);
	}
}

void coefficientSearch(void)
{
	int i, j, u, v, z, flag, temp, temp2;
	int poly[zb_Ysize][zb_Xsize];
	int max_a,u_temp, v_temp, index_x, index_y, choosen_index, lod_temp;
	int LM[2][zb_alpha];
	int effTable[n][tableSize_alpha][tableSize_a];

	//Initialisazion
	for(i=0;i<n;i++)
		for(u=0;u<tableSize_alpha;u++)
			for(v=0;v<tableSize_a;v++)
			{
				effTable[i][u][v]=0;
			}

	for(u=0;u<2;u++)
		for(v=0;v<zb_alpha;v++)
		{
			LM[u][v]=0;
		}

	max_a = cal_max_a();
	if( max_a >= tableSize_a )
	{
		printf("\n\nmax_a = %d, tableSize_a should be set larger\n\n", max_a);
	}


	//calculate coeff_table
	for(i=0;i<n;i++)
	{
		//initialization poly
		for(u=0;u<zb_Ysize;u++)
			for(v=0;v<zb_Xsize;v++)
			{
				poly[u][v]=0;
			}

		zerobasis(i,LM);	//calculate all the zb which index under tableSize_alpha and its leading monomial

		for(j=max_a;j>=0;j--)
		{
			//initialize coefficient
			for(z=0;z<(max_alpha+1);z++)
			{
				effTable[i][z][j]=0;
			}
			
			//search the choosen zb
			index_x = tg_order[j][0];
			index_y = tg_order[j][1];

			//cal the first coeffTable value
			for(z=zb_alpha-1;z>=0;z--)
				if( LM[0][z]==index_x && LM[1][z]==index_y )
				{
					choosen_index = z;
					if( z<tableSize_alpha )
					{
						effTable[i][z][j]=1;
					}
					else
					{
						printf("\n\n0.tableSize_alpha is not large enough\n\n");
					}
					
					break;
				}

			//initialize the first poly 
			for(u=0;u<zb_Ysize;u++)
				for(v=0;v<zb_Xsize;v++)
					if( zb[choosen_index][u][v]!=0 )
					{
						poly[u][v]=zb[choosen_index][u][v];
					}

			//set the item of poly as LM[choosen_index] to be 0
			poly[index_y][index_x]=0;

			//cal the coeffTable
			flag=0;	//judge if the poly is empty 
			for(u=0;u<zb_Ysize;u++)
			{
				for(v=0;v<zb_Xsize;v++)
					if( poly[u][v]!=0 )
					{
						flag=1;
						break;
					}

				if( flag==1 )
				{
					break;
				}
			}

			//start to search coefficient
			while( flag!=0 )
			{
				//find the poly's the second largest item
				lod_temp = -1;
				index_x = -1;
				index_y = -1;
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( poly[u][v]!=0 && lod_temp<mono_order[0][u][v] )
						{
							lod_temp=mono_order[0][u][v];
							index_x = v;
							index_y = u;
						}

				//check the correctess of the progress of search second largest item
				if( index_x==-1 || index_y==-1 )
				{
					printf("\n\nsearch progress has error\n\n");
				}
				
				//find the choose_index of zb
				for(z=zb_alpha-1;z>=0;z--)
					if( LM[0][z]==index_x && LM[1][z]==index_y )
					{
						choosen_index = z;
						if( z<tableSize_alpha )
						{
							effTable[i][z][j] = poly[index_y][index_x];
						}
						else
						{
							printf("\n\n1.tableSize_alpha is not large enough\n\n");
						}
						break;
					}

				//update the poly
				temp = poly[index_y][index_x];
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb[choosen_index][u][v]!=0 )
						{
							temp2 =	mul( temp,zb[choosen_index][u][v] );
							poly[u][v] = add( temp2,poly[u][v] );
						}
				
				//check if the poly is empty
				flag=0;	//judge if the poly is empty 
				for(u=0;u<zb_Ysize;u++)
				{
					for(v=0;v<zb_Xsize;v++)
						if( poly[u][v]!=0 )
						{
							flag=1;
							break;
						}

					if( flag==1 )
					{
						break;
					}
				}
				//**************************

			}
			//************************
		}
		//**********************
	}

	for(i=0;i<n;i++)
		for(u=0;u<(max_alpha+1);u++)
			for(v=0;v<(max_a+1);v++)
			{
				coeff_table[i][u][v] = effTable[i][u][v];
			}

#ifdef printCoeffTable
	//*********debug: printf the coeffTable*********
	printf("\ncoeffTable for Herm(8,4)");
	for(i=0;i<n;i++)
	{
		printf("\n**************************");
		printf("\npoint_%d(%d,%d):\n", i, point[i][0], point[i][1]);
		
		for(v=0;v<(max_a+1);v++)
		{
			printf("\t%d",v);
		}
		
		for(u=0;u<(max_alpha+1);u++)
		{
			printf("\n\n");
			for(v=0;v<(max_a+1);v++)
			{
				printf("\t%d", coeff_table[i][u][v]);
			}
		}

		printf("\n\n");

	}
	printf("\n");
	/**********************************/
#endif

}	
//********************************


#ifdef cheatingDecoding
void decoding(void)
{
	int i, j, value;
	unsigned int mask=1;
	float proba_temp;


	//cal the deg(1,wz) corresponding to Q_min
	Deg_iterNum = -1;
	int	flag_iterNum = 1;	//make the search part more effciency
		for(int u=0; u<(lm+1) && flag_iterNum ;u++)
			for(j=0; j<monoTable_Ysize && flag_iterNum ;j++)
				for(i=0; i<monoTable_Xsize && flag_iterNum ;i++)
					if ( mono_order[u][j][i] == iterNum )
					{
						Deg_iterNum = i*w + j*(w+1) + u*weiz;
						flag_iterNum = 0;
					}

	if(Deg_iterNum < 0)
	{
		printf("\n\nDeg_poly is error\n\n");
	}

	//decoding
	if( codewordScore > Deg_poly )	//decdoing correctly
	{
		//codeword
		for(i=0;i<n;i++)
			dec_codeword[i] = codeword[i];

		//bi_codeword
		for(i=0;i<n*p;i++)
			dec_bicodeword[i] = bi_codeword[i];
	}
	else if(codewordScore <= Deg_poly)	//decoding incorrectly
	{
		//search dec_codeword[n] from RM
		for(i=0;i<n;i++)
		{
			proba_temp = 0.0;
			for(j=0;j<q;j++)
				if( proba_temp<RM[j][i] )
				{
					proba_temp = RM[j][i];
					dec_codeword[i] = root[j];
				}
		}

		//bi_codeword[n*p]
		//nonbinary --> binary
		for(i=0;i<n;i++)
		{	
			value=dec_codeword[i];
			mask=1;
			for(j=0;j<p;j++)
			{
				if((value & mask)>0)
					dec_bicodeword[p*i+j]=1;
				else
					dec_bicodeword[p*i+j]=0;
				mask=mask<<1;
			}
		}	

	}

}
#endif







