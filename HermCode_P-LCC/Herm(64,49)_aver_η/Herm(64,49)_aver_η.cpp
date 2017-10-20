#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "main.h"
#include "FiniteFieldBasisGF(16).h"
#include "FiniteFieldBasisCal.h"



//main()
float N0;
float pi=3.141593;	// pai

//*************
int bi_message[k*p], message[k];	//transmitted messge[k]
int codeword[n], bi_codeword[n*p]; //codewords
float tx_symbol[p*n][2], rx_symbol[p*n][2], sgm;
float RM[q][n];
//findpoint()
int point[n][2];	//rational points[2][w^3], [0]->x, [1]->y
//tgorder()
int tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
//mono_table()
int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
int gmatrix[k][n];	//generator matrix[k][n]
int dec_message[k],dec_codeword[n],dec_bicodeword[n*p]; //decoding result
//test_vec_construction
int eta, test_vec_num;
int large_vec[choose_num][n];
float test_set_ordered[2][n];
//int test_vec_com[n];
int test_vec_com[n];
int x_ordered[2][n];
//interpolation
int bi_rxword[p*n], rxword[n]; //received wordss
int epcount1, testCount1[test_vec_num_max], testCount2, testCount1_com;
int Q_interpoly[test_vec_num_max][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];	//chosen interpolated polynomial [z-deg+1][max(deg_y)+1][w+1]
int uu;	//factorisation step index
int l, listNum[test_vec_num_max];	//candidate output index
//factorization() and rcs()
int Q[k][facpoly_Zsize][facpoly_Ysize][facpoly_Xsize];	//sequtial deduction polynomial [number of fac steps=k][rs][y_size][w+1], y_size> maxdeg_y]+rs*(deg_¦µ(k-1))
//int rootlist[k][lm+1];	//the list of roots [number of fac steps=k][expected number of roots in each solution, 5>rs]
int output[lm+1][k], outputList[test_vec_num_max][lm+1][k];	//the list of candidate message [expeced number of candidate message, >rs][length of message, k]
int expoly[2][expoly_Ysize][expoly_Xsize];	//expanded polynomial in [z+f_k-1-u*pb_k-1-u]^rs, expoly[rs][3>(max(deg_y) in encoding functions)*(rs-1)][3>(max(deg_x) in encoding functions)*(rs-1)]
//**********************


int Q_uncom_elem[test_vec_num_max][init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	//****debug**********
	int degree_test_temp[test_vec_num_max][init_polyNum];  //store the degree of leading monomial of polynomial
	int degree_test[test_vec_num_max];		//store the degree of the choosen polynomial

	unsigned long int seq_num_Now;	//record the seq_num at the moment
	unsigned long int DecSucc_SeqNum, ChosenWrong_SeqNum;	//record the seq_num should be decoding success in genius mode
	double CWR;	//record the rate of chosen uncorrectly

	//some var about computation complexity
	unsigned long int addNum, mulNum, addNum_temp, mulNum_temp;
	int flag_addNum=0, flag_mulNum=0;


	//r_thread
	float r_thread_min = 1.0;
	
	//GammaValue
	float gamma[eta_max][3];	//gamma_i(min, max, average)
	int eta_count[eta_max];
	double outputlist_count;
	unsigned long int outputlistNum;

	//CountAverEta
	int index_eta;
	double aver_eta;

	//*******************

void findpoints(void);
void mono_table(void);
void polyexp1(int, int, int, int, int poly[][expoly_Ysize][expoly_Xsize]);
//void zerobasis(void);
//void polebasis(void);
void tgorder(void);
//void coefficient(void);
void test_vec_contruction(void);
void interpolation(void);
void com_elem_interpolation(int g[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int interpoint[][3]);
void uncom_elem_interpolation(int interpoint[4],int inGroup[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int outGroup1[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int outGroup2[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize]);
//void factorisation(void);
void factorisation(int Q_input[test_vec_num_max][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int output_list[test_vec_num_max][lm + 1][k], int list_num[test_vec_num_max]);
void rcs(int);
//void choose(void);
void choose(int output_codeword[n], int output_bicodeword[n*p], int output_list[test_vec_num_max][lm + 1][k], int list_num[test_vec_num_max], int flag_decoding_alg);
void generator(void);
void encoder(int message_temp[], int codeword_temp[]);
void modulation(void);
void channel(void);
void demodulation(void);
float PDF(int,int);
int result_compare(int output_list[test_vec_num_max][lm + 1][k], int output_codeword[n], int output_list_BF[test_vec_num_max][lm + 1][k], int output_codeword_BF[n]);

//cheating
void cheatDec(void);
void cheatFac(void);


void main()
{
	int i, u, m, num, value;
	unsigned int mask = 1;
	long int error, ferror;
	unsigned long int seq_num;	//number of input binary sequences
	unsigned long int j,v;
	float start, finish, SNR;
	double BER, FER;
	double progress;
	//double channelError_count, successError_count;
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
		printf("affine point[%d]=(%d,%d)\n",i,point[i][0],point[i][1]);
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
	FILE * fp_input = fopen("basic_input.txt","r");
	if(fp_input!=NULL) {
		fscanf(fp_input,"start SNR:%f\n", &start);
		fscanf(fp_input,"finish SNR:%f\n", &finish);
		fscanf(fp_input,"seq_num:%d", &seq_num);

		fclose(fp_input);
	}else{
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



#ifdef _PrintGammaValue_
	FILE * fp_ch = fopen("Channel_GammaValue.txt", "w");
#endif

	for(SNR=start; SNR<=finish; SNR=SNR+interval)
	{
		N0=(1.0/((float)k/(float)n))/pow(10.0, SNR/10.0);
		sgm=sqrt(N0/2);
		
		error=0;
		ferror=0;

		//channelError_count=0.0;
		//successError_count=0.0;
		v=0;

		DecSucc_SeqNum=0;
		ChosenWrong_SeqNum=0;

		outputlist_count = 0.0;
		addNum_count = 0.0;
		mulNum_count = 0.0;
		totalNum_count = 0.0;
		
#ifdef CountAverEta
		aver_eta = 0.0;
#endif
		
		r_thread_min = 1.0;

#ifdef _PrintGammaValue_
		fprintf(fp_ch, "\nSNR=%0.2f", SNR);
#endif
		
		//*****debug******
		for (i = 0; i < eta_max; i++)
		{
			gamma[i][0] = 1.0;
			gamma[i][1] = 0;
			gamma[i][2] = 0;
		}
		for (i = 0; i < eta_max; i++)
			eta_count[i] = 0;

		//*****************

		for(j=1;j<=seq_num;j++)
		{
			addNum=0;
			mulNum=0;

			//*****debug*****
			seq_num_Now=j;
			//***************
#ifndef cheatingEncoding
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
//			message[0]=1;	message[1]=2;	message[2]=1;	message[3]=3;
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
#endif

			//modulation
			modulation();

			//channel
			channel();

			//demodulation
			demodulation();

			//LIST DECODER
			//test vector construction
			test_vec_contruction();

#ifdef _PrintGammaValue_
			//fprintf(fp_ch, "\n");
			//for (u = 3*n/4-1; u < n; u++)
			//	fprintf(fp_ch, "%0.2f\t", test_set_ordered[0][u]);
#endif

#ifndef CountAverEta
#ifdef cheatingDec
			//********LCC interpolation*****************
			//cheatingDecoding
			cheatDec();

			//frame error rate calculation
			int temp = ferror;
			for(u=0;u<n;u++)
				if(dec_codeword[u]!=codeword[u])
				{
					ferror++;
					break;
				}

			if(temp==ferror)
				v++;

#else
			//interpolation
			interpolation();

#ifndef cheatingFac
			//factorisation
			factorisation(Q_interpoly, outputList, listNum);
			//choose
			choose(dec_codeword, dec_bicodeword, outputList, listNum, 0);

			//bit error rate calculation
			int temp = error;
			for (u = 0; u<n*p; u++)	//n*4
			{
				if (dec_bicodeword[u] != bi_codeword[u])
					error++;
			}

			//frame error rate calculation - LCC
			if (error>temp)
				ferror++;
			else if (error == temp) //decoding success
			{
				//calculate the aver error of decoding success situation
				v++;
				//successError_count = successError_count + (epcount1 - successError_count) / (double)(v);
			}

#else
			//cheatingFac
			cheatFac();

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

#endif

#endif

#endif
			//channelError_count = channelError_count + (epcount1 - channelError_count) / (double)(j);
			addNum_count = addNum_count + (addNum - addNum_count) / (double)(j);
			mulNum_count = mulNum_count + (mulNum - mulNum_count) / (double)(j);
			totalNum_count = addNum_count + mulNum_count;

#ifdef	CountAverEta
			aver_eta = aver_eta + (index_eta - aver_eta) / (double)(j);
#endif

			progress = (double)(j * 100) / (double)seq_num;
			BER = (double)(error) / (double)(n*p*j);
			FER = (double)(ferror) / (double)(j);
			//***********************

			if (j % 400 == 0)
				printf("Progress=%0.1f, seqNow=%d, SNR=%2.2f, BitErrors=%2.1d, BER=%E, FrameErrors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\r", progress, seq_num_Now, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);

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

		printf("Progress=%0.1f, seqNow=%d, SNR=%2.2f, BitErrors=%2.1d, BER=%E, FrameErrors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n", progress, seq_num_Now, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);

		OpenFile;
		fprintf(fp, "Progress=%0.1f, seqNow=%d, SNR=%2.2f, BitErrors=%2.1d, BER=%E, FrameErrors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n", progress, seq_num_Now, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);
		fclose(fp);


		//****debug******
#ifdef CountAverEta
		FILE *fp_eta = fopen("Average.txt", "a");
		fprintf(fp_eta, "SNR=%0.1f\tAverEta=%E\n", SNR, aver_eta);
		fclose(fp_eta);
#endif

#ifdef _PrintRThread_
		FILE *fp_r_thread = fopen("Record_r_thread.txt", "a");
		fprintf(fp_r_thread, "η=%d, SNR = %0.1f:\tr_min in n-η = %0.3f\n", eta_max, SNR, r_thread_min);
		fclose(fp_r_thread);
#endif

#ifdef _PrintGammaValue_
		fprintf(fp_ch,"\n");
		for (i = 0; i < eta_max; i++)
		{
			fprintf(fp_ch, "location[%d]:\tmin(%0.3f)\tmax(%0.3f)\taver(%0.3f)\n", i, gamma[i][0], gamma[i][1], gamma[i][2]);
		}
		fprintf(fp_ch, "\neta count:\n");
		int eta_temp = 0;
		for (i = 0; i < eta_max; i++)
		{
			fprintf(fp_ch, "eta[%d]=%d\t", i+1, eta_count[i]);
			eta_temp += eta_count[i];
		}
		fprintf(fp_ch, "\neta count rate:\n");
		for (i = 0; i < eta_max; i++)
		{
			fprintf(fp_ch, "rate[%d]=%0.2f\t", i+1, (float)eta_count[i]/eta_temp*100);
		}
		fprintf(fp_ch, "\n");

#endif

		//****************

	}
	//*****debug********
#ifdef _PrintGammaValue_
	fclose(fp_ch);
#endif
	//****************
	
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
			gmatrix[i][j]=mul(power(point[j][0], tg_order[i][0]), power(point[j][1], tg_order[i][1]));

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

void modulation(void)
{
	int i;

	//BPSK
	for(i=0;i<n*p;i++)	
	{
		tx_symbol[i][0]=(double)(-2*bi_codeword[i]+1);
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
		for(int x=0;x<pointNum;x++)
		for(int y=0;y<pointNum;y++)
//		for(int h=0;h<pointNum;h++)
		{
			RM[j][i]=(Pr[0][y]*Pr[1][x]*Pr[2][v]*Pr[3][u]);  
			j++;
		}
	
	}


/*	//*****debug**********
	printf("\n\ncodeword");
	for(v=0;v<n;v++)
		printf("\t%d",codeword[v]);

	printf("\n\nrxword\t");
	for(v=0;v<n;v++)
		printf("\t%d",rxword[v]);
	
	printf("\n\n");
	for(v=0;v<n*p/4;v++)
		printf("\t(%f,%f,%f,%f)",rx_symbol[v*2+0][0],rx_symbol[v*2+0][1],rx_symbol[v*2+1][0],rx_symbol[v*2+1][0]);

	printf("\n\n");
	for(u=0;u<q;u++)
	{
		printf("\n");
		printf("\t%d",root[u]);
		for(v=0;v<n;v++)
			printf("\t%f",RM[u][v]);
	}
	printf("\n\n");
*/	//*******************

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


void test_vec_contruction()
{
	//construct the test_vec[test_vec_num_max]
	int i,j,temp_index;
	float temp;
	float test_set[2][n]; // [0] stores value, [1] stores index
	float large_vec_proba[choose_num][n]; // largest index is 0, second largest index is 1

	//cal test_set
	for(i=0;i<n;i++)
	{
		//finde the largest vec
		large_vec_proba[0][i]=0;
		for(j=0;j<q;j++)
			if( large_vec_proba[0][i]<RM[j][i] )
			{
				large_vec_proba[0][i]=RM[j][i];
				large_vec[0][i]=root[j];
			}
		//find the second largest vec
		large_vec_proba[1][i]=0;
		for(j=0;j<q;j++)
			if( (large_vec_proba[1][i]<RM[j][i]) && (root[j]!=large_vec[0][i]) )
			{
				large_vec_proba[1][i]=RM[j][i];
				large_vec[1][i]=root[j];
			}
		
		test_set[0][i]=large_vec_proba[1][i]/large_vec_proba[0][i];
		test_set[1][i]=i;
	}

	
	epcount1=0;
	for(i=0;i<n;i++)	//n
		if(codeword[i]!=large_vec[0][i])
			epcount1++;

/*	//***********debug*************
	printf("test_set:\n");
	for(i=0;i<n;i++)
		printf("\t%6.3f",test_set[0][i]);

	printf("\n\n");
	for(i=0;i<n;i++)
		printf("\t%d",(int)test_set[1][i]);
	printf("\n\n");
*/	//*************************************

	//order test_vec with bubble sorting order
	//key part!!
	for(i=0;i<n;i++)
	{
		test_set_ordered[0][i]=test_set[0][i];
		test_set_ordered[1][i]=test_set[1][i];
	}

	for(i=0;i<(n-1);i++)
		for(j=0;j<(n-(i+1));j++)
			if( test_set_ordered[0][j]>test_set_ordered[0][j+1] )
			{
				temp=test_set_ordered[0][j];
				test_set_ordered[0][j]=test_set_ordered[0][j+1];
				test_set_ordered[0][j+1]=temp;

				temp_index=test_set_ordered[1][j];
				test_set_ordered[1][j]=test_set_ordered[1][j+1];
				test_set_ordered[1][j+1]=temp_index;
			}

/*	//**************debug********************
	printf("test_set_order:\n");
	for(i=0;i<n;i++)
		printf("\t%6.3f",test_set_ordered[0][i]);

	printf("\n\n");
	for(i=0;i<n;i++)
		printf("\t%d",(int)test_set_ordered[1][i]);
	printf("\n\n");

*/	//***************************************

	//**********debug************
	if (r_thread_min > test_set_ordered[0][n - eta_max+4])
		r_thread_min = test_set_ordered[0][n - eta_max+4];

	//***************************

	//set the value of eta
	eta = 0;
	for (i = n - 1; i >= 0; i--)
		if (test_set_ordered[0][i] > r_thread)
		{
			eta++;
			if (eta > eta_max)
			{
				eta--;
				break;
			}
		}
		else 
			break;
		
	if (eta < eta_min)	//make sure η>=η_min
		eta = eta_min;

	if (eta == 0 && eta > eta_max)
		printf("\neta > eta_max is error\n");

	test_vec_num = (int)pow(2.0, eta);

#ifdef CountAverEta
	index_eta = 0;
	index_eta = eta;
#endif

	//*******debug************
	//caculate the value of gamma value
	for (i = 0; i <= eta_max-1; i++)
	{
		//min_value
		if (gamma[i][0] > test_set_ordered[0][n - 1 - i])
			gamma[i][0] = test_set_ordered[0][n - 1 - i];

		//max_value
		if (gamma[i][1] < test_set_ordered[0][n - 1 - i])
			gamma[i][1] = test_set_ordered[0][n - 1 - i];

		//averge_value, need to be initialize
		gamma[i][2] = gamma[i][2] + (test_set_ordered[0][n - 1 - i] - gamma[i][2]) / (float)seq_num_Now;
	}
	//count the num of difference eta
	eta_count[eta-1]++;
	//************************

	//contruction test_vec
	for (i = 0; i < n; i++)
		if (i < n - eta)
			test_vec_com[i] = large_vec[0][(int)test_set_ordered[1][i]];
		else
			test_vec_com[i] = -1;

	//construct x_ordered
	for(i=0;i<n;i++)
		for(j=0;j<2;j++)
			x_ordered[j][i]=point[(int)test_set_ordered[1][i]][j];

/*	//*****debug*********
	//caculate the num of diffs between codeword and test_vec[i] 
	testCount1[0]=0;
	for(i=0;i<n-eta_max;i++)
		if(codeword[(int)test_set_ordered[1][i]]!=test_vec_com[i])
			testCount1[0]++;

//	printf("\n\nfirst n-eta_max error num :%d\n\n",testCount1[0]);
	
	testCount1_com=testCount1[0];
	for(i=1;i<test_vec_num_max;i++)
		testCount1[i]=testCount1[0];

	//the num of loop equals to eta_max
	unsigned long int mask;
	for(i=0;i<test_vec_num_max;i++)
	{
		mask=1;
		for(int j=n-1;j>=n-eta_max-1;j--)
		{
			//bit caculate
			if( (i & mask)>0 )
				temp_index=1;
			else 
				temp_index=0;

			if( codeword[(int)test_set_ordered[1][j]] != large_vec[temp_index][(int)test_set_ordered[1][j]])
				testCount1[i]++;
		
			mask=mask<<1;
		}	
	}
*/	//*******************


/*	//************debug***************
	printf("\n\n");
	printf("large_vec[]\t");
	for(j=0;j<n;j++)
		printf("\t%d",large_vec[0][j]);

	printf("\n\n");
	printf("large_vec2[]\t");
	for(j=0;j<n;j++)
		printf("\t%d",large_vec[1][j]);

	printf("\n\n");
	printf("x_order[]\t");
	for(i=0;i<n;i++)
		printf("\t(%d,%d)",x_ordered[0][i],x_ordered[1][i]);

	
	printf("\n\n");
	printf("large_vec_order[]");
	for(j=0;j<n;j++)
		printf("\t%d",large_vec[0][(int)test_set_ordered[1][j]]);
	
	printf("\n\n");
	printf("large_vec2_order[]");
	for(j=0;j<n;j++)
		printf("\t%d",large_vec[1][(int)test_set_ordered[1][j]]);
*/
	//printf("\n\n");
	//for(int u=0;u<q;u++)
	//{
	//	printf("\n");
	//	printf("\t%d",root[u]);
	//	for(int v=0;v<n;v++)
	//		printf("\t%f",RM[u][v]);
	//}

	//********************************		



}

void interpolation()
{
	int i,j,u,v,z,num,temp,temp_index;
	int com_elem_interpoint[n][3];
	int Q_com_elem[init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];
//	int uncom_elem_interpoint[eta_max][4];
	int degree_temp[test_vec_num_max][init_polyNum];
	//int Q_uncom_elem[test_vec_num_max][init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	//common element interpolation
	//set common element interpoint (xi,ri)
	for(i=0;i<n;i++)
		if (i<n-eta)
		{
			com_elem_interpoint[i][0] = x_ordered[0][i];
			com_elem_interpoint[i][1] = x_ordered[1][i];
			com_elem_interpoint[i][2] = test_vec_com[i];
		}
		else
		{
			com_elem_interpoint[i][0] = -1;
			com_elem_interpoint[i][1] = -1;
			com_elem_interpoint[i][2] = -1;
		}

/*	//*********debug***********
	printf("\n\ncom_elem_interpoint\n");
	for(i=0;i<n-eta_max;i++)
		printf("\t%d",com_elem_interpoint[i][0]);
	printf("\n");
	for(i=0;i<n-eta_max;i++)
		printf("\t%d",com_elem_interpoint[i][1]);
*/	//***************************

	for(i=0;i<init_polyNum;i++)	//num_polys
		for(j=0;j<interpoly_Zsize;j++)	//z-deg+1
			for(u=0;u<interpoly_Ysize;u++)	//max(deg_y)+1
				for(v=0;v<interpoly_Xsize;v++)	//w+1
					Q_com_elem[i][j][u][v]=0;

#ifdef _PrintLCCLod_
	printf("\nLCC interpolation Lod:");
#endif

	com_elem_interpolation(Q_com_elem,com_elem_interpoint);
//	NewInter(Q_com_elem, com_elem_interpoint);


	//com_elem interpolation finish

//conditional compile
#ifndef _GS_Normal_
	//uncommon element interpolation
	int uncom_elem_interpoint[eta_max][4];
	for (i = 0; i < eta_max; i++)
		for (u = 0; u < 4; u++)
			uncom_elem_interpoint[i][u] = -1;
	//set uncomon element interpolation
	for(i=n-eta;i<n;i++)
	{
		uncom_elem_interpoint[i-n+eta][0]=x_ordered[0][i];
		uncom_elem_interpoint[i-n+eta][1]=x_ordered[1][i];
		uncom_elem_interpoint[i-n+eta][2]=large_vec[0][(int)test_set_ordered[1][i]];
		uncom_elem_interpoint[i-n+eta][3]=large_vec[1][(int)test_set_ordered[1][i]];
	}

/*	//******debug********
	printf("\n\nuncom_elem_interpoint\n");
	for(i=0;i<eta_max;i++)
		printf("\t%d",uncom_elem_interpoint[i][0]);
	printf("\n");
	for(i=0;i<eta_max;i++)
		printf("\t%d",uncom_elem_interpoint[i][1]);
	printf("\n");
	for(i=0;i<eta_max;i++)
		printf("\t%d",uncom_elem_interpoint[i][2]);
*/	//*********************

	//initialize
	for(z=0;z<test_vec_num_max;z++)
		for(i=0;i<init_polyNum;i++)
			for(j=0;j<interpoly_Zsize;j++)
				for(u=0;u<interpoly_Ysize;u++)
					for(v=0;v<interpoly_Xsize;v++)
						Q_uncom_elem[z][i][j][u][v]=0;

	for(i=0;i<init_polyNum;i++)
		for(j=0;j<interpoly_Zsize;j++)
			for(u=0;u<interpoly_Ysize;u++)
				for(v=0;v<interpoly_Xsize;v++)
					Q_uncom_elem[0][i][j][u][v]=Q_com_elem[i][j][u][v];	

	//interpolation
	for(i=0;i<eta;i++)
	{
		num=(int)pow(2.0,i);
		for(u=num-1;u>=0;u--)
			uncom_elem_interpolation(uncom_elem_interpoint[i],Q_uncom_elem[u],Q_uncom_elem[2*u+0],Q_uncom_elem[2*u+1]);
	}		

	//******debug******

	//****************

	//choose the poly to factorization

	//initialize
	for(i=0;i<test_vec_num_max;i++)
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					Q_interpoly[i][u][v][z]=0;

	//eta_max>0
#ifdef _PrintGroupLod_
	printf("\n");
#endif
	for(i=0;i<test_vec_num;i++)
	{
		//calculate the degree of poly
		for(j=0;j<init_polyNum;j++)
		{
			temp=-1;
			degree_temp[i][j]=0;
			//*******debug*******
			degree_test_temp[i][j]=0;
			//******************
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						if( Q_uncom_elem[i][j][u][v][z]!=0 )
							if( temp < mono_order[u][v][z] )
							{
								temp=mono_order[u][v][z];

								//***debug**************
								degree_test_temp[i][j]=z*w+v*(w+1)+u*weiz;
								//***********************
							}

			degree_temp[i][j]=temp;
		}

		//choose the min degree
		temp_index=0;
		temp=degree_temp[i][0];
		for(j=1;j<init_polyNum;j++)
			if( temp>degree_temp[i][j] && degree_temp[i][j]<=iterNum )
			{
				temp=degree_temp[i][j];
				temp_index=j;

			}
		
		//****debug*******
		degree_test[i]=degree_test_temp[i][temp_index];
		
#ifdef _PrintGroupLod_
			printf("\ntest_vector[%d]'s lod = [%d,\t%d,\t%d,\t%d]", i, degree_temp[i][0], degree_temp[i][1], degree_temp[i][2], degree_temp[i][3]);
#endif
			//******************
		
		//assignment
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					Q_interpoly[i][u][v][z]=Q_uncom_elem[i][temp_index][u][v][z];
	}
#ifdef _PrintGroupLod_
	printf("\n");
#endif
	//******************************
#else
	//eta_max=0
	for(i=0;i<test_vec_num_max;i++)
	{
		//calculate the degree of poly
		for(j=0;j<init_polyNum;j++)
		{
			temp=-1;
			degree_temp[i][j]=0;
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						if( Q_com_elem[j][u][v][z]!=0 )
							if( temp < mono_order[u][v][z] )
								temp=mono_order[u][v][z];

			degree_temp[i][j]=temp;
		}

		//choose the min degree
		temp_index=0;
		temp=degree_temp[i][0];
		for(j=1;j<init_polyNum;j++)
			if( temp>degree_temp[i][j] )
			{
				temp=degree_temp[i][j];
				temp_index=j;
			}

		degree_test[i] = degree_temp[i][temp_index];

		//assignment
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					Q_interpoly[i][u][v][z]=Q_com_elem[temp_index][u][v][z];
	}
	//********************
#endif
	
/*	//********debug*************
	//used to prove the polynomials choosen for fac equals to 0 over all the interpoint
	//prove the efficiencies of the interpolation
	//printf("\ninterpolation check procesing result:\n");
	int temp_x, temp_y, temp_z, temp_flag;
	for(i=0;i<test_vec_num;i++)
	{
		temp_flag = 0;
		for (j = 0; j < n - eta; j++)
		{
			temp = 0;
			for (u = 0; u < interpoly_Zsize; u++)	//rs
				for (v = 0; v < interpoly_Ysize; v++)	//max(deg_y)+1
					for (z = 0; z < interpoly_Xsize; z++)	//w+1
						if (Q_interpoly[i][u][v][z] != 0)
						{
							temp_x = power(com_elem_interpoint[j][0], z);
							temp_y = power(com_elem_interpoint[j][1], v);
							temp_z = power(com_elem_interpoint[j][2], u);
							temp = add(temp, mul(Q_interpoly[i][u][v][z], mul(temp_z, mul(temp_y, temp_x))));
						}
			if (temp != 0)
			{
				printf("\nQ_interpoly[%d] in com_point[%d](%d,%d,%d) = %d", i, j, com_elem_interpoint[j][0], com_elem_interpoint[j][1], com_elem_interpoint[j][2], temp);
				temp_flag = 1;
			}
		}
		unsigned int mask = (int)pow(2.0, eta-1);
		int value = 0;

		for (j = 0; j < eta; j++)
		{
			value = (i&mask) > 0 ? 1 : 0;
			mask = mask >> 1;

			temp = 0;
			for (u = 0; u < interpoly_Zsize; u++)	//rs
				for (v = 0; v < interpoly_Ysize; v++)	//max(deg_y)+1
					for (z = 0; z < interpoly_Xsize; z++)	//w+1
						if (Q_interpoly[i][u][v][z] != 0)
						{
							temp_x = power(uncom_elem_interpoint[j][0], z);
							temp_y = power(uncom_elem_interpoint[j][1], v);
							temp_z = power(uncom_elem_interpoint[j][value + 2], u);
							temp = add(temp, mul(Q_interpoly[i][u][v][z], mul(temp_z, mul(temp_y, temp_x))));
						}
			if (temp != 0)
			{
				printf("\nQ_interpoly[%d] in uncom_point[%d](%d,%d,%d) = %d", i, j, uncom_elem_interpoint[j][0], uncom_elem_interpoint[j][1], uncom_elem_interpoint[j][value + 2], temp);
				temp_flag = 1;
			}
		}
		
		if (temp_flag)
			printf("\n");
	}
*/	//**************************

}


void com_elem_interpolation(int g[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int interpoint[][3])
{
	int i, j, u, v, z, delta_temp, delta_temp1, lod_temp, temp1, temp2; 
	int delta[init_polyNum],  J[init_polyNum], act[init_polyNum], lod[init_polyNum], lod_min, j_min;	//(delta, J, act, lod)[num of polys]
	int f[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];  //g[num_polys][z-deg+1][max(deg_y)+1][w+1], g1[z-deg+1][max(deg_y)+1][w+2], (g2, f)[z-deg+1][max(deg_y)+1][w+1]	
	int g1[interpoly_Zsize][interpoly_Ysize+w][interpoly_Xsize+1], g2[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];
	int interpoly_Ysize_1, interpoly_Xsize_1;

	interpoly_Ysize_1 = interpoly_Ysize + w;	//used to initiliaze g1[]
	interpoly_Xsize_1 = interpoly_Xsize + 1;	//used to initiliaze g1[]

	//Initialisation
	for(i=0;i<init_polyNum;i++)	//num_polys
		for(j=0;j<interpoly_Zsize;j++)	//z-deg+1
			for(u=0;u<interpoly_Ysize;u++)	//max(deg_y)+1
				for(v=0;v<interpoly_Xsize;v++)	//w+1
					g[i][j][u][v]=0;
	
	for(i=0;i<(lm+1);i++)	//rs
		for(j=0;j<w;j++)	//w
			g[w*i+j][i][j][0]=1;	//j+w*i
	

	//Interpolation
	for(i=0;i<n-eta;i++)	//wrt each point
	{
		//Calculate each polynomial's leading order
		for(j=0;j<init_polyNum;j++)	//num_poly
		{
			lod_temp=0;
			lod[j]=0;
			for(u=0;u<interpoly_Zsize;u++)	//rs
			{
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				{
					for(z=0;z<interpoly_Xsize;z++)	//w+1
					{
						if(g[j][u][v][z]!=0)
						{
							lod_temp=mono_order[u][v][z];
							if(lod_temp>lod[j])
								lod[j]=lod_temp;
						}
					}
				}
			}
		}
			
#ifndef _NoReductionCom_
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
			
		flag_addNum=1;
		flag_mulNum=1;



		//Calculate the hasse derivative mapping of each polynomials
		j_min=-1;
		for(j=0;j<init_polyNum;j++)	//wrt each polynomial
		{
			J[j]=0;
			delta[j]=0;
			if(act[j]==1)	//for polynomials with leading order less of equal to C
			{
				//Hasse derivative mapping
				for(u=0;u<interpoly_Zsize;u++)	//deg_z
				{
					delta_temp1=0;

					for(v=0;v<interpoly_Ysize;v++)	//deg_y
					{
						for(z=0;z<interpoly_Xsize;z++)	//deg_x
							if(g[j][u][v][z]!=0)	//coputation num concerning coefficients
							{
								delta_temp = mul( power(interpoint[i][0],z),power(interpoint[i][1],v) );
								delta_temp=mul(delta_temp, g[j][u][v][z]);
								//Hasse derivative mapping
								delta_temp1=add( delta_temp1,delta_temp );
							}
					}

					if(u==0)	//deg_z==0
					{
						delta[j] = delta_temp1;
					}
					else if(u>0)	//deg_z>0
					{
						delta_temp1 = mul( delta_temp1,power(interpoint[i][2],u) );
						delta[j] = add( delta[j],delta_temp1 );
					}

				}

				if(delta[j]!=0)
				{
					J[j]=1;	//record those polynomial with a nonzero hasse derivative mapping
					lod_min=lod[j];
					j_min=j;
				}
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

#ifdef _PrintLCCLod_
		printf("\nthe %d iter lod = {%d\t%d\t%d\t%d}", i, lod[0], lod[1], lod[2], lod[3]);
		printf("\nj_min = %d", j_min);
		printf("\nJ =\t");
		for (u = 0; u < init_polyNum; u++)
			if (J[u] != 0)
				printf("%d\t", u);
		printf("\n");
#endif

		if(j_min!=-1)
		{
			//f=g[j_min]
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						f[u][v][z]=g[j_min][u][v][z];
	
			//Modify nonzero polynomials
			for(j=0;j<init_polyNum;j++)	//num of polys
			{
				if(J[j]==1)
				{
					if(j!=j_min)
					{
						//delta*g_k+delta_k*f
						for(u=0;u<interpoly_Zsize;u++)	//rs
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
								for(z=0;z<interpoly_Xsize;z++)	//w+1
								{
									if(	g[j][u][v][z]!=0 )
									{
										temp1 = mul(delta[j_min],g[j][u][v][z]);
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
										g[j][u][v][z]=add(temp1,temp2);	
									}
									else
										g[j][u][v][z] = 0;
								}
					}
					else if(j==j_min)
					{
						for(u=0;u<interpoly_Zsize;u++)	//rs
						{
							for(v=0;v<interpoly_Ysize_1;v++)	//max(deg_y)+1
								for(z=0;z<interpoly_Xsize_1;z++)	//w+2
								{
									g1[u][v][z]=0;
								}

							for(v=0;v<interpoly_Ysize;v++)
								for(z=0;z<interpoly_Xsize;z++)	//w+1
								{
									g2[u][v][z]=0;
								}
						}
						
						//g1=x*f
						for(u=0;u<interpoly_Zsize;u++)	//rs
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
								for(z=0;z<interpoly_Xsize;z++)	//w+1
									if(g[j][u][v][z]!=0)
									{
										g1[u][v][z+1]=g[j][u][v][z];
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
						for(u=0;u<interpoly_Zsize;u++)	//rs
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
								for(z=0;z<interpoly_Xsize;z++)	//w+1
									if(g[j][u][v][z]!=0)
									{
										g2[u][v][z] = mul(interpoint[i][0],g[j][u][v][z]);
									}
						//g=g1+g2
						for(u=0;u<interpoly_Zsize;u++)	//rs
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
								for(z=0;z<interpoly_Xsize;z++)	//w+1	
									if( g1[u][v][z]!=0 || g2[u][v][z]!=0 )
									{
										g[j][u][v][z] = add(g1[u][v][z],g2[u][v][z]);
									}
									else
										g[j][u][v][z] = 0;
					}
				}
			}
		}



		flag_addNum = 0;
		flag_mulNum = 0;

		//debug: dectect the Q_com_elem[i] who has factor (x+a_i)
//		DectectIfFactorInPoly(g, init_polyNum, interpoint[0][i]);


/*	//****debug**********
	//Calculate each polynomial's leading order
	for(j=0;j<init_polyNum;j++)	//num_poly
	{
		lod_temp=0;
		lod[j]=0;
		for(u=0;u<interpoly_Zsize;u++)	//rs
		{
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
			{
				for(z=0;z<interpoly_Xsize;z++)	//w+1
				{
					if(g[j][u][v][z]!=0)
					{
						lod_temp=mono_order[u][v][z];
						if(lod_temp>lod[j])
							lod[j]=lod_temp;
					}
				}
			}
		}
	}
*/

	//int temp_x, temp_y, temp_z,temp;

	//for(j=0;j<init_polyNum;j++)   
	//{
	//	temp=0;
	//	for(u=0;u<interpoly_Zsize;u++)	//rs
	//		for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
	//			for(z=0;z<interpoly_Xsize;z++)	//w+1
	//				if(g[j][u][v][z]!=0)
	//				{
	//					temp_x = power(interpoint[i][0], z);
	//					temp_y = power(interpoint[i][1], v);
	//					temp_z = power(interpoint[i][2], u);
	//					temp= add( temp, mul( g[j][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
	//				}
	//	
	//	if( (lod[j]<=iterNum) && temp!=0 )
	//		printf("\ng[%d] with lod[%d]=%d in point[%d](%d,%d,%d) = %d", j, j, lod[j], i, interpoint[i][0], interpoint[i][1], interpoint[i][2], temp);

	//}
	//printf("\n");
	//*************************

	}


}

void uncom_elem_interpolation(int interpoint[4],int inGroup[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int outGroup1[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int outGroup2[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize])
{
	int i, j, u, v, z, J[2][init_polyNum], act[init_polyNum], lod_temp, lod[init_polyNum], lod_min[2], j_min[2];	//(delta, J, act, lod)[num of polys]
	int g1[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize+1], g2[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];
	int inGroup_temp[init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];
	int c[init_polyNum][lm+1],result_temp[init_polyNum][2];
	int index_min,poly_temp1,poly_temp2;

	//initialization
	for(j=0;j<init_polyNum;j++)
		for(i=0;i<2;i++)
		{	
			result_temp[j][i]=0;
		}

	for(j=0;j<init_polyNum;j++)
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1	
				{
					inGroup_temp[j][u][v][z]=inGroup[j][u][v][z];	
				}

	//Interpolation

	//Calculate each polynomial's leading order
	for(j=0;j<init_polyNum;j++)	//num_poly
	{
		lod_temp=0;
		lod[j]=0;
		for(u=0;u<interpoly_Zsize;u++)	//rs
		{
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
			{
				for(z=0;z<interpoly_Xsize;z++)	//w+1
				{
					if(inGroup_temp[j][u][v][z]!=0)
					{
						lod_temp=mono_order[u][v][z];
						if(lod_temp>lod[j])
							lod[j]=lod_temp;
					}
				}
			}
		}
	}
		
		
	//Initialise the eliminator array act[num_poly]
#ifndef _NoReductionUncom_
	for(j=0;j<init_polyNum;j++)	//num_poly
	{
		if(lod[j]<=iterNum)	//C=n when multiplicity = 1
			act[j]=1;
		else
			act[j]=0;
	}
#else
	for(j=0;j<init_polyNum;j++)
		act[j]=1;
#endif
	
	flag_addNum = 1;
	flag_mulNum = 1;

   

	//Calculate the hasse derivative mapping of each polynomials
	for(j=0;j<init_polyNum;j++)
		if(act[j]==1)
			for(u=0;u<interpoly_Zsize;u++)
			{
				c[j][u]=0;
	
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
					if(inGroup_temp[j][u][v][z]!=0)
					{
						poly_temp1 = mul( power(interpoint[0],z),power(interpoint[1],v) );
						poly_temp1 = mul( inGroup_temp[j][u][v][z],poly_temp1 );
						c[j][u] = add( c[j][u],poly_temp1 );
					//	c[j][u] = add( c[j][u],mul( inGroup_temp[j][u][v][z],mul( power(interpoint[0],z),power(interpoint[1],v) ) ) );
					}
			}
	
	for(j=0;j<init_polyNum;j++)
		if(act[j]==1)
			for(i=0;i<2;i++)
			{
				result_temp[j][i]=0;
		
				for(u=0;u<interpoly_Zsize;u++)
				{
					poly_temp1 = power(interpoint[i+2],u );
					poly_temp1 = mul( c[j][u],poly_temp1 );
					result_temp[j][i] = add( result_temp[j][i],poly_temp1 );
				//	result_temp[j][i] = add( result_temp[j][i],mul( c[j][u],power(interpoint[i+2],u) ) );
				}						
			}
	
	for(i=0;i<2;i++)
	{
		j_min[i]=-1;
		for(j=0;j<init_polyNum;j++)
		{
			J[i][j]=0;
			if( act[j]==1 && result_temp[j][i]!=0 )
			{
				J[i][j]=1;	//record those polynomial with a nonzero hasse derivative mapping
				lod_min[i]=lod[j];
				j_min[i]=j;
			}
		}
	}
	//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
	for(i=0;i<2;i++)
		for(j=0;j<init_polyNum;j++)	//num_polys
		{
			if(J[i][j]==1 && lod[j]<lod_min[i])
			{
				lod_min[i]=lod[j];
				j_min[i]=j;
			}
		}
	//printf("\nj_min=%d\n", j_min);

	//initialize outGroup1 and outGroup2
	for(j=0;j<init_polyNum;j++)
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1	
				{
					outGroup1[j][u][v][z]=inGroup_temp[j][u][v][z];
					outGroup2[j][u][v][z]=inGroup_temp[j][u][v][z];
				}

	//update the poly of outGroup1 
	if(j_min[0]!=-1)
	{
		index_min=j_min[0];	//index_min = j'

		//Modify nonzero polynomials
		for(j=0;j<init_polyNum;j++)	//num of polys
		{
			if(J[0][j]==1)
			{
				if(j!=j_min[0])
				{
					//delta*g_k+delta_k*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
							{
								if(inGroup_temp[j][u][v][z]!=0)
								{
									poly_temp1 = mul( result_temp[index_min][0],inGroup_temp[j][u][v][z] );
								}
								else
									poly_temp1 = 0;

								if(inGroup_temp[index_min][u][v][z]!=0)
								{
									poly_temp2 = mul( result_temp[j][0],inGroup_temp[index_min][u][v][z] ) ;
								}
								else
									poly_temp2 = 0;

								if( poly_temp1!=0 || poly_temp2!=0 )
								{
									outGroup1[j][u][v][z] = add(poly_temp1,poly_temp2);
								}
								else
									outGroup1[j][u][v][z] = 0;

							//	poly_temp1 = mul( result_temp[index_min][0],inGroup_temp[j][u][v][z] );
							//	poly_temp2 = mul( result_temp[j][0],inGroup_temp[index_min][u][v][z] ) ;
							//	outGroup1[j][u][v][z]=add(poly_temp1,poly_temp2);
							//	outGroup1[j][u][v][z]=add(2,2);
							}
				}
				else if(j==j_min[0])
				{
					for(u=0;u<interpoly_Zsize;u++)	//rs
					{
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
						{
							for(z=0;z<(interpoly_Xsize+1);z++)	//w+2
								g1[u][v][z]=0;
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								g2[u][v][z]=0;
						}
					}
					
					//g1=x*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								g1[u][v][z+1]=inGroup_temp[index_min][u][v][z];

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
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								if(inGroup_temp[index_min][u][v][z]!=0)
								{
									g2[u][v][z]=mul( interpoint[0],inGroup_temp[index_min][u][v][z]);
								}
					//g=g1+g2
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								if( g1[u][v][z]!=0 || g2[u][v][z]!=0 )
								{
									outGroup1[index_min][u][v][z]=add(g1[u][v][z],g2[u][v][z]);
								}
								else
									outGroup1[index_min][u][v][z] = 0;
				}
			}
		}
	}

	//update the poly of outGroup2 
	if(j_min[1]!=-1)
	{
		index_min = j_min[1];

		//Modify nonzero polynomials
		for(j=0;j<init_polyNum;j++)	//num of polys
		{
			if(J[1][j]==1)
			{
				if(j!=j_min[1])
				{
					//delta*g_k+delta_k*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
							{
								if(inGroup_temp[j][u][v][z]!=0)
								{
									poly_temp1 = mul( result_temp[index_min][1],inGroup_temp[j][u][v][z] );
								}
								else 
									poly_temp1 = 0;

								if(inGroup_temp[index_min][u][v][z]!=0)
								{
									poly_temp2 = mul( result_temp[j][1],inGroup_temp[index_min][u][v][z] );
								}
								else
									poly_temp2 = 0;

								if( poly_temp1!=0 || poly_temp2!=0 )
								{
									outGroup2[j][u][v][z]=add( poly_temp1,poly_temp2 );	
								}
								else
									outGroup2[j][u][v][z] = 0;

							//	poly_temp1 = mul( result_temp[index_min][1],inGroup_temp[j][u][v][z] );
							//	poly_temp2 = mul( result_temp[j][1],inGroup_temp[index_min][u][v][z] );
							//	outGroup2[j][u][v][z]=add( poly_temp1,poly_temp2 );	
							}
				}
				else if(j==j_min[1])
				{
					for(u=0;u<interpoly_Zsize;u++)	//rs
					{
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
						{
							for(z=0;z<(interpoly_Xsize+1);z++)	//w+2
								g1[u][v][z]=0;
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								g2[u][v][z]=0;
						}
					}
					
					//g1=x*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								g1[u][v][z+1]=inGroup_temp[index_min][u][v][z];

					//convert x^w+1=y^w+y, difference with diff code
					for(u=0;u<interpoly_Zsize;u++)	//rs
					{
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
						{
							if(g1[u][v][w+1]!=0)
							{
								g1[u][v+1][0]=add(g1[u][v+1][0],g1[u][v][w+1]);
								g1[u][v+w][0]=add(g1[u][v+w][0],g1[u][v][w+1]);
								g1[u][v][w+1]=0;
							}
						}
					}

					//g2=xi*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								if(inGroup_temp[index_min][u][v][z]!=0)
								{
									g2[u][v][z]=mul( interpoint[0],inGroup_temp[index_min][u][v][z]);
								}

					//g=g1+g2
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								if( g1[u][v][z]!=0 || g2[u][v][z]!=0 )
								{
									outGroup2[index_min][u][v][z]=add(g1[u][v][z],g2[u][v][z]);
								}
								else
									outGroup2[index_min][u][v][z] = 0;
				}
			}
		}
	}

  

	flag_addNum = 0;
	flag_mulNum = 0;
	
	//*******debug*********
	//outGroup1
	//for(j=0;j<init_polyNum;j++)
	//{
	//	lod_temp=0;
	//	lod[j]=0;
	//	for(u=0;u<interpoly_Zsize;u++)	//rs
	//	{
	//		for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
	//		{
	//			for(z=0;z<interpoly_Xsize;z++)	//w+1
	//			{
	//				if(outGroup1[j][u][v][z]!=0)
	//				{
	//					lod_temp=mono_order[u][v][z];
	//					if(lod_temp>lod[j])
	//						lod[j]=lod_temp;
	//				}
	//			}
	//		}
	//	}
	//}

	//int temp_x,temp_y,temp_z,temp;

	//for(j=0;j<init_polyNum;j++)
	//{
	//	temp=0;
	//	for(u=0;u<interpoly_Zsize;u++)	//rs
	//		for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
	//			for(z=0;z<interpoly_Xsize;z++)	//w+1
	//				if(outGroup1[j][u][v][z]!=0)
	//				{
	//					temp_x=power( interpoint[0],z );
	//					temp_y=power( interpoint[1],v );
	//					temp_z=power( interpoint[2],u );
	//					temp= add( temp, mul( outGroup1[j][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
	//				}
	//	if (temp)
	//		printf("\noutGroup1[%d] with lod[%d]=%d in point[7](%d,%d,%d) = %d",j,j,lod[j],interpoint[0],interpoint[1],interpoint[2],temp);
	//}
	//printf("\n");

	//outGroup2
	//for(j=0;j<init_polyNum;j++)
	//{
	//	lod_temp=0;
	//	lod[j]=0;
	//	for(u=0;u<interpoly_Zsize;u++)	//rs
	//	{
	//		for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
	//		{
	//			for(z=0;z<interpoly_Xsize;z++)	//w+1
	//			{
	//				if(outGroup2[j][u][v][z]!=0)
	//				{
	//					lod_temp=mono_order[u][v][z];
	//					if(lod_temp>lod[j])
	//						lod[j]=lod_temp;
	//				}
	//			}
	//		}
	//	}
	//}
	
	//for(j=0;j<init_polyNum;j++)
	//{
	//	temp=0;
	//	for(u=0;u<interpoly_Zsize;u++)	//rs
	//		for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
	//			for(z=0;z<interpoly_Xsize;z++)	//w+1
	//				if(outGroup2[j][u][v][z]!=0)
	//				{
	//					temp_x=power( interpoint[0],z );
	//					temp_y=power( interpoint[1],v );
	//					temp_z=power( interpoint[3],u );
	//					temp= add( temp, mul( outGroup2[j][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
	//				}
	//	if (temp)
	//		printf("\noutGroup2[%d] with lod[%d]=%d in point[7](%d,%d,%d) = %d",j,j,lod[j],interpoint[0],interpoint[1],interpoint[3],temp);
	//}
	//printf("\n");
	//******************

}

void factorisation(int Q_input[test_vec_num_max][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int output_list[test_vec_num_max][lm + 1][k], int list_num[test_vec_num_max])
{
	int i, j, u, v, z;

	//initialization
	for (i = 0; i < test_vec_num_max; i++)
		for (u = 0; u < lm + 1; u++)
			for (v = 0; v < k; v++)
				output_list[i][u][v] = -1;
	
	for(i=0;i<test_vec_num;i++)
	{
		//Initialisation
		for(j=0;j<k;j++)	//number of fac steps=k
			for(u=0;u<facpoly_Zsize;u++)	//rs
				for(v=0;v<facpoly_Ysize;v++)	//y_size
					for(z=0;z<facpoly_Xsize;z++)	//w+1
						Q[j][u][v][z]=0;

//		for(u=0;u<k;u++)	//number of fac steps=k
//			for(v=0;v<lm+1;v++)	//5>rs=expected number of roots
//				rootlist[u][v]=-1;	

		for(u=0;u<lm+1;u++)	//5>(rs-1)=expected number of output lists
			for(v=0;v<k;v++)	//k
				output[u][v]=-1;

		//Initialisation of factorisation
		uu=0;	//recursive deduction index
		l=0;	//candidate output index
		//q_0(z)=itp(z)
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					Q[uu][u][v][z]=Q_input[i][u][v][z];


		flag_addNum = 1;
		flag_mulNum = 1;



		//recursive coefficient search
		rcs(uu);

  
			
		flag_addNum = 0;
		flag_mulNum = 0;


		//store the necessary data
		for(u=0;u<lm+1;u++)	//5>(rs-1)=expected number of output lists
			for(v=0;v<k;v++)	//k
				output_list[i][u][v]=output[u][v];

		list_num[i]=l;

	}

}

//this rcs() only can be used in lm=1
void rcs(int uu)
{
	int i, j, u, v, m, z, t, r, i_1, j_1, i_2, j_2, a, b, leadMono, leadMono_temp, alpha, act, temp;
	int lc[lm+1], q_temp[lm+1][rcspoly_Ysize+w][rcspoly_Xsize];	//q_temp[z-deg+1][y_size>max(deg_y)+1+(rs-1)*(max(deg_y) in encoding functions)][14>w+(rs-1)*w], lc[rs]--leading coefficient polynomial
	int index_flag;
	int rcspoly_Ysize_1,rcspoly_Xsize_1;	//the q_temp size of the first step of factorization
	int rootlist[lm+1];	//the list of roots [number of fac steps=k][expected number of roots in each solution, 5>rs]

	//array size initialization
	rcspoly_Ysize_1= faiMax_Ysize + facpoly_Ysize + 1;	
	rcspoly_Xsize_1= faiMax_Xsize + lm*facpoly_Xsize + 1;	//ensure the Xsize can let the q_temp mod curve H_w normally

	if( (rcspoly_Ysize_1>=rcspoly_Ysize) || (rcspoly_Xsize_1>=rcspoly_Xsize) )
	{
		printf("\n fac size.0 has error\n");
	}

	leadMono=0; leadMono_temp=0;	//leading monomial index
	act=0;	//judge value for recursive search of each f_k-1-u


	//initialization
	for(v=0;v<lm+1;v++)	//5>rs=expected number of roots
		rootlist[v]=-1;	

	//find pb_k-1-uu
	j_1 = tg_order[k-1-uu][1];
	i_1 = tg_order[k-1-uu][0];


	//initialise q_temp
	//for(i=0;i<lm+1;i++)	//rs
	//	for(j=0;j<rcspoly_Ysize;j++)	
	//		for(u=0;u<rcspoly_Xsize;u++)	
	//			q_temp[i][j][u]=0;
	memset(q_temp,0,sizeof(q_temp));

	//Calculate q_temp[pb_k-1-u]=q[u][pb_k-1-u]
	for(i=0;i<facpoly_Zsize;i++)	//rs
	{
		for(j=0;j<facpoly_Ysize;j++)	//y_size=max(deg_y) in encoding function*(rs-1), and 31>=max(deg_y)+1
		{
			for(u=0;u<facpoly_Xsize;u++)
				if( Q[uu][i][j][u]!=0 )
				{
					q_temp[i][j+i*j_1][u+i*i_1] = Q[uu][i][j][u];	//this time, q_temp size is [rs+1][maxY(j_1)+maxY(Q)][maxX(i_1)+maxX(Q)]

					if( (j+i*j_1) > rcspoly_Ysize )
					{
						printf("\n fac size.1 has error__Y\n");
					}
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
	//printf("The leading monomial in q[u][pb_k-1-u]");

	//find the leading coefficient polynomial lc[rs]
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

	//printf("leading coefficient polynomial");
	
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
		if(b==0)
		{
			act=1;
			rootlist[u]=root[i];
			u++;
		}
	}		

	//For each distinct root of rootlist[u]
	if(act==1)	//act==1 means there is at least one root in rootlist[u];
	{
		for(i=0;i<lm+1;i++)	//2>rs
		{
			if(rootlist[i]!=-1)
			{	
				alpha=rootlist[i];
				
				output[l][k-1-uu]=alpha;	//output[l][k-1-uu]
					
				//printf("current outputs");

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
							//for(u=0;u<lm+1;u++)	//rs
							//	for(m=0;m<rcspoly_Ysize+w;m++)	//y_size
							//		for(z=0;z<rcspoly_Xsize;z++)	//w+1
							//			q_temp[u][m][z]=0;
							memset(q_temp,0,sizeof(q_temp));

							//calculate (z+f_k-1-u*pb_k-1-u)^j
							polyexp1(alpha, i_1, j_1, j, expoly);
							
							//calculate q_temp
							for(u=0;u<j+1;u++)	//rs
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

														q_temp[u][v+m][t+z] = add( q_temp[u][v+m][t+z],temp ); 													}
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
											q_temp[u][m+w][temp] = add( q_temp[u][m+w][temp],q_temp[u][m][index_flag] );	//y^(w+1), may be overflow!!
											q_temp[u][m][index_flag] = 0;
										}
									index_flag = index_flag-1;
								}

							}

							//Q[r] = add( Q[r],q_temp )
							for(u=0;u<lm+1;u++)	//rs
								for(m=0;m<facpoly_Ysize;m++)	//y_size
									for(z=0;z<facpoly_Xsize;z++)	//w+1
										if( Q[r][u][m][z]!=0 || q_temp[u][m][z]!=0 )
										{
											Q[r][u][m][z]=add( Q[r][u][m][z],q_temp[u][m][z] );
										}

						}
				
					//next coefficient searching
					
					rcs(r);
					//make Q[r] set 0, let other branch use it
					//for(u=0;u<facpoly_Zsize;u++)	//rs
					//	for(m=0;m<facpoly_Ysize;m++)	//y_size
					//		for(z=0;z<facpoly_Xsize;z++)	//w+1
					//		{
					//			Q[r][u][m][z]=0;
					//		}
					memset(Q[r],0,sizeof(Q[r]));

				}
			}
		}
	}

}

void polyexp1(int c, int i, int j, int deg_z, int poly[][expoly_Ysize][expoly_Xsize])
{
	int u, m, z;	//p1[rs+1][18>(max(deg_y) in encoding functions)*(rs-1)][10>(max(deg_x) in encoding functions)*(rs-1)], p2[rs][26=18+max(deg_y) in encoding functions][14=10+max(deg_x) in encoding functions]
	int poly_temp[lm+1+1][expoly_Ysize+faiMax_Ysize+w][expoly_Xsize+faiMax_Xsize];
	int temp, temp_Ysize, temp_Xsize;
	int index_flag;

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
			for(m=0;m<temp_Ysize+w;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
				for(z=0;z<temp_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
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
		{
			for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
				for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
					if( poly[u][m][z]!=0 )
					{
						temp = mul( c,poly[u][m][z] );
						poly_temp[u][j+m][i+z] = add( poly_temp[u][j+m][i+z],temp );
					}

			//convert x^w+1=y^w+y, difference with diff code
			index_flag = temp_Xsize-1;
			while( index_flag>w )
			{
				temp = index_flag-(w+1);
				for(m=0;m<temp_Ysize;m++)
					if( poly_temp[u][m][index_flag]!=0 )
					{
						poly_temp[u][m+1][temp] = add( poly_temp[u][m+1][temp],poly_temp[u][m][index_flag] );
						poly_temp[u][m+w][temp] = add( poly_temp[u][m+w][temp],poly_temp[u][m][index_flag] );
						poly_temp[u][m][index_flag] = 0;
					}
				index_flag = index_flag-1;
			}

		}

		for(u=0;u<lm+1;u++)	//rs
			for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
				for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
				{
					poly[u][m][z] = poly_temp[u][m][z];
				}
	}
}

void choose(int output_codeword[n], int output_bicodeword[n*p], int output_list[test_vec_num_max][lm + 1][k], int list_num[test_vec_num_max], int flag_decoding_alg)
{
	int i, j, u, v, z, value, flag, min_index_1, min_index_2;
	unsigned int mask=1;
//	float proba[lm*test_vec_num_max], proba_temp, temp;
	double proba[lm*test_vec_num_max], proba_temp, temp;	
	int codeword_temp[n];

	flag_addNum = 1;
	flag_mulNum = 1;
	//normal mode
	//Initialise hamming distance counter
	for(u=0;u<lm*test_vec_num_max;u++)
		proba[u]=-1.0;

	proba_temp=0.0;
	flag = 0;
	min_index_1 = -1;
	min_index_2 = -1;

	for(i=0;i<test_vec_num;i++)
		if( list_num[i]!=0 )
		{
			for(j=0;j<list_num[i];j++)
			{
				//reencoding
				encoder(output_list[i][j],codeword_temp);
				//calculate the posteriori probablity
				temp = 1.0;
				for(u=0;u<n;u++)
				{
					for(v=0;v<q;v++)
						if(codeword_temp[u]==root[v])
							temp = temp*RM[v][u];
				}

				if( proba_temp < temp )	//< or <=
				{
					proba_temp = temp;
					min_index_1 = i;
					min_index_2 = j;

					for(v=0;v<n;v++)
						output_codeword[v] = codeword_temp[v];

					flag = 1;	//exist at less one valid solution
				}
			}
		}

	//output the decoding result 
	if(flag==0)	// not exist a valid solution
	{
		for(u=0;u<n;u++)
			output_codeword[u]=large_vec[0][u];
		
		//nonbinary --> binary
		for(u=0;u<n;u++)
		{	
			value=output_codeword[u];
			mask=1;
			for(v=0;v<p;v++)
			{
				if((value & mask)>0)
					output_bicodeword[p*u+v]=1;
				else
					output_bicodeword[p*u + v] = 0;
				mask=mask<<1;
			}
		}			
	}
	else if(flag==1)	//exist a valid solution
	{

		//encode message[min_index]
//		encoder(output_list[min_index_1][min_index_2],output_codeword);
		//nonbinary --> binary
			for(u=0;u<n;u++)
			{	
				value=output_codeword[u];
				mask=1;
				for(v=0;v<p;v++)
				{
					if((value & mask)>0)
						output_bicodeword[p*u + v] = 1;
					else
						output_bicodeword[p*u + v] = 0;
					mask=mask<<1;
				}
			}
	}
	flag_addNum = 0;
	flag_mulNum = 0;
	//*********************



	//jugde CWR
	int count_temp;
	//calculate ChosenWrong_SeqNum
	if(flag)	//flag: 1--> valid solution, 0--> invalid solution
	{
		count_temp=0;
		for(u=0;u<k;u++)
			if(output_list[min_index_1][min_index_2][u]!=message[u])
			{
				++count_temp;	
			}

		if(count_temp)	//count_temp: 0-->correct, >0-->incorrect
		{
			++ChosenWrong_SeqNum;
		}
	}

	//calculate DecSucc_SeqNum
	flag=0;
	for(i=0;i<test_vec_num;i++)
	{
		for(j=0;j<lm+1;j++)
		{
			count_temp=0;
			for(u=0;u<k;u++)
				if(output_list[i][j][u]==message[u])
				{
					count_temp++;
				}

			if(count_temp==k)
			{
				flag=1;	//there is a correct answer at least
				DecSucc_SeqNum++;	//this frame has decoded successful
				break;
			}
		}
	
		if(flag==1)
		{	
			break;
		}
		else if(flag==0)
		{
			flag=0;
		}
	}
	//*****************************
	
	//**********debug*************
	int epcount2=0;
	for(u=0;u<n;u++)
		if(codeword[u]!=output_codeword[u])
			epcount2++;

	if (epcount1 <= able_correct && epcount2 != 0)	//this seq_num has chosen the wrong one
	{
		ChosenWrong_SeqNum++;
		//		printf("\n\nChoice error\n\n");
	}

	//******debug*******
	//caculate the outputlistNum
	outputlistNum = 0;
	for (i = 0; i < test_vec_num; i++)
		outputlistNum += list_num[i];

	//caculate the degree_test[i]
	for (i = 0; i < test_vec_num; ++i)
	{
		int mono_temp = -1;
		degree_test[i] = 0;
		for (u = 0; u<interpoly_Zsize; u++)	//rs
			for (v = 0; v<interpoly_Ysize; v++)	//max(deg_y)+1
				for (z = 0; z<interpoly_Xsize; z++)	//w+1
					if (Q_interpoly[i][u][v][z] != 0 && mono_temp < mono_order[u][v][z])
					{
						mono_temp = mono_order[u][v][z];
						degree_test[i] = z*w + v*(w + 1) + u*weiz;
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
		for (int j = n - 1; j >= n - eta; j--)
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

	if (flag_decoding_alg==0)
	{
		for(i=0;i<test_vec_num;i++)
			if( (n-testCount1[i])>degree_test[i] )
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

				if (flag == 0)
				{


					printf("\n\nseq_num_Now=%d, No.%d test vector factorization has failed!!, testCount1=%d, testCount2=%d, degreeOfPolynomail=%d, list_num=%d\n", seq_num_Now, i, testCount1[i], testCount2, degree_test[i], list_num[i]);
					//					printf("\nNo.%d test vector has error!!, testCount1=%d, testCount2=%d, testCount1_uncom=%d, list_num=%d\n\n", i, testCount1[i], testCount2[i], (testCount1[i]-testCount1_com), list_num[i]); 
					//					printf("\nseq_num_Now=%d, This sequence is decoding failed， and error num is %d!\n\n", seq_num_Now, epcount1);
					//********debug*************
					int temp_x, temp_y, temp_z;
					int temp, h, z;

					//**************
					printf("\n");
					for (u = 0; u < interpoly_Zsize; u++)	//rs
						for (v = 0; v < interpoly_Ysize; v++)	//max(deg_y)+1
							for (z = 0; z < interpoly_Xsize; z++)	//w+1
								if (Q_interpoly[i][u][v][z] != 0)
								{
									printf("Q_interpoly[%d][%d][%d][%d]=%d\n", i, u, v, z, Q_interpoly[i][u][v][z]);
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
					printf("x[]\t\t");
					for (i = 0; i < n; ++i)
						printf("\t(%d,%d)", point[i][0], point[i][1]);

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
						printf("\t(%d,%d)", x_ordered[0][i], x_ordered[1][i]);


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

			//*******************
	}

}

void cheatFac(void)
{
	int i, u, v, value;
	unsigned int mask = 1;
	int flag_judge;	//used to judge if there is correct codeword 

	//start
	flag_judge = 0;
	for (i = 0; i<test_vec_num; i++)
		if ((n - testCount1[i])>degree_test[i])
		{
			flag_judge = 1;
			for (u = 0; u<n; u++)
			{
				dec_codeword[u] = codeword[u];
			}

			////nonbinary --> binary
			//for (u = 0; u<n; u++)
			//{
			//	value = dec_codeword[u];
			//	mask = 1;
			//	for (v = 0; v<p; v++)
			//	{
			//		if ((value & mask)>0)
			//			dec_bicodeword[p*u + v] = 1;
			//		else
			//			dec_bicodeword[p*u + v] = 0;
			//		mask = mask << 1;
			//	}
			//}

			break;
		}

	if (!flag_judge)	//flag_judge==0 means that there is not correct codeword
	{
		for (u = 0; u<n; u++)
		{
			dec_codeword[u] = large_vec[0][u];
		}

		////nonbinary --> binary
		//for (u = 0; u<n; u++)
		//{
		//	value = dec_codeword[u];
		//	mask = 1;
		//	for (v = 0; v<p; v++)
		//	{
		//		if ((value & mask)>0)
		//			dec_bicodeword[p*u + v] = 1;
		//		else
		//			dec_bicodeword[p*u + v] = 0;
		//		mask = mask << 1;
		//	}
		//}
	}

}

void cheatDec(void)
{
	int i, u, v, value, flag_judge;
	unsigned int mask = 1;
	int genius, errorCorrectionNum;

	//calculate the error correction ability
	genius = w*(w - 1) / 2;
	errorCorrectionNum = (n - k - genius) / 2;

	//start judgment
	flag_judge = 0;		//used to judge if there is correct codeword, 0-->No, 1-->Yes
	for (i = 0; i<test_vec_num; i++)
		if (testCount1[i] <= errorCorrectionNum)
		{
			flag_judge = 1;
			for (u = 0; u<n; u++)
			{
				dec_codeword[u] = codeword[u];
			}

			//nonbinary --> binary
			//for (u = 0; u<n; u++)
			//{
			//	value = dec_codeword[u];
			//	mask = 1;
			//	for (v = 0; v<p; v++)
			//	{
			//		if ((value & mask)>0)
			//			dec_bicodeword[p*u + v] = 1;
			//		else
			//			dec_bicodeword[p*u + v] = 0;
			//		mask = mask << 1;
			//	}
			//}

			break;
		}

	if (!flag_judge)	//flag_judge==0 means that there is not correct codeword
	{
		for (u = 0; u<n; u++)
		{
			dec_codeword[u] = large_vec[0][u];
		}

		//nonbinary --> binary
		//for (u = 0; u<n; u++)
		//{
		//	value = dec_codeword[u];
		//	mask = 1;
		//	for (v = 0; v<p; v++)
		//	{
		//		if ((value & mask)>0)
		//			dec_bicodeword[p*u + v] = 1;
		//		else
		//			dec_bicodeword[p*u + v] = 0;
		//		mask = mask << 1;
		//	}
		//}
	}
}


void mono_table()
{
	int i, j, u, v, z, weight[weight_Zsize][weight_XYsize], mono_order_1[weight_Zsize][weight_XYsize];	//26=(100/weight(z))+1, 100=pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][max(deg_y)+1][w+1]

	//Initialisation
	for (i = 0; i<weight_Zsize; i++)
		for (j = 0; j<weight_XYsize; j++)
		{
			weight[i][j] = -1;
			mono_order_1[i][j] = -1;
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


	//construction 2-dimensional mono table
	z = 0;
	for (v = 0; v<mono_ordinSize; v++)	//for each possible weight until weight(term[rs-1][max(deg_y)+1][w-1])+1, note term[rs-1][max(deg_y)+1][w-1] is the term next to term[rs-1][max(deg_y)][w]
	{
		for (i = 0; i<weight_Zsize; i++)	//230/weight(z)
		{
			for (j = 0; j<weight_XYsize; j++)	//230>the index of the largest term with deg(z)=0, and weight=weight(term[rs-1][max(deg_y)+1][w-1])
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

#ifdef PrintfMonoTable
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
	printf("\n\n");

	//******************
#endif
}


void findpoints()
{
	int i, j, u, a1, a2, a3, x, y;

	//Initialisation
	for(i=0;i<2;i++)
		for(j=0;j<n;j++)	//w^3
			point[j][i]=0;

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
				point[u][0]=x;
				point[u][1]=y;
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
	index_temp = tg_size;	//according to the last formualtion to calculate out

	j = 0;
	for (i = 0; i<index_temp; i++)
		for (u = i; u >= 0; u--)
			if (u <= w)
			{
				tg_order[j][0] = u;	// 0<deg_x<=w
				tg_order[j][1] = i - u;	// 0<deg_y
				j++;
			}

}
