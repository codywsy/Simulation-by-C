#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "FiniteFieldBasisGF(64).h"
#include "FiniteFieldBasisCal.h"

//#define _GS_Normal_
//#define _Complexity_
//#define _NoReduction_
//#define _PolyCoeffNumUncom_
//#define _PolyCoeffNumFac_


#define w 8	//genius = w*(w-1)/2 = 28
#define weiz 180
#define lm 1
#define iterNum 512	//when m=1, C is equal to n
#define able_correct 164
#define pointNum 2
#define interval 0.3

//conditional compile
#ifndef _GS_Normal_
	#define eta 5
	#define test_vec_num 32 //the num of test_vec is 2^eta
#else
	#define eta 0
	#define test_vec_num 1 //the num of test_vec is 2^eta
#endif

//test_vec_construction
#define choose_num 2	//fixed

//tgorder()
#define tg_size 4447	//more than k, tg_size represent the probably used pole basis num, (w+1)*w/2 + (w+1) * ( (w+interpoly_Ysize)-w+1 ) 
//mono_table()
#define weight_Zsize 27
#define weight_XYsize 4447	//equals to the tg_size
#define mono_ordinSize (2*4480)	//choose the value with debugging
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 490	//large than interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize 4447 //>index of term[max(deg_y)][w] in the pole basis + 1, nearly equals to tg_size
//polebasis()


//interpolation()
#define init_polyNum (w*(lm+1))	//change with diff lm, the poly num of the init polyGroup
#define interpoly_Zsize (lm+1)	//maxValue of z is lm=1, so the Zsize should add 1 more.
#define interpoly_Ysize 460	//maxdeg of y is (w-1) + w*(n/(w+1)), so the Ysize should add 1 more.
#define interpoly_Xsize	(w+1)	//maxdeg of x is w, so the Xsize should add 1 more.
//factorization()
#define facpoly_Zsize (lm+1)// same as the interpoly_Zsize
#define facpoly_Ysize 460		// same as the interpoly_Ysize
#define facpoly_Xsize (w+1)	// same as the interpoly_Xsize
//rcs()
#define rcspoly_Ysize 489	//degY+1, degY = interpoly_Ysize + max(j_1) + w*( (max(i_1)+w)/(w+1) )
#define rcspoly_Xsize 17	//degX+1, degX = max(i_1) + w
#define faiMax_Ysize 20	//change with diff w and k, the max degY of probably used polebasis(polebasis_k-1)
#define faiMax_Xsize 8	//change with diff w and k, the max degX of probably used polebasis(polebasis_k-1)
//expoly()-->expanded polynomial
#define expoly_Ysize (faiMax_Ysize+1)	
#define expoly_Xsize (faiMax_Xsize+1)


//main()
unsigned long int seq_num;	//number of input binary sequences
float SNR,N0;
double BER, FER;
float pi=3.141593;	// pai
//GF elements
//extern int mularray[];
//extern int root[];
//extern int logNum[];
//*************
int bi_message[k*p], message[k];	//transmitted messge[k]
int codeword[n], bi_codeword[n*p]; //codewords
float tx_symbol[p*n][2], rx_symbol[p*n][2], sgm;
float RM[q][n];
//polebasis()


//findpoint()
int point[2][n];	//rational points[2][w^3], [0]->x, [1]->y
//tgorder()
int tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
//mono_table()
int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
int gmatrix[k][n];	//generator matrix[k][n]
int dec_message[k],dec_codeword[n],dec_bicodeword[n*p]; //decoding result
//test_vec_construction
int large_vec[choose_num][n];
float test_set_ordered[2][n];
//int test_vec_com[n];
int test_vec_com[n-eta];
int x_ordered[2][n];
//interpolation
int bi_rxword[p*n], rxword[n]; //received wordss
int epcount1, epcount2, testCount1[test_vec_num], testCount2[test_vec_num], testCount1_com;
int Q_interpoly[test_vec_num][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];	//chosen interpolated polynomial [z-deg+1][max(deg_y)+1][w+1]
int uu;	//factorisation step index
int l, listNum[test_vec_num];	//candidate output index
//factorization() and rcs()
int Q[k][facpoly_Zsize][facpoly_Ysize][facpoly_Xsize];	//sequtial deduction polynomial [number of fac steps=k][rs][y_size][w+1], y_size> maxdeg_y]+rs*(deg_¦µ(k-1))
int rootlist[k][lm+1];	//the list of roots [number of fac steps=k][expected number of roots in each solution, 5>rs]
int output[lm+1][k], outputList[test_vec_num][lm+1][k];	//the list of candidate message [expeced number of candidate message, >rs][length of message, k]
int expoly[2][expoly_Ysize][expoly_Xsize];	//expanded polynomial in [z+f_k-1-u*pb_k-1-u]^rs, expoly[rs][3>(max(deg_y) in encoding functions)*(rs-1)][3>(max(deg_x) in encoding functions)*(rs-1)]

int Q_uncom_elem[test_vec_num][init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	//****debug**********
	int degree_test_temp[test_vec_num][init_polyNum];  //store the degree of leading monomial of polynomial
	int degree_test[test_vec_num];		//store the degree of the choosen polynomial

	//some var about choosing incorrectly
	unsigned long int seq_num_Now;	//record the seq_num at the moment
	unsigned long int DecSucc_SeqNum, ChosenWrong_SeqNum;	//record the seq_num should be decoding success in genius mode
	double CWR;	//record the rate of chosen uncorrectly

	//some var about computation complexity
	unsigned long int addNum, mulNum, addNum_temp, mulNum_temp;
	int flag_addNum=0, flag_mulNum=0;

	//*******************

void findpoints(void);
void mono_table(void);
void polyexp(int, int, int);
void polyexp1(int, int, int, int);
//void zerobasis(void);
//void polebasis(void);
void tgorder(void);
//void coefficient(void);
void test_vec_contruction(void);
void interpolation(void);
void com_elem_interpolation(int g[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int interpoint[][n-eta]);
void uncom_elem_interpolation(int interpoint[4],int inGroup[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int outGroup1[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int outGroup2[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize]);
void factorisation(void);
void rcs(int);
void choose(void);
void generator(void);
void encoder(int message_temp[], int codeword_temp[]);
void modulation(void);
void channel(void);
void demodulation(void);
float PDF(int,int);


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
	if((fp=fopen("LCC_Herm(512,153)_eta=5.txt","a"))==NULL)
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

			//LIST DECODER
			//test vector construction
			test_vec_contruction();
			//interpolation
			interpolation();

			//factorisation
			factorisation();

			//choose
			choose();

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

			if(ferror>209)
				break;
		}

		if(ferror>209)
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

		fp=fopen("LCC_Herm(512,153)_eta=5.txt","a");
		fprintf(fp,"Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f, Choice Errors=%2.1d, CWR=%E\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count, ChosenWrong_SeqNum, CWR);
		fclose(fp);

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
		for(int z=0;z<pointNum;z++)
		for(u=0;u<pointNum;u++)
		for(v=0;v<pointNum;v++)
		for(int x=0;x<pointNum;x++)
		for(int y=0;y<pointNum;y++)
		for(int h=0;h<pointNum;h++)
		{
			RM[j][i]=(Pr[0][h]*Pr[1][y]*Pr[2][x]*Pr[3][v]*Pr[4][u]*Pr[5][z]);  
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
	//construct the test_vec[test_vec_num]
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

//order test_vec
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

	//contruction test_vec
	for(i=0;i<n-eta;i++)
		test_vec_com[i]=large_vec[0][(int)test_set_ordered[1][i]];

	//construct x_ordered
	for(i=0;i<n;i++)
		for(j=0;j<2;j++)
			x_ordered[j][i]=point[j][(int)test_set_ordered[1][i]];

	//*****debug*********
	//caculate the num of diffs between codeword and test_vec[i] 
	testCount1[0]=0;
	for(i=0;i<n-eta;i++)
		if(codeword[(int)test_set_ordered[1][i]]!=test_vec_com[i])
			testCount1[0]++;

//	printf("\n\nfirst n-eta error num :%d\n\n",testCount1[0]);
	
	testCount1_com=testCount1[0];
	for(i=1;i<test_vec_num;i++)
		testCount1[i]=testCount1[0];

	//the num of loop equals to eta
	unsigned long int mask;
	for(i=0;i<test_vec_num;i++)
	{
		mask=1;
		for(int j=n-1;j>=n-eta-1;j--)
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
	//*******************


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
	
*/	//********************************		



}

void interpolation()
{
	int i,j,u,v,z,num,temp,temp_index;
	int com_elem_interpoint[3][(n-eta)];
	int Q_com_elem[init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];
//	int uncom_elem_interpoint[eta][4];
	int degree_temp[test_vec_num][init_polyNum];
//	int Q_uncom_elem[test_vec_num][init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	//common element interpolation

	//set common element interpoint (xi,ri)
	for(i=0;i<n-eta;i++)
	{
		com_elem_interpoint[0][i]=x_ordered[0][i];
		com_elem_interpoint[1][i]=x_ordered[1][i];
		com_elem_interpoint[2][i]=test_vec_com[i];
	}

/*	//*********debug***********
	printf("\n\ncom_elem_interpoint\n");
	for(i=0;i<n-eta;i++)
		printf("\t%d",com_elem_interpoint[i][0]);
	printf("\n");
	for(i=0;i<n-eta;i++)
		printf("\t%d",com_elem_interpoint[i][1]);
*/	//***************************

	for(i=0;i<init_polyNum;i++)	//num_polys
		for(j=0;j<interpoly_Zsize;j++)	//z-deg+1
			for(u=0;u<interpoly_Ysize;u++)	//max(deg_y)+1
				for(v=0;v<interpoly_Xsize;v++)	//w+1
					Q_com_elem[i][j][u][v]=0;

	com_elem_interpolation(Q_com_elem,com_elem_interpoint);
	
	//com_elem interpolation finish

#ifdef _Complexity_
	//****debug********
	printf("\n\nCint_1=%d\tmul=%d\tadd=%d\n\n", (addNum+mulNum), mulNum, addNum);
	addNum_temp = addNum;
	mulNum_temp = mulNum;
	//*****************
#endif

//conditional compile
#ifndef _GS_Normal_
	//uncommon element interpolation
	int uncom_elem_interpoint[eta][4];
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
	for(z=0;z<test_vec_num;z++)
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
		{
			uncom_elem_interpolation(uncom_elem_interpoint[i],Q_uncom_elem[u],Q_uncom_elem[2*u+0],Q_uncom_elem[2*u+1]);
		}
	}		

#ifdef _Complexity_
	//****debug********
	printf("\n\nCint_2=%d\tmul=%d\tadd=%d\n\n", (addNum+mulNum-addNum_temp-mulNum_temp), (mulNum-mulNum_temp), (addNum-addNum_temp));
	addNum_temp = addNum;
	mulNum_temp = mulNum;
	//*****************
#endif

	//choose the poly to factorization

	//initialize
	for(i=0;i<test_vec_num;i++)
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					Q_interpoly[i][u][v][z]=0;

	//eta>0
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
		//******************

		//assignment
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					Q_interpoly[i][u][v][z]=Q_uncom_elem[i][temp_index][u][v][z];
	}
	//******************************
#else
	//eta=0
	for(i=0;i<test_vec_num;i++)
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
	int temp_x, temp_y, temp_z;
	for(i=0;i<test_vec_num;i++)
	{
		for(j=0;j<n-eta;j++)
		{
			temp=0;
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						if(Q_interpoly[i][u][v][z]!=0)
						{
							temp_x=power( com_elem_interpoint[0][j],z );
							temp_y=power( com_elem_interpoint[1][j],v );
							temp_z=power( com_elem_interpoint[2][j],u );
							temp= add( temp, mul( Q_interpoly[i][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
						}
			printf("\nQ_interpoly[%d] in point[%d](%d,%d,%d) = %d",i,j,com_elem_interpoint[0][j],com_elem_interpoint[1][j],com_elem_interpoint[2][j],temp);
		}

		for(j=0;j<eta;j++)
		{	
			temp=0;
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						if(Q_interpoly[i][u][v][z]!=0)
						{
							temp_x=power( uncom_elem_interpoint[j][0],z );
							temp_y=power( uncom_elem_interpoint[j][1],v );
							temp_z=power( uncom_elem_interpoint[j][i+2],u );
							temp= add( temp, mul( Q_interpoly[i][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
						}
			printf("\nQ_interpoly[%d] in point[%d](%d,%d,%d) = %d",i,j,uncom_elem_interpoint[j][0],uncom_elem_interpoint[j][1],uncom_elem_interpoint[j][i+2],temp);
		}
		
	}
*/	//**************************

}


void com_elem_interpolation(int g[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int interpoint[][n-eta])
{
	int i, j, u, v, z, delta_temp, lod_temp, temp1, temp2; 
	int delta[init_polyNum],  J[init_polyNum], act[init_polyNum], lod[init_polyNum], lod_min, j_min;	//(delta, J, act, lod)[num of polys]
	int f[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];  //g[num_polys][z-deg+1][max(deg_y)+1][w+1], g1[z-deg+1][max(deg_y)+1][w+2], (g2, f)[z-deg+1][max(deg_y)+1][w+1]	
	int g1[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize+1], g2[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

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
		flag_addNum=1;
		flag_mulNum=1;
		j_min=-1;
		for(j=0;j<init_polyNum;j++)	//wrt each polynomial
		{
			J[j]=0;
			delta[j]=0;
			if(act[j]==1)	//for polynomials with leading order less of equal to C
			{
				//Hasse derivative mapping
				for(u=0;u<interpoly_Zsize;u++)	//rs
				{
					for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					{
						for(z=0;z<interpoly_Xsize;z++)	//w+1
							if(g[j][u][v][z]!=0)	//coputation num concerning coefficients
							{
								flag_mulNum=0;
								delta_temp=mul( mul( power(interpoint[0][i],z),power(interpoint[1][i],v) ), power(interpoint[2][i],u) );
								flag_mulNum=1;

								delta_temp=mul(delta_temp, g[j][u][v][z]);
								//Hasse derivative mapping
								delta[j]=add(delta[j],delta_temp);
							}
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
										g2[u][v][z] = mul(interpoint[0][i],g[j][u][v][z]);
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


	int temp_x, temp_y, temp_z,temp;

	for(j=0;j<init_polyNum;j++)   
	{
		temp=0;
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					if(g[j][u][v][z]!=0)
					{
						temp_x=power( interpoint[0][i],z );
						temp_y=power( interpoint[1][i],v );
						temp_z=power( interpoint[2][i],u );
						temp= add( temp, mul( g[j][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
					}
		
		if( (lod[j]<=iterNum) && temp!=0 )
			printf("\ng[%d] with lod[%d]=%d in point[%d](%d,%d,%d) = %d",j,j,lod[j],i,interpoint[0][i],interpoint[1][i],interpoint[2][i],temp);

	}
//	printf("\n");
*/	//*************************

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

#ifdef _PolyCoeffNum_
	int NumTemp=0;
	for(j=0;j<init_polyNum;j++)
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1	
					if(inGroup_temp[j][u][v][z]!=0)
						NumTemp++;

	NumTemp = NumTemp/init_polyNum;
	//*********************
#endif
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
	flag_addNum=1;
	flag_mulNum=1;		
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
							flag_mulNum = 0;
							poly_temp1 = mul( power(interpoint[0],z),power(interpoint[1],v) );
							flag_mulNum = 1;

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
	
/*	//*******debug*********
	//outGroup1
	for(j=0;j<init_polyNum;j++)
	{
		lod_temp=0;
		lod[j]=0;
		for(u=0;u<interpoly_Zsize;u++)	//rs
		{
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
			{
				for(z=0;z<interpoly_Xsize;z++)	//w+1
				{
					if(outGroup1[j][u][v][z]!=0)
					{
						lod_temp=mono_order[u][v][z];
						if(lod_temp>lod[j])
							lod[j]=lod_temp;
					}
				}
			}
		}
	}

	int temp_x,temp_y,temp_z,temp;

	for(j=0;j<init_polyNum;j++)
	{
		temp=0;
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					if(outGroup1[j][u][v][z]!=0)
					{
						temp_x=power( interpoint[0],z );
						temp_y=power( interpoint[1],v );
						temp_z=power( interpoint[2],u );
						temp= add( temp, mul( outGroup1[j][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
					}
		printf("\noutGroup1[%d] with lod[%d]=%d in point[7](%d,%d,%d) = %d",j,j,lod[j],interpoint[0],interpoint[1],interpoint[2],temp);
	}
	printf("\n");

	//outGroup2
	for(j=0;j<init_polyNum;j++)
	{
		lod_temp=0;
		lod[j]=0;
		for(u=0;u<interpoly_Zsize;u++)	//rs
		{
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
			{
				for(z=0;z<interpoly_Xsize;z++)	//w+1
				{
					if(outGroup2[j][u][v][z]!=0)
					{
						lod_temp=mono_order[u][v][z];
						if(lod_temp>lod[j])
							lod[j]=lod_temp;
					}
				}
			}
		}
	}
	
	for(j=0;j<init_polyNum;j++)
	{
		temp=0;
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					if(outGroup2[j][u][v][z]!=0)
					{
						temp_x=power( interpoint[0],z );
						temp_y=power( interpoint[1],v );
						temp_z=power( interpoint[3],u );
						temp= add( temp, mul( outGroup2[j][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
					}
		printf("\noutGroup2[%d] with lod[%d]=%d in point[7](%d,%d,%d) = %d",j,j,lod[j],interpoint[0],interpoint[1],interpoint[3],temp);
	}
	printf("\n");
*/	//******************

}

void factorisation(void)
{
	int i, j, u, v, z;
	
	for(i=0;i<test_vec_num;i++)
	{
		//Initialisation
		for(j=0;j<k;j++)	//number of fac steps=k
			for(u=0;u<facpoly_Zsize;u++)	//rs
				for(v=0;v<facpoly_Ysize;v++)	//y_size
					for(z=0;z<facpoly_Xsize;z++)	//w+1
						Q[j][u][v][z]=0;

		for(u=0;u<k;u++)	//number of fac steps=k
			for(v=0;v<lm+1;v++)	//5>rs=expected number of roots
				rootlist[u][v]=-1;	

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
					Q[uu][u][v][z]=Q_interpoly[i][u][v][z];

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

		//store the necessary data
		for(u=0;u<lm+1;u++)	//5>(rs-1)=expected number of output lists
			for(v=0;v<k;v++)	//k
				outputList[i][u][v]=output[u][v];

		listNum[i]=l;

	}

#ifdef _Complexity_
	//****debug********
	printf("\n\nCfac=%d\tmul=%d\tadd=%d\n\n", (addNum+mulNum-addNum_temp-mulNum_temp), (mulNum-mulNum_temp), (addNum-addNum_temp));
	addNum_temp = addNum;
	mulNum_temp = mulNum;
	//*****************
#endif

}

void rcs(int uu)
{
	int i, j, u, v, m, z, t, r, i_1, j_1, i_2, j_2, a, b, leadMono, leadMono_temp, alpha, act, temp;
	int lc[lm+1], q_temp[lm+1][rcspoly_Ysize][rcspoly_Xsize];	//q_temp[z-deg+1][y_size>max(deg_y)+1+(rs-1)*(max(deg_y) in encoding functions)][14>w+(rs-1)*w], lc[rs]--leading coefficient polynomial
	int d1,d2;
	int rcspoly_Ysize_1,rcspoly_Xsize_1;	//the q_temp size of the first step of factorization
//	int pbY,pbX;

	//array size initialization
	rcspoly_Ysize_1= faiMax_Ysize + interpoly_Ysize + 1;
	rcspoly_Xsize_1= faiMax_Xsize + w + 1;
//	pbY = rcspoly_Ysize+1;
//	pbX = w+1;

	//printf("\nWhen u=%d\n", u);
	leadMono=0; leadMono_temp=0;	//leading monomial index
	act=0;	//judge value for recursive search of each f_k-1-u

	flag_addNum = 1;
	flag_mulNum = 1;
	//find pb_k-1-uu
	j_1 = tg_order[k-1-uu][1];
	i_1 = tg_order[k-1-uu][0];

	/*
	for(j=0;j<pbY;j++)	//y_size, efficiently enough
	{	
		for(i=0;i<pbX;i++)	//w+1
		{
			if(pb[k-1-uu][j][i]==1)	//k-1-uu
			{
				j_1=j;	//record pb_k-1-u's y degree
				i_1=i;	//record pb_k-1-u's x degree
				break;
			}
		}
		if((j_1!=-1) && (i_1!=-1))
			break;
	}
	*/
	//printf("pb_k-1-u");

	//initialise q_temp
	for(i=0;i<lm+1;i++)	//rs
		for(j=0;j<rcspoly_Ysize;j++)	//y_size
			for(u=0;u<rcspoly_Xsize;u++)	//5>w*(rs-1)+w
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

	//convert x^w+1=y^w+y, difference with diff code
	d1=0;
	d2=0;

	//*****debug*****
	if( rcspoly_Xsize_1 < (w+1) )
		printf("\n\nrcs() has error\n\n");
	//**************

	for(i=0;i<lm+1;i++)	//rs
		{	
			for(j=0;j<rcspoly_Ysize_1;j++)	//y_size=degY+1
			{
				for(u=(w+1);u<rcspoly_Xsize_1;u++)	//w+1, 5>w*(rs-1)+w
				{
					if(q_temp[i][j][u]!=0)	//refering to two finitefield
					{
						d1=u/(w+1);	//d1 in (8,4) is less than 2;
						d2=u%(w+1);	//u%(w+1)

						//should add a funciton to calculate the poly (z+hixiyi)^j
						q_temp[i][j+d1*w][d2]= add(q_temp[i][j+d1*w][d2], q_temp[i][j][u]);	//this time, q_temp size is [rs+1][1+last_maxY(q_temp)+w*(last_maxX(q_temp)/(w+1))][w+1]
						q_temp[i][j+d1][d2]	 = add(q_temp[i][j+d1][d2], q_temp[i][j][u]);
						q_temp[i][j][u] = 0;
					}
					//q_temp[i][j][u] = 0;
					
				}
			}
		}
	//printf("q[u][pb_k-1-u]");

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

/*	//find the leading coefficient polynomial lc[rs]
	for(i=0;i<lm+1;i++)	//rs
	{
		flag=1;
		lc[i]=-1;
		for(j=0;j<rcspoly_Ysize;j++)	//y_size
		{
			for(u=0;u<w+1;u++)	//w+1
			{
				if(q_temp[i][j][u]!=0)
				{
					flag=0;
					if(j==j_2 && u==i_2)
					{
						lc[i]=q_temp[i][j][u];
						break;
					}
				}
			}
			if(lc[i]!=-1)
				break;
		}
		if(flag==1)
			lc[i]=0;
	}

	//if q_temp[i][j][u]==0, then setting lc[i]==0, and flag==2
	flag=0;
	for(i=0;i<lm+1;i++)
		if(lc[i]==0)
			flag++;

	//printf("leading coefficient polynomial");
	
	u=0;	//root index
	//find the roots of leading coefficient polynomial lc[rs]
	for(i=0;i<q;i++)	//number of elements in GF(4)
		if(flag!=(lm+1))	//rs
		{
			b=0;
			for(j=0;j<(lm+1);j++)	//rs
			{
				a=mul(lc[j], power(root[i], j));
				b=add(a, b);
			}
			if(b==0)
			{
				act=1;
				rootlist[uu][u]=root[i];
				u++;
			}
		}
		else if(flag==(lm+1))
			{
				rootlist[uu][0]=0;
				act=1;
			}
*/	//printf("roots");

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
			rootlist[uu][u]=root[i];
			u++;
#ifndef _NoReduction_
			break;	//jump out of the loop
#endif
		}
	}

	//For each distinct root of rootlist[u]
	if(act==1)
	{
		for(i=0;i<lm+1;i++)	//2>rs
		{
			if(rootlist[uu][i]!=-1)
			{	
				alpha=rootlist[uu][i];
				
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
					//Initialise q_temp
					for(j=0;j<lm+1;j++)	//rs
						for(u=0;u<rcspoly_Ysize;u++)	//y_size
							for(m=0;m<rcspoly_Xsize;m++)	//3>w*(rs-1)+w
								q_temp[j][u][m]=0;	
			
					//q_temp=q[uu][z+f_k-1-u*pb_k-1-u]
					for(j=0;j<lm+1;j++)	//rs
					{
						//calculate (z+f_k-1-u*pb_k-1-u)^j
						if(j==0)
						{
							for(u=0;u<lm+1;u++)
								for(m=0;m<expoly_Ysize;m++)
									for(z=0;z<expoly_Xsize;z++)
										expoly[u][m][z]=0;
							expoly[0][0][0]=1;
						}
						else if(j>0)
								polyexp1(alpha, i_1, j_1, j);
						
						//calculate q[uu][z+f_k-1-u*pb_k-1-u]
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
													q_temp[u][v+m][t+z] = add( q_temp[u][v+m][t+z],temp );
												}
												// caution overflow problem of q_temp
												// here q_temp, ysize=exploy_Ysize+facpoly_Ysize-2+1=9, xsize=expoly_Xsize+facpoly_Xsize-2+1=6
											}
										}
									}
								}
							}
						}
					}
				
					//convert x^w+1=y^w+y, difference with diff code
					d1=0;
					d2=0;
					for(j=0;j<lm+1;j++)	//rs
					{
						for(u=0;u<rcspoly_Ysize;u++)	//10 y_size
						{
							for(m=w+1;m<rcspoly_Xsize;m++)	//w+1, 14>w*(rs-1)+w
							{
								if(q_temp[j][u][m]!=0)	//refering to two finitefield add
								{
									d1=m/(w+1);	//d1 in (8,4) is less than 2;
									d2=m%(w+1);	//u%(w+1)
									q_temp[j][u+d1*w][d2] = add(q_temp[j][u+d1*w][d2], q_temp[j][u][m]);
									q_temp[j][u+d1][d2]	  = add(q_temp[j][u+d1][d2], q_temp[j][u][m]);
									q_temp[j][u][m] = 0;
								}
								
							}
						}
					}
				
					r=uu+1;

					//q[u+1]=q_temp
					for(j=0;j<lm+1;j++)	//rs
						for(u=0;u<facpoly_Ysize;u++)	//y_size
							for(m=0;m<facpoly_Xsize;m++)	//w+1
								Q[r][j][u][m]=q_temp[j][u][m];

					//printf("q[u+1]");

					//next coefficient searching
					
					rcs(r);
				}
			}
		}
	}

	flag_addNum = 0;
	flag_mulNum = 0; 
	
}

void choose()
{
	int i, j, u, v, value, flag, min_index_1, min_index_2;
	unsigned int mask=1;
	float proba[lm*test_vec_num], proba_temp, temp;
	int codeword_temp[n];

	flag_addNum=1;
	flag_mulNum=1;
	//normal mode
	//Initialise hamming distance counter
	for(u=0;u<lm*test_vec_num;u++)
		proba[u]=-1.0;

	proba_temp=0.0;
	flag = 0;
	min_index_1 = -1;
	min_index_2 = -1;

	for(i=0;i<test_vec_num;i++)
		if( listNum[i]!=0 )
		{
			for(j=0;j<listNum[i];j++)
			{
				//reencoding
				encoder(outputList[i][j],codeword_temp);
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
					min_index_1 = i;
					min_index_2 = j;

					for(v=0;v<n;v++)
						dec_codeword[v] = codeword_temp[v];

					flag = 1;	//exist at less one valid solution
				}
			}
		}

	//output the decoding result 
	if(flag==0)	// not exist a valid solution
	{
		for(u=0;u<n;u++)
			dec_codeword[u]=large_vec[0][u];
		
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

		//encode message[min_index]
//		encoder(outputList[min_index_1][min_index_2],dec_codeword);
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
	flag_addNum=0;
	flag_mulNum=0;

#ifdef _Complexity_
	//****debug********
	printf("\n\nCcho=%d\tmul=%d\tadd=%d\n\n", (addNum+mulNum-addNum_temp-mulNum_temp), (mulNum-mulNum_temp), (addNum-addNum_temp));
	addNum_temp = addNum;
	mulNum_temp = mulNum;
	//*****************
#endif

	//*********************

	//genius mode
	flag=0;
	int count_temp;
	for(i=0;i<test_vec_num;i++)
	{
		for(j=0;j<lm+1;j++)
		{
			count_temp=0;
			for(u=0;u<k;u++)
				if(outputList[i][j][u]==message[u])
				{
					count_temp++;
				}

			if(count_temp==k)
			{
				flag=1;
				DecSucc_SeqNum++;

				if( min_index_1>=0 || min_index_2>=0 )	//min_index==-1 represent there is no valid solution
					if( min_index_1!=i || min_index_2!=j )
					{
						ChosenWrong_SeqNum++;
					}

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
	
	//Check Herm(64, 39)'s working perperty
	epcount2=0;
	for(u=0;u<n;u++)
		if(codeword[u]!=dec_codeword[u])
			epcount2++;
	//printf("\nepcount2=%d\n", epcount2);

	//******debug*******
	//caculate the testCount2
	for(i=0;i<test_vec_num;i++)
	{
		testCount2[i]=0;
		for(j=0;j<k;j++)
			if( message[j]!=outputList[i][0][j] )
				testCount2[i]++;
	}

/*	int m=0;
	for(i=0;i<test_vec_num;i++)
		if( (n-testCount1[i])>degree_test[i] && testCount2[i]!=0 )
		{
			printf("\n\nseq_num_Now=%d, No.%d test vector factorization has failed!!, testCount1=%d, testCount2=%d, degreeOfPolynomail=%d, listNum=%d\n", seq_num_Now, i, testCount1[i], testCount2[i], degree_test[i], listNum[i]);
//			printf("\nNo.%d test vector has error!!, testCount1=%d, testCount2=%d, testCount1_uncom=%d, listNum=%d\n\n", i, testCount1[i], testCount2[i], (testCount1[i]-testCount1_com), listNum[i]); 
//			printf("\nseq_num_Now=%d, This sequence is decoding failed， and error num is %d!\n\n", seq_num_Now, epcount1);
			//********debug*************
//			int temp_x, temp_y, temp_z, temp, h;
			int z;

			//**************
			printf("\n");
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						if(Q_interpoly[i][u][v][z]!=0)
						{
							printf("Q_interpoly[%d][%d][%d][%d]=%d\n",i,u,v,z,Q_interpoly[i][u][v][z]);
						}
			printf("\n");

			//************
			printf("\nTx message is:\n");
			for(u=0;u<k;u++)	//k
				printf("%d\t", message[u]);

			for(m=0;m<listNum[i];m++)
			{
				printf("\nlist_%d:\n",m);
				for(u=0;u<k;u++)
					printf("%d\t",outputList[i][m][u]);
			}
			printf("\n\n");
			//***************
			printf("\nTx message is:\n");
			for(u=0;u<k;u++)	//k
				printf("%d ", message[u]);
			printf("\nTx code word is:\n");
			for(u=0;u<n;u++)	//n
				printf("%d ", codeword[u]);

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

			printf("\n\n");
			printf("codeword_order[]");
			for(j=0;j<n;j++)
				printf("\t%d",codeword[(int)test_set_ordered[1][j]]);

			printf("\n%d errors in received!", epcount1);
			printf("\noutputlist:");
			
			for(i=0;i<test_vec_num;i++)
			{
				printf("\n\nlistNum[%d]=%d",i,listNum[i]);
				for(m=0;m<listNum[i];m++)
				{
					printf("\n[%d][%d]",i,m);
					for(u=0;u<k;u++)
						printf("\t%d",outputList[i][m][u]);
					printf("\tproba[%d][%d]=%f",i,m,proba[i*(lm+1)+m]);
					
				}
			}
			printf("\nDecoded code word is:\n");
			for(u=0;u<n;u++)	//n
				printf("%d ", dec_codeword[u]);
			printf("\n%d errors after decoding!\n\n", epcount2);
			//***********
		}

*/			//*******************

	if( epcount1<=able_correct && epcount2!=0 )	//this seq_num has chosen the wrong one
	{
		printf("\n\nDecoding is error, epcount1 = %d, epcount2 = %d\n\n", epcount1, epcount2);
	}
	
}


void mono_table()
{
	int i, j, u, v, z, weight[weight_Zsize][weight_XYsize], mono_order_1[weight_Zsize][weight_XYsize];	//26=(100/weight(z))+1, 100=pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][max(deg_y)+1][w+1]

	//Initialisation
	for(i=0;i<weight_Zsize;i++)	//(230/weight(z))+1
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
		for(i=0;i<weight_Zsize;i++)	//230/weight(z)
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

/*	//**********debug********
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
	for(i=0;i<weight_Zsize;i++)
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

	printf("\n\nweight degree for rcs():");
	for(j=0;j<interpoly_Ysize;j++)
	{
		printf("\n%d",j);
		for(i=0;i<interpoly_Xsize;i++)
			printf("\t%d",weiDeg_rcs[j][i]);
	}

	printf("\n\nweight degree index for rcs():");
	for(j=0;j<interpoly_Ysize;j++)
	{
		printf("\n%d",j);
		for(i=0;i<interpoly_Xsize;i++)
			printf("\t%d",weiDeg_rcs_index[j][i]);
	}
*/
	//******************

}


void polyexp1(int c, int i, int j, int deg_z)
{
	int u, m, z;	//p1[rs+1][18>(max(deg_y) in encoding functions)*(rs-1)][10>(max(deg_x) in encoding functions)*(rs-1)], p2[rs][26=18+max(deg_y) in encoding functions][14=10+max(deg_x) in encoding functions]

	//Initialisations
	for(u=0;u<lm+1;u++)	//rs
		for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
			for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
				expoly[u][m][z]=0;

	if(deg_z==0)
		expoly[0][0][0]=1;
	else if(deg_z>0)	//because deg_z is equal to or less than 1, so deg_z>0 <-> deg_z==1
	{
		expoly[1][0][0]=1;
		expoly[0][j][i]=c;

	}

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

	//judge the index's scale of coresponding tgsize
	/*pesudo code
	if(  0 <= index < w ) tgsize = (index+2) * (index+1) / 2
	else if( w <= index ) tgsize = w*(w+1)/2 + (w+1)*(index-w+1)
	*/
	index_temp = 498;	//according to the last formualtion to calculate out

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