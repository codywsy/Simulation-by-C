#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _GS_Normal_
//#define _Complexity_
#define _NoReduction_
//#define _PolyCoeffNumUncom_
//#define _PolyCoeffNumFac_
#define PrintfMonoTable

#define OpenFile fp=fopen("LCC_Herm(8,4)GS.txt","a")

#define q 4
#define n 8
#define k 4
#define p 2
#define w 2
#define weiz 4
#define iterNum 8	//when m=1, C is equal to n
#define lm 1
#define able_correct 1
#define pointNum 2
#define interval 1

//conditional compile
#ifndef _GS_Normal_
	#define eta 3
	#define test_vec_num 8 //the num of test_vec is 2^eta
#else
	#define eta 0
	#define test_vec_num 1 //the num of test_vec is 2^eta
#endif

//test_vec_construction
#define choose_num 2	//fixed

//tgorder()
#define tg_size 186	//more than k, tg_size represent the probably used pole basis num, (w+1)*w/2 + (w+1) * ( (w+interpoly_Ysize)-w+1 ) 
//mono_table()
#define weight_Zsize 4
#define weight_XYsize 186	//equals to the tg_size
#define mono_ordinSize 200	//choose the value with debugging, equals to max weigtdegree in weight[monoTable_Zsize][weight_XYsize]
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 60	//large than interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize 186 //>index of term[max(deg_y)][w] in the pole basis + 1, nearly equals to tg_size
//interpolation()
#define init_polyNum (w*(lm+1))	//change with diff lm, the poly num of the init polyGroup
#define interpoly_Zsize (lm+1)	//maxValue of z is lm=1, so the Zsize should add 1 more.
#define interpoly_Ysize 6	//maxdeg of y is (w-1) + w*(n/(w+1)), so the Ysize should add 1 more.
#define interpoly_Xsize	(w+1)	//maxdeg of x is w, so the Xsize should add 1 more.
//factorization()
#define facpoly_Zsize (lm+1)// same as the interpoly_Zsize
#define facpoly_Ysize 6		// same as the interpoly_Ysize
#define facpoly_Xsize (w+1)	// same as the interpoly_Xsize
//rcs()
#define rcspoly_Ysize 11	//degY+1, degY = interpoly_Ysize + max(j_1) + w*( (max(i_1)+w)/(w+1) )
#define rcspoly_Xsize 6	//degX+1, degX = max(i_1) + w
#define faiMax_Ysize 2	//change with diff w and k, the max degY of probably used polebasis
#define faiMax_Xsize 3	//change with diff w and k, the max degX of probably used polebasis
//expoly()-->expanded polynomial
#define expoly_Ysize (faiMax_Ysize+1)	
#define expoly_Xsize (faiMax_Xsize+1)


//main()
unsigned long int seq_num;	//number of input binary sequences
float SNR,N0;
double BER, FER;
//GF elements
int mularray[]={1,2,3};
int root[]={0, 1, 2, 3};	//elements in GF(4)
//*************
int bi_message[k*p], message[k];	//transmitted messge[k]
int codeword[n], bi_codeword[n*p]; //codewords
float tx_symbol[p*n][2], rx_symbol[p*n][2], sgm;
float RM[q][n];
//findpoint()
int point[2][n];	//rational points[2][w^3], [0]->x, [1]->y
//tgorder()
int tg_order_1[tg_size][2],tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
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
float pi=3.141593;	// pai

//int Q_uncom_elem[test_vec_num][init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

	//****debug**********
//	int com_elem_interpoint[3][(n-eta)];
//	int uncom_elem_interpoint[eta][4];
//	int degree_temp[test_vec_num][init_polyNum];
	int degree_test_temp[test_vec_num][init_polyNum];  //store the degree of leading monomial of polynomial
	int degree_test[test_vec_num];		//store the degree of the choosen polynomial

	unsigned long int seq_num_Now;	//record the seq_num at the moment
	unsigned long int DecSucc_SeqNum, ChosenWrong_SeqNum;	//record the seq_num should be decoding success in genius mode
	double CWR;	//record the rate of chosen uncorrectly

	//some var about computation complexity
	unsigned long int addNum, mulNum, addNum_temp, mulNum_temp;
	int flag_addNum=0, flag_mulNum=0;

	//*******************


int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
void findpoints(void);
void mono_table(void);
void polyexp1(int, int, int, int, int poly[][expoly_Ysize][expoly_Xsize]);
//void zerobasis(void);
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

		for(j=1;j<=seq_num;j++)
		{
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
				successError_count = successError_count + (epcount1-successError_count)/(double)(v);
			}

			channelError_count = channelError_count + (epcount1-channelError_count)/(double)(j);

			progress=(double)(j*100)/(double)seq_num;
			BER=(double)(error)/(double)(n*p*j);
			FER=(double)(ferror)/(double)(j);
			
			CWR=(double)(ChosenWrong_SeqNum)/(double)(DecSucc_SeqNum); 

			printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, Choice Errors=%2.1d, CWR=%E, channelError_count=%0.2f, successError_count=%0.2f\r", progress, SNR, error, BER, ferror, FER, ChosenWrong_SeqNum, CWR, channelError_count, successError_count);

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

		printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, Choice Errors=%2.1d, CWR=%E, channelError_count=%0.2f, successError_count=%0.2f\n", progress, SNR, error, BER, ferror, FER, ChosenWrong_SeqNum, CWR, channelError_count, successError_count);

		OpenFile;
		fprintf(fp,"Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, Choice Errors=%2.1d, CWR=%E, channelError_count=%0.2f, successError_count=%0.2f\n", progress, SNR, error, BER, ferror, FER, ChosenWrong_SeqNum, CWR, channelError_count, successError_count);
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
//		for(int x=0;x<pointNum;x++)
//		for(int y=0;y<pointNum;y++)
//		for(int h=0;h<pointNum;h++)
		{
			RM[j][i]=(Pr[0][v]*Pr[1][u]);  
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
	int Q_uncom_elem[test_vec_num][init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];

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
			uncom_elem_interpolation(uncom_elem_interpoint[i],Q_uncom_elem[u],Q_uncom_elem[2*u+0],Q_uncom_elem[2*u+1]);
	}		


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
	int i, j, u, v, z, delta[init_polyNum], delta_temp, J[init_polyNum], act[init_polyNum], lod_temp, lod[init_polyNum], lod_min, j_min;	//(delta, J, act, lod)[num of polys]
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
			
			
		//Initialise the eliminator array act[num_poly]
		for(j=0;j<init_polyNum;j++)	//num_poly
		{
			if(lod[j]<=iterNum)	//C=n when multiplicity = 1
				act[j]=1;
			else
				act[j]=0;
		}
			
			
		//Calculate the hasse derivative mapping of each polynomials
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
						{
							if(g[j][u][v][z]!=0)
							{
								delta_temp=mul( mul( power(interpoint[0][i],z),power(interpoint[1][i],v) ), power(interpoint[2][i], u) );
								delta_temp=mul(delta_temp, g[j][u][v][z]);
								//Hasse derivative mapping
								delta[j]=add(delta[j],delta_temp);
							}
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
									g[j][u][v][z]=add(mul(delta[j_min],g[j][u][v][z]),mul(delta[j],f[u][v][z]));	
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
									g1[u][v][z+1]=g[j][u][v][z];
						//convert x^5=y^4+y, difference with diff code
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
									g2[u][v][z]=mul(interpoint[0][i],g[j][u][v][z]);
						//g=g1+g2
						for(u=0;u<interpoly_Zsize;u++)	//rs
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
								for(z=0;z<interpoly_Xsize;z++)	//w+1
									g[j][u][v][z]=add(g1[u][v][z],g2[u][v][z]);
					}
				}
			}
		}


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
	int f[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];  //g[num_polys][z-deg+1][max(deg_y)+1][w+1], g1[z-deg+1][max(deg_y)+1][w+2], (g2, f)[z-deg+1][max(deg_y)+1][w+1]	
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
	for(j=0;j<init_polyNum;j++)	//num_poly
	{
		if(lod[j]<=iterNum)	//C=n when multiplicity = 1
			act[j]=1;
		else
			act[j]=0;
	}
		
		
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
								poly_temp1 = mul( result_temp[index_min][0],inGroup_temp[j][u][v][z] );
								poly_temp2 = mul( result_temp[j][0],inGroup_temp[index_min][u][v][z] ) ;
								outGroup1[j][u][v][z]=add(poly_temp1,poly_temp2);
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
					//convert x^5=y^4+y, difference with diff code
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
								g2[u][v][z]=mul( interpoint[0],inGroup_temp[index_min][u][v][z]);
					//g=g1+g2
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								outGroup1[index_min][u][v][z]=add(g1[u][v][z],g2[u][v][z]);
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
								poly_temp1 = mul( result_temp[index_min][1],inGroup_temp[j][u][v][z] );
								poly_temp2 = mul( result_temp[j][1],inGroup_temp[index_min][u][v][z] );
								outGroup2[j][u][v][z]=add( poly_temp1,poly_temp2 );	
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
					//convert x^3=y^2+y, difference with diff code
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
								g2[u][v][z]=mul( interpoint[0],inGroup_temp[index_min][u][v][z]);
					//g=g1+g2
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								outGroup2[index_min][u][v][z]=add(g1[u][v][z],g2[u][v][z]);
				}
			}
		}
	}
	
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

		//recursive coefficient search
		rcs(uu);

		//store the necessary data
		for(u=0;u<lm+1;u++)	//5>(rs-1)=expected number of output lists
			for(v=0;v<k;v++)	//k
				outputList[i][u][v]=output[u][v];

		listNum[i]=l;

	}

}

//this rcs() only can be used in lm=1
void rcs(int uu)
{
	int i, j, u, v, m, z, t, r, i_1, j_1, i_2, j_2, a, b, leadMono, leadMono_temp, alpha, act, temp;
	int lc[lm+1], q_temp[lm+1][rcspoly_Ysize][rcspoly_Xsize];	//q_temp[z-deg+1][y_size>max(deg_y)+1+(rs-1)*(max(deg_y) in encoding functions)][14>w+(rs-1)*w], lc[rs]--leading coefficient polynomial
	int d1,d2,index_flag;
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
					r=uu+1;

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
							for(m=0;m<facpoly_Ysize;m++)
								for(z=0;z<facpoly_Xsize;z++)
									q_temp[j][m][z] = Q[uu][j][m][z];
						}
						else if(j>0)
						{
							polyexp1(alpha, i_1, j_1, j, expoly);
						
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

								//convert x^w+1=y^w+y, difference with diff code
								int Ysize_temp = (facpoly_Ysize-1) + (expoly_Ysize-1) + 1;	//keep universality
								index_flag = (facpoly_Xsize-1) + (expoly_Xsize-1) ;	//the max deg_x for q_temp, rcspoly_Xsize = (expoly_Xsize-1) + (facpoly_Xsize-1) + 1
								while( index_flag>w )
								{
									temp = index_flag-(w+1);	// deg_x - (w+1)
									for(m=0;m<Ysize_temp;m++)
										if( q_temp[u][m][index_flag]!=0 )
										{
											q_temp[u][m+1][temp] = add( q_temp[u][m+1][temp],q_temp[u][m][index_flag] );	//y^1
											q_temp[u][m+w][temp] = add( q_temp[u][m+w][temp],q_temp[u][m][index_flag] );	//y^(w+1)
											q_temp[u][m][index_flag] = 0;
										}
									index_flag = index_flag-1;
								}

							}
						}

					}
								
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
		for(u=0;u<lm+1;u++)	//rs
			for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
				for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
					poly[u][m][z]=0;

		poly[1][0][0]=1;
		poly[0][j][i]=c;
	}
	else if(deg_z>1)
	{
		printf("\n\n expoly has error!!");
	}
/*	else if(deg_z>1)	//when lm=1, deg_z<=1
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
	*/
}

void choose()
{
	int i, j, u, m, v, min_num, value, flag, min_index;
	unsigned int mask=1;
	float proba[lm*test_vec_num],temp_x,temp_y,proba_temp;
	int bi_codeword_temp[n*p],codeword_temp[n];
	int tx_symbol_temp[n*p][2];

	//normal mode
	//Initialise hamming distance counter
	for(u=0;u<lm*test_vec_num;u++)
		proba[u]=-1.0;

	for(i=0;i<test_vec_num;i++)
		if( listNum[i]!=0 )
		{
			for(j=0;j<listNum[i];j++)
			{
				//reencoding
				encoder(outputList[i][j],codeword_temp);
				//calculate the posteriori probablity
				proba_temp=1.0;
				for(u=0;u<n;u++)
				{
					for(v=0;v<q;v++)
						if(codeword_temp[u]==root[v])
							proba_temp=proba_temp*RM[v][u];
				}

				proba[i*lm+j]=proba_temp;
			}
		}

	//choose the min
	proba_temp=0.0;
	min_index=-1;
	flag=0;
	for(i=0;i<lm*test_vec_num;i++)
		if(proba[i]!=-1.0)
		{
			flag=1;	//exist at less one valid solution

			if(proba_temp<proba[i])
			{
				proba_temp=proba[i];
				min_index=i;
			}
		}

	//encode message[min_index]
	encoder(outputList[min_index/lm][min_index%lm],codeword_temp);

	//nonbinary --> binary
		for(u=0;u<n;u++)
		{	
			value=codeword_temp[u];
			mask=1;
			for(v=0;v<p;v++)
			{
				if((value & mask)>0)
					bi_codeword_temp[p*u+v]=1;
				else
					bi_codeword_temp[p*u+v]=0;
				mask=mask<<1;
			}
		}
		//******************


	//output the correct 
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
		for(u=0;u<n;u++)
			dec_codeword[u]=codeword_temp[u];
		for(u=0;u<n*p;u++)
			dec_bicodeword[u]=bi_codeword_temp[u];
	}
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

/*	if(flag==0)	// not exist a valid solution
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
	else if(flag==1)
	{
		for(u=0;u<n;u++)
			dec_codeword[u]=codeword[u];		

		for(u=0;u<n*p;u++)
			dec_bicodeword[u]=bi_codeword[u];
	}
*/
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

	for(i=0;i<test_vec_num;i++)
		if( (n-testCount1[i])>degree_test[i] && testCount2[i]!=0 )
		{
			printf("\n\nseq_num_Now=%d, No.%d test vector factorization has failed!!, testCount1=%d, testCount2=%d, degreeOfPolynomail=%d, listNum=%d\n", seq_num_Now, i, testCount1[i], testCount2[i], degree_test[i], listNum[i]);
//			printf("\nNo.%d test vector has error!!, testCount1=%d, testCount2=%d, testCount1_uncom=%d, listNum=%d\n\n", i, testCount1[i], testCount2[i], (testCount1[i]-testCount1_com), listNum[i]); 
//			printf("\nseq_num_Now=%d, This sequence is decoding failed， and error num is %d!\n\n", seq_num_Now, epcount1);
			//********debug*************
			int temp_x, temp_y, temp_z;
			int temp,h,z;

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
			printf("x[]\t\t");
			for (i = 0; i < n; ++i)
				printf("\t(%d,%d)", point[0][i], point[1][i]);

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

			//*******************

	if( epcount1<=able_correct && epcount2!=0)	//this seq_num has chosen the wrong one
	{
		ChosenWrong_SeqNum++;
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

void mono_table()
{
	int i, j, u, v, z, weight[weight_Zsize][weight_XYsize], mono_order_1[monoTable_Zsize][weight_XYsize];	//26=(100/weight(z))+1, 100=pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][max(deg_y)+1][w+1]

	//Initialisation
	for (i = 0; i<weight_Zsize; i++)
		for (j = 0; j<weight_XYsize; j++)
		{
			weight[i][j] = -1;
		}

	for (i = 0; i<monoTable_Zsize; i++)
		for (j = 0; j<weight_XYsize; j++)
		{
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

int mul(int fac1,int fac2)
{
	int mulresult=0;

	if(fac1==0||fac2==0)
		mulresult=0;
	else
	{
		for(int i=0;i<(q-1);i++)
			{
			for(int j=0;j<(q-1);j++)
				{
				if(fac1==mularray[i]&&fac2==mularray[j])
					mulresult=mularray[(i+j)%(q-1)];
				}
			}
	}

	return mulresult;
}

int add(int fac1,int fac2)
{
	int i,c[p],d[p],e[p],f;
	unsigned int mask=1;

	for(i=0;i<p;i++){
		if((fac1 & mask)>0)
			c[i]=1;
		else
			c[i]=0;

		mask=mask<<1; //shift 1 bit left
	}

	mask=1;
	for(i=0;i<p;i++){
		if((fac2 & mask)>0)
			d[i]=1;
		else
			d[i]=0;

		mask=mask<<1;
	}

	for(i=0;i<p;i++)  //p=6
		e[i]=c[i]^d[i];
	
	int num=1;
	f=0;
	for(i=0;i<p;i++)
	{
		f+=(e[i]*num);
		num*=2;
	}

	return f;
}


int power(int a, int b)
{
	int i,temp,pow_result=-1;

	if(b==0)
		pow_result=1;
	else if(a==0 && b!=0)
		pow_result=0;
	else if(a>0 && b!=0)
	{
		for(i=0;i<q-1;i++)
			if(a==mularray[i])
			{
				temp=(i*b)%(q-1);
				pow_result=mularray[temp];
			}
	}

	return pow_result;

}

int inv(int fac)
{
	int i,invresult;
	if(fac==0)
		invresult=0;
	else
		for(i=0;i<q-1;i++){	//size of mularray
			if(fac==mularray[i])
			invresult=mularray[((q-1)-i)%(q-1)];	//size of mularray
	}

	return invresult;
}


