/********change element with different code************
1.n,k,q,p			m-->iter_num,tm,lm,mono_size,fa_size
2.a[n],mul_array[n]
3.demodulation()--> RM[][] contruction
*****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//mode_EuclideanDisdance	mode = 1
//mode_PostProba			mode = 2
//mode_MLrule				mode = 3
//mode_genius				mode = 4
#define mode 2
#define finiteField 16

#define checkDecoding
//#define checkPrint
//#define checkInterpolation	

//************variable define***************
#define k 7  //length of message
#define eta 2 // chosen num
#define test_vec_num 4 //the num of test_vec is 2^eta
#define able_correct_num ((n-k+1)/2)
#define m 1 // the multiplicity of LCC algorithm is 1
#define lm 1 // LCC's lm is equal to 1, solid!!
#define V 1 //projected basic function BPSK=1, QPSK=2
#define pointNum 2 // the number of constell point == 2^V
#define choose_num 2 //the num of codeword choice, choose 2 in BASIC LCC algorithm

#define interval 1
 
//****modify variable**********
#define monoTable_Ysize 2	//designed to cover fully, larger than monoOrder_Ysize
#define monoTable_Xsize 70	// designed to cover fully, larger than monoOrder_Xsize and (n+k+1)
#define interpoly_Ysize (lm+1)	//value+1, fix value is lm 
#define interpoly_Xsize (n+1)	//value+1, max value is iterNum in this algorithm

//*********************************************

//************************
#if finiteField == 16
//(15,X)
#define n 15	//length of codeword
#define q 16  //GF(q) q=n+1
#define p 4 //GF(2^p) = GF(q)
int mularray[]={1,2,4,8,3,6,12,11,5,10,7,14,15,13,9}; 	//this array is used to the infinite field element mul
int root[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};			//n+1		be used to the factorization
int logNum[] = {-1,0,1,4,2,8,5,10,3,14,9,7,6,13,11,12};	//used to locate the degree of finitefield element through the value of finitefield element

#elif finiteField == 32
//(31,X)
#define n 31	//length of codeword
#define q 32  //GF(q) q=n+1
#define p 5	//GF(2^p) = GF(q)
int mularray[]={1,2,4,8,16,5,10,20,13,26,17,7,14,28,29,31,27,19,3,6,12,24,21,15,30,25,23,11,22,9,18}; 
				//this array is used to the infinite field element mul
int root[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};	
			//n+1		be used to the factorization
int logNum[] = {-1,0,1,18,2,5,19,11,3,29,6,27,20,8,12,23,
				4, 10,30,17,7,7,22,28,26,21,25,9,16,13,14,
				15};	//used to locate the degree of finitefield element through the value of finitefield element

#elif finiteField == 64
//(63,X)
int mularray[]={1,2,4,8,16,32,3,6,12,24,48,35,5,10,20,40,
				19,38,15,30,60,59,53,41,17,34,7,14,28,56,51,37,
				9,18,36,11,22,44,27,54,47,29,58,55,45,25,50,39,
				13,26,52,43,21,42,23,46,31,62,63,61,57,49,33};  //this array is used to the infinite field element mul

int root[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,    
			17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
			33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,
			49,50,51,52,53,54,55,56,57,58,59,60,61,62,63};	//n+1		be used to the factorization

int logNum[] = {-1,0,1,6,2,12,7,26,3,32,13,35,8,48,27,18,
				   4,24,33,16,14,52,36,54,9,45,49,38,28,41,19,
				   56,5,62,25,11,34,31,17,47,15,23,53,51,37,44,
				   55,40,10,61,46,30,50,22,39,43,29,60,42,21,20,
				   59,57,58};	//used to locate the degree of finitefield element through the value of finitefield element
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
float pi=3.141593;
float tx_symbol[n*p][2];	//modulated tranmitted signal
float rx_symbol[n*p][2];	//received destorted signal

int mono_order[monoTable_Ysize][monoTable_Xsize];	//multiplexity matrix
float N0, sgm;
int gmatrix[k][n];
//***************

//************demodulation var*************
float RM[q][n];   // reliability matrix
//int test_vec[test_vec_num][n]; //store test vectors
int test_vec_com[n-eta];	//store first n-eta symbols
float test_set_ordered[2][n]; // [0] stores value, [1] stores index
int large_vec[choose_num][n]; // largest index is 0, second largest index is 1
int x_ordered[n];
int T_poly[k]; // store T(x)
#ifdef checkDecoding
	int test_vec_check[test_vec_num][eta];
#endif


int psi_codeword_temp[n]; 	// store psi 
int reen_codeword_ordered[n-eta];
int v_poly[k+1]; //store v(x)
int Q_interpoly[test_vec_num][lm+1][n+1];
int message_temp[test_vec_num][k];
int valid_flag[test_vec_num]; //store if the message_temp[j] of fac is vaild


	//*******debug*************
	unsigned long int addNum, mulNum, addNum_temp, mulNum_temp;
	int flag_addNum=0, flag_mulNum=0;
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
float PDF(int,int);
void test_set_construction(void);
void decoding(void);

void main()
{
	int i, x,value,num;
	float start,finish;
	unsigned long int j;
	unsigned int mask=1;
	long int error, ferror;
	double progress, b;
	double addNum_count, mulNum_count, totalNum_count;

	FILE *fp;
	if((fp=fopen("LCC(15,7)¦Ç=2.txt","a"))==NULL)
	{
		printf("Can't open the data.txt\n");
		exit(0);
	}
	fclose(fp);


	srand(1977);

	flag_addNum=1;
	flag_mulNum=1;

	mono_table();

	generator();

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

	for(SNR=start;SNR<=finish;SNR+=interval)
	{	
		N0=(1.0/(float(k)/float(n)))/pow(10.0, SNR/10.0);
		sgm=sqrt(N0/2);
		b=1.0;
		error=0;
		ferror=0; 

		addNum_count=0.0;
		mulNum_count=0.0;
		totalNum_count=0.0;

		for(j=1;j<=seq_num;j++)
		{
			addNum=0;
			mulNum=0;
			//printf("\n\n*************For the %dth frame*************:\n", j);
			//Generate binary message
			for(i=0;i<k*p;i++)
				bi_message[i]=rand()%2;

			//Convert to nonbinary
			
			for(i=0;i<k;i++)
			{
				num=1;
				message[i]=0;
				for(x=0;x<p;x++)
				{
					message[i]=message[i]+(bi_message[p*i+x]*num);
					num=num*2;
				}
			}
					

			encoder(message,codeword);

			//Convert the codeword into binary
			for(i=0;i<n;i++)
			{
				value=codeword[i];
				mask=1;
				for(x=0;x<p;x++) //for(m=p-1;m>=0;m--)
				{
					if((value & mask)>0)
						bi_codeword[p*i+x]=1;
					else
						bi_codeword[p*i+x]=0;
					mask=mask<<1;
				}
			}
			

			modulation();
			
			channel();
			
			/////////To be done?////////////
			demodulation(); // output soft information

			test_set_construction();
			
			decoding();	

			////////////////////////////////
			
			//bit error rate calculation
			int temp = error;
			for(i=0;i<p*n;i++)
			{
				if(bi_codeword[i]!=dec_bicodeword[i])
				{
					error++;
				}
			}

			//frame error rate calculation
			if( error>temp )
				ferror++;
/*			
			if(f_error_temp==1)
				printf("\nerror=%d\n",error);
			f_error_temp=0;
*/			/*
			if(error!=0)
				ferror++;
			*/
			addNum_count = addNum_count + (addNum-addNum_count)/(double)(j);
			mulNum_count = mulNum_count + (mulNum-mulNum_count)/(double)(j);
			totalNum_count = addNum_count + mulNum_count;

			progress=(double)(j*100)/(double)seq_num;

			BER=(double)(error)/(double)(n*p*j);
			FER=(double)(ferror)/(double)(j);
			printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\r", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);
			if(ferror>309)
				break;

		}

		if(ferror>309)
		{
			BER=(double)error/(double)(n*p*j);
			FER=(double)(ferror)/(double)(j);
		}
		else
		{
			BER=(double)error/(double)(n*p*seq_num);
			FER=(double)(ferror)/(double)(seq_num);
		}

		printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);

		fp=fopen("LCC(15,7)¦Ç=2.txt","a");
		fprintf(fp,"Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n",progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);
		fclose(fp);

	}
	
}

void mono_table()
{
	int i;
	
	// 1st row
	for(i=0;i<monoTable_Xsize;i++)
		mono_order[0][i]=2*i+1;
	
	// 2nd row
	mono_order[1][0]=0;
	for(i=1;i<monoTable_Xsize;i++)
		mono_order[1][i]=2*i;

}

int power(int a, int b)
{
	int temp,pow_result=-1;

	if(b==0)
	{
		pow_result=1;
		return pow_result;
	}
	else if(a==0 && b!=0)
	{
		pow_result=0;
		return pow_result;
	}
	else if(a>0 && b!=0)
	{
		if(flag_mulNum==1)
		{
			mulNum += b;
		}		
		temp = (logNum[a]*b) % (q-1);
		pow_result=mularray[temp];
		return pow_result;
	}

	if( a<0 )
	{	
		printf("\n\n power has error!!");
	}

}

int mul(int fac1,int fac2)
{
	int mulresult=0,temp;

	//mulNum increasing one
	if(flag_mulNum==1)
	{
		mulNum++;
	}

	if(fac1==0||fac2==0)
	{
		mulresult = 0;
		return mulresult;
	}
	else 
	{
		temp = (logNum[fac1]+logNum[fac2]) % (q-1);
		mulresult = mularray[temp];
		return mulresult;
	}

}

int add(int fac1, int fac2)
{
	//addNum increasing one
	if( flag_addNum==1 )
	{
		addNum++;
	}
	return (fac1 ^ fac2);
}

int inv(int fac)
{
	int i, invresult;

	if(fac==0)
		invresult=0;
	else
	{
/*		for(i=0;i<n;i++)
		{
			if(fac==mularray[i])
				invresult=mularray[(n-i)%n];
		}
*/
		i = logNum[fac]%n;
		invresult = mularray[(n-i)%n];
	}

	return invresult;
}

void generator(void)
{
	int i, j;
	for(i=0;i<k;i++)
		for(j=0;j<n;j++)	//n
			gmatrix[i][j]=0;

	//generator matrix
	for(i=0;i<k;i++)	//k
		for(j=0;j<n;j++)	//n
			gmatrix[i][j] = power( mularray[j],i );

}

void encoder(int message_temp[], int codeword_temp[])
{
	int i, j, temp;

	//Encoding
	for(i=0;i<n;i++)	//n
	{
		codeword_temp[i]=0;
		for(j=0;j<k;j++)	//k
		{
			temp = mul( message_temp[j],gmatrix[j][i] );
			codeword_temp[i] = add( codeword_temp[i],temp );
		}
	}

}

void modulation(void)
{
	int i;

	//BPSK
	for(i=0;i<n*p;i++)
	{
		tx_symbol[i][0]=-2*bi_codeword[i]+1;	//0-->(1,0)
		tx_symbol[i][1]=0;  					//1-->(-1,0)
	}
}

void channel(void)
{
	int i, j;
	float u, r, g;

	//Add AWGN 
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

void demodulation(void)  //reconed constellation diagram
{	//soft decision
	int i,j,u,v;
	float sum,proba[pointNum];  
	float Pr[p][2];

	//initialize
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
//		for(u=0;u<pointNum;u++)
		for(v=0;v<pointNum;v++)
		for(int x=0;x<pointNum;x++)
		for(int y=0;y<pointNum;y++)
		for(int h=0;h<pointNum;h++)
		{
//			RM[j][i]=(Pr[4][u]*Pr[3][v]*Pr[2][x]*Pr[1][y]*Pr[0][h]); 
			RM[j][i]=(Pr[0][h]*Pr[1][y]*Pr[2][x]*Pr[3][v]);
			j++;
		}
	

	}

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

void test_set_construction(void)
{	//construct the test_vec[test_vec_num]
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
	
	//contruction test_vec
	for(i=0;i<n-eta;i++)
		test_vec_com[i]=large_vec[0][(int)test_set_ordered[1][i]];

	//construct x_ordered
	for(i=0;i<n;i++)
		x_ordered[i]=mularray[(int)test_set_ordered[1][i]];

#ifdef checkPrint	
	//************debug***************
	printf("\n\n");
	printf("xi[]\t");
	for(j=0;j<n;j++)
		printf("\t%d",mularray[j]);

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
	for(j=0;j<n;j++)
		printf("\t%d",x_ordered[j]);

	
	printf("\n\n");
	printf("large_vec_order[]");
	for(j=0;j<n;j++)
		printf("\t%d",large_vec[0][(int)test_set_ordered[1][j]]);
	
	printf("\n\n");
	printf("large_vec2_order[]");
	for(j=0;j<n;j++)
		printf("\t%d",large_vec[1][(int)test_set_ordered[1][j]]);
	
	printf("\n\n");
	printf("codeword[]\t");
	for(j=0;j<n;j++)
		printf("\t%d",codeword[(int)test_set_ordered[1][j]]);
	//********************************
#endif
}

void decoding(void)
{
	int i, j, v, u, value;
	unsigned int mask=1;
	int temp, min_num=-1;
	int codeword_temp[test_vec_num][n], bicodeword_temp[test_vec_num][n*p], bicodeword_coor[n*p];

	int flag=0, index_temp, choosen_temp;
	int epcount1[test_vec_num];	// the num of uncorrect codeword

	//initialize
	for(i=0;i<test_vec_num;i++)
	{
		for(j=0;j<eta;j++)
			test_vec_check[i][j] = -1;

		epcount1[i] = 0;
	}

	//start
	for(i=0;i<test_vec_num;i++)
	{
		mask = 1;
		for(j=0;j<eta;j++)
		{
			if( (i&mask) > 0 )
				choosen_temp = 1;
			else
				choosen_temp = 0;

			index_temp = n-j-1;
			test_vec_check[i][eta-1-j] = large_vec[choosen_temp][(int)test_set_ordered[1][index_temp]];
			mask = mask << 1;
		}
	}

	//error count
	for(j=0;j<n-eta;j++)
		if( codeword[(int)test_set_ordered[1][j]] != test_vec_com[j] )
			epcount1[0]++;

	for(i=1;i<test_vec_num;i++)
		epcount1[i] = epcount1[0];

	for(i=0;i<test_vec_num;i++)
	{
		for(j=0;j<eta;j++)
		{
			if( codeword[(int)test_set_ordered[1][n-eta+j]] != test_vec_check[i][j] )
				epcount1[i]++;
		}
	}

	//************strange part****************
	for(j=0;j<test_vec_num;j++)
	{
		flag=0;
		for(i=0;i<k;i++)
			if( message[i]==message_temp[j][i] )
				flag++;

		if(flag==k)
			break;
	}
	//****************************
	flag = 0;
	for(j=0;j<test_vec_num;j++)
	{
		if(epcount1[j] <= able_correct_num)
			{
				for(u=0;u<n;u++)
				{
					dec_codeword[u] = codeword[u]; 
				}
				for(u=0;u<p*n;u++)
				{
					dec_bicodeword[u]=bi_codeword[u];
				}

				flag = 1;

				break;
			}
	}

	if( flag == 0)	//there is not codeword can be correct
	{
		for(u=0;u<n;u++)
		{
			dec_codeword[u] = large_vec[0][u];
		}

		//Convert the codeword into binary
		for(i=0;i<n;i++)
			{
				value=dec_codeword[i];
				mask=1;
				for(v=0;v<p;v++)
				{
					if((value & mask)>0)
						dec_bicodeword[p*i+v]=1;
					else
						dec_bicodeword[p*i+v]=0;
					mask=mask<<1;
				}
			}
	}


}
