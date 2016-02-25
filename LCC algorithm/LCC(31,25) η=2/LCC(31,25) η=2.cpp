/**************************************
	LCC algorithm simulation programme
	RS(31,25) η=2  

****************************************/

/********change element with different code************
n,k,q,p			m-->iter_num,tm,lm,mono_size,fa_size
a[n],mul_array[n]
demodulation()
add(int a, int b)
*****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//************variable define***************
#define n 31	//length of codeword
#define k 25 //length of message
#define q 32 //GF(q) q=n+1
#define p 5 //GF(2^p) = GF(q)
#define m 1 // the multiplicity of LCC algorithm is 1
#define lm 1 // LCC's lm is equal to 1
#define y_size 2	// mono_table size, equal to lm+1
#define x_size 40	// mono_table size, equal to x+5

#define V 1 //projected basic function, BPSK=1, QPSK=2
#define pointNum 2 // the number of constell point == 2^V
#define eta 2 // chosen num
#define test_vec_num 4 //the num of test_vec is 2^eta
#define choose_num 2 //the num of codeword choice, choose 2 in BASIC LCC algorithm


//*********************************************

//***************basic function statement**************
int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
void mono_table(void);
void encoder(void);
void modulation(void);
void channel(void);
//**********************************************
//***************my function***************** 
void demodulation(void);
float PDF(int,int);
void test_set_construction(void);
void reencoding(void);
void cal_T_poly(void);
void cal_v_poly(void);
void interpolation(void);
void com_elem_interpolation(int Q_temp[][lm+1][n+1],int interpoint[][2]);
void uncom_elem_interpolation(int,int,int,int inGroup[][lm+1][n+1],int outGroup1[][lm+1][n+1],int outGroup2[][lm+1][n+1]);
void factorization(void);
void choose(void);


//**************basic var*****************
unsigned long int seq_num;	//number of input binary sequences
float SNR;
double BER, FER;
int mularray[]={1,2,4,8,16,5,10,20,13,26,17,7,14,28,29,31,27,19,3,6,12,24,21,15,30,25,23,11,22,9,18}; 
				//this array is used to the infinite field element mul
int a[]={1,2,4,8,16,5,10,20,13,26,17,7,14,28,29,31,27,19,3,6,12,24,21,15,30,25,23,11,22,9,18};   
		//caution: this array must change with different sequence of infinite element
int root[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};	
			//n+1		be used to the factorization

int bi_message[k*p];	//binary input message
int message[k];			//decimal input message
int codeword[n];		//encoded code word
int bi_codeword[n*p];   //encoded binary codeword
int dec_bicodeword[n*p];  //decoded binary codeword
int dec_codeword[n];    //decoded codeword
float pi=3.141593;
float tx_symbol[n*p][2];	//modulated tranmitted signal
float rx_symbol[n*p][2];	//received destorted signal

int mono_order[y_size][x_size];	//multiplexity matrix
int channel_type, kv_signal;
float N0, sgm;
unsigned long int add_num,mul_num;
//***************

//************demodulation var*************
float RM[q][n];   // reliability matrix
int test_vec[test_vec_num][n]; //store test vectors
float test_set[2][n]; // [0] stores value, [1] stores index
float test_set_ordered[2][n]; // [0] stores value, [1] stores index
int large_vec[choose_num][n]; // largest index is 0, second largest index is 1
//int reencoding_poly[k];
int x_ordered[n];
int T_poly[k]; // store T(x)

int rx_codeword_ordered[test_vec_num][n];
int psi_codeword_temp[n]; 	// store psi 
int reen_codeword_ordered[test_vec_num][n];
int v_poly[k+1]; //store v(x)
int Q_interpoly[test_vec_num][lm+1][n+1];
int message_temp[test_vec_num][k];
int valid_flag[test_vec_num]; //store if the message_temp[j] of fac is vaild

int ferror_temp;



void main()
{
	int i, x,value;
	float start,finish;
	unsigned long int j;
	unsigned int mask=1;
	long int error, ferror;
	double progress, b;
	unsigned long int aver_num;

	FILE *fp;
	if((fp=fopen("LCC_rs(31,25)η=2.txt","a"))==NULL)
	{
		printf("Can't open the data.txt\n");
		exit(0);
	}
	fclose(fp);


	srand(1977);

	mono_table();


	printf("\nEnter start SNR: ");
	scanf("%f", &start);
	printf("\nEnter finish SNR: ");
	scanf("%f", &finish);
	printf("\nseq_num: ");
	scanf("%d", &seq_num);

	for(SNR=start;SNR<=finish;SNR++)
	{	
		mul_num=0;
		add_num=0;
		aver_num=0;

		N0=(1.0/(float(k)/float(n)))/pow(10.0, SNR/10.0);
		sgm=sqrt(N0/2);
		b=1.0;
		error=0;
		ferror=0; 
		for(j=1;j<=seq_num;j++)
		{

		//printf("\n\n*************For the %dth frame*************:\n", j);
		//Generate binary message
			for(i=0;i<k*p;i++)
				bi_message[i]=rand()%2;

			//Convert to nonbinary
			for(i=0;i<k;i++)
				message[i]=bi_message[p*i]+2*bi_message[p*i+1]+4*bi_message[p*i+2]+8*bi_message[p*i+3]+16*bi_message[p*i+4];

/*		//*************debug**************
		message[0]=8;	message[1]=11;	message[2]=14;	message[3]=15;	message[4]=7;	message[5]=3;
		message[6]=5;	message[7]=3;	message[8]=1;

*/		//*******************************


			encoder();

/*		//*************debug**********************
				codeword[0]=1;	codeword[1]=0;	codeword[2]=6;	codeword[3]=12;	codeword[4]=15;	codeword[5]=4;
			codeword[6]=2;	codeword[7]=12;	codeword[8]=0;	codeword[9]=8;	codeword[10]=7;	codeword[11]=2;
			codeword[12]=12;	codeword[13]=6;	codeword[14]=1;


*/		//***************************************


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
			
			reencoding();
			
			interpolation();
			
			factorization();
			
			choose();
			////////////////////////////////

			
			//Number of errors calculation
			//error=0;
			for(i=0;i<p*n;i++)
				if(bi_codeword[i]!=dec_bicodeword[i])
					error++;

			//ferror
			int temp=0;
			ferror_temp=1;
			for(i=0;i<n;i++)
				if(codeword[i]==dec_codeword[i])
					temp++;

			if(temp==n)
				ferror_temp=0;
			else if(temp<n)
				ferror_temp=1;
			else if(temp>n)
				printf("\n\ncount ferror error\n\n");

			if(ferror_temp==1)
				ferror++;
/*			
			if(f_error_temp==1)
				printf("\nerror=%d\n",error);
			f_error_temp=0;
*/			/*
			if(error!=0)
				ferror++;
			*/
			aver_num= mul_num/(long int)(j);

			progress=(double)(j*100)/(double)seq_num;
			
			BER=(double)(error)/(double)(n*p*j);
			FER=(double)(ferror)/(double)(j);
			printf("Progress=%0.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, add_num=%ld, mul_num=%ld, aver_num=%ld\r", progress, SNR, error, BER, ferror, FER, add_num, mul_num, aver_num);
//			printf(" Progress=%0.1f,dectect=%d\r",progress,dectec);
			
		}
		BER=(double)error/(double)(n*p*seq_num);
		FER=(double)(ferror)/(double)(seq_num);
		printf("Progress=%0.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, add_num=%ld, mul_num=%ld, aver_num=%ld\n", progress, SNR, error, BER, ferror, FER, add_num, mul_num, aver_num);
//		printf(" Progress=%0.1f,dectect=%d\r",progress,dectec);
		fp=fopen("LCC_rs(31,25)η=2.txt","a");
		fprintf(fp,"Progress=%0.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, add_num=%ld, mul_num=%ld, aver_num=%ld\n",progress, SNR, error, BER, ferror, FER, add_num, mul_num, aver_num);
		fclose(fp);

	}


	getchar();
	getchar();
		
}

void mono_table()
{
	int i;
	
	// 1st row
	for(i=0;i<x_size;i++)
		mono_order[0][i]=2*(i+1);
	
	// 2nd row
	mono_order[1][0]=0;
	for(i=1;i<x_size;i++)
		mono_order[1][i]=2*i-1;

}


int power(int a, int b)
{
	int i, pow_result;

	pow_result=1;
	for(i=0;i<b;i++)
		pow_result=mul(pow_result,a);

	return pow_result;
}

int mul(int fac1,int fac2)
{
	int mulresult=0;

	mul_num++;

	if(fac1==0||fac2==0)
		mulresult=0;
	else
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				if(fac1==mularray[i]&&fac2==mularray[j])
					mulresult=mularray[(i+j)%n];
			}
		}

	return mulresult;
}

int add(int fac1,int fac2)
{
	int i,c[p],d[p],e[p],f;
	unsigned int mask=1;

	add_num++;

	for(i=(p-1);i>=0;i--){
		if((fac1 & mask)>0)
			c[i]=1;
		else
			c[i]=0;

		mask=mask<<1; //shift 1 bit left
	}

	mask=1;
	for(i=(p-1);i>=0;i--){
		if((fac2 & mask)>0)
			d[i]=1;
		else
			d[i]=0;

		mask=mask<<1;
	}

	for(i=0;i<p;i++)  //p=6
		e[i]=c[i]^d[i];
	
	f=e[0]*16+e[1]*8+e[2]*4+e[3]*2+e[4]*1; //conver to decimal
//	f=e[0]*32+e[1]*16+e[2]*8+e[3]*4+e[4]*2+e[5]*1; //conver to decimal

	return f;
}


int inv(int fac)
{
	int i,invresult;
	if(fac==0)
		invresult=0;
	else
		for(i=0;i<n;i++){
			if(fac==mularray[i])
			invresult=mularray[(n-i)%n];
	}

	return invresult;
}



void encoder(void)
{
	int i, j, temp;

	//Encoding
	for(i=0;i<n;i++)	
	{
		codeword[i]=0;
		for(j=0;j<k;j++)	
		{
			temp=mul(message[j], power(a[i], j));
			codeword[i]=add(codeword[i], temp);

		}
	}
}


void modulation(void)
{
	int i;

	//BPSK
	for(i=0;i<n*p;i++)
	{
		tx_symbol[i][0]=-2*bi_codeword[i]+1;	//0-->1  1-->-1
		tx_symbol[i][1]=0;
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

			Pr[j][0]=proba[0]/sum;	//proba-->0
			Pr[j][1]=proba[1]/sum;	//proba-->1 

		}

		j=0;
		for(u=0;u<pointNum;u++)
			for(v=0;v<pointNum;v++)
				for(int x=0;x<pointNum;x++)
					for(int y=0;y<pointNum;y++)
						for(int h=0;h<pointNum;h++)
						{
							RM[j][i]=(Pr[4][u]*Pr[3][v]*Pr[2][x]*Pr[1][y]*Pr[0][h]);  
							j++;
						}
	
	}

/*	//*****debug**********
	printf("\n\ncodeword");
	for(v=0;v<n;v++)
		printf("\t%d",codeword[v]);
	

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

void test_set_construction(void)
{	//construct the test_vec[test_vec_num]
	int i,j,temp_index;
	float temp;
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
	for(j=0;j<test_vec_num;j++)
		for(i=0;i<n-eta;i++)
			test_vec[j][i]=large_vec[0][(int)test_set_ordered[1][i]];

	j=0;
	for(i=0;i<2;i++)
		for(int u=0;u<2;u++)
			{
				test_vec[j][n-eta]	=large_vec[i][(int)test_set_ordered[1][n-eta]];
				test_vec[j][n-eta+1]=large_vec[u][(int)test_set_ordered[1][n-eta+1]];
				j++;
			}

	//construct x_ordered
	for(i=0;i<n;i++)
		x_ordered[i]=a[(int)test_set_ordered[1][i]];

	
/*	//************debug***************
	printf("\n\n");
	printf("large_vec[]");
	for(j=0;j<n;j++)
		printf("\t%d",large_vec[0][j]);

	printf("\n\n");
	printf("large_vec2[]");
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
	

	for(j=0;j<test_vec_num;j++)
	{
		printf("\n\ntest_vec[%d]\t",j);
		for(i=0;i<n;i++)
			printf("\t%d",test_vec[j][i]);
	}	
*/	//********************************

}


void reencoding(void)
{
	int i,j,u,temp1;


	for(j=0;j<test_vec_num;j++)
		for(i=0;i<n;i++)
			rx_codeword_ordered[j][i]=test_vec[j][i];

	cal_T_poly();

	//cal psi, the codeword_temp
	for(i=0;i<n;i++)
	{	
		psi_codeword_temp[i]=0;
		for(u=0;u<k;u++)
		{
			temp1=mul( T_poly[u],power( x_ordered[i],u ) );
			psi_codeword_temp[i]=add( psi_codeword_temp[i],temp1 );
		}
	}			


	//cal y'=y+psi
	for(j=0;j<test_vec_num;j++)
		for(i=0;i<n;i++)
//			if(i<k)
//				reen_codeword_ordered[j][i]=0;
//			else
				reen_codeword_ordered[j][i]=add(rx_codeword_ordered[j][i],psi_codeword_temp[i]);

/*	//*******debug*****************
	printf("\n\npsi_codeword_temp\t");
	for(i=0;i<n;i++)
		printf("\t%d",psi_codeword_temp[i]);

	for(j=0;j<test_vec_num;j++)
	{
		printf("\n\nreen_codeword_ordered[%d]",j);
		for(i=0;i<n;i++)
			printf("\t%d",reen_codeword_ordered[j][i]);
	}

	printf("\n\n");
*/	//*********************************

}


void cal_T_poly(void)
{
	int i,u,v;
	int temp1,temp2;
	int poly[k],poly_temp1[k+1],poly_temp2[k];

	//initialize
	for(i=0;i<k;i++)
		T_poly[i]=0;

	for(u=0;u<k;u++)
	{
		//initialize
		for(i=0;i<k;i++)
			poly[i]=poly_temp2[i]=0;
		for(i=0;i<k+1;i++)
			poly_temp1[i]=0;

		poly[0]=1;
		temp1=1;
		for(v=0;v<k;v++)
			if(v!=u)
			{
				for(i=0;i<k+1;i++)
					poly_temp1[i]=0;
				for(i=0;i<k;i++)
					poly_temp2[i]=0;

				//calculate poly
				for(i=0;i<k;i++)
					if(poly[i]!=0)
					{
						poly_temp1[i+1]=poly[i];
						poly_temp2[i]=mul(x_ordered[v],poly[i]);
					}

				for(i=0;i<k;i++)
					poly[i]=add( poly_temp1[i],poly_temp2[i] );

				//calculate the factor of denominator
				temp1=mul( temp1,add(x_ordered[u],x_ordered[v]) );

			}

		temp2=mul( rx_codeword_ordered[0][u],inv(temp1) );
		for(i=0;i<k;i++)
			poly[i]=mul( temp2,poly[i] );

		for(i=0;i<k;i++)
			T_poly[i]=add( T_poly[i],poly[i] );

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

		temp1=add( temp1,rx_codeword_ordered[0][i] );

		if( temp1!=0 )
			printf("\n\n\tT_poly is error");
	}
*/	//***********************************

}

void cal_v_poly(void)
{
	int i,u;
	int poly_temp1[k+2],poly_temp2[k+1];

	//initialize
	for(u=0;u<k+1;u++)
		v_poly[u]=0;	


	v_poly[0]=1;
	for(i=0;i<k;i++)
	{
		for(u=0;u<k+1;u++)
		{
			poly_temp1[u]=0;
			poly_temp2[u]=0;
		}
		poly_temp1[k+1]=0;
		
		for(u=0;u<k+1;u++)
			if( v_poly[u]!=0 )
			{
				poly_temp1[u+1]=v_poly[u];
				poly_temp2[u]=mul( x_ordered[i],v_poly[u] );
			}

		for(u=0;u<k+1;u++)
			v_poly[u]=add( poly_temp1[u],poly_temp2[u] );
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
	int i,j,u,v,temp,temp1,temp2;
	int com_elem_interpoint[n-k-eta][2];
	int Q_com_elem[lm+1][lm+1][n+1];
	int uncom_elem_interpoint[eta][3];
	int Q_uncom_elem[test_vec_num][lm+1][lm+1][n+1];
	int degree_temp[test_vec_num][2];

	//common element interpolation

	//initialize
	for(j=0;j<lm+1;j++)
		for(u=0;u<lm+1;u++)
			for(v=0;v<n+1;v++)
				Q_com_elem[j][u][v]=0;

	//set common element interpoint (xi,ri)
	cal_v_poly();

	for(i=k;i<n-eta;i++)
	{
		temp1=0;
		//cal v(ai)
		for(u=0;u<k+1;u++)
		{
			temp=mul( v_poly[u], power( x_ordered[i],u ) );
			temp1=add( temp1,temp );
		}

		//****debug***
//		printf("\nv(a%d)=%d",i,temp1);
		//***********

		com_elem_interpoint[i-k][0]=x_ordered[i];
		com_elem_interpoint[i-k][1]=mul( reen_codeword_ordered[0][i],inv(temp1) );
	}

/*	//*********debug***********
	printf("\n\ncom_elem_interpoint\n");
	for(i=0;i<n-k-eta;i++)
		printf("\t%d",com_elem_interpoint[i][0]);
	printf("\n");
	for(i=0;i<n-k-eta;i++)
		printf("\t%d",com_elem_interpoint[i][1]);
	
*/	//***************************


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


	//uncommon element interpolation
	//set uncomon element interpolation
	for(i=n-eta;i<n;i++)
	{
		temp1=0;
		//cal v(ai)
		for(u=0;u<k+1;u++)
		{
			temp=mul( v_poly[u], power( x_ordered[i],u ) );
			temp1=add( temp1,temp );
		}

		//****debug***
//		printf("\nv(a%d)=%d",i,temp1);
		//***********

		uncom_elem_interpoint[i-n+eta][0]=x_ordered[i];
		uncom_elem_interpoint[i-n+eta][1]=mul( add( large_vec[0][(int)test_set_ordered[1][i]],psi_codeword_temp[i] ),inv(temp1) );
		uncom_elem_interpoint[i-n+eta][2]=mul( add( large_vec[1][(int)test_set_ordered[1][i]],psi_codeword_temp[i] ),inv(temp1) );
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
	for(j=0;j<test_vec_num;j++)
		for(i=0;i<lm+1;i++)
			for(u=0;u<lm+1;u++)
				for(v=0;v<n+1;v++)
					Q_uncom_elem[j][i][u][v]=0;

	//interpolation
	//i=n-eta
	i=0;
	uncom_elem_interpolation(uncom_elem_interpoint[i][0],uncom_elem_interpoint[i][1],uncom_elem_interpoint[i][2],Q_com_elem,Q_uncom_elem[0],Q_uncom_elem[2]);
	
	//i=n-eta+1
	i=1;
	uncom_elem_interpolation(uncom_elem_interpoint[i][0],uncom_elem_interpoint[i][1],uncom_elem_interpoint[i][2],Q_uncom_elem[0],Q_uncom_elem[0],Q_uncom_elem[1]);
	uncom_elem_interpolation(uncom_elem_interpoint[i][0],uncom_elem_interpoint[i][1],uncom_elem_interpoint[i][2],Q_uncom_elem[2],Q_uncom_elem[2],Q_uncom_elem[3]);

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



/*	//***********debug****************
	for(j=0;j<test_vec_num;j++)
		for(i=0;i<lm+1;i++)
		{
			printf("\n\nQ_uncom_elem[%d][%d]",j,i);
			for(u=0;u<lm+1;u++)
			{
				printf("\n");
				for(v=0;v<n+1;v++)
					printf("\t%d",Q_uncom_elem[j][i][u][v]);
			}
		}
*/	//********************************

	//initialize
	for(j=0;j<test_vec_num;j++)
		for(u=0;u<lm+1;u++)
			for(v=0;v<n+1;v++)
				Q_interpoly[j][u][v]=0;

	//choose the poly to factorization
	for(j=0;j<test_vec_num;j++)
	{
		//calculate the degree of poly
		for(i=0;i<2;i++)
		{	
			temp=-1;
			degree_temp[j][i]=0;
			for(u=0;u<lm+1;u++)
				for(v=0;v<n+1;v++)
					if( Q_uncom_elem[j][i][u][v]!=0 )
						if(temp < mono_order[u][v] )
						{
							temp=mono_order[u][v];
							degree_temp[j][i]=v-u;
						}
		}
		//choose the poly with the min degree
		if(degree_temp[j][0]<degree_temp[j][1])
		{
			for(u=0;u<lm+1;u++)
				for(v=0;v<n+1;v++)
					if(Q_uncom_elem[j][0][u][v]!=0)
						Q_interpoly[j][u][v]=Q_uncom_elem[j][0][u][v];
		}
		else if(degree_temp[j][0]>degree_temp[j][1])
		{
			for(u=0;u<lm+1;u++)
				for(v=0;v<n+1;v++)
					if(Q_uncom_elem[j][1][u][v]!=0)
						Q_interpoly[j][u][v]=Q_uncom_elem[j][1][u][v];			
		}
		else if(degree_temp[j][0]==degree_temp[j][1])
		{
			for(u=0;u<lm+1;u++)
				for(v=0;v<n+1;v++)
					if(Q_uncom_elem[j][1][u][v]!=0)
						Q_interpoly[j][u][v]=Q_uncom_elem[j][1][u][v];

			printf("in choose, deg[%d][0]==deg[%d][1]",j,j);
		}


	}


/*	//*************debug*****************
	for(j=0;j<test_vec_num;j++)
	{
		printf("\n\nQ_interpoly[%d]",j);
		for(u=0;u<lm+1;u++)
		{
			printf("\n");
			for(v=0;v<n+1;v++)
				printf("\t%d",Q_interpoly[j][u][v]);
		}
	}
*/	//***********************************




}


void com_elem_interpolation(int Q_temp[][lm+1][n+1],int interpoint[][2])
{
	int i,j,u,v;
	int temp1,temp2,result_temp[2],deg_temp[2],update_poly_index,index_temp[2];
	int poly_temp1[lm+1][n-k-eta+2],poly_temp2[lm+1][n-k-eta+2];


	Q_temp[0][0][0]=1;
	Q_temp[1][1][0]=1;

	//interpolation
	for(i=0;i<n-k-eta;i++)
	{
		//compute
		for(j=0;j<2;j++)
		{
			temp2=0;
			for(u=0;u<lm+1;u++)
				for(v=0;v<n-k-eta+1;v++)
					if(Q_temp[j][u][v]!=0)
					{
						temp1=mul( power(interpoint[i][0],v),power(interpoint[i][1],u) );
						temp1=mul( Q_temp[j][u][v],temp1 );
						temp2=add( temp2, temp1 );
					}

			result_temp[j]=temp2;
		}

		//choose
		update_poly_index=-1;

		for(j=0;j<2;j++)
			{
				index_temp[0]=-1;
				index_temp[1]=-1;
				temp1=-1;
				for(u=0;u<lm+1;u++)
					for(v=0;v<n-k-eta+1;v++)
						if(Q_temp[j][u][v]!=0)
							if(temp1 < mono_order[u][v] )
							{
								temp1=mono_order[u][v];
								index_temp[1]=u;
								index_temp[0]=v;
							}

				deg_temp[j]= index_temp[0]-index_temp[1];

			}

		if( result_temp[0]!=0 && result_temp[1]!=0 )
		{
			if( deg_temp[0]<deg_temp[1] )
				update_poly_index=0;
			else if( deg_temp[0]>deg_temp[1] )
				update_poly_index=1;
			else if( deg_temp[0]==deg_temp[1] )
				update_poly_index=0;  // arbitrarily choose
				
		}
		else if( result_temp[0]!=0 && result_temp[1]==0 )
				update_poly_index=0;
		else if( result_temp[0]==0 && result_temp[1]!=0 )
				update_poly_index=1;
		else if( result_temp[0]==0 && result_temp[1]==0 )
				update_poly_index=-1;

		//update
		if( update_poly_index!=-1 )
		{
			//g(x,y)
			for(j=0;j<2;j++)
				if(j!=update_poly_index)
				{

					temp1=mul( result_temp[j],inv(result_temp[update_poly_index]) );

					if(temp1!=0)
					{	
						for(u=0;u<lm+1;u++)
							for(v=0;v<(n-k-eta+1);v++)
								Q_temp[j][u][v]=add( Q_temp[j][u][v],mul(Q_temp[update_poly_index][u][v],temp1) );
					}
				}

			//f(x,y)
			for(u=0;u<(lm+1);u++)
				for(v=0;v<(n-k-eta+1);v++)
					poly_temp1[u][v]=mul( interpoint[i][0],Q_temp[update_poly_index][u][v] );

			poly_temp2[0][0]=0;
			poly_temp2[1][0]=0;

			for(u=0;u<lm+1;u++)
				for(v=0;v<(n-k-eta+1);v++)
					poly_temp2[u][v+1]=Q_temp[update_poly_index][u][v];

			for(u=0;u<lm+1;u++)
				for(v=0;v<(n-k-eta+1);v++)
					Q_temp[update_poly_index][u][v]=add( poly_temp1[u][v],poly_temp2[u][v] );

		}

/*		//********debug***************
		for(int h=0;h<lm+1;h++)
		{
			temp1=0;
			for(u=0;u<lm+1;u++)
				for(v=0;v<n+1;v++)
					if( Q_temp[h][u][v]!=0 )
					{
						temp2=mul( power(interpoint[i][0],v),power(interpoint[i][1],u) );
						temp2=mul( temp2,Q_temp[h][u][v] );
						temp1=add( temp1,temp2);
					}

			if( temp1!=0 )
				printf("\ninterpoint[%d], poly[%d] is error",i,h);
//			else if(temp1==0)
//				printf("\nsucceed");
		
		}


*/		//***************************





	}	

/*	//***********debug******************
	for(j=0;j<n-k-eta;j++)
	{
		for(i=0;i<lm+1;i++)
		{
			temp1=0;
			for(u=0;u<lm+1;u++)
				for(v=0;v<n+1;v++)
					if( Q_temp[i][u][v]!=0 )
					{
						temp2=mul( power(interpoint[j][0],v),power(interpoint[j][1],u) );
						temp2=mul( temp2,Q_temp[i][u][v] );
						temp1=add( temp1,temp2);
					}

			if( temp1!=0 )
				printf("\ninterpoint[%d], poly[%d] is error",j,i);
//			else if(temp1==0)
//				printf("\nsucceed");
		
		}
	}
*/	//*********************************

}

void uncom_elem_interpolation(int a,int y1_HD,int y2_HD,int inGroup[][lm+1][n+1],int outGroup1[][lm+1][n+1],int outGroup2[][lm+1][n+1])
{
	int i,j,u,v;
	int c0_0,c0_1,c1_0,c1_1;
	int result_temp[2][2],index_temp[2],deg_temp[2],update_poly_index[2];
	int temp1,temp2;
	int poly_temp1[lm+1][n+2],poly_temp2[lm+1][n+2];
	int inGroup_temp[lm+1][lm+1][n+1];



	//initialize
	for(j=0;j<lm+1;j++)
		for(u=0;u<lm+1;u++)
			for(v=0;v<n+1;v++)
				inGroup_temp[j][u][v]=inGroup[j][u][v];

	//compute
	c0_0=0;
	c0_1=0;
	c1_0=0;
	c1_1=0;
	for(v=0;v<n+1;v++)
	{
		c0_0=add( c0_0,mul( inGroup[0][0][v],power(a,v) ) );
		c0_1=add( c0_1,mul( inGroup[0][1][v],power(a,v) ) );

		c1_0=add( c1_0,mul( inGroup[1][0][v],power(a,v) ) );
		c1_1=add( c1_1,mul( inGroup[1][1][v],power(a,v) ) );
	}

	result_temp[0][0]=add( c0_0,mul(c0_1,y1_HD) );
	result_temp[0][1]=add( c1_0,mul(c1_1,y1_HD) );
	result_temp[1][0]=add( c0_0,mul(c0_1,y2_HD) );
	result_temp[1][1]=add( c1_0,mul(c1_1,y2_HD) );

	//choose 
	update_poly_index[0]=-1;
	update_poly_index[1]=-1;

	for(j=0;j<2;j++)
	{
		index_temp[0]=-1;
		index_temp[1]=-1;
		temp1=-1;	
		for(u=0;u<lm+1;u++)
			for(v=0;v<n-k-eta+1;v++)
				if(inGroup[j][u][v]!=0)
					if(temp1 < mono_order[u][v] )
					{
						temp1=mono_order[u][v];
						index_temp[1]=u;
						index_temp[0]=v;
					}

		deg_temp[j]= index_temp[0]-index_temp[1];

	}	

	for(j=0;j<2;j++)
	{
		if( result_temp[j][0]!=0 && result_temp[j][1]!=0 )
		{
			if( deg_temp[0]<deg_temp[1] )
				update_poly_index[j]=0;
			else if( deg_temp[0]>deg_temp[1] )
				update_poly_index[j]=1;
			else if( deg_temp[0]==deg_temp[1] )
				update_poly_index[j]=0;
		}
		else if( result_temp[j][0]!=0 && result_temp[j][1]==0 )
				update_poly_index[j]=0;
		else if( result_temp[j][0]==0 && result_temp[j][1]!=0 )
				update_poly_index[j]=1;
		else if( result_temp[j][0]==0 && result_temp[j][1]==0 )
				update_poly_index[j]=-1;
	}
	
	//update
	if( update_poly_index[0]!=-1 )	
	{
		//update g_1(x,y)
		for(i=0;i<2;i++)
			if(i!=update_poly_index[0])
			{
				temp1=mul( result_temp[0][i],inv( result_temp[0][update_poly_index[0]] ) );
				if(temp1!=0)
				{
					for(u=0;u<lm+1;u++)
						for(v=0;v<n+1;v++)
							outGroup1[i][u][v]=add( inGroup_temp[i][u][v],mul(temp1,inGroup_temp[ update_poly_index[0] ][u][v] ) );
				}
				else if( temp1==0 )
				{
					for(u=0;u<lm+1;u++)
						for(v=0;v<n+1;v++)
							outGroup1[i][u][v]=inGroup_temp[i][u][v];
				}
			}

		for(u=0;u<lm+1;u++)
			for(v=0;v<n+1;v++)
				poly_temp1[u][v]=poly_temp2[u][v]=0;

		//update f_1(x,y),update_poly_index
		for(u=0;u<(lm+1);u++)
			for(v=0;v<(n+1);v++)
				if( inGroup_temp[update_poly_index[0]][u][v]!=0 )
					poly_temp1[u][v]=mul( a,inGroup_temp[update_poly_index[0]][u][v] );

		poly_temp2[0][0]=0;
		poly_temp2[1][0]=0;

		for(u=0;u<lm+1;u++)
			for(v=0;v<(n+1);v++)
				poly_temp2[u][v+1]=inGroup_temp[update_poly_index[0]][u][v];

		for(u=0;u<lm+1;u++)
			for(v=0;v<n+1;v++)
				outGroup1[update_poly_index[0]][u][v]=add( poly_temp1[u][v],poly_temp2[u][v]);
	}

	if( update_poly_index[1]!=-1 )
	{
		//update g_2(x,y)
		for(i=0;i<2;i++)
			if(i!=update_poly_index[1])
			{
				temp1=mul( result_temp[1][i],inv( result_temp[1][update_poly_index[1]] ) );
				if( temp1!=0 )
				{
					for(u=0;u<lm+1;u++)
						for(v=0;v<n+1;v++)
							outGroup2[i][u][v]=add( inGroup_temp[i][u][v],mul(temp1,inGroup_temp[ update_poly_index[1] ][u][v] ) );
				}
				else if( temp1==0 )
				{
					for(u=0;u<lm+1;u++)
						for(v=0;v<n+1;v++)
							outGroup2[i][u][v]=inGroup_temp[i][u][v];
				}
			}

		for(u=0;u<lm+1;u++)
			for(v=0;v<n+1;v++)
				poly_temp1[u][v]=poly_temp2[u][v]=0;

		//update f_2(x,y)
		for(u=0;u<(lm+1);u++)
			for(v=0;v<(n+1);v++)
				if( inGroup_temp[update_poly_index[1]][u][v]!=0 )
				poly_temp1[u][v]=mul( a,inGroup_temp[update_poly_index[1]][u][v] );

		poly_temp2[0][0]=0;
		poly_temp2[1][0]=0;

		for(u=0;u<lm+1;u++)
			for(v=0;v<(n+1);v++)
				poly_temp2[u][v+1]=inGroup_temp[update_poly_index[1]][u][v];

		for(u=0;u<lm+1;u++)
			for(v=0;v<n+1;v++)
				outGroup2[update_poly_index[1]][u][v]=add( poly_temp1[u][v],poly_temp2[u][v] );
	}


/*	//*******debug************
	for(j=0;j<lm+1;j++)
	{
		temp1=0;
		for(u=0;u<lm+1;u++)
			for(v=0;v<n+1;v++)
				if( outGroup1[j][u][v]!=0 )
				{
					temp2=mul( power(a,v),power(y1_HD,u) );
					temp2=mul( outGroup1[j][u][v],temp2 );
					temp1=add( temp1,temp2 );
				}

		if(temp1!=0)
			printf("\ninterpoint(%d,%d),Group1 is error");
//		else if(temp1==0)
//			printf("\nsucceed");
	}


	for(j=0;j<lm+1;j++)
	{
		temp1=0;
		for(u=0;u<lm+1;u++)
			for(v=0;v<n+1;v++)
				if( outGroup2[j][u][v]!=0 )
				{
					temp2=mul( power(a,v),power(y2_HD,u) );
					temp2=mul( outGroup2[j][u][v],temp2 );
					temp1=add( temp1,temp2 );
				}

		if(temp1!=0)
			printf("\ninterpoint(%d,%d), Group2 is error");
//		else if(temp1==0)
//			printf("\nsucceed");
	}
*/	//*********************

}

	
void factorization(void)
{
	int i,j,u,v;
	int X_poly_size,X_temp_size;
//	int Vpoly[k+1],Vpoly_temp1[k+2],Vpoly_temp2[k+1];
	int X_poly[test_vec_num][n+k+1],X_temp[n+k+1+k];
	int divisor[2],dividend[2],quo[2];
	int div_poly[n+k+1];

	int test_poly[test_vec_num][lm+1][n+k+1];


	X_poly_size=n+k+1;
	X_temp_size=n+k+1+k;

	//initialize
	for(j=0;j<test_vec_num;j++)
		for(v=0;v<k;v++)
			message_temp[j][v]=0;

	//calculate message_temp poly
	for(j=0;j<test_vec_num;j++)
	{
		for(v=0;v<X_poly_size;v++)
			X_poly[j][v]=0;


		//calculate the real q0(x,y) * v(x)
		for(i=0;i<k+1;i++)
			if(v_poly[i]!=0)
			{
				for(v=0;v<X_temp_size;v++)
					X_temp[v]=0;

				for(v=0;v<n+1;v++)
					if(Q_interpoly[j][0][v]!=0)
						X_temp[v+i]=mul( v_poly[i],Q_interpoly[j][0][v] );
	
				for(v=0;v<X_poly_size;v++)
					X_poly[j][v]=add( X_poly[j][v],X_temp[v] );
			}

/*		//**********debug***********
		printf("\n\nX_poly[%d]",j);
		for(v=0;v<X_poly_size;v++)
			printf("\t%d",X_poly[j][v]);
*/		//*****************************

		//*********debug*************
		for(u=0;u<lm+1;u++)
			for(v=0;v<X_poly_size;v++)
				test_poly[j][u][v]=0;

		for(v=0;v<X_poly_size;v++)
			test_poly[j][0][v]=X_poly[j][v];
		for(v=0;v<n+1;v++)
			test_poly[j][1][v]=Q_interpoly[j][1][v];
		//***************************



		//calculate the div of poly
		//find the divisor[2]
		divisor[0]=divisor[1]=0;
		for(v=n;v>=0;v--)
			if(Q_interpoly[j][1][v]!=0)
			{
				divisor[0]=Q_interpoly[j][1][v];
				divisor[1]=v;
				break;
			}

		//processing
		quo[0]=0;
		quo[1]=0;
		valid_flag[j]=-1;
		for(i=0;i<k+1;i++)
		{
			//initialize
			for(v=0;v<X_poly_size;v++)
				div_poly[v]=0;


			dividend[0]=0;
			dividend[1]=0;
			//find the dividend[2]
			for(v=X_poly_size-1;v>=0;v--)
				if(X_poly[j][v]!=0)
				{
					dividend[0]=X_poly[j][v];
					dividend[1]=v;
					break;
				}
			
			//effective solution, no remainder
			if( v==-1 && X_poly[j][0]==0 )
			{
				valid_flag[j]=1;
				break;
			}
			//uneffective solution, remainder
			if( dividend[1]<divisor[1] )
			{	valid_flag[j]=0;
				break;
			}
			//find the quotient
			quo[0]=mul( dividend[0],inv(divisor[0]) );
			quo[1]=dividend[1]-divisor[1];

			if( (0<=quo[1]) && (quo[1]<=(k-1)) )
			{
				message_temp[j][quo[1]]=quo[0];

				// div_poly=divisor*quotient
				for(v=0;v<n+1;v++)
					if(Q_interpoly[j][1][v]!=0)
						div_poly[v+quo[1]]=mul( quo[0],Q_interpoly[j][1][v] );

				//dividend-div_poly
				for(v=0;v<X_poly_size;v++)
					X_poly[j][v]=add( X_poly[j][v],div_poly[v] );
			}
			else if( (quo[1]<0) || (quo[1]>(k-1)) )
			{
				valid_flag[j]=0;
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
	int i,j,v,u,value;
	unsigned int mask=1;
	int temp,min_num=-1;
	int codeword_temp[test_vec_num][n],bicodeword_temp[test_vec_num][n*p],bicodeword_coor[n*p];
	int tx_symbol_temp[test_vec_num][(n*p)/2][2],coordin_point[(n*p)/2][2];


	//Encoding
		//calculate the real message_temp
		for(j=0;j<test_vec_num;j++)
			for(v=0;v<k;v++)
				if(valid_flag[j]==1)
					message_temp[j][v]=add(message_temp[j][v],T_poly[v]);

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
		



	for(j=0;j<test_vec_num;j++)
		if(valid_flag[j]==1)
			for(i=0;i<n;i++)
			{
				codeword_temp[j][i]=0;
				for(v=0;v<k;v++)
					{
						temp=mul(message_temp[j][v],power(a[i],v));
						codeword_temp[j][i]=add(codeword_temp[j][i],temp);
					}

			}


	//Convert the codeword into binary
	for(j=0;j<test_vec_num;j++)
		if(valid_flag[j]==1)
			for(i=0;i<n;i++)
				{
					value=codeword_temp[j][i];
					mask=1;
					for(v=0;v<p;v++)
					{
						if((value & mask)>0)
							bicodeword_temp[j][p*i+v]=1;
						else
							bicodeword_temp[j][p*i+v]=0;
						mask=mask<<1;
					}
				}
	
	//rx_bicodeword, bicodeword_coor
	for(j=0;j<test_vec_num;j++)
		if(valid_flag[j]==1)
			for(i=0;i<n;i++)
				{
					value=large_vec[0][i];
					mask=1;
					for(v=0;v<p;v++)
					{
						if((value & mask)>0)
							bicodeword_coor[p*i+v]=1;
						else
							bicodeword_coor[p*i+v]=0;
						mask=mask<<1;
					}
				}


/*	//***********debug***************
	printf("\n\ntx_codeword\t\t");
	for(i=0;i<n;i++)
		printf("\t%d",codeword[i]);


	printf("\n\nrx_codeword\t\t");
	for(i=0;i<n;i++)
		printf("\t%d",large_vec[0][i]);

	for(j=0;j<test_vec_num;j++)
	{
		printf("\n\nvalid_flag=%d  codeword_temp[%d]",valid_flag[j],j);
		for(i=0;i<n;i++)
			printf("\t%d",codeword_temp[j][i]);
	}

	printf("\n");
*/	//********************************

/*	//mode1: calculate Euclidean distance
		//cal symbol_temp
		for(j=0;j<test_vec_num;j++)
			if(valid_flag[j]==1)
				for(i=0;i<(n*p)/2;i++)
				{
					tx_symbol_temp[j][i][0]=-2*bicodeword_temp[j][i*2+0]+1;	//0-->1
					tx_symbol_temp[j][i][1]=-2*bicodeword_temp[j][i*2+1]+1;  	//1-->-1
				}	


		//calculate the distance
		float temp_x,temp_y;
		float distance_temp[test_vec_num];
		min_num=-1;
		float dis_temp=1000.0;
		for(j=0;j<test_vec_num;j++)
			if(valid_flag[j]==1)
			{
				distance_temp[j]=0.0;
				for(i=0;i<(n*p)/2;i++)
				{
					temp_x=tx_symbol_temp[j][i][0]-rx_symbol[i][0];
					temp_y=tx_symbol_temp[j][i][1]-rx_symbol[i][1];
					distance_temp[j]+=sqrt( (temp_x*temp_x)+(temp_y*temp_y) ) ;				
				}

				if(dis_temp>distance_temp[j])
				{
					dis_temp=distance_temp[j];
					min_num=j;
				}

			}
*/		//***********************

/*		//mode2: posteriori probablility
		float proba[test_vec_num];
		float proba_temp=0.0;
		min_num=-1;
		for(j=0;j<test_vec_num;j++)
			if(valid_flag[j]==1)
			{
				proba[j]=1.0;
				for(i=0;i<n;i++)
				{
					for(u=0;u<q;u++)
						if(codeword_temp[j][i]==root[u])
							proba[j]=proba[j]*RM[u][i];
				}

				if(proba_temp<proba[j])
				{
					proba_temp=proba[j];
					min_num=j;
				}

			}
*/
/*		//*****debug*************
		printf("\n\nmin_num=%d\n\n",min_num);
		for(j=0;j<test_vec_num;j++)
			if(valid_flag[j]==1)
				printf("\nproba[%d]=%f\n",j,proba[j]);
*/		//***********************

		//***********************************

/*		//mode3: ML rule
		int index_c=0;
		int index_r=0;
		float total_temp[test_vec_num];
		for(j=0;j<test_vec_num;j++)
			total_temp[j]=10000.0;

		for(j=0;j<test_vec_num;j++)
			if(valid_flag[j]==1)
			{
				float temp_c=0.0;
				float temp_r=0.0;
				for(i=0;i<n;i++)
					if(codeword_temp[j][i]!=large_vec[0][i])
					{
						for(u=0;u<q;u++)
							if(codeword_temp[j][i]==root[u])
								index_c=u;

						temp_c+=log(RM[index_c][i]);

						for(u=0;u<q;u++)
							if(large_vec[0][i]==root[u])
								index_r=u;

						temp_r+=log(RM[index_r][i]);
					}

				total_temp[j]=temp_r-temp_c;
			}

		min_num=-1;
		float temp1=1000.0;
		for(j=0;j<test_vec_num;j++)
			if(valid_flag[j]==1)
				if(temp1>total_temp[j])
				{
					temp1=total_temp[j];
					min_num=j;
				}
*/		//*****************************************

/*		//*****debug*************
		printf("\n\nmin_num=%d\n\n",min_num);
		for(j=0;j<test_vec_num;j++)
			if(valid_flag[j]==1)
				printf("\ntotal_temp[%d]=%f\n",j,total_temp[j]);
*/		//***********************

		//**********genius: optimal choose***********
		ferror_temp=1;
		for(j=0;j<test_vec_num;j++)
		if(valid_flag[j]==1)
			{	
				int temp=0;
				ferror_temp=1;
				for(i=0;i<n;i++)
					if(codeword[i]==codeword_temp[j][i])
						temp++;

				if(temp==n)
					{
						ferror_temp=0;
						min_num=j;
						break;
					}

				else if(temp!=n)
					{
						ferror_temp=1;
						min_num=-1;
					}
			}

	//***************************






					

	//judge
	if(min_num==-1)
	{	
//		printf("\ndecoding failed\n");

		for(i=0;i<n*p;i++)
			dec_bicodeword[i]=bicodeword_coor[i];

		for(i=0;i<n;i++)
			dec_codeword[i]=large_vec[0][i];

	}
	else if( (min_num>=0) && (min_num<test_vec_num) )
	{
		for(i=0;i<n*p;i++)
			dec_bicodeword[i]=bicodeword_temp[min_num][i];

		for(i=0;i<n;i++)
			dec_codeword[i]=codeword_temp[min_num][i];

/*		//*****debug*************
		printf("\n\nmin_num=%d\n\n",min_num);
		for(j=0;j<test_vec_num;j++)
			if(valid_flag[j]==1)
				printf("\ndistance_temp[%d]=%f\n",j,distance_temp[j]);
*/		//***********************
	}
	else printf("\n\nchoose error, min_num=%d\n",min_num);


	//**********debug*******************
	int flag=0;
	temp=0;
	for(i=0;i<n;i++)
		if(codeword[i]!=large_vec[0][i])
			temp++;

	for(j=0;j<test_vec_num;j++)
	{
		flag=0;
		for(i=0;i<k;i++)
			if( message[i]==message_temp[j][i] )
				flag++;

		if(flag==k)
			break;
	}



//	if(temp>2)
//		printf("\n\nsymbol_error: %d\n\n",temp);

	if( (temp<=2) && (flag!=k) )
		printf("\n\ndecoding failed, temp=%d\n\n",temp);
//		printf("\n\nsymbol_error: %d\n\n",temp);
	//*************************************




}