/****************************************
RS(15,9)code   GS algorithm  cheating bound
*****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//*****************************************
#define n 15	//length of codeword
#define k 9 //length of message
#define q 16 //GF(q) q=n+1
#define p 4 //GF(2^p) = GF(q)
#define tm 3  // error correction capability

//**************************


int mul(int, int);
int add(int, int);
int power(int, int);
void mono_table(void);
void encoder(void);
void modulation(void);
void channel(void);
void demodulation(void);
void decoding(void);


//**************************
unsigned long int seq_num;	//number of input binary sequences
int SNR;
double BER, FER;
int mularray[]={1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};  //this array is used to the infinite field element mul
int a[15]={1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};      //caution: this array must change with different sequence of infinite element

int root[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};	//n+1		be used to the factorization

int bi_message[k*p];	//binary input message
int message[k];			//decimal input message
int codeword[n];		//encoded code word
int bi_codeword[n*p];   //encoded binary codeword
int dec_codeword[n*p];  //decoded binary codeword
float pi=3.141593;
float tx_symbol[n*p/2][2];	//modulated tranmitted signal
float rx_symbol[n*p/2][2];	//received destorted signal

int mono_order[40][200], channel_type, kv_signal;	//multiplexity matrix
float N0, sgm;

//*******************************************
int rx_bicodeword[n*p];   //demodulated binary codeword
int rx_codeword[n];   //demodulated codeword to be decoded




void main()
{
	int i, x, s, start, finish, value, count;
	unsigned long int j;
	unsigned int mask=1;
	long int error, ferror;
	double progress, b;

	srand(1977);

	printf("\nEnter start SNR: ");
	scanf("%d", &start);
	printf("\nEnter finish SNR: ");
	scanf("%d", &finish);
	printf("\nseq_num: ");
	scanf("%d", &seq_num);

	for(SNR=start;SNR<=finish;SNR++)
	{	
		N0=(1.0/(float(k)/float(n)))/pow(10.0, float(SNR)/10.0);
		sgm=sqrt(N0/2);
		b=1.0;
		error=0;
		//ferror=0; 
		for(j=1;j<=seq_num;j++)
		{

			for(i=0;i<k*p;i++)
				bi_message[i]=rand()%2;

			//Convert to nonbinary
			for(i=0;i<k;i++)
				message[i]=bi_message[p*i]+2*bi_message[p*i+1]+4*bi_message[p*i+2]+8*bi_message[p*i+3];

			encoder();
	
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
			demodulation();

			decoding();
			//Number of errors calculation
			//error=0;
			for(i=0;i<p*n;i++)
				if(bi_codeword[i]!=dec_codeword[i])
					error++;
				
			/*
			if(error!=0)
				ferror++;
			*/

			progress=(double)(j*100)/(double)seq_num;
			
			BER=(double)(error)/(double)(n*p*j);
			//FER=(double)(ferror)/(double)(j),
			printf("Progress=%0.1f, SNR=%2.1d, Bit Errors=%2.1d, BER=%E\r", progress, SNR, error, BER);
			
		}
		BER=(double)error/(double)(n*p*seq_num);
		//FER=(double)(ferror)/(double)(seq_num);
		printf("Progress=%0.1f, SNR=%2.1d, Bit Errors=%2.1d, BER=%E\n",progress, SNR, error, BER);

	}

	getchar();
	getchar();
		
}


void encoder(void)
{
	int i, j, codeword_temp;

	//Encoding
	for(i=0;i<n;i++)	
	{
		codeword[i]=0;
		for(j=0;j<k;j++)	
		{
			codeword_temp=mul(message[j], power(a[i], j));
			codeword[i]=add(codeword[i], codeword_temp);

		}
	}
}


void modulation(void)
{
	int i;

	//QPSK
	for(i=0;i<(n*p)/2;i++)
	{
		tx_symbol[i][0]=2*bi_codeword[i*2+0]-1;	//0-->-1
		tx_symbol[i][1]=2*bi_codeword[i*2+1]-1;	//1-->1
	}
}

void channel(void)
{
	int i, j;
	float u, r, g;

	//Add AWGN 
	for(i=0;i<(n*p)/2;i++)	
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
{	//hard decision
	int i,j;
	float d[4];
	float temp;
	int point[4][2]={ {1,1} , {0,1} , {0,0} , {1,0} };

	for(i=0;i<(n*p)/2;i++)	
	{
		d[0]=(rx_symbol[i][0]-1)*(rx_symbol[i][0]-1)+(rx_symbol[i][1]-1)*(rx_symbol[i][1]-1);   //(1,1)-->11
		d[1]=(rx_symbol[i][0]+1)*(rx_symbol[i][0]+1)+(rx_symbol[i][1]-1)*(rx_symbol[i][1]-1);	  //(-1,1)-->01
		d[2]=(rx_symbol[i][0]+1)*(rx_symbol[i][0]+1)+(rx_symbol[i][1]+1)*(rx_symbol[i][1]+1);	  //(-1,-1)-->00
		d[3]=(rx_symbol[i][0]-1)*(rx_symbol[i][0]-1)+(rx_symbol[i][1]+1)*(rx_symbol[i][1]+1);	  //(1,-1)-->10

		temp=d[0];
		for(j=1;j<4;j++)
			if( temp>d[j] )
				temp=d[j];

		for(j=0;j<4;j++)
			if( temp==d[j] )
			{
				rx_bicodeword[i*2+0]=point[j][0];
				rx_bicodeword[i*2+1]=point[j][1];
			}

	}

	//Convert rx_bicodeword to nonbinary's rx_codeword
	for(i=0;i<n;i++)
		rx_codeword[i] = rx_bicodeword[p*i]+2*rx_bicodeword[p*i+1]+4*rx_bicodeword[p*i+2]+8*rx_bicodeword[p*i+3];	

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
	int i,c[4],d[4],e[4],f;
	unsigned int mask=1;

	for(i=3;i>=0;i--){
		if((fac1 & mask)>0)
			c[i]=1;
		else
			c[i]=0;

		mask=mask<<1; //shift 1 bit left
	}
	mask=1;
	for(i=3;i>=0;i--){
		if((fac2 & mask)>0)
			d[i]=1;
		else
			d[i]=0;

		mask=mask<<1;
	}
	for(i=0;i<4;i++)
		e[i]=c[i]^d[i];


	f=e[0]*8+e[1]*4+e[2]*2+e[3]*1; //conver to decimal

	return f;
}



void decoding(void)
{
	int i,j,flag,value;
	int flag_codeword[n];
	unsigned int mask=1;

	for(i=0,flag=0;i<n;i++)
		if( rx_codeword[i]!=codeword[i])
			flag++;

	if(flag<=tm)
		for(i=0;i<n;i++)
			flag_codeword[i]=codeword[i];
	else
		for(i=0;i<n;i++)
			flag_codeword[i]=rx_codeword[i];	


	//Convert the codeword into binary
			for(i=0;i<n;i++)
			{
				value=flag_codeword[i];
				mask=1;
				for(j=0;j<p;j++) //for(m=p-1;m>=0;m--)
				{
					if((value & mask)>0)
						dec_codeword[p*i+j]=1;
					else
						dec_codeword[p*i+j]=0;
					mask=mask<<1;
				}
			}

}