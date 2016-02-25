#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 7	//length of codeword
#define k 5 //length of message
#define q 8 //GF(q)
#define p 3 //GF(2^p) = GF(q)
#define lm 2 //Desgined factorization output list size
#define maxcost 300 //Designed maximal interpolation cost

int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
void mono_table(void);
void encoder(void);
void modulation(void);
void channel(void);



unsigned long int seq_num;	//number of input binary sequences
int SNR;
double BER, FER;
int mularray[]={1,2,4,3,6,7,5};
int root[]={0,1,2,3,4,5,6,7};	//n+1	
int a[7]={1,2,4,3,6,7,5};

int bi_message[k*p];	//binary input message
int message[k];			//decimal input message
int codeword[n];		//encoded code word
int bi_codeword[n*p];   //encoded binary codeword
int rx_codeword[n*p];   //demodulated binary codeword
int dec_codeword[n*p];  //decoded binary codeword
float pi=3.141593;
float tx_symbol[n*p][2];	//modulated tranmitted signal
float rx_symbol[n*p][2];	//received destorted signal

int mono_order[40][200], channel_type, kv_signal;	//multiplexity matrix
float N0, sgm;


void main()
{
	int i, m, s, start, finish, value, count;
	unsigned long int j;
	unsigned int mask=1;
	long int error, ferror;
	double progress, b;

	srand(1977);

	mono_table();

	printf("\nEnter start SNR: ");
	scanf("%d", &start);
	printf("\nEnter finish SNR: ");
	scanf("%d", &finish);
	printf("\nPlease input the number of k*p-bit sequence: ");
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
			//printf("\n\n*************For the %dth frame*************:\n", j);
			//Generate binary message
			for(i=0;i<k*p;i++)
				bi_message[i]=rand()%2;

			//Convert to nonbinary
			for(i=0;i<k;i++)
				message[i]=bi_message[p*i]+2*bi_message[p*i+1]+4*bi_message[p*i+2];

			encoder();


			//Convert the codeword into binary
			for(i=0;i<n;i++)
			{
				value=codeword[i];
				mask=1;
				for(m=0;m<p;m++) //for(m=p-1;m>=0;m--)
				{
					if((value & mask)>0)
						bi_codeword[p*i+m]=1;
					else
						bi_codeword[p*i+m]=0;
					mask=mask<<1;
				}
			}

			modulation();
			
			channel();
			
			/////////To be done?////////////
			demodulation..?

			decoding..?

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
}

void mono_table(void)
{
	int i, j, v, l;

	for(j=0;j<40;j++) 				
		for(i=0;i<200;i++)
			mono_order[j][i]=-1; 

	j=0;	//represent row
	v=0;	//increasing counter for monomial order
	for(i=0;i<200;i++)
	{
		if(i<(k-1))	
		{
			mono_order[j][i]=v;
			v++;
		}
		else
		{
			l=floor((float)i/(float)(k-1));	
			for(j=0;j<=l;j++)
			{
				mono_order[j][i-(k-1)*j]=v;	//i-v*j
				v++;
			}
		}
	}
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
	int i, j;

	//BPSK
	for(i=0;i<n*p;i++)	
	{
		tx_symbol[i][0]=-1*(2*bi_codeword[i]-1);
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
	int mulresult;

	if(fac1==0||fac2==0)
		mulresult=0;
	else
		for(int i=0;i<7;i++){
			for(int j=0;j<7;j++){
				if(fac1==mularray[i]&&fac2==mularray[j])
					mulresult=mularray[(i+j)%7];
			}
		}

	return mulresult;
}

int add(int fac1,int fac2)
{
	int i,c[3],d[3],e[3],f;
	unsigned int mask=1;

	for(i=2;i>=0;i--){
		if((fac1 & mask)>0)
			c[i]=1;
		else
			c[i]=0;

		mask=mask<<1; //shift 1 bit left
	}
	mask=1;
	for(i=2;i>=0;i--){
		if((fac2 & mask)>0)
			d[i]=1;
		else
			d[i]=0;

		mask=mask<<1;
	}
	for(i=0;i<3;i++)
		e[i]=c[i]^d[i];


	f=e[0]*4+e[1]*2+e[2]*1; //conver to decimal

	return f;
}



int inv(int fac)
{
	int i,invresult;
	if(fac==0)
		invresult=0;
	else
		for(i=0;i<7;i++){
			if(fac==mularray[i])
			invresult=mularray[(7-i)%7];
	}

	return invresult;
}

	
void demodulation(void)
{
	int i;
	float d1,d2;
	for(i=0;i<n*p;i++)	
	{
		d1=(rx_symbol[i][0]-1)*(rx_symbol[i][0]-1)+rx_symbol[i][1]*rx_symbol[i][1];   // a*a= the square of a 
		d2=(rx_symbol[i][0]+1)*(rx_symbol[i][0]+1)+rx_symbol[i][1]*rx_symbol[i][1];
		if(d1<d2)
			rx_codeword[i]=0;
		else
			rx_codeword[i]=1;

	}
}



	
	
	
