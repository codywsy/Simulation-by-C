/****************************************
QPSK uncoded
*****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//*******************
#define n 15	//length of codeword
#define k 9 //length of message
#define q 16 //GF(q) q=n+1
#define p 4 //GF(2^p) = GF(q)
//**********************************

void modulation(void);
void channel(void);
void demodulation(void);

//******************************
unsigned long int seq_num;	//number of input binary sequences
int SNR;
double BER, FER;

//*******************
int bi_codeword[n*p];   //encoded binary codeword
float pi=3.141593;
float tx_symbol[n*p/2][2];	//modulated tranmitted signal
float rx_symbol[n*p/2][2];	//received destorted signal
float N0, sgm;
int rx_bicodeword[n*p];   //demodulated binary codeword



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
		N0=(1.0/(float(n)/float(n)))/pow(10.0, float(SNR)/10.0);
		sgm=sqrt(N0/2);
		b=1.0;
		error=0;
		//ferror=0; 
		for(j=1;j<=seq_num;j++)
		{

		//printf("\n\n*************For the %dth frame*************:\n", j);
		//Generate binary message
			for(i=0;i<n*p;i++)
				bi_codeword[i]=rand()%2;
	
			modulation();
			
			channel();
			
			/////////To be done?////////////
			demodulation();

			//Number of errors calculation
			//error=0;
			for(i=0;i<n*p;i++)
				if(bi_codeword[i]!=rx_bicodeword[i])
					error++;
				
			/*
			if(error!=0)
				ferror++;
			*/

			progress=(double)(j*100)/(double)seq_num;
			
			BER=(double)(error)/(double)(n*p*j);
			//FER=(double)(ferror)/(double)(j),
			printf("Progress=%0.1f, SNR=%2.1d, Bit Errors=%2.1d, BER=%E\r", progress, SNR, error, BER);
//			printf(" Progress=%0.1f,dectect=%d\r",progress,dectec);
			
		}
		BER=(double)error/(double)(n*p*seq_num);
		//FER=(double)(ferror)/(double)(seq_num);
		printf("Progress=%0.1f, SNR=%2.1d, Bit Errors=%2.1d, BER=%E\n",progress, SNR, error, BER);
//		printf(" Progress=%0.1f,dectect=%d\r",progress,dectec);

	}

	getchar();
	getchar();
		
}

void modulation(void)
{
	int i;

	//QPSK
	for(i=0;i<(n*p)/2;i++)
	{
		tx_symbol[i][0]=2*bi_codeword[i*2+0]-1;
		tx_symbol[i][1]=2*bi_codeword[i*2+1]-1;
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
		d[0]=(rx_symbol[i][0]-1)*(rx_symbol[i][0]-1)+(rx_symbol[i][1]-1)*(rx_symbol[i][1]-1);   //(1,1)-->(1,1)
		d[1]=(rx_symbol[i][0]+1)*(rx_symbol[i][0]+1)+(rx_symbol[i][1]-1)*(rx_symbol[i][1]-1);	  //(-1,1)-->(0,1)
		d[2]=(rx_symbol[i][0]+1)*(rx_symbol[i][0]+1)+(rx_symbol[i][1]+1)*(rx_symbol[i][1]+1);	  //(-1,-1)-->(0,0)
		d[3]=(rx_symbol[i][0]-1)*(rx_symbol[i][0]-1)+(rx_symbol[i][1]+1)*(rx_symbol[i][1]+1);	  //(1,-1)-->(1,0)

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

}