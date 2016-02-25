/*******************************
(63,15) m=2 simulation programme

*******************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/********change element with different code************
n,k,q,p			m-->iter_num,tm,lm,mono_size,fa_size
a[n]
*****************************************************/

//***********var define*****************
#define n 63	//length of codeword
#define k 15 //length of message
#define q 64 //GF(q) q=n+1
#define p 6 //GF(2^p) = GF(q)
#define m 4 // zero of multiplicity
#define iter_num 630  // C
#define tm 31
#define lm 8 //Desgined factorization output list size   (should be mended)
//********************************

//***************basic function statement**************
int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
void encoder(void);
void modulation(void);
void channel(void);
void demodulation(void);
void decoding(void);
//**********************************************



//**************basic var*****************
unsigned long int seq_num;	//number of input binary sequences
float SNR;
double BER, FER;
int mularray[]={1,2,4,8,16,32,3,6,12,24,48,35,5,10,20,40,
				19,38,15,30,60,59,53,41,17,34,7,14,28,56,51,37,
				9,18,36,11,22,44,27,54,47,29,58,55,45,25,50,39,
				13,26,52,43,21,42,23,46,31,62,63,61,57,49,33};  //this array is used to the infinite field element mul
int a[]={1,2,4,8,16,32,3,6,12,24,48,35,5,10,20,40,
		19,38,15,30,60,59,53,41,17,34,7,14,28,56,51,37,
		9,18,36,11,22,44,27,54,47,29,58,55,45,25,50,39,
		13,26,52,43,21,42,23,46,31,62,63,61,57,49,33};     //caution: this array must change with different sequence of infinite element

int root[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,    
			17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
			33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,
			49,50,51,52,53,54,55,56,57,58,59,60,61,62,63};	//n+1		be used to the factorization

float N0, sgm;

int bi_message[k*p];	//binary input message
int message[k];			//decimal input message
int codeword[n];		//encoded code word
int bi_codeword[n*p];   //encoded binary codeword
float tx_symbol[n*p/2][2];	//modulated tranmitted signal
float rx_symbol[n*p/2][2];	//received destorted signal
int rx_bicodeword[n*p];   //demodulated binary codeword
int rx_codeword[n];   //demodulated codeword to be decoded
int dec_codeword[n*p];  //decoded binary codeword
float pi=3.141593;

//***********************


void main()
{
	int i, x, s,value, count;
	float start,finish;
	unsigned long int j;
	unsigned int mask=1;
	long int error, ferror;
	double progress, b;

	FILE *fp;
	if((fp=fopen("data_(63,15)m=2.txt","a"))==NULL)
	{
		printf("Can't open the data.txt\n");
		exit(0);
	}
	fclose(fp);


	srand(1977);


	printf("\nEnter start SNR: ");
	scanf("%f", &start);
	printf("\nEnter finish SNR: ");
	scanf("%f", &finish);
	printf("\nseq_num: ");
	scanf("%d", &seq_num);


	for(SNR=start;SNR<=finish;SNR++)
	{	
		N0=(1.0/(float(k)/float(n)))/pow(10.0, SNR/10.0);
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
				message[i]=bi_message[p*i]+2*bi_message[p*i+1]+4*bi_message[p*i+2]+8*bi_message[p*i+3]
					   +16*bi_message[p*i+4]+32*bi_message[p*i+5];

/*			//******debug*******
			message[0]=63;	message[1]=15;	message[2]=45;	message[3]=45;	message[4]=20;	message[5]=22;	message[6]=34;	
			message[7]=44;	message[8]=8;	message[9]=51;	message[10]=29;	message[11]=23;	message[12]=5;	message[13]=41;
			message[14]=20;
	
			//*******************
*/
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
			if(f_error_temp==1)
				printf("\nerror=%d\n",error);
			f_error_temp=0;
*/			/*
			if(error!=0)
				ferror++;
			*/

			progress=(double)(j*100)/(double)seq_num;
			
			BER=(double)(error)/(double)(n*p*j);
			//FER=(double)(ferror)/(double)(j),
			printf("Progress=%0.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\r", progress, SNR, error, BER);
//			printf(" Progress=%0.1f,dectect=%d\r",progress,dectec);
			
		}
		BER=(double)error/(double)(n*p*seq_num);
		//FER=(double)(ferror)/(double)(seq_num);
		printf("Progress=%0.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\n",progress, SNR, error, BER);
//		printf(" Progress=%0.1f,dectect=%d\r",progress,dectec);
		fp=fopen("data_(63,15)m=2.txt","a");
		fprintf(fp,"Progress=%0.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\n",progress, SNR, error, BER);
		fclose(fp);
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
		tx_symbol[i][1]=2*bi_codeword[i*2+1]-1;  //1-->1
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

void demodulation(void)  //reconed constellation diagram
{	//hard decision
	int i,j;
	float d[4],temp;
	int point[4][2]={ {1,1} , {0,1} , {0,0} , {1,0} };

	for(i=0;i<(n*p)/2;i++)	
	{
		d[0]=pow( (rx_symbol[i][0]-1),2 ) + pow( (rx_symbol[i][1]-1),2 );   //(1,1)-->11
		d[1]=pow( (rx_symbol[i][0]+1),2 ) + pow( (rx_symbol[i][1]-1),2 ); 	//(-1,1)-->01
		d[2]=pow( (rx_symbol[i][0]+1),2 ) + pow( (rx_symbol[i][1]+1),2 );   //(-1,-1)-->00
		d[3]=pow( (rx_symbol[i][0]-1),2 ) + pow( (rx_symbol[i][1]+1),2 );   //(1,-1)-->10

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
		rx_codeword[i] = rx_bicodeword[p*i]+2*rx_bicodeword[p*i+1]+4*rx_bicodeword[p*i+2]+8*rx_bicodeword[p*i+3]
					 	+16*rx_bicodeword[p*i+4]+32*rx_bicodeword[p*i+5];

}

void decoding(void)
{
	int i,temp;

	temp=0;
	for(i=0;i<n;i++)
		if( rx_codeword[i]!=codeword[i])
			temp++;

	if(temp<=tm)
	{
		for(i=0;i<n*p;i++)
			dec_codeword[i]=bi_codeword[i];
	}
	else if(temp>tm)
	{
		for(i=0;i<n*p;i++)
			dec_codeword[i]=rx_bicodeword[i];
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
	int i,c[p],d[p],e[p],f;
	unsigned int mask=1;

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


	f=e[0]*32+e[1]*16+e[2]*8+e[3]*4+e[4]*2+e[5]*1; //conver to decimal

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

