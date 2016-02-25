#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//************variable define***************
#define n 15	//length of codeword
#define k 9 //length of message
#define q 16 //GF(q) q=n+1
#define p 4 //GF(2^p) = GF(q)
#define m 1 // the multiplicity of LCC algorithm is 1

#define V 2 //projected basic function
#define pointNum 4 // the number of constell point == 2^V

//*******************************************

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

//**************basic var*****************
unsigned long int seq_num;	//number of input binary sequences
float SNR;
double BER, FER;
int mularray[]={1,2,4,8,3,6,12,11,5,10,7,14,15,13,9}; 	//this array is used to the infinite field element mul
int a[15]={1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};		//caution: this array must change with different sequence of infinite element
int root[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};			//n+1		be used to the factorization

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
//***************

//************demodulation var*************
float RM[q][n];   // reliability matrix



void main()
{
	int i,j,u;
    for(i=0;i<n*2;i++)
	{
		rx_symbol[i][0]=-1.0;
		rx_symbol[i][1]=1.0;
	}

	N0=1.0;

	demodulation();

    for(i=0;i<q;i++)
    {
    	for(j=0;j<n;j++)
    		printf("%f\t",RM[i][j]);

    	printf("\n");
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
	
	f=e[0]*8+e[1]*4+e[2]*2+e[3]*1; //conver to decimal
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
		tx_symbol[i][0]=-2*bi_codeword[i*2+0]+1;	//0-->1
		tx_symbol[i][1]=-2*bi_codeword[i*2+1]+1;  //1-->-1
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
{	//soft decision
	int i,j,u,v;
	float sum,proba[pointNum];  
	float Pr[2][4];

	for(i=0;i<n;i++)
	{	
		for(j=0;j<V;j++)
		{
			sum=0;
			for(u=0;u<pointNum;u++)
			{	
				proba[u]=PDF(i*V+j,u);
				sum+=proba[u];
			}

			for(u=0;u<pointNum;u++)
				Pr[j][u]=proba[u]/sum;

		}

		for(u=0;u<pointNum;u++)
			for(v=0;v<pointNum;v++)
				RM[u*4+v][i]=(Pr[1][u]*Pr[0][v]);
	}
}


float PDF(int i,int u)
{
	float temp1,temp2,temp3;
	float constell_point[4][2]={ {1,1},{-1,1},{-1,-1},{1,-1} }; // S1{1,1}-->00,S2{-1,1}-->01,S3{-1,-1}-->11,S4{1,-1}-->10

	temp3= 1.0/(pi*N0);
	temp1= pow((rx_symbol[i][0]-constell_point[u][0]),2) + pow((rx_symbol[i][1]-constell_point[u][1]),2);
	temp2= temp3*exp(-(0.5*temp1)/N0);  // 0.5 is the normalization coefficient
	
	return temp2;
}