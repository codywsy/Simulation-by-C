#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define OpenFile fp=fopen("KVoptimal_Herm(64,39).txt","a")
#define FrameError 209
#define CheatEncoding

//**********define*****************
#define q 16
#define n 64
#define k 39
#define p 4
#define w 4
#define weiz 44
#define iterNum 64	//when m=1, C is equal to n
#define lm 1
#define able_correct 9
#define pointNum 2
//#define eta 6
//#define test_vec_num 64 //the num of test_vec is 2^eta
#define interval 0.25

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
//#define pbNum 80	//key elment, less than tg_size
//#define pb_Ysize 20	//large than max(degY) + 1
//#define pb_Xsize (w+1)	//w+1
//***********************************


//*************global var**************************
//main()
unsigned long int seq_num;	//number of input binary sequences
float SNR,N0;
double BER, FER;
float pi=3.141593;	// pai
//GF elements
int mularray[] = {1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};  //this array is used to the infinite field element mul
int root[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};	//n+1		be used to the factorization
int logNum[] = {-1,0,1,4,2,8,5,10,3,14,9,7,6,13,11,12};	//used to locate the degree of finitefield element through the value of finitefield element

//*************
int bi_message[k*p], message[k];	//transmitted messge[k]
int codeword[n], bi_codeword[n*p]; //codewords
float tx_symbol[p*n][2], rx_symbol[p*n][2], sgm;
float RM[q][n];
//polebasis()
//int pb[pbNum][pb_Ysize][pb_Xsize];	//pole basis[num_pb][y_size][x_size]
//findpoint()
int point[2][n];	//rational points[2][w^3]
//tgorder()
int tg_order_1[tg_size][2],tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
//mono_table()
int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
int gmatrix[k][n];	//generator matrix[k][n]
int dec_message[k],dec_codeword[n],dec_bicodeword[n*p]; //decoding result
//******************************************************


int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
void findpoints(void);
void mono_table(void);
void tgorder(void);
void generator(void);
void polebasis(void);
void encoder(int message_temp[], int codeword_temp[]);
void modulation(void);
void channel(void);
void demodulation(void);
float PDF(int,int);





void main()
{
	int i, u, m, s, num, value, counter;	
	float start, finish;
	unsigned long int j;
	unsigned int mask=1;
	long int error,ferror;
	double progress, b;

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

	
//	polebasis();

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

#ifdef CheatEncoding
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


	for(SNR=start; SNR<=finish; SNR=SNR+interval)
	{
		N0=(1.0/((float)k/(float)n))/pow(10.0, SNR/10.0);
		sgm=sqrt(N0/2);
		
		error=0;
		ferror=0;

		for(j=1;j<=seq_num;j++)
		{

#ifndef CheatEncoding
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
#endif


			//modulation
			modulation();

			//channel
			channel();

			//demodulation
			demodulation();

			//LIST DECODER

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

			progress=(double)(j*100)/(double)seq_num;
			BER=(double)(error)/(double)(n*p*j);
			FER=(double)(ferror)/(double)(j);

			printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E\r", progress, SNR, error, BER, ferror, FER);

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

		printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E\n", progress, SNR, error, BER, ferror, FER);

		OpenFile;
		fprintf(fp,"Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E\n", progress, SNR, error, BER, ferror, FER);
		fclose(fp);

	}

	getchar();
	
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
	int i, j, u, v, index, value, rxword[n], bi_rxword[n*p];
	unsigned long int mask;
	float sum, proba[pointNum], maxValue;
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
//		for(u=0;u<pointNum;u++)
		for(v=0;v<pointNum;v++)
		for(int x=0;x<pointNum;x++)
		for(int y=0;y<pointNum;y++)
		for(int h=0;h<pointNum;h++)
		{
			RM[j][i]=(Pr[0][h]*Pr[1][y]*Pr[2][x]*Pr[3][v]);  
			j++;
		}
	
	}

	//find the rxword and bi_rxword
	for(j=0;j<n;j++)
	{
		maxValue=RM[0][j];
		index=0;
		for(i=1;i<q;i++)
			if(maxValue<RM[i][j])
			{
				maxValue=RM[i][j];
				index=i;
			}
			rxword[j]=root[index];
	}

	//convert to binary 
	for(u=0;u<n*p;u++)	//n*4
		bi_rxword[u]=0;

	for(u=0;u<n;u++)	//n
	{
		value=rxword[u];
		mask=1;
		for(v=0;v<p;v++)
		{
			if((value & mask)>0)
				bi_rxword[p*u+v]=1;
			else
				bi_rxword[p*u+v]=0;
			mask=mask<<1;
		}
	}
	
	//*********KV optimal rule**************
	float score_1=0;
	float score_2=0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<q;j++)
		{
			if( codeword[i]==root[j] )
				score_1+=RM[j][i];
			
			score_2+=(RM[j][i]*RM[j][i]);
		}
	}

	score_2=weiz*score_2;
	score_2 = pow( score_2,0.5f );

	if( score_1>score_2 )
	{
		for(u=0;u<n;u++)
			dec_codeword[u]=codeword[u];
		
		for(u=0;u<n*p;u++)
			dec_bicodeword[u]=bi_codeword[u];
	}
	else
	{
		for(u=0;u<n;u++)
			dec_codeword[u]=rxword[u];
		for(u=0;u<n*p;u++)
			dec_bicodeword[u]=bi_rxword[u];
	}
	//*************************************

/*
	for(j=0;j<n;j++)
	{
		maxValue=RM[0][j];
		index=0;
		for(i=1;i<q;i++)
			if(maxValue<RM[i][j])
			{
				maxValue=RM[i][j];
				index=i;
			}
			rxword[j]=root[index];
	}

	epcount1=0;
	for(i=0;i<n;i++)	//n
		if(rxword[i]!=codeword[i])
			epcount1++;
*/	
//	printf("\n\nepcount1=%d\n", epcount1);



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
		for(j=0;j<monoTable_Ysize;j++)
		{
			printf("\n");
			for(z=0;z<monoTable_Xsize;z++)
				printf("\t%d ", mono_order[i][j][z]);
		}
	}
*/	//******************

}

/*
void polebasis()
{
	int v,i,j;

	//Initialisation
	for(v=0;v<pbNum;v++)	//num of pb
		for(i=0;i<pb_Ysize;i++)	//y_size
			for(j=0;j<pb_Xsize;j++)	//x_size
				pb[v][i][j]=0;

	for(v=0;v<pbNum;v++)	//num of pb
		pb[v][tg_order[v][1]][tg_order[v][0]]=1;
}
*/