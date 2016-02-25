#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FiniteFieldSize 16
#define mode 1
#define CheatingDecoding
#define myWay
#define OpenFile fp=fopen("[cheating]Herm(64,49)_l=1.txt","a")
#define FrameError 310

//#define checkDecoding

//#define printCoeffTable
//#define printDemodulation
//#define printfMonotable
//#define printMatrixConvert
//****global define***********
#define k 49
#define weiz 54
#define lm 1	//design length
#define pointNum 2
#define interval 0.5

#if FiniteFieldSize == 4
#define q 4
#define n 8
#define p 2
#define w 2

#elif FiniteFieldSize == 16
#define q 16
#define n 64
#define p 4
#define w 4

#elif FiniteFieldSize == 64
#define q 64
#define n 512
#define p 6
#define w 8

#endif

#define genius (w*(w-1)/2)
#define able_correct ((n-k-genius)/2)

//tgorder()
#define tg_size 297	//more than k, tg_size represent the probably used pole basis num, (w+1)*w/2 + (w+1) * ( (w+interpoly_Ysize)-w+1 ) 
//mono_table()
#define weight_Zsize 75
#define weight_XYsize 297	//equals to the tg_size
#define mono_ordinSize 700	//choose the value with debugging
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 100	//large than interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize 297 //>index of term[max(deg_y)][w] in the pole basis + 1, nearly equals to tg_size
//polebasis()


//finitfieldSize
#if FiniteFieldSize == 4
int mularray[] = {1,2,3};  //this array is used to the infinite field element mul
int root[] = {0, 1, 2, 3};	//n+1		be used to the factorization
int logNum[] = {-1, 0, 1, 2};	//used to locate the degree of finitefield element through the value of finitefield element

#elif FiniteFieldSize == 16
int mularray[] = {1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};  //this array is used to the infinite field element mul
int root[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};	//n+1		be used to the factorization
int logNum[] = {-1,0,1,4,2,8,5,10,3,14,9,7,6,13,11,12};	//used to locate the degree of finitefield element through the value of finitefield element

#elif FiniteFieldSize == 64
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

//main()
unsigned long int seq_num;	//number of input binary sequences
float SNR,N0;
double BER, FER;
float pi=3.141593;	// pai
int iterNum=0; //iteration num
int bi_message[k*p], message[k];	//transmitted messge[k]
int codeword[n], bi_codeword[n*p]; //codewords
float tx_symbol[p*n][2], rx_symbol[p*n][2], sgm;
float RM[q][n];	//Reliability Matrix
int MM[q][n];	//Multiplicity Matrix

//findpoint()
int point[2][n];	//rational points[2][w^3], [0]->x, [1]->y
//tgorder()
int tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
//mono_table()
int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
int gmatrix[k][n];	//generator matrix[k][n]
int dec_message[k],dec_codeword[n],dec_bicodeword[n*p]; //decoding result
//MatrixConvert()


	
	//****debug*********
#ifdef CheatingDecoding
	int Deg_iterNum, codewordScore;
#endif


	//some var about computation complexity
	unsigned long int addNum, mulNum, addNum_temp, mulNum_temp;
	int flag_addNum=0, flag_mulNum=0;

	//*******************


//**********function**********
int mul(int, int);
int add(int, int);
int power(int, int);
void findpoints(void);
void mono_table(void);
void tgorder(void);
void generator(void);
void encoder(int message_temp[], int codeword_temp[]);
void modulation(void);
void channel(void);
void demodulation(void);
float PDF(int,int);
void decoding(void);
//*****************
void MatrixConvert(void);
//****************



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

		addNum_count=0.0;
		mulNum_count=0.0;
		totalNum_count=0.0;

		flag_addNum=1;
		flag_mulNum=1;

		for(j=1;j<=seq_num;j++)
		{
			addNum=0;
			mulNum=0;

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

			//MatrixConvert
			iterNum = 0;
			MatrixConvert();

			//LIST DECODER
			//test vector construction
			decoding();

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
						
			printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\r", progress, SNR, error, BER, ferror, FER,	addNum_count, mulNum_count, totalNum_count);

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

		printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);

		OpenFile;
		fprintf(fp,"Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E, addNum=%0.2f, mulNum=%0.2f, total_num=%0.2f\n", progress, SNR, error, BER, ferror, FER, addNum_count, mulNum_count, totalNum_count);
		fclose(fp);

	}

	flag_addNum=0;
	flag_mulNum=0;
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

#if FiniteFieldSize == 4
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
#elif FiniteFieldSize == 16
		for(int y=0;y<pointNum;y++)
		for(int h=0;h<pointNum;h++)
		for(u=0;u<pointNum;u++)
		for(v=0;v<pointNum;v++)
		{
			RM[j][i]=(Pr[0][v]*Pr[1][u]*Pr[2][h]*Pr[3][y]);  
			j++;
		}
#elif FiniteFieldSize == 64
		for(int z=0;z<pointNum;z++)
		for(int x=0;x<pointNum;x++)
		for(int y=0;y<pointNum;y++)
		for(int h=0;h<pointNum;h++)
		for(u=0;u<pointNum;u++)
		for(v=0;v<pointNum;v++)
		{
			RM[j][i]=(Pr[0][v]*Pr[1][u]*Pr[2][h]*Pr[3][y]*Pr[4][x]*Pr[5][z]);  
			j++;
		}
#endif
	}

#ifdef printDemodulation
	//*****debug**********
	printf("\n\ncodeword");
	for(v=0;v<n;v++)
		printf("\t%d",codeword[v]);

/*	printf("\n\nrxword\t");
	for(v=0;v<n;v++)
		printf("\t%d",rxword[v]);
	
	printf("\n\n");
	for(v=0;v<n*p/4;v++)
		printf("\t(%f,%f,%f,%f)",rx_symbol[v*2+0][0],rx_symbol[v*2+0][1],rx_symbol[v*2+1][0],rx_symbol[v*2+1][0]);
*/
	printf("\n\n");
	for(u=0;u<q;u++)
	{
		printf("\n");
		printf("\t%d",root[u]);
		for(v=0;v<n;v++)
			printf("\t%f",RM[u][v]);
	}
	printf("\n\n");
	//*******************
#endif

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

void MatrixConvert(void)
{
	int i, j, u, v, index_i, index_j, lM, flag_iterNum;
	unsigned long int degree_iterNum, iterNum_temp;
	float RM_temp[q][n], tempRM, max_temp;
	int MulpMatrix[q][n], tempMM;
	int s=0;	//large enough to make the circle of MatrixConvert() keep moving

#ifdef printMatrixConvert
	FILE *fp;

	fp=fopen("TestResult.txt","a");
	fprintf(fp, "\n\nRM_temp:\n");
	for(j=0;j<q;j++)
	{
		for(i=0;i<n;i++)
			fprintf(fp, "%f\t",RM[j][i]);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
	fclose(fp);
#endif

	//initialization
	for(i=0;i<n;i++)
		for(j=0;j<q;j++)
		{
			RM_temp[j][i] = RM[j][i];
			MulpMatrix[j][i] = 0; 
			MM[j][i] = 0;	//Further, this global var 
		}

	s = 3000;	//inportant part, must be check
	lM = 0;

	//start
	while( s>0 && lM<=lm )
	{
		max_temp = 0.0;
		index_j = -1;
		index_i = -1;


		//find the maximal entry factor in RM
		for(i=0;i<n;i++)
			for(j=0;j<q;j++)
			{
				if( max_temp < RM_temp[j][i] )
				{
//					sumNum1 += 4;
					max_temp = RM_temp[j][i];
					index_j = j;
					index_i = i;

#ifdef myWay
//					sumNum2 += 5;
					//because the sum of a cloumn of RM is equal to 1. so 
					//for a cloumn, if there is a factor larger than 0.5, that must
					//be the largest one of the cloumn
					if( max_temp > 0.5 )
					{
//						sumNum2 += 4;
						max_temp = RM_temp[j][i];
						index_j = j;
						index_i = i;
						break; //jump out the cloumn circle
					}
#endif
				}
			}
		if(index_j==-1 && index_i==-1 )
			printf("\n\nMatrixConvert() error\n\n");


		//update
		RM_temp[index_j][index_i] = RM_temp[index_j][index_i] / (MulpMatrix[index_j][index_i]+2);
		MulpMatrix[index_j][index_i] += 1;
		s -= 1; 

#ifdef printMatrixConvert
	fp=fopen("TestResult.txt","a");
	fprintf(fp, "\nindex_j = %d\tindex_i = %d\ts=%d", index_j, index_i, s);

	fprintf(fp, "\n\nMM:\n");
	for(j=0;j<q;j++)
	{
		for(i=0;i<n;i++)
			fprintf(fp, "%d\t",MulpMatrix[j][i]);
		fprintf(fp,"\n");
	}

	fprintf(fp, "\n\nRM_temp:\n");
	for(j=0;j<q;j++)
	{
		for(i=0;i<n;i++)
			fprintf(fp, "%f\t",RM_temp[j][i]);
		fprintf(fp,"\n");
	}

	fprintf(fp,"\n");
//	printf("\n\nsumNum1=%d\tsumNum2=%d", sumNum1, sumNum2);
	fclose(fp);
#endif
		//cal CM
		iterNum_temp = 0;
		for(j=0;j<q;j++)
			for(i=0;i<n;i++)
				if(MulpMatrix[j][i]!=0)
				{
					iterNum_temp +=( MulpMatrix[j][i] * (MulpMatrix[j][i]+1));
				}
		iterNum_temp = iterNum_temp/2;

		//cal the deg(1,wz) corresponding to CM
		degree_iterNum = 0;
		flag_iterNum = 1;	//make the search part more effciency
		for(u=0; u<(lm+1) && flag_iterNum ;u++)
			for(j=0; j<monoTable_Ysize && flag_iterNum ;j++)
				for(i=0; i<monoTable_Xsize && flag_iterNum ;i++)
					if ( mono_order[u][j][i] == iterNum_temp )
					{
						degree_iterNum = i*w + j*(w+1) + u*weiz;
						flag_iterNum = 0;
						break;
					}

		//cal lM
		lM = degree_iterNum/weiz;
		
	}

	//transform to the MM[q][n]
	for(j=0;j<q;j++)
		for(i=0;i<n;i++)
		{
			MM[j][i] = MulpMatrix[j][i];
		}

	iterNum = iterNum_temp;

#ifdef CheatingDecoding
	//cal degree_iterNum
	Deg_iterNum = degree_iterNum;
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

void tgorder()
{
	int i, j, u, index_temp;

	for(j=0;j<tg_size;j++)
		for(i=0;i<2;i++)
			tg_order[j][i]=-1;

	//judge the index's scale of coresponding tgsize
	/*pesudo code
	if(  0 <= index < w ) tgsize = (index+2) * (index+1) / 2
	else if( w <= index ) tgsize = w*(w+1)/2 + (w+1)*(index-w+1)
	*/
	index_temp = 99;	//according to the last formualtion to calculate out

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

void mono_table()
{
	int i, j, u, v, z, weight[weight_Zsize][weight_XYsize], mono_order_1[monoTable_Zsize][weight_XYsize];	//26=(100/weight(z))+1, 100=pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][max(deg_y)+1][w+1]

	//Initialisation
	for(i=0;i<monoTable_Zsize;i++)	//(230/weight(z))+1
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
		for(i=0;i<monoTable_Zsize;i++)	//230/weight(z)
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

#ifdef printfMonotable
	//**********debug********
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
	for(i=0;i<monoTable_Zsize;i++)
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

	//******************
#endif
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

void decoding(void)
{
	int i, j, value, epcount1, score_temp;
	float proba_temp, temp;
	unsigned int mask=1;
	int large_vec[n];

#ifdef printDemodulation
	printf("\n\n");
	for(i=0;i<q;i++)
	{
		printf("\n");
		printf("\t%d",root[i]);
		for(j=0;j<n;j++)
			printf("\t%d",MM[i][j]);
	}
	printf("\n\n");
#endif 

#if mode == 1	

	//cal codewordScore
	codewordScore = 0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<q;j++)
			if( codeword[i]==root[j] )
				{
					codewordScore += MM[j][i];
				}
	}

	//decoding
	if( codewordScore > Deg_iterNum )	//decdoing correctly
	{
		//codeword
		for(i=0;i<n;i++)
			dec_codeword[i] = codeword[i];

		//bi_codeword
		for(i=0;i<n*p;i++)
			dec_bicodeword[i] = bi_codeword[i];
	}
	else if(codewordScore <= Deg_iterNum)	//decoding incorrectly
	{
		//search dec_codeword[n] from RM
		for(i=0;i<n;i++)
		{
			proba_temp = 0.0;
			for(j=0;j<q;j++)
				if( proba_temp<RM[j][i] )
				{
					proba_temp = RM[j][i];
					dec_codeword[i] = root[j];
				}
		}

		//bi_codeword[n*p]
		//nonbinary --> binary
		for(i=0;i<n;i++)
		{	
			value=dec_codeword[i];
			mask=1;
			for(j=0;j<p;j++)
			{
				if((value & mask)>0)
					dec_bicodeword[p*i+j]=1;
				else
					dec_bicodeword[p*i+j]=0;
				mask=mask<<1;
			}
		}	

	}

#elif mode == 2	//using R Matrix, KV optimal
	proba_temp = 0.0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<q;j++)
			if( codeword[i]==root[j] )
				{
					proba_temp += RM[j][i];
				}
	}
	proba_temp = proba_temp * proba_temp;

	temp = 0.0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<q;j++)
			if(RM[j][i]!=0)
			{
				temp += (RM[j][i]*RM[j][i]);
			}
	}
	
	if( proba_temp > (temp*weiz) )
	{
		//codeword
		for(i=0;i<n;i++)
			dec_codeword[i] = codeword[i];

		//bi_codeword
		for(i=0;i<n*p;i++)
			dec_bicodeword[i] = bi_codeword[i];
	}
	else if(proba_temp <= (temp*weiz))
	{
		//search dec_codeword[n] from RM
		for(i=0;i<n;i++)
		{
			proba_temp = 0.0;
			for(j=0;j<q;j++)
				if( proba_temp<RM[j][i] )
				{
					proba_temp = RM[j][i];
					dec_codeword[i] = root[j];
				}
		}

		//bi_codeword
		//nonbinary --> binary
		for(i=0;i<n;i++)
		{	
			value=dec_codeword[i];
			mask=1;
			for(j=0;j<p;j++)
			{
				if((value & mask)>0)
					dec_bicodeword[p*i+j]=1;
				else
					dec_bicodeword[p*i+j]=0;
				mask=mask<<1;
			}
		}	
	}

#endif


#ifdef checkDecoding
	//search dec_codeword[n] from RM
	for(i=0;i<n;i++)
	{
		proba_temp = 0.0;
		for(j=0;j<q;j++)
			if( proba_temp<RM[j][i] )
			{
				proba_temp = RM[j][i];
				large_vec[i] = root[j];
			}
	}

	epcount1 = 0;
	for(i=0;i<n;i++)
		if( codeword[i]!=large_vec[i] )
			epcount1++;

	if( (epcount1<=able_correct) && (codewordScore<=Deg_iterNum) )	//error num can be correct
		{
			printf("\n\ndecdoing is error!!, epcount1=%d\tcodwordScore=%d\tDeg_iterNum=%d\n\n", epcount1, codewordScore, Deg_iterNum);

			printf("\n\ncodeword:\t");
			for(i=0;i<n;i++)
				printf("%d\t",codeword[i]);
			printf("\n\n");

			printf("\n\n");
			for(i=0;i<q;i++)
			{
				printf("\n");
				printf("\t%d",root[i]);
				for(j=0;j<n;j++)
					printf("\t%d",MM[i][j]);
			}
			printf("\n\n");

			printf("\n\n");
			for(i=0;i<q;i++)
			{
				printf("\n");
				printf("\t%d",root[i]);
				for(j=0;j<n;j++)
					printf("\t%f",RM[i][j]);
			}
			printf("\n\n");
		}
#endif

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
		temp = (logNum[a]*b) % (q-1);
		pow_result=mularray[temp];
		return pow_result;
	}

	if( a<0 )
	{	
		printf("\n\n power has error!!");
	}

}

//********************************