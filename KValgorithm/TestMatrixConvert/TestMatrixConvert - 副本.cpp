//**********************************************
//作用：将reliability matrix 转化成 Multiplicity matrix, 并输出iterNum
//输入：
//		1. 
//		2.
//输出：Multiplicity matrix MM[q][n]
//***********************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "FiniteFieldBasisGF(16).h"
#include "FiniteFieldBasisCal.h"

#define myWay
#define OpenFile fp=fopen("Herm(64,49)_KV_l=5.txt","a")
//#define printfMonotable
//#define checkMulpMatrix



#define w 4
#define weiz 54
#define lm 5
#define pointNum 2
#define interval 1


//*****need to be modify_second**************
//tgorder()
#define tg_size 815	//more than k, tg_size represent the probably used pole basis num, (w+1)*w/2 + (w+1) * ( (w+interpoly_Ysize)-w+1 ) 
//mono_table()
#define weight_Zsize 35
#define weight_XYsize 815	//equals to the tg_size
#define mono_ordinSize 855	//choose the value with debugging
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 166	//large than interpoly_Ysize, but we should avoid to the array crossing the border
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize 815 //>index of term[max(deg_y)][w] in the pole basis + 1, nearly equals to tg_size
//*******************

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
double aver_iterNum;
//findpoint()
int point[n][2];	//rational points[2][w^3], [0]->x, [1]->y
//tgorder()
int tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
//mono_table()
int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
int gmatrix[k][n];	//generator matrix[k][n]
int dec_message[k],dec_codeword[n],dec_bicodeword[n*p]; //decoding result
//MatrixConvert()
int s=0;	//large enough to make the circle of MatrixConvert() keep moving

//********degbug***************
	//some var about computation complexity
	unsigned long int addNum, mulNum, addNum_temp, mulNum_temp;
	int flag_addNum=0, flag_mulNum=0;

	//some var about Max_m in table
	int Max_m_intable;
	int Max_m;
//*****************************

void findpoints(void);
void mono_table(void);
void tgorder(void);
void generator(void);
void encoder(int message_temp[], int codeword_temp[]);
void modulation(void);
void channel(void);
void demodulation(void);
float PDF(int,int);
void MatrixConvert(void);

void main(void)
{
	int i, u, m, num, value;	
	float start, finish;
	unsigned long int j,v;
	unsigned int mask=1;
	double progress;

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
		printf("affine point[%d]=(%d,%d)\n",i,point[i][0],point[i][1]);
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
		
		aver_iterNum=0.0;
		Max_m=0;

		for(j=1;j<=seq_num;j++)
		{

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

			progress=(double)(j*100)/(double)seq_num;
			aver_iterNum = aver_iterNum + (iterNum-aver_iterNum)/(double)j;
			
			if(Max_m < Max_m_intable)
			{
				Max_m = Max_m_intable;
			}

			
			printf("Progress=%0.1f, SNR=%2.2f, aver_iterNum=%0.2f, Max_m=%d\r", progress, SNR, aver_iterNum, Max_m);

		}

		printf("Progress=%0.1f, SNR=%2.2f, aver_iterNum=%0.2f, Max_m=%d\n", progress, SNR, aver_iterNum, Max_m);

		OpenFile;
		fprintf(fp,"Progress=%0.1f, SNR=%2.2f, aver_iterNum=%0.2f, Max_m=%d\n", progress, SNR, aver_iterNum, Max_m);
		fclose(fp);

	}

}


void MatrixConvert()
{
	int i, j, u, v, index_i, index_j, lM, flag_iterNum;
	unsigned long int degree_iterNum, iterNum_temp;
	float RM_temp[q][n], tempRM, max_temp;
	int MulpMatrix[q][n], tempMM;

	//initialization
	for(i=0;i<n;i++)
		for(j=0;j<q;j++)
		{
			RM_temp[j][i] = RM[j][i];
			MulpMatrix[j][i] = 0; 
			MM[j][i] = 0;
		}

	s = INT_MAX;	//important part, must be check
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

		//update
		RM_temp[index_j][index_i] = RM_temp[index_j][index_i] / (MulpMatrix[index_j][index_i]+2);
		MulpMatrix[index_j][index_i] += 1;
		s -= 1; 

#ifdef checkMulpMatrix	
	fp=fopen("TestResult.txt","a");
	fprintf(fp, "\nindex_j = %d\tindex_i = %d\ts=%d", index_j, index_i, s);

	fprintf(fp, "\n\nMM:\n");
	for(j=0;j<q;j++)
	{
		for(i=0;i<n;i++)
		{
			fprintf(fp, "%d\t",MulpMatrix[j][i]);
		}
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
					if ( mono_order[u][j][i] == iterNum_temp )	//here refering to the size of mono_table 
					{
						degree_iterNum = i*w + j*(w+1) + u*weiz;
						flag_iterNum = 0;
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

	iterNum = iterNum_temp;	//CM

#ifdef checkFac
	//cal degree_iterNum
	Deg_iterNum = 10000;
	Deg_iterNum = degree_iterNum;

	//cal codewordScore
	codewordScore = 0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<q;j++)
			if( codeword[i]==root[j] )
				{
					codewordScore += MulpMatrix[j][i];
				}
	}
#endif

	Max_m_intable=0;

	for(j=0;j<q;j++)
		for(i=0;i<n;i++)
		{
			if(Max_m_intable<MulpMatrix[j][i])
			{
				Max_m_intable=MulpMatrix[j][i];
			}
		}


	
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




	//weight of monomial over function x^(w+1)-y^w-y=0
	for(i=0;i<weight_Zsize;i++)	
	{
		if(i==0)
			for(j=0;j<weight_XYsize;j++)
				weight[i][j]=tg_order[j][0]*w + tg_order[j][1]*(w+1);
		else if(i>0)
			for(j=0;j<weight_XYsize;j++)
				weight[i][j]= weiz*i + weight[0][j];
	}


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

void tgorder()
{
	int i, j, u, index_temp;

	for(j=0;j<tg_size;j++)
		for(i=0;i<2;i++)
			tg_order[j][i]=-1;

	//judge the index's scale of coresponding tgsize
	/*pesudo code (index = deg_x + deg_y)
	if(  0 <= index < w ) tgsize = (index+2) * (index+1) / 2
	else if( w <= index ) tgsize = w*(w+1)/2 + (w+1)*(index-w+1)
	*/
	index_temp = 164;	//according to the last formualtion to calculate out

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

void findpoints()
{
	int i, j, u, a1, a2, a3, x, y;

	//Initialisation
	for(j=0;j<n;j++)	//w^3
		for(i=0;i<2;i++)
			point[j][i]=0;

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
				point[u][0]=x;
				point[u][1]=y;
				u++;
			}
		}
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
			gmatrix[i][j]=mul(power(point[j][0], tg_order[i][0]), power(point[j][1], tg_order[i][1]));

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


