#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define q 4
#define n 8
#define k 4
#define p 2
#define w 2
#define weiz 4
#define iterNum 8	//when m=1, C is equal to n
#define lm 1
#define able_correct 1
#define pointNum 2
#define eta 7
#define test_vec_num 128 //the num of test_vec is 2^eta
#define interval 1

//test_vec_construction
#define choose_num 2

//tgorder()
#define tg_size 100	//tg_size represent the probably used pole basis num
//mono_table()
#define weight_Zsize 4
#define weight_XYsize 100
#define mono_ordinSize 100
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 8	//large than the interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize 22 //>index of term[max(deg_y)][w] in the pole basis + 1
//polebasis()
#define pbNum 80	//num of pb
#define pb_Ysize 20	//large than max(degY) + 1
#define pb_Xsize (w+1)	//w+1
//interpolation()
#define init_polyNum 4	//the poly num of the init polyGroup
#define interpoly_Zsize (lm+1)	//maxValue of z is lm=1, so the Zsize should add 1 more.
#define interpoly_Ysize 6	//maxdeg of y is (w-1) + w*(n/(w+1)), so the Ysize should add 1 more.
#define interpoly_Xsize	(w+1)	//maxdeg of x is w, so the Xsize should add 1 more.
//factorization()
#define facpoly_Zsize (lm+1)
#define facpoly_Ysize 6	
#define facpoly_Xsize (w+1)	
//rcs()
#define rcspoly_Ysize 10	//degY+1, degY = interpoly_Ysize + max(j_1) + w*( (max(i_1)+w)/(w+1) )
#define rcspoly_Xsize 5	//degX+1, degX = max(i_1) + w
#define faiMax_Ysize 1	//the max degY of probably used polebasis
#define faiMax_Xsize 2	//the max degX of probably used polebasis
//expoly()-->expanded polynomial
#define expoly_Ysize (faiMax_Ysize+1)	
#define expoly_Xsize (faiMax_Xsize+1)


//main()
unsigned long int seq_num;	//number of input binary sequences
float SNR,N0;
double BER, FER;
int mularray[]={1,2,3};
int root[]={0, 1, 2, 3};	//elements in GF(4)
int bi_message[k*p], message[k];	//transmitted messge[k]
int codeword[n], bi_codeword[n*p]; //codewords
float tx_symbol[p*n][2], rx_symbol[p*n][2], sgm;
float RM[q][n];
//polebasis()
int pb[pbNum][pb_Ysize][pb_Xsize];	//pole basis[num_pb][y_size][x_size]
//findpoint()
int point[2][n];	//rational points[2][w^3]
//tgorder()
int tg_order_1[tg_size][2],tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
//mono_table()
int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
int gmatrix[k][n];	//generator matrix[k][n]
int dec_message[k],dec_codeword[n],dec_bicodeword[n*p]; //decoding result
//test_vec_construction
int large_vec[choose_num][n];
float test_set_ordered[2][n];
//nt test_vec_com[n];
int test_vec_com[n-eta];
int x_ordered[2][n];
//interpolation
int bi_rxword[p*n], rxword[n]; //received wordss
int epcount1, epcount2;
int Q_interpoly[test_vec_num][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];	//chosen interpolated polynomial [z-deg+1][max(deg_y)+1][w+1]
int uu;	//factorisation step index
int l, listNum[test_vec_num];	//candidate output index
//factorization() and rcs()
int Q[k][facpoly_Zsize][facpoly_Ysize][facpoly_Xsize];	//sequtial deduction polynomial [number of fac steps=k][rs][y_size][w+1], y_size> maxdeg_y]+rs*(deg_¦µ(k-1))
int rootlist[k][lm+1];	//the list of roots [number of fac steps=k][expected number of roots in each solution, 5>rs]
int output[lm+1][k], outputList[test_vec_num][lm+1][k];	//the list of candidate message [expeced number of candidate message, >rs][length of message, k]
int expoly[2][expoly_Ysize][expoly_Xsize];	//expanded polynomial in [z+f_k-1-u*pb_k-1-u]^rs, expoly[rs][3>(max(deg_y) in encoding functions)*(rs-1)][3>(max(deg_x) in encoding functions)*(rs-1)]
float pi=3.141593;	// pai

//int part_1[25][25],part_2[25][25],part_com[25][25];	//max(deg_y)+1
//int root_order[2][3]={{0,0,1},{0,1,0}}; //multiplicity m=2



int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
int permutation(int, int);
void findpoints(void);
void mono_table(void);
//void polyexp(int, int, int);
void polyexp1(int, int, int, int);
//void zerobasis(void);
void polebasis(void);
void tgorder(void);
//void coefficient(void);
void test_vec_contruction(void);
void interpolation(void);
void com_elem_interpolation(int g[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int interpoint[][n-eta]);
void uncom_elem_interpolation(int interpoint[4],int inGroup[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int outGroup1[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int outGroup2[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize]);
void factorisation(void);
void rcs(int);
void choose(void);
void generator(void);
void encoder(int message_temp[], int codeword_temp[]);
void modulation(void);
void channel(void);
void demodulation(void);
float PDF(int,int);

void main()
{
	int i, u, m, s, num, value, counter;	//hd[5>rs]
	float start, finish;
	unsigned long int j;
	unsigned int mask=1;
	long int error,ferror;
	double progress, b;

	FILE *fp;
	if((fp=fopen("LCC_Herm(8,4)eta=7.txt","a"))==NULL)
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
		printf("fai_%d=x^%d*y^%d\t",i,tg_order[i][0],tg_order[i][1]);
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

	
	polebasis();

	srand(1977);

	
	printf("\nPlease enter start SNR: ");
	scanf("%f", &start);
	printf("\nPlease enter finish SNR: ");
	scanf("%f", &finish);
	printf("\nseq_num: ");
	scanf("%d", &seq_num);

	for(SNR=start; SNR<=finish; SNR=SNR+interval)
	{
		N0=(1.0/((float)k/(float)n))/pow(10.0, SNR/10.0);
		sgm=sqrt(N0/2);
		
		error=0;
		ferror=0;

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

			//LIST DECODER
			//test vector construction
			test_vec_contruction();
			//interpolation
			interpolation();

			//factorisation
			factorisation();

			//choose
			choose();

/*
			//Convert the decoded codeword into binary
			for(u=0;u<n*p;u++)	//n*4
				dec_bicodeword[u]=0;
			
			for(u=0;u<n;u++)	//n
			{
				value=dec_codeword[u];
				mask=1;
				for(m=0;m<p;m++)
				{
					if((value & mask)>0)
						dec_bicodeword[p*u+m]=1;
					else
						dec_bicodeword[p*u+m]=0;
					mask=mask<<1;
				}
			}
*/
			//bit error rate calculation
			for(u=0;u<n*p;u++)	//n*4
			{
				if(dec_bicodeword[u]!=bi_codeword[u])
					error++;
			}

			//frame error rate calculation
			int temp=0;
			int ferror_temp=-1;
			for(u=0;u<n;u++)
				if(codeword[u]==dec_codeword[u])
					temp++;

			if(temp==n)
				ferror_temp=0;
			else if(temp<n)
				ferror_temp=1;

			if(ferror_temp==1)
				ferror++;

			progress=(double)(j*100)/(double)seq_num;
			BER=(double)(error)/(double)(n*p*j);
			FER=(double)(ferror)/(double)(j);

			printf("Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E\r", progress, SNR, error, BER, ferror, FER);

			if(ferror>200)
				break;
		}

		if(ferror>200)
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

		fp=fopen("LCC_Herm(8,4)eta=7.txt","a");
		fprintf(fp,"Progress=%0.1f, SNR=%2.2f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E\n", progress, SNR, error, BER, ferror, FER);
		fclose(fp);

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
	int i,j,u,v,index;
	float sum,proba[pointNum],maxValue;
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
		for(u=0;u<pointNum;u++)
		for(v=0;v<pointNum;v++)
//		for(int x=0;x<pointNum;x++)
//		for(int y=0;y<pointNum;y++)
//		for(int h=0;h<pointNum;h++)
		{
			RM[j][i]=(Pr[1][u]*Pr[0][v]);  
			j++;
		}
	
	}
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

/*
void demodulation()
{
	int i;

	for(i=0;i<n*p;i++)
	{
		if(rx_symbol[i][0]>=0)
			bi_rxword[i]=0;
		else
			bi_rxword[i]=1;
	}
	
	for(i=0;i<n;i++)	
		rxword[i]=1*bi_rxword[p*i]+2*bi_rxword[p*i+1];
	

	epcount1=0;
	for(i=0;i<n;i++)	//n
		if(rxword[i]!=codeword[i])
			epcount1++;
	
//	printf("\n\nepcount1=%d\n", epcount1);
}
*/

void test_vec_contruction()
{
	//construct the test_vec[test_vec_num]
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

	
	epcount1=0;
	for(i=0;i<n;i++)	//n
		if(codeword[i]!=large_vec[0][i])
			epcount1++;

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
	for(i=0;i<n-eta;i++)
		test_vec_com[i]=large_vec[0][(int)test_set_ordered[1][i]];

	//construct x_ordered
	for(i=0;i<n;i++)
		for(j=0;j<2;j++)
			x_ordered[j][i]=point[j][(int)test_set_ordered[1][i]];

/*	//************debug***************
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
	for(i=0;i<n;i++)
		printf("\t(%d,%d)",x_ordered[0][i],x_ordered[1][i]);

	
	printf("\n\n");
	printf("large_vec_order[]");
	for(j=0;j<n;j++)
		printf("\t%d",large_vec[0][(int)test_set_ordered[1][j]]);
	
	printf("\n\n");
	printf("large_vec2_order[]");
	for(j=0;j<n;j++)
		printf("\t%d",large_vec[1][(int)test_set_ordered[1][j]]);
	
*/	//********************************		



}

void interpolation()
{
	int i,j,u,v,z,num,temp,temp1,temp_index;
	int com_elem_interpoint[3][(n-eta)];
	int Q_com_elem[init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];
	int uncom_elem_interpoint[eta][4];
	int Q_uncom_elem[test_vec_num][init_polyNum][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];
	int degree_temp[test_vec_num][init_polyNum];

	//common element interpolation

	//set common element interpoint (xi,ri)
	for(i=0;i<n-eta;i++)
	{
		com_elem_interpoint[0][i]=x_ordered[0][i];
		com_elem_interpoint[1][i]=x_ordered[1][i];
		com_elem_interpoint[2][i]=test_vec_com[i];
	}

/*	//*********debug***********
	printf("\n\ncom_elem_interpoint\n");
	for(i=0;i<n-eta;i++)
		printf("\t%d",com_elem_interpoint[i][0]);
	printf("\n");
	for(i=0;i<n-eta;i++)
		printf("\t%d",com_elem_interpoint[i][1]);
*/	//***************************

	for(i=0;i<init_polyNum;i++)	//num_polys
		for(j=0;j<interpoly_Zsize;j++)	//z-deg+1
			for(u=0;u<interpoly_Ysize;u++)	//max(deg_y)+1
				for(v=0;v<interpoly_Xsize;v++)	//w+1
					Q_com_elem[i][j][u][v]=0;

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
		uncom_elem_interpoint[i-n+eta][0]=x_ordered[0][i];
		uncom_elem_interpoint[i-n+eta][1]=x_ordered[1][i];
		uncom_elem_interpoint[i-n+eta][2]=large_vec[0][(int)test_set_ordered[1][i]];
		uncom_elem_interpoint[i-n+eta][3]=large_vec[1][(int)test_set_ordered[1][i]];
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
	for(z=0;z<test_vec_num;z++)
		for(i=0;i<init_polyNum;i++)
			for(j=0;j<interpoly_Zsize;j++)
				for(u=0;u<interpoly_Ysize;u++)
					for(v=0;v<interpoly_Xsize;v++)
						Q_uncom_elem[z][i][j][u][v]=0;

	for(i=0;i<init_polyNum;i++)
		for(j=0;j<interpoly_Zsize;j++)
			for(u=0;u<interpoly_Ysize;u++)
				for(v=0;v<interpoly_Xsize;v++)
					Q_uncom_elem[0][i][j][u][v]=Q_com_elem[i][j][u][v];	


	//interpolation
	for(i=0;i<eta;i++)
	{
		num=(int)pow(2.0,i);
		for(u=num-1;u>=0;u--)
			uncom_elem_interpolation(uncom_elem_interpoint[i],Q_uncom_elem[u],Q_uncom_elem[2*u+0],Q_uncom_elem[2*u+1]);
	}		


	//choose the poly to factorization

	//initialize
	for(i=0;i<test_vec_num;i++)
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					Q_interpoly[i][u][v][z]=0;
	//eta>0
	for(i=0;i<test_vec_num;i++)
	{
		//calculate the degree of poly
		for(j=0;j<init_polyNum;j++)
		{
			temp=-1;
			degree_temp[i][j]=0;
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						if( Q_uncom_elem[i][j][u][v][z]!=0 )
							if( temp < mono_order[u][v][z] )
								temp=mono_order[u][v][z];

			degree_temp[i][j]=temp;
		}

		//choose the min degree
		temp_index=0;
		temp=degree_temp[i][0];
		for(j=1;j<init_polyNum;j++)
			if( temp>degree_temp[i][j] )
			{
				temp=degree_temp[i][j];
				temp_index=j;
			}

		//assignment
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					Q_interpoly[i][u][v][z]=Q_uncom_elem[i][temp_index][u][v][z];

	}

/*	//eta=0
	for(i=0;i<test_vec_num;i++)
	{
		//calculate the degree of poly
		for(j=0;j<init_polyNum;j++)
		{
			temp=-1;
			degree_temp[i][j]=0;
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						if( Q_com_elem[j][u][v][z]!=0 )
							if( temp < mono_order[u][v][z] )
								temp=mono_order[u][v][z];

			degree_temp[i][j]=temp;
		}

		//choose the min degree
		temp_index=0;
		temp=degree_temp[i][0];
		for(j=1;j<init_polyNum;j++)
			if( temp>degree_temp[i][j] )
			{
				temp=degree_temp[i][j];
				temp_index=j;
			}

		//assignment
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					Q_interpoly[i][u][v][z]=Q_com_elem[temp_index][u][v][z];
	}
*/

/*	
	//********debug*************
	int temp_x, temp_y, temp_z;
	for(i=0;i<test_vec_num;i++)
	{
		for(j=0;j<n-eta;j++)
		{
			temp=0;
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						if(Q_interpoly[i][u][v][z]!=0)
						{
							temp_x=power( com_elem_interpoint[0][j],z );
							temp_y=power( com_elem_interpoint[1][j],v );
							temp_z=power( com_elem_interpoint[2][j],u );
							temp= add( temp, mul( Q_interpoly[i][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
						}
			printf("\nQ_interpoly[%d] in point[%d](%d,%d,%d) = %d",i,j,com_elem_interpoint[0][j],com_elem_interpoint[1][j],com_elem_interpoint[2][j],temp);
		}

		for(j=0;j<eta;j++)
		{	
			temp=0;
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						if(Q_interpoly[i][u][v][z]!=0)
						{
							temp_x=power( uncom_elem_interpoint[j][0],z );
							temp_y=power( uncom_elem_interpoint[j][1],v );
							temp_z=power( uncom_elem_interpoint[j][i+2],u );
							temp= add( temp, mul( Q_interpoly[i][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
						}
			printf("\nQ_interpoly[%d] in point[%d](%d,%d,%d) = %d",i,j,uncom_elem_interpoint[j][0],uncom_elem_interpoint[j][1],uncom_elem_interpoint[j][i+2],temp);
		}
		
	}
	//**************************
*/


}


void com_elem_interpolation(int g[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int interpoint[][n-eta])
{
	int i, j, u, v, z, delta[init_polyNum], delta_temp, J[init_polyNum], act[init_polyNum], lod_temp, lod[init_polyNum], lod_min, j_min;	//(delta, J, act, lod)[num of polys]
	int f[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];  //g[num_polys][z-deg+1][max(deg_y)+1][w+1], g1[z-deg+1][max(deg_y)+1][w+2], (g2, f)[z-deg+1][max(deg_y)+1][w+1]	
	int g1[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize+1], g2[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];
	int poly_index;

	//Initialisation
	for(i=0;i<init_polyNum;i++)	//num_polys
		for(j=0;j<interpoly_Zsize;j++)	//z-deg+1
			for(u=0;u<interpoly_Ysize;u++)	//max(deg_y)+1
				for(v=0;v<interpoly_Xsize;v++)	//w+1
					g[i][j][u][v]=0;
	
	for(i=0;i<(lm+1);i++)	//rs
		for(j=0;j<w;j++)	//w
			g[w*i+j][i][j][0]=1;	//j+w*i
	

	//Interpolation
	for(i=0;i<n-eta;i++)	//wrt each point
	{
		//Calculate each polynomial's leading order
		for(j=0;j<init_polyNum;j++)	//num_poly
		{
			lod_temp=0;
			lod[j]=0;
			for(u=0;u<interpoly_Zsize;u++)	//rs
			{
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				{
					for(z=0;z<interpoly_Xsize;z++)	//w+1
					{
						if(g[j][u][v][z]!=0)
						{
							lod_temp=mono_order[u][v][z];
							if(lod_temp>lod[j])
								lod[j]=lod_temp;
						}
					}
				}
			}
		}
			
			
		//Initialise the eliminator array act[num_poly]
		for(j=0;j<init_polyNum;j++)	//num_poly
		{
			if(lod[j]<=iterNum)	//C=n when multiplicity = 1
				act[j]=1;
			else
				act[j]=0;
		}
			
			
		//Calculate the hasse derivative mapping of each polynomials
		j_min=-1;
		for(j=0;j<init_polyNum;j++)	//wrt each polynomial
		{
			J[j]=0;
			delta[j]=0;
			if(act[j]==1)	//for polynomials with leading order less of equal to C
			{
				//Hasse derivative mapping
				for(u=0;u<interpoly_Zsize;u++)	//rs
				{
					for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					{
						for(z=0;z<interpoly_Xsize;z++)	//w+1
						{
							if(g[j][u][v][z]!=0)
							{
								delta_temp=mul( mul( power(interpoint[0][i],z),power(interpoint[1][i],v) ), power(interpoint[2][i], u) );
								delta_temp=mul(delta_temp, g[j][u][v][z]);
								//Hasse derivative mapping
								delta[j]=add(delta[j],delta_temp);
							}
						}
					}
				}
				if(delta[j]!=0)
				{
					J[j]=1;	//record those polynomial with a nonzero hasse derivative mapping
					lod_min=lod[j];
					j_min=j;
				}
			}
		}
			
		//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
		for(j=0;j<init_polyNum;j++)	//num_polys
		{
			if(J[j]==1 && lod[j]<lod_min)
			{
				lod_min=lod[j];
				j_min=j;
			}
		}
			//printf("\nj_min=%d\n", j_min);

		if(j_min!=-1)
		{
			//f=g[j_min]
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						f[u][v][z]=g[j_min][u][v][z];
	
			//Modify nonzero polynomials
			for(j=0;j<init_polyNum;j++)	//num of polys
			{
				if(J[j]==1)
				{
					if(j!=j_min)
					{
						//delta*g_k+delta_k*f
						for(u=0;u<interpoly_Zsize;u++)	//rs
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
								for(z=0;z<interpoly_Xsize;z++)	//w+1
									g[j][u][v][z]=add(mul(delta[j_min],g[j][u][v][z]),mul(delta[j],f[u][v][z]));	
					}
					else if(j==j_min)
					{
						for(u=0;u<interpoly_Zsize;u++)	//rs
						{
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							{
								for(z=0;z<(interpoly_Xsize+1);z++)	//w+2
									g1[u][v][z]=0;
								for(z=0;z<interpoly_Xsize;z++)	//w+1
									g2[u][v][z]=0;
							}
						}
						
						//g1=x*f
						for(u=0;u<interpoly_Zsize;u++)	//rs
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
								for(z=0;z<interpoly_Xsize;z++)	//w+1
									g1[u][v][z+1]=g[j][u][v][z];
						//convert x^3=y^2+y, difference with diff code
						for(u=0;u<interpoly_Zsize;u++)	//rs
						{
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							{
								if(g1[u][v][w+1]!=0)
								{
									g1[u][v+1][0]=add(g1[u][v+1][0],g1[u][v][w+1]);
									g1[u][v+w][0]=add(g1[u][v+w][0],g1[u][v][w+1]);
									g1[u][v][w+1]=0;
								}
							}
						}
						//g2=xi*f
						for(u=0;u<interpoly_Zsize;u++)	//rs
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
								for(z=0;z<interpoly_Xsize;z++)	//w+1
									g2[u][v][z]=mul(interpoint[0][i],g[j][u][v][z]);
						//g=g1+g2
						for(u=0;u<interpoly_Zsize;u++)	//rs
							for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
								for(z=0;z<interpoly_Xsize;z++)	//w+1
									g[j][u][v][z]=add(g1[u][v][z],g2[u][v][z]);
					}
				}
			}
		}
	}
/*
	//****debug**********
	int temp_x, temp_y, temp_z,temp;
	for(i=0;i<init_polyNum;i++)
	{
		for(j=0;j<n-eta;j++)
		{
			temp=0;
			for(u=0;u<interpoly_Zsize;u++)	//rs
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
						if(g[i][u][v][z]!=0)
						{
							temp_x=power( interpoint[0][j],z );
							temp_y=power( interpoint[1][j],v );
							temp_z=power( interpoint[2][j],u );
							temp= add( temp, mul( g[i][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
						}
			printf("\ng[%d] with lod[%d]=%d in point[%d](%d,%d,%d) = %d",i,i,lod[i],j,interpoint[0][j],interpoint[1][j],interpoint[2][j],temp);
		}
		printf("\n");
	}
*/	//*************************
}

void uncom_elem_interpolation(int interpoint[4],int inGroup[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int outGroup1[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize],int outGroup2[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize])
{

	int i, j, u, v, z, J[2][init_polyNum], act[init_polyNum], lod_temp, lod[init_polyNum], lod_min[2], j_min[2];	//(delta, J, act, lod)[num of polys]
	int f[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];  //g[num_polys][z-deg+1][max(deg_y)+1][w+1], g1[z-deg+1][max(deg_y)+1][w+2], (g2, f)[z-deg+1][max(deg_y)+1][w+1]	
	int g1[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize+1], g2[interpoly_Zsize][interpoly_Ysize][interpoly_Xsize];
	int c[init_polyNum][lm+1],result_temp[init_polyNum][2];

	int poly_index;

	for(j=0;j<init_polyNum;j++)
		for(i=0;i<2;i++)
			result_temp[j][i]=0;

	//Interpolation

	//Calculate each polynomial's leading order
	for(j=0;j<init_polyNum;j++)	//num_poly
	{
		lod_temp=0;
		lod[j]=0;
		for(u=0;u<interpoly_Zsize;u++)	//rs
		{
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
			{
				for(z=0;z<interpoly_Xsize;z++)	//w+1
				{
					if(inGroup[j][u][v][z]!=0)
					{
						lod_temp=mono_order[u][v][z];
						if(lod_temp>lod[j])
							lod[j]=lod_temp;
					}
				}
			}
		}
	}
		
		
	//Initialise the eliminator array act[num_poly]
	for(j=0;j<init_polyNum;j++)	//num_poly
	{
		if(lod[j]<=iterNum)	//C=n when multiplicity = 1
			act[j]=1;
		else
			act[j]=0;
	}
		
		
	//Calculate the hasse derivative mapping of each polynomials
	for(j=0;j<init_polyNum;j++)
		if(act[j]==1)
			for(u=0;u<interpoly_Zsize;u++)
			{
				c[j][u]=0;
	
				for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
					for(z=0;z<interpoly_Xsize;z++)	//w+1
					if(inGroup[j][u][v][z]!=0)
					{
						c[j][u] = add( c[j][u],mul( inGroup[j][u][v][z],mul( power(interpoint[0],z),power(interpoint[1],v) ) ) );
					}
			}
	
	for(j=0;j<init_polyNum;j++)
		if(act[j]==1)
			for(i=0;i<2;i++)
			{
				result_temp[j][i]=0;
		
				for(u=0;u<interpoly_Zsize;u++)
				{
					result_temp[j][i] = add( result_temp[j][i],mul( c[j][u],power(interpoint[i+2],u) ) );
				}						
			}
	
	for(i=0;i<2;i++)
	{
		j_min[i]=-1;
		for(j=0;j<init_polyNum;j++)
		{
			J[i][j]=0;
			if( act[j]==1 && result_temp[j][i]!=0 )
			{
				J[i][j]=1;	//record those polynomial with a nonzero hasse derivative mapping
				lod_min[i]=lod[j];
				j_min[i]=j;
			}
		}
	}
	//Identify the minimal polynomial with a nonzero Hasse derivative evaluation
	for(i=0;i<2;i++)
		for(j=0;j<init_polyNum;j++)	//num_polys
		{
			if(J[i][j]==1 && lod[j]<lod_min[i])
			{
				lod_min[i]=lod[j];
				j_min[i]=j;
			}
		}
	//printf("\nj_min=%d\n", j_min);

	//initialize outGroup1 and outGroup2
	for(j=0;j<init_polyNum;j++)
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1	
					outGroup1[j][u][v][z]=inGroup[j][u][v][z];

	for(j=0;j<init_polyNum;j++)
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1	
					outGroup2[j][u][v][z]=inGroup[j][u][v][z];

	//update the poly of outGroup1 
	if(j_min[0]!=-1)
	{
		//f=g[j_min]
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					f[u][v][z]=inGroup[j_min[0]][u][v][z];

		//Modify nonzero polynomials
		for(j=0;j<init_polyNum;j++)	//num of polys
		{
			if(J[0][j]==1)
			{
				if(j!=j_min[0])
				{
					//delta*g_k+delta_k*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								outGroup1[j][u][v][z]=add( mul( result_temp[j_min[0]][0],inGroup[j][u][v][z] ),mul( result_temp[j][0],f[u][v][z] ) );	
				}
				else if(j==j_min[0])
				{
					for(u=0;u<interpoly_Zsize;u++)	//rs
					{
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
						{
							for(z=0;z<(interpoly_Xsize+1);z++)	//w+2
								g1[u][v][z]=0;
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								g2[u][v][z]=0;
						}
					}
					
					//g1=x*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								g1[u][v][z+1]=inGroup[j][u][v][z];
					//convert x^3=y^2+y, difference with diff code
					for(u=0;u<interpoly_Zsize;u++)	//rs
					{
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
						{
							if(g1[u][v][w+1]!=0)
							{
								g1[u][v+1][0]=add(g1[u][v+1][0],g1[u][v][w+1]);
								g1[u][v+w][0]=add(g1[u][v+w][0],g1[u][v][w+1]);
								g1[u][v][w+1]=0;
							}
						}
					}
					//g2=xi*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								g2[u][v][z]=mul( interpoint[0],inGroup[j][u][v][z]);
					//g=g1+g2
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								outGroup1[j][u][v][z]=add(g1[u][v][z],g2[u][v][z]);
				}
			}
		}
	}

	//update the poly of outGroup2 
	if(j_min[1]!=-1)
	{
		//f=g[j_min]
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					f[u][v][z]=inGroup[j_min[1]][u][v][z];

		//Modify nonzero polynomials
		for(j=0;j<init_polyNum;j++)	//num of polys
		{
			if(J[1][j]==1)
			{
				if(j!=j_min[1])
				{
					//delta*g_k+delta_k*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								outGroup2[j][u][v][z]=add( mul( result_temp[j_min[1]][1],inGroup[j][u][v][z] ),mul( result_temp[j][1],f[u][v][z] ) );	
				}
				else if(j==j_min[1])
				{
					for(u=0;u<interpoly_Zsize;u++)	//rs
					{
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
						{
							for(z=0;z<(interpoly_Xsize+1);z++)	//w+2
								g1[u][v][z]=0;
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								g2[u][v][z]=0;
						}
					}
					
					//g1=x*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								g1[u][v][z+1]=inGroup[j][u][v][z];
					//convert x^3=y^2+y, difference with diff code
					for(u=0;u<interpoly_Zsize;u++)	//rs
					{
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
						{
							if(g1[u][v][w+1]!=0)
							{
								g1[u][v+1][0]=add(g1[u][v+1][0],g1[u][v][w+1]);
								g1[u][v+w][0]=add(g1[u][v+w][0],g1[u][v][w+1]);
								g1[u][v][w+1]=0;
							}
						}
					}
					//g2=xi*f
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								g2[u][v][z]=mul( interpoint[0],inGroup[j][u][v][z]);
					//g=g1+g2
					for(u=0;u<interpoly_Zsize;u++)	//rs
						for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
							for(z=0;z<interpoly_Xsize;z++)	//w+1
								outGroup2[j][u][v][z]=add(g1[u][v][z],g2[u][v][z]);
				}
			}
		}
	}
	
/*	//*******debug*********
	int temp_x,temp_y,temp_z,temp;

	
	for(j=0;j<init_polyNum;j++)
	{
		temp=0;
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					if(outGroup1[j][u][v][z]!=0)
					{
						temp_x=power( interpoint[0],z );
						temp_y=power( interpoint[1],v );
						temp_z=power( interpoint[2],u );
						temp= add( temp, mul( outGroup1[j][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
					}
		printf("\noutGroup1[%d] with lod[%d]=%d in point[7](%d,%d,%d) = %d",j,j,lod[j],interpoint[0],interpoint[1],interpoint[2],temp);
	}
	printf("\n");
	
	for(j=0;j<init_polyNum;j++)
	{
		temp=0;
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					if(outGroup2[j][u][v][z]!=0)
					{
						temp_x=power( interpoint[0],z );
						temp_y=power( interpoint[1],v );
						temp_z=power( interpoint[3],u );
						temp= add( temp, mul( outGroup2[j][u][v][z],mul( temp_z,mul( temp_y,temp_x ) ) ) );	
					}
		printf("\noutGroup2[%d] with lod[%d]=%d in point[7](%d,%d,%d) = %d",j,j,lod[j],interpoint[0],interpoint[1],interpoint[3],temp);
	}
	printf("\n");
*/	//******************

}

void factorisation(void)
{
	int i, j, u, v, z;
	
	for(i=0;i<test_vec_num;i++)
	{
		//Initialisation
		for(j=0;j<k;j++)	//number of fac steps=k
			for(u=0;u<facpoly_Zsize;u++)	//rs
				for(v=0;v<facpoly_Ysize;v++)	//y_size
					for(z=0;z<facpoly_Xsize;z++)	//w+1
						Q[j][u][v][z]=0;

		for(u=0;u<k;u++)	//number of fac steps=k
			for(v=0;v<lm+1;v++)	//5>rs=expected number of roots
				rootlist[u][v]=-1;	

		for(u=0;u<lm+1;u++)	//5>(rs-1)=expected number of output lists
			for(v=0;v<k;v++)	//k
				output[u][v]=-1;

		//Initialisation of factorisation
		uu=0;	//recursive deduction index
		l=0;	//candidate output index
		//q_0(z)=itp(z)
		for(u=0;u<interpoly_Zsize;u++)	//rs
			for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
				for(z=0;z<interpoly_Xsize;z++)	//w+1
					Q[uu][u][v][z]=Q_interpoly[i][u][v][z];

		//recursive coefficient search
		rcs(uu);

		//store the necessary data
		for(u=0;u<lm+1;u++)	//5>(rs-1)=expected number of output lists
			for(v=0;v<k;v++)	//k
				outputList[i][u][v]=output[u][v];

		listNum[i]=l;

	}

}

void rcs(int uu)
{
	int i, j, u, v, m, z, t, r, i_1, j_1, i_2, j_2, a, b, d, leadMono, leadMono_temp, alpha, act, lc[lm+1], q_temp[lm+1][rcspoly_Ysize][rcspoly_Xsize];	//q_temp[z-deg+1][y_size>max(deg_y)+1+(rs-1)*(max(deg_y) in encoding functions)][14>w+(rs-1)*w], lc[rs]--leading coefficient polynomial
	int d1,d2,flag;
	int rcspoly_Ysize_1,rcspoly_Xsize_1;	//the q_temp size of the first step of factorization
	int pbY,pbX;

	//array size initialization
	rcspoly_Ysize_1= faiMax_Ysize + interpoly_Ysize + 1;
	rcspoly_Xsize_1= faiMax_Xsize + w + 1;
	pbY = rcspoly_Ysize+1;
	pbX = w+1;

	//printf("\nWhen u=%d\n", u);
	leadMono=0; leadMono_temp=0;	//leading monomial index
	act=0;	//judge value for recursive search of each f_k-1-u
	//find pb_k-1-uu
	j_1=-1;
	i_1=-1;
	for(j=0;j<pbY;j++)	//y_size, efficiently enough
	{	
		for(i=0;i<pbX;i++)	//w+1
		{
			if(pb[k-1-uu][j][i]==1)	//k-1-uu
			{
				j_1=j;	//record pb_k-1-u's y degree
				i_1=i;	//record pb_k-1-u's x degree
				break;
			}
		}
		if((j_1!=-1) && (i_1!=-1))
			break;
	}

	//printf("pb_k-1-u");

	//initialise q_temp
	for(i=0;i<lm+1;i++)	//rs
		for(j=0;j<rcspoly_Ysize;j++)	//y_size
			for(u=0;u<rcspoly_Xsize;u++)	//5>w*(rs-1)+w
				q_temp[i][j][u]=0;


	//Calculate q_temp[pb_k-1-u]=q[u][pb_k-1-u]
	for(i=0;i<facpoly_Zsize;i++)	//rs
	{
		for(j=0;j<facpoly_Ysize;j++)	//y_size=max(deg_y) in encoding function*(rs-1), and 31>=max(deg_y)+1
		{
			for(u=0;u<facpoly_Xsize;u++)	//
			{
				q_temp[i][j+i*j_1][u+i*i_1]=Q[uu][i][j][u];	//this time, q_temp size is [rs+1][maxY(j_1)+maxY(Q)][maxX(i_1)+maxX(Q)]
			}
		}
	}


	//convert x^3=y^2+y, difference with diff code
	d1=0;
	d2=0;
	for(i=0;i<lm+1;i++)	//rs
		{	
			for(j=0;j<rcspoly_Ysize_1;j++)	//y_size=degY+1
			{
				for(u=(w+1);u<rcspoly_Xsize_1;u++)	//w+1, 5>w*(rs-1)+w
				{
					if(q_temp[i][j][u]!=0)
					{
						d1=u/(w+1);	//d1 in (8,4) is less than 2;
						d2=u%(w+1);	//u%(w+1)

						//should add a funciton to calculate the poly (z+hixiyi)^j

						q_temp[i][j+d1*w][d2]=add(q_temp[i][j+d1*w][d2], q_temp[i][j][u]);	//this time, q_temp size is [rs+1][1+last_maxY(q_temp)+w*(last_maxX(q_temp)/(w+1))][w+1]
						q_temp[i][j+d1][d2]	 =add(q_temp[i][j+d1][d2], q_temp[i][j][u]);
					}
					q_temp[i][j][u]=0;
				}
			}
		}


	//printf("q[u][pb_k-1-u]");

	//find the leading monomial of q_temp[pb_k-1-u]
	for(i=0;i<lm+1;i++)	//rs
	{
		for(j=0;j<rcspoly_Ysize;j++)	//y_size
		{
			for(u=0;u<w+1;u++)	//w+1
			{
				if(q_temp[i][j][u]!=0)
				{
					for(v=0;v<pbNum;v++)	//number of pb, *take a note*--pole basis vialation might happen, as 95 pole basis functions can not promise to term[y_size-1][w]
					{
						for(m=0;m<rcspoly_Ysize+1;m++)	//y_size+1
						{	
							for(z=0;z<w+1;z++)	//w+1
							{
								if(pb[v][m][z]==1 && j==m && u==z)
								{
									leadMono_temp=v;
									if(leadMono_temp>=leadMono)
									{
										leadMono=leadMono_temp;	//record the leading monomial's index
										i_2=u;	//record the leading monomial's x degree
										j_2=j;	//record the leading monomial's y degree
									}	
								}
							}
						}
					}
				}
			}
		}
	}

	//printf("The leading monomial in q[u][pb_k-1-u]");

/*	//find the leading coefficient polynomial lc[rs]
	for(i=0;i<lm+1;i++)	//rs
	{
		flag=1;
		lc[i]=-1;
		for(j=0;j<rcspoly_Ysize;j++)	//y_size
		{
			for(u=0;u<w+1;u++)	//w+1
			{
				if(q_temp[i][j][u]!=0)
				{
					flag=0;
					if(j==j_2 && u==i_2)
					{
						lc[i]=q_temp[i][j][u];
						break;
					}
				}
			}
			if(lc[i]!=-1)
				break;
		}
		if(flag==1)
			lc[i]=0;
	}

	//if q_temp[i][j][u]==0, then setting lc[i]==0, and flag==2
	flag=0;
	for(i=0;i<lm+1;i++)
		if(lc[i]==0)
			flag++;

	//printf("leading coefficient polynomial");
	
	u=0;	//root index
	//find the roots of leading coefficient polynomial lc[rs]
	for(i=0;i<q;i++)	//number of elements in GF(4)
		if(flag!=(lm+1))	//rs
		{
			b=0;
			for(j=0;j<(lm+1);j++)	//rs
			{
				a=mul(lc[j], power(root[i], j));
				b=add(a, b);
			}
			if(b==0)
			{
				act=1;
				rootlist[uu][u]=root[i];
				u++;
			}
		}
		else if(flag==(lm+1))
			{
				rootlist[uu][0]=0;
				act=1;
			}
*/	//printf("roots");

	//find the leading coefficient polynomial lc[rs]
	for(i=0;i<lm+1;i++)	//rs
	{
		lc[i]=0;
		for(j=0;j<rcspoly_Ysize;j++)	//y_size
		{
			for(u=0;u<w+1;u++)	//w+1
			{
				if(q_temp[i][j][u]!=0 && j==j_2 && u==i_2)
					lc[i]=q_temp[i][j][u];
			}
		}
	}

	//printf("leading coefficient polynomial");
	
	u=0;	//root index
	//find the roots of leading coefficient polynomial lc[rs]
	for(i=0;i<q;i++)	//number of elements in GF(4)
	{
		b=0;
		for(j=0;j<lm+1;j++)	//rs
		{
			a=mul(lc[j], power(root[i], j));
			b=add(a, b);
		}
		if(b==0)
		{
			act=1;
			rootlist[uu][u]=root[i];
			u++;
		}
	}

	//For each distinct root of rootlist[u]
	if(act==1)
	{
	for(i=0;i<lm+1;i++)	//2>rs
	{
		if(rootlist[uu][i]!=-1)
		{	
			alpha=rootlist[uu][i];
			
			output[l][k-1-uu]=alpha;	//output[l][k-1-uu]
				
			//printf("current outputs");

			if(uu==(k-1))	//when u=k-1, output candidate polynomial
			{
				for(j=0;j<k;j++)	//k
				{
					if(output[l][j]==-1)
						output[l][j]=output[l-1][j];
				}

				l++;	//locate next candidate output
			}
			else  //update the q[uu+1]
			{
				//Initialise q_temp
				for(j=0;j<lm+1;j++)	//rs
					for(u=0;u<rcspoly_Ysize;u++)	//y_size
						for(m=0;m<rcspoly_Xsize;m++)	//3>w*(rs-1)+w
							q_temp[j][u][m]=0;	
				
				//q_temp=q[uu][z+f_k-1-u*pb_k-1-u]
				for(j=0;j<lm+1;j++)	//rs
				{
					//calculate (z+f_k-1-u*pb_k-1-u)^j
					if(j==0)
					{
						for(u=0;u<lm+1;u++)
							for(m=0;m<expoly_Ysize;m++)
								for(z=0;z<expoly_Xsize;z++)
									expoly[u][m][z]=0;
						expoly[0][0][0]=1;
					}
					else if(j>0)
							polyexp1(alpha, i_1, j_1, j);
					
					//calculate q[uu][z+f_k-1-u*pb_k-1-u]
					for(u=0;u<lm+1;u++)	//rs
					{
						for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
						{
							for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
							{
								if(expoly[u][m][z]!=0)
								{
									for(v=0;v<facpoly_Ysize;v++)	//y_size
									{
										for(t=0;t<facpoly_Xsize;t++)	//w+1
										{
											if(Q[uu][j][v][t]!=0)
												q_temp[u][v+m][t+z]=add(q_temp[u][v+m][t+z], mul(expoly[u][m][z], Q[uu][j][v][t]));
											// caution overflow problem of q_temp
											// here q_temp, ysize=exploy_Ysize+facpoly_Ysize-2+1=9, xsize=expoly_Xsize+facpoly_Xsize-2+1=6
										}
									}
								}
							}
						}
					}
				}
				
				//convert x^3=y^2+y, difference with diff code
				d1=0;
				d2=0;
				for(j=0;j<lm+1;j++)	//rs
				{
					for(u=0;u<rcspoly_Ysize;u++)	//10 y_size
					{
						for(m=w+1;m<rcspoly_Xsize;m++)	//w+1, 14>w*(rs-1)+w
						{
							if(q_temp[j][u][m]!=0)
							{
								d1=m/(w+1);	//d1 in (8,4) is less than 2;
								d2=m%(w+1);	//u%(w+1)
								q_temp[j][u+d1*w][d2]=add(q_temp[j][u+d1*w][d2], q_temp[j][u][m]);
								q_temp[j][u+d1][d2]	 =add(q_temp[j][u+d1][d2], q_temp[j][u][m]);
							}
							q_temp[j][u][m]=0;
						}
					}
				}
			
				r=uu+1;

				//q[u+1]=q_temp
				for(j=0;j<lm+1;j++)	//rs
					for(u=0;u<facpoly_Ysize;u++)	//y_size
						for(m=0;m<facpoly_Xsize;m++)	//w+1
							Q[r][j][u][m]=q_temp[j][u][m];

				//printf("q[u+1]");

				//next coefficient searching
				
				rcs(r);
			}
		}
	}
	}
	
}

void choose()
{
	int i, j, u, m, v, min_num, value, flag, min_index;
	unsigned int mask=1;
	float hd[(lm+1)*test_vec_num+1],temp_x,temp_y,dis_temp;
	int bi_codeword_temp[n*p],codeword_temp[n];
	int tx_symbol_temp[n*p][2];

	//Initialise hamming distance counter
	for(u=0;u<(lm+1)*test_vec_num;u++)
		hd[u]=-1.0;

	for(i=0;i<test_vec_num;i++)
		if( listNum[i]==0 )
		{
			for(j=0;j<lm+1;j++)
				hd[i*(lm+1)+j]=-1.0;
		}
		else if( listNum[i]!=0 )
		{
			for(j=0;j<listNum[i];j++)
			{
				//reencoding
				encoder(outputList[i][j],codeword_temp);
				//nonbinary --> binary
				for(u=0;u<n;u++)
				{	
					value=codeword_temp[u];
					mask=1;
					for(v=0;v<p;v++)
					{
						if((value & mask)>0)
							bi_codeword_temp[p*u+v]=1;
						else
							bi_codeword_temp[p*u+v]=0;
						mask=mask<<1;
					}
				}
				//transfer to the symbol_temp
				for(u=0;u<n*p;u++)
				{
					tx_symbol_temp[u][0]=-2*bi_codeword_temp[u]+1;//0-->(1,0)
					tx_symbol_temp[u][1]=0;						  //1-->(-1,0)
				}
				//calculate the Euclidean distance
				dis_temp=0.0;
				for(u=0;u<n*p;u++)
				{
					temp_x=tx_symbol_temp[u][0]-rx_symbol[u][0];
					temp_y=tx_symbol_temp[u][1]-rx_symbol[u][1];
					dis_temp+=sqrt( (temp_x*temp_x) + (temp_y*temp_y) );
				}
				hd[i*(lm+1)+j]=dis_temp;
			}
		}


		//choose the min
		dis_temp=1000.0;
		min_index=-1;
		flag=0;
		for(i=0;i<(lm+1)*test_vec_num;i++)
			if(hd[i]!=-1.00)
			{
				flag=1;	//exist at less one valid solution

				if(dis_temp>hd[i])
				{
					dis_temp=hd[i];
					min_index=i;
				}
			}

		//encode message[min_index]
		encoder(outputList[min_index/(lm+1)][min_index%(lm+1)],codeword_temp);

		//nonbinary --> binary
			for(u=0;u<n;u++)
			{	
				value=codeword_temp[u];
				mask=1;
				for(v=0;v<p;v++)
				{
					if((value & mask)>0)
						bi_codeword_temp[p*u+v]=1;
					else
						bi_codeword_temp[p*u+v]=0;
					mask=mask<<1;
				}
			}


		//output the correct 
		if(flag==0)	// not exist a valid solution
		{
			for(u=0;u<n;u++)
				dec_codeword[u]=large_vec[0][u];
			
			//nonbinary --> binary
			for(u=0;u<n;u++)
			{	
				value=dec_codeword[u];
				mask=1;
				for(v=0;v<p;v++)
				{
					if((value & mask)>0)
						dec_bicodeword[p*u+v]=1;
					else
						dec_bicodeword[p*u+v]=0;
					mask=mask<<1;
				}
			}			
		}
		else if(flag==1)	//exist a valid solution
		{
			for(u=0;u<n;u++)
				dec_codeword[u]=codeword_temp[u];
			for(u=0;u<n*p;u++)
				dec_bicodeword[u]=bi_codeword_temp[u];
		}



			//Check Herm(64, 39)'s working perperty
			epcount2=0;
			for(u=0;u<n;u++)
				if(codeword[u]!=dec_codeword[u])
					epcount2++;
			//printf("\nepcount2=%d\n", epcount2);
/*
			if( (epcount1<=able_correct) && (epcount2!=0) )
			{
//				printf("\nFrame %d errors!", j);
				printf("\nTx message is:\n");
				for(u=0;u<k;u++)	//k
					printf("%d ", message[u]);
				printf("\nTx code word is:\n");
				for(u=0;u<n;u++)	//n
					printf("%d ", codeword[u]);

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


				printf("large_vec[]\t");
				for(j=0;j<n;j++)
					printf("\t%d",large_vec[0][j]);
			
				printf("\n\n");
				printf("large_vec2[]\t");
				for(j=0;j<n;j++)
					printf("\t%d",large_vec[1][j]);
			
				printf("\n\n");
				printf("x_order[]\t");
				for(i=0;i<n;i++)
					printf("\t(%d,%d)",x_ordered[0][i],x_ordered[1][i]);
			
				
				printf("\n\n");
				printf("large_vec_order[]");
				for(j=0;j<n;j++)
					printf("\t%d",large_vec[0][(int)test_set_ordered[1][j]]);
				
				printf("\n\n");
				printf("large_vec2_order[]");
				for(j=0;j<n;j++)
					printf("\t%d",large_vec[1][(int)test_set_ordered[1][j]]);

				printf("\n%d errors in received!", epcount1);
				printf("\noutputlist:");
				
				for(i=0;i<test_vec_num;i++)
				{
					printf("\n\nlistNum[%d]=%d",i,listNum[i]);
					for(m=0;m<listNum[i];m++)
					{
						printf("\n[%d][%d]",i,m);
						for(u=0;u<k;u++)
							printf("\t%d",outputList[i][m][u]);
					}
				}
//				printf("\nDecoded message is:\n");
//				for(u=0;u<k;u++)	//n
//					printf("%d ", dec_message[u]);
				printf("\nDecoded code word is:\n");
				for(u=0;u<n;u++)	//n
					printf("%d ", dec_codeword[u]);
				printf("\n%d errors after decoding!\n\n", epcount2);
			}
*/		
}


void tgorder()
{
	int i,j;

	for(i=0;i<tg_size;i++)
	{
		tg_order_1[i][0]=-1;
		tg_order_1[i][1]=-1;
		tg_order[i][0]=-1;
		tg_order[i][1]=-1;
	}

	//total graduate order
	tg_order_1[0][0]=0;
	tg_order_1[0][1]=0;

	for(i=1;i<tg_size*3;i++)
	{
		if(tg_order_1[i-1][0]==0)
		{
			tg_order_1[i][0]=tg_order_1[i-1][1]+1;
			tg_order_1[i][1]=0;
		}
		else
		{
			tg_order_1[i][0]=tg_order_1[i-1][0]-1;
			tg_order_1[i][1]=tg_order_1[i-1][1]+1;
		}
	}

	//eliminate terms containing x^3
	j=0;
	for(i=0;i<tg_size;i++)
	{
		if(tg_order_1[i][0]<w+1)	//w+1
		{
			tg_order[j][0]=tg_order_1[i][0];
			tg_order[j][1]=tg_order_1[i][1];
			j++;
		}
	}
	/*
	for(i=0;i<100;i++)
		printf("\n(%d, %d)", tg_order[i][0], tg_order[i][1]);
	*/
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

/*	for(i=0;i<weight_Zsize;i++)	
	{
		if(i==0)
			for(j=0;j<weight_XYsize;j++)
				weight[i][j]=tg_order[j][0]*w + tg_order[j][1]*(w+1);
		else if(i>0)
			for(j=0;j<weight_XYsize;j++)
				weight[i][j]= weiz*i + weight[0][j];
	}
*/
	for(i=0;i<weight_Zsize;i++)	
	{
		weight[i][0]=weiz*i;
		for(j=1;j<weight_XYsize;j++)
			weight[i][j]=weiz*i+j+1;
	}


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

/*	
	printf("\nMonmial basis is:\n");
	for(i=0;i<3;i++)
	{
		for(j=0;j<25;j++)
		{
			for(z=0;z<5;z++)
				printf("%d ", mono_order[i][j][z]);
			printf("\n");
		}
		printf("\n");
	}
*/	
}

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

void polyexp1(int c, int i, int j, int deg_z)
{
	int u, m, z;	//p1[rs+1][18>(max(deg_y) in encoding functions)*(rs-1)][10>(max(deg_x) in encoding functions)*(rs-1)], p2[rs][26=18+max(deg_y) in encoding functions][14=10+max(deg_x) in encoding functions]

	//Initialisations
	for(u=0;u<lm+1;u++)	//rs
		for(m=0;m<expoly_Ysize;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
			for(z=0;z<expoly_Xsize;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
				expoly[u][m][z]=0;

	if(deg_z==0)
		expoly[0][0][0]=1;
	else if(deg_z>0)	//because deg_z is equal to or less than 1, so deg_z>0 <-> deg_z==1
	{
		expoly[1][0][0]=1;
		expoly[0][j][i]=c;

	}

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

int inv(int fac)
{
	int i,invresult;
	if(fac==0)
		invresult=0;
	else
		for(i=0;i<q-1;i++){	//size of mularray
			if(fac==mularray[i])
			invresult=mularray[((q-1)-i)%(q-1)];	//size of mularray
	}

	return invresult;
}

int permutation(int a, int b)
{
	int i;
	double fac1, fac2, perm_result;

	fac1=1;
	for(i=a;i>(a-b);i--)
		fac1=fac1*(float)i;

	fac2=1;
	for(i=b;i>0;i--)
		fac2=fac2*(float)i;

	perm_result=fac1/fac2;

	return perm_result;
}
