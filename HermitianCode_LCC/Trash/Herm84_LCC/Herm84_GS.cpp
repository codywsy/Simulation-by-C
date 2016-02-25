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
//#define lm 1
#define able_correct 1

#define interpoly_Ysize 8
#define interpoly_Xsize	3	//w+1
#define facpoly_Ysize 8
#define facpoly_Xsize 3	//w+1
#define rcspoly_Ysize 12
#define rcspoly_Xsize 7

unsigned long int seq_num;	//number of input binary sequences
float SNR;
double BER, FER;
int mularray[]={1,2,3};
int root[q]={0, 1, 2, 3};	//elements in GF(4)
int pb[80][20][w+1];	//pole basis[num_pb][y_size][x_size]
int point[2][n];	//rational points[2][w^3]
int tg_order_1[100][2],tg_order[100][2];	//[promise the tg order beyond (w, max(deg_y))][2]
int mono_order[2][8][w+1];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
int gmatrix[k][n];	//generator matrix[k][n]
int bi_message[k*p], message[k];	//transmitted messge[k]
int dec_message[k],dec_codeword[n],dec_bicodeword[n*p]; //decoding result
int codeword[n], bi_codeword[n*p]; //codewords
float tx_symbol[p*n][2], rx_symbol[p*n][2], sgm;	//modulated/received symbols
int bi_rxword[p*n], rxword[n]; //received wordss
int epcount1, epcount2;
int itp[2][8][w+1];	//chosen interpolated polynomial [z-deg+1][max(deg_y)+1][w+1]
int uu;	//factorisation step index
int l;	//candidate output index
int Q[k][2][11][5];	//sequtial deduction polynomial [number of fac steps=k][rs][y_size][w+1], y_size> maxdeg_y]+rs*(deg_¦µ(k-1))
int rootlist[k][2];	//the list of roots [number of fac steps=k][expected number of roots in each solution, 5>rs]
int output[2][k];	//the list of candidate message [expeced number of candidate message, >rs][length of message, k]
int expoly[2][4][4];	//expanded polynomial in [z+f_k-1-u*pb_k-1-u]^rs, expoly[rs][3>(max(deg_y) in encoding functions)*(rs-1)][3>(max(deg_x) in encoding functions)*(rs-1)]
float pi=3.141593;	// pai



int part_1[25][25],part_2[25][25],part_com[25][25];	//max(deg_y)+1


int root_order[2][3]={{0,0,1},{0,1,0}}; //multiplicity m=2



int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
int permutation(int, int);
void findpoints(void);
void mono_table(void);
//void polyexp(int, int, int);
void polyexp1(int, int, int, int);
void zerobasis(void);
void polebasis(void);
void tgorder(void);
void coefficient(void);
void interpolation(void);
void factorisation(void);
void rcs(int);
void generator(void);
void encoder(int message_temp[], int codeword_temp[]);
void modulation(void);
void channel(void);
void demodulation(void);

void main()
{
	int i, u, m, s, num, value, hd[5], counter, d_min, list_num;	//hd[5>rs]
	float N0, start, finish;
	unsigned long int j;
	unsigned int mask=1;
	long int error,ferror;
	double progress, b;

	
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

	for(SNR=start; SNR<=finish; SNR++)
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

			encoder(message,codeword);

			//convert to binary 
			for(u=0;u<n*p;u++)	//n*4
				bi_codeword[u]=0;

			for(u=0;u<n;u++)	//n
			{
				value=codeword[u];
				mask=1;
				for(m=p-1;m>=0;m--)
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
			//interpolation
			interpolation();

			//factorisation
			factorisation();

			//choose the decoded message from the output list
			//Initialise hamming distance counter
			for(u=0;u<2;u++)	//2>rs
				hd[2]=-1;
			
			if(l==0)	//if there is no output
			{
				for(u=0;u<n;u++)	//n
					dec_codeword[u]=rxword[u];
				
				for(u=0;u<k;u++)
					dec_message[u]=-1;

			}
			else	//if there is output
			{
				for(u=0;u<l;u++)
				{
					hd[u]=0;

					//re-encoded, and calculate the hamming distance from the received word
					int output_codeword[n];
					encoder(output[u],output_codeword);
					for(m=0;m<n;m++)	//n
						if(output_codeword[m]!=rxword[m])
							hd[u]++;
				}
				//determine the minimum distance and the respective polynomial
				d_min=hd[0];
				list_num=0;
				for(u=0;u<l;u++)
				{
					if(hd[u]<d_min && hd[u]!=-1)
					{
						d_min=hd[u];
						list_num=u;
					}
				}
				//transmitted message judgement priority
				for(u=0;u<l;u++)
				{
					if(hd[u]==d_min)
					{
						counter=0;
						for(m=0;m<k;m++)	//k
							if(output[u][m]==message[m])
								counter++;
						if(counter==k)	//k
							list_num=u;
					}
				}

				//choose the result from the candidate message list
				for(u=0;u<k;u++)
					dec_message[u]=output[list_num][u];
				
				encoder(dec_message,dec_codeword);

			}

			//Check Herm(64, 39)'s working perperty
			epcount2=0;
			for(u=0;u<n;u++)
				if(codeword[u]!=dec_codeword[u])
					epcount2++;
			//printf("\nepcount2=%d\n", epcount2);

			if( (epcount1<=able_correct) && (epcount2!=0) )
			{
//				printf("\nFrame %d errors!", j);
				printf("\nTx message is:\n");
				for(u=0;u<k;u++)	//k
					printf("%d ", message[u]);
				printf("\nTx code word is:\n");
				for(u=0;u<n;u++)	//n
					printf("%d ", codeword[u]);
				printf("\nRx code word is:\n");
				for(u=0;u<n;u++)	//n
					printf("%d ", rxword[u]);
				printf("\n%d errors in received!", epcount1);
				printf("\noutputlist:");
				for(m=0;m<2;m++)
				{
					printf("\n[%d]");
					for(u=0;u<k;u++)
						printf("\t%d",output[m][u]);
				}
				printf("\nDecoded message is:\n");
				for(u=0;u<k;u++)	//n
					printf("%d ", dec_message[u]);
				printf("\nDecoded code word is:\n");
				for(u=0;u<n;u++)	//n
					printf("%d ", dec_codeword[u]);
				printf("\n%d errors after decoding!\n\n", epcount2);
			}
			
			//Convert the decoded codeword into binary
			for(u=0;u<n*p;u++)	//n*4
				dec_bicodeword[u]=0;
			
			for(u=0;u<n;u++)	//n
			{
				value=dec_codeword[u];
				mask=1;
				for(m=p-1;m>=0;m--)
				{
					if((value & mask)>0)
						dec_bicodeword[p*u+m]=1;
					else
						dec_bicodeword[p*u+m]=0;
					mask=mask<<1;
				}
			}

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

			printf("Progress=%0.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E\r", progress, SNR, error, BER, ferror, FER);
/*
			if(progress==b)
			{
				BER=(double)error/(double)(n*j);	//n*4
				printf("Progress=%0.1f, SNR=%2.1d, Errors=%2.1d, BER=%E\r", progress, SNR, error, BER);
				b++;
			}
*/
		}

		BER=(double)error/(double)(n*p*(j-1));
		FER=(double)(ferror)/(double)(j-1);
		printf("Progress=%0.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E, Frame Errors=%2.1d, FER=%E\n", progress, SNR, error, BER, ferror, FER);
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
		rxword[i]=2*bi_rxword[p*i]+1*bi_rxword[p*i+1];
	

/*	//******debug********
	printf("\n\nmessage:\n");
	for(i=0;i<k;i++)
		printf("%d\t",message[i]);
	printf("\n\ncodeword:\n");
	for(i=0;i<n;i++)
		printf("%d\t",codeword[i]);
	printf("\n\nrxword:\n");
	for(i=0;i<n;i++)
		printf("%d\t",rxword[i]);
*/	//******************

	epcount1=0;
	for(i=0;i<n;i++)	//n
		if(rxword[i]!=codeword[i])
			epcount1++;
	
//	printf("\n\nepcount1=%d\n", epcount1);
}

void interpolation()
{
	int  i, j, u, v, z, g[4][2][8][w+1], poly, delta[4], delta_temp, J[4], act[4], lod_temp, lod[4], lod_min, j_min, g1[2][8][w+2], g2[2][8][w+1], f[2][8][w+1];	//(delta, J, act, lod)[num of polys], g[num_polys][z-deg+1][max(deg_y)+1][w+1], g1[z-deg+1][max(deg_y)+1][w+2], (g2, f)[z-deg+1][max(deg_y)+1][w+1]	
	
	
	//Initialisation
	for(i=0;i<4;i++)	//num_polys
		for(j=0;j<2;j++)	//z-deg+1
			for(u=0;u<8;u++)	//max(deg_y)+1
				for(v=0;v<w+1;v++)	//w+1
					g[i][j][u][v]=0;
	
	for(i=0;i<2;i++)	//rs
		for(j=0;j<w;j++)	//w
			g[w*i+j][i][j][0]=1;	//j+w*i
	

	//Interpolation
	for(i=0;i<n;i++)	//wrt each point
	{
		//Calculate each polynomial's leading order
		for(j=0;j<4;j++)	//num_poly
		{
			lod_temp=0;
			lod[j]=0;
			for(u=0;u<2;u++)	//rs
			{
				for(v=0;v<8;v++)	//max(deg_y)+1
				{
					for(z=0;z<w+1;z++)	//w+1
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
		for(j=0;j<4;j++)	//num_poly
		{
			if(lod[j]<=iterNum)	//C=n when multiplicity = 1
				act[j]=1;
			else
				act[j]=0;
		}
			
			
		//Calculate the hasse derivative mapping of each polynomials
		j_min=-1;
		for(j=0;j<4;j++)	//wrt each polynomial
		{
			J[j]=0;
			delta[j]=0;
			if(act[j]==1)	//for polynomials with leading order less of equal to C
			{
				//Hasse derivative mapping
				for(u=0;u<2;u++)	//rs
				{
					for(v=0;v<8;v++)	//max(deg_y)+1
					{
						for(z=0;z<w+1;z++)	//w+1
						{
							if(g[j][u][v][z]!=0)
							{
								delta_temp=mul(mul(power(point[0][i], z), power(point[1][i], v)), power(rxword[i], u));
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
		for(j=0;j<4;j++)	//num_polys
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
			for(u=0;u<2;u++)	//rs
				for(v=0;v<8;v++)	//max(deg_y)+1
					for(z=0;z<w+1;z++)	//w+1
						f[u][v][z]=g[j_min][u][v][z];
	
			//Modify nonzero polynomials
			for(j=0;j<4;j++)	//num of polys
			{
				if(J[j]==1)
				{
					if(j!=j_min)
					{
						//delta*g_k+delta_k*f
						for(u=0;u<2;u++)	//rs
							for(v=0;v<8;v++)	//max(deg_y)+1
								for(z=0;z<w+1;z++)	//w+1
									g[j][u][v][z]=add(mul(delta[j_min],g[j][u][v][z]),mul(delta[j],f[u][v][z]));	
					}
					else if(j==j_min)
					{
						for(u=0;u<2;u++)	//rs
						{
							for(v=0;v<8;v++)	//max(deg_y)+1
							{
								for(z=0;z<w+2;z++)	//w+2
									g1[u][v][z]=0;
								for(z=0;z<w+1;z++)	//w+1
									g2[u][v][z]=0;
							}
						}
						
						//g1=x*f
						for(u=0;u<2;u++)	//rs
							for(v=0;v<8;v++)	//max(deg_y)+1
								for(z=0;z<w+1;z++)	//w+1
									g1[u][v][z+1]=g[j][u][v][z];
						//convert x^5=y^4+y
						for(u=0;u<2;u++)	//rs
						{
							for(v=0;v<8;v++)	//max(deg_y)+1
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
						for(u=0;u<2;u++)	//rs
							for(v=0;v<8;v++)	//max(deg_y)+1
								for(z=0;z<w+1;z++)	//w+1
									g2[u][v][z]=mul(point[0][i],g[j][u][v][z]);
						//g=g1+g2
						for(u=0;u<2;u++)	//rs
							for(v=0;v<8;v++)	//max(deg_y)+1
								for(z=0;z<w+1;z++)	//w+1
									g[j][u][v][z]=add(g1[u][v][z],g2[u][v][z]);
					}
				}
			}
		}
	}

	//Find the minimal polynomial
	for(j=0;j<4;j++)	//num of polys
	{
		lod[j]=0;
		lod_temp=0;
		for(u=0;u<2;u++)	//rs
		{
			for(v=0;v<8;v++)	//max(deg_y)+1
			{
				for(z=0;z<w+1;z++)	//w+1
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
	lod_min=lod[0];
	poly=0;
	for(j=0;j<4;j++)	//num of polys
	{
		if(lod[j]<lod_min)
		{
			lod_min=lod[j];
			poly=j;
		}
	}
	for(u=0;u<2;u++)	//rs
		for(v=0;v<8;v++)	//max(deg_y)+1
			for(z=0;z<w+1;z++)	//w+1
					itp[u][v][z]=g[poly][u][v][z];

/*	//*********debug*******
	printf("inter_porly:\n");
	for(u=0;u<2;u++)	//rs
	{	
		for(v=0;v<8;v++)	//max(deg_y)+1
		{
			for(z=0;z<w+1;z++)	//w+1
				printf("%d\t",g[poly][u][v][z]);
			printf("\n");
		}
		printf("\n");
	}
*/	//*********************

}

void factorisation()
{
	int i, j, u, v, z;
	
	//Initialisation
	for(i=0;i<k;i++)	//number of fac steps=k
		for(u=0;u<2;u++)	//rs
			for(v=0;v<facpoly_Ysize;v++)	//y_size
				for(z=0;z<facpoly_Xsize;z++)	//w+1
					Q[i][u][v][z]=0;

	for(i=0;i<k;i++)	//number of fac steps=k
		for(j=0;j<2;j++)	//5>rs=expected number of roots
			rootlist[i][j]=-1;	

	for(i=0;i<2;i++)	//5>(rs-1)=expected number of output lists
		for(j=0;j<k;j++)	//k
			output[i][j]=-1;

	//Initialisation of factorisation
	uu=0;	//recursive deduction index
	l=0;	//candidate output index
	//q_0(z)=itp(z)
	for(u=0;u<2;u++)	//rs
		for(v=0;v<interpoly_Ysize;v++)	//max(deg_y)+1
			for(z=0;z<interpoly_Xsize;z++)	//w+1
				Q[uu][u][v][z]=itp[u][v][z];

	//recursive coefficient search
	rcs(uu);
}

void rcs(int uu)
{
	int i, j, u, v, m, z, t, r, i_1, j_1, i_2, j_2, a, b, d, lm, lm_temp, alpha, act, lc[3], q_temp[2][rcspoly_Ysize][rcspoly_Xsize];	//q_temp[z-deg+1][y_size>max(deg_y)+1+(rs-1)*(max(deg_y) in encoding functions)][14>w+(rs-1)*w], lc[rs]--leading coefficient polynomial
	int d1,d2,flag;
	//printf("\nWhen u=%d\n", u);
	lm=0; lm_temp=0;	//leading monomial index
	act=0;	//judge value for recursive search of each f_k-1-u
	//find pb_k-1-uu
	j_1=-1;
	i_1=-1;
	for(j=0;j<40;j++)	//y_size, efficiently enough
	{	
		for(i=0;i<w+1;i++)	//w+1
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
	for(i=0;i<2;i++)	//rs
		for(j=0;j<rcspoly_Ysize;j++)	//y_size
			for(u=0;u<rcspoly_Xsize;u++)	//5>w*(rs-1)+w
				q_temp[i][j][u]=0;


	//Calculate q_temp[pb_k-1-u]=q[u][pb_k-1-u]
	for(i=0;i<2;i++)	//rs
	{
		for(j=0;j<facpoly_Ysize;j++)	//y_size=max(deg_y) in encoding function*(rs-1), and 31>=max(deg_y)+1
		{
			for(u=0;u<facpoly_Xsize;u++)	//
			{
				q_temp[i][j+i*j_1][u+i*i_1]=Q[uu][i][j][u];	//this time, q_temp size is [rs+1][maxY(j_1)+maxY(Q)][maxX(i_1)+maxX(Q)]
			}
		}
	}


	//convert x^3=y^2+y;
	d1=0;
	d2=0;
	for(i=0;i<2;i++)	//rs
		{	
			for(j=0;j<9;j++)	//y_size=degY+1
			{
				for(u=3;u<5;u++)	//w+1, 5>w*(rs-1)+w
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
	for(i=0;i<2;i++)	//rs
	{
		for(j=0;j<11;j++)	//y_size
		{
			for(u=0;u<w+1;u++)	//w+1
			{
				if(q_temp[i][j][u]!=0)
				{
					for(v=0;v<20;v++)	//number of pb, *take a note*--pole basis vialation might happen, as 95 pole basis functions can not promise to term[y_size-1][w]
					{
						for(m=0;m<12;m++)	//y_size+1
						{	
							for(z=0;z<w+1;z++)	//w+1
							{
								if(pb[v][m][z]==1 && j==m && u==z)
								{
									lm_temp=v;
									if(lm_temp>=lm)
									{
										lm=lm_temp;	//record the leading monomial's index
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

	//find the leading coefficient polynomial lc[rs]
	for(i=0;i<2;i++)	//rs
	{
		flag=1;
		lc[i]=-1;
		for(j=0;j<11;j++)	//y_size
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
	for(i=0;i<2;i++)
		if(lc[i]==0)
			flag++;

	//printf("leading coefficient polynomial");
	
	u=0;	//root index
	//find the roots of leading coefficient polynomial lc[rs]
	for(i=0;i<4;i++)	//number of elements in GF(4)
		if(flag!=2)	//rs
		{
			b=0;
			for(j=0;j<2;j++)	//rs
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
		else if(flag==2)
			{
				rootlist[uu][0]=0;
				act=1;
			}
	//printf("roots");

	//For each distinct root of rootlist[u]
	if(act==1)
	{
	for(i=0;i<2;i++)	//2>rs
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
				for(j=0;j<2;j++)	//rs
					for(u=0;u<rcspoly_Ysize;u++)	//y_size
						for(m=0;m<rcspoly_Xsize;m++)	//3>w*(rs-1)+w
							q_temp[j][u][m]=0;	
				
				//q_temp=q[uu][z+f_k-1-u*pb_k-1-u]
				for(j=0;j<2;j++)	//rs
				{
					//calculate (z+f_k-1-u*pb_k-1-u)^j
					if(j==0)
					{
						for(u=0;u<2;u++)
							for(m=0;m<4;m++)
								for(z=0;z<4;z++)
									expoly[u][m][z]=0;
						expoly[0][0][0]=1;
					}
					else if(j>0)
							polyexp1(alpha, i_1, j_1, j);
					
					//calculate q[uu][z+f_k-1-u*pb_k-1-u]
					for(u=0;u<2;u++)	//rs
					{
						for(m=0;m<4;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
						{
							for(z=0;z<4;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
							{
								if(expoly[u][m][z]!=0)
								{
									for(v=0;v<facpoly_Ysize;v++)	//y_size
									{
										for(t=0;t<facpoly_Xsize;t++)	//w+1
										{
											if(Q[uu][j][v][t]!=0)
												q_temp[u][v+m][t+z]=add(q_temp[u][v+m][t+z], mul(expoly[u][m][z], Q[uu][j][v][t]));
										}
									}
								}
							}
						}
					}
				}
				
				//convert x^3=y^2+y
				d1=0;
				d2=0;
				for(j=0;j<2;j++)	//rs
				{
					for(u=0;u<10;u++)	//y_size
					{
						for(m=w+1;m<4;m++)	//w+1, 14>w*(rs-1)+w
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
				for(j=0;j<3;j++)	//rs
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

void tgorder()
{
	int i,u;

	for(i=0;i<100;i++)
	{
		tg_order_1[i][0]=-1;
		tg_order_1[i][1]=-1;
		tg_order[i][0]=-1;
		tg_order[i][1]=-1;
	}

	//total graduate order
	tg_order_1[0][0]=0;
	tg_order_1[0][1]=0;

	for(i=1;i<100;i++)
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
	u=0;
	for(i=0;i<100;i++)
	{
		if(tg_order_1[i][0]<w+1)	//w+1
		{
			tg_order[u][0]=tg_order_1[i][0];
			tg_order[u][1]=tg_order_1[i][1];
			u++;
		}
	}
	/*
	for(i=0;i<100;i++)
		printf("\n(%d, %d)", tg_order[i][0], tg_order[i][1]);
	*/
}

void mono_table()
{
	int i, j, u, v, l, weight[26][100], mono_order_1[26][100];	//26=(100/weight(z))+1, 100=pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][max(deg_y)+1][w+1]

	//Initialisation
	for(i=0;i<26;i++)	//(230/weight(z))+1
	{
		for(j=0;j<100;j++)	//pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][y_size][w+1].
		{
			weight[i][j]=-1;
			mono_order_1[i][j]=-1;
		}
	}
	for(i=0;i<2;i++)	// rs
		for(j=0;j<8;j++)	//max(deg_y)+1
			for(u=0;u<w+1;u++)	//w+1
				mono_order[i][j][u]=-1;

	//weight of monomial over function x^3-y^2-y=0
	for(i=0;i<26;i++)	
	{
		weight[i][0]=weiz*i;
		for(j=1;j<100;j++)
			weight[i][j]=weiz*i+j+1;
	}

	l=0;
	for(v=0;v<100;v++)	//for each possible weight until weight(term[rs-1][max(deg_y)+1][w-1])+1, note term[rs-1][max(deg_y)+1][w-1] is the term next to term[rs-1][max(deg_y)][w]
	{
		for(i=0;i<26;i++)	//230/weight(z)
		{
			for(j=0;j<100;j++)	//230, >the index of the largest term with deg(z)=0, and weight=weight(term[rs-1][max(deg_y)+1][w-1])
			{
				if(weight[i][j]==v)
				{
					mono_order_1[i][j]=l;
					l++;
				}
			}
		}
	}

	// 2-dimensional table transform into 3-dimensional table
	for(u=0;u<2;u++)	//rs
		for(l=0;l<40;l++)	//>index of term[max(deg_y)][w] in the pole basis + 1
			mono_order[u][tg_order[l][1]][tg_order[l][0]]=mono_order_1[u][l];

	/*
	printf("\nMonmial basis is:\n");
	for(i=0;i<3;i++)
	{
		for(j=0;j<25;j++)
		{
			for(k=0;k<5;k++)
				printf("%d ", mono_order[i][j][k]);
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
	for(v=0;v<80;v++)	//num of pb
		for(i=0;i<20;i++)	//y_size
			for(j=0;j<w+1;j++)	//x_size
				pb[v][i][j]=0;


	for(v=0;v<80;v++)	//num of pb
		pb[v][tg_order[v][1]][tg_order[v][0]]=1;
}
/*
void polyexp(int xi, int yi, int b)
{
	int i, j, k, c1, p1[26][25], p2[25][26], p3[25][25];	//p1[max(deg_y)+2][max(deg_y)+1], p2[max(deg_y)+1][max(deg_y)+2], p3[max(deg_y)+1][max(deg_y)+1]
	for(i=0;i<25;i++)	//max(deg_y)+1
	{
		for(j=0;j<25;j++)	//max(deg_y)+1
		{
			part_com[i][j]=0;
			p3[i][j]=0;
		}
	}
	for(i=0;i<26;i++)	//max(deg_y)+2
		for(j=0;j<25;j++)	//max(deg_y)+1
			p1[i][j]=0;
	for(i=0;i<25;i++)	//max(deg_y)+1
		for(j=0;j<26;j++)	//max(deg_y)+2
			p2[i][j]=0;

	if(b==0)
		part_com[0][0]=1;

	else
	{
	part_com[1][0]=1;
	part_com[0][1]=power(xi, 4);	//xi^w
	part_com[0][0]=add(yi, power(xi, 5));	//xi^(w+1)+yi
	c1=part_com[0][0];

	for(k=1;k<b;k++)
	{
		for(i=0;i<25;i++)	//max(deg_y)+1
		{
			for(j=0;j<25;j++)	//max(deg_y)+1
			{
				p1[i+1][j]=part_com[i][j];
				p2[i][j+1]=mul(power(xi, 4), part_com[i][j]);	//xi^w
				p3[i][j]=mul(c1,part_com[i][j]);
			}
		}

		for(i=0;i<25;i++)	//max(deg_y)+1
			for(j=0;j<25;j++)	//max(deg_y)+1
				part_com[i][j]=add(add(p1[i][j],p2[i][j]),p3[i][j]);

	}
	}
}
*/
void polyexp1(int c, int i, int j, int deg_z)
{
	int u, m, z, v, p1[2][4][4], p2[1][6][6];	//p1[rs+1][18>(max(deg_y) in encoding functions)*(rs-1)][10>(max(deg_x) in encoding functions)*(rs-1)], p2[rs][26=18+max(deg_y) in encoding functions][14=10+max(deg_x) in encoding functions]

	//Initialisations
	for(u=0;u<2;u++)	//rs
		for(m=0;m<4;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
			for(z=0;z<4;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
				expoly[u][m][z]=0;
	for(u=0;u<2;u++)	//rs+1
		for(m=0;m<4;m++)	//18>(max(deg_y) in encoding functions)*(rs-1)
			for(z=0;z<4;z++)	//10>(max(deg_x) in encoding functions)*(rs-1)
				p1[u][m][z]=0;
	for(u=0;u<1;u++)	//rs
		for(m=0;m<6;m++)	//26=18+max(deg_y) in encoding functions
			for(z=0;z<6;z++)	//14=10+max(deg_x) in encoding functions
				p2[u][m][z]=0;

	if(deg_z==0)
		expoly[0][0][0]=1;
	else
	{
		expoly[1][0][0]=1;
		expoly[0][j][i]=c;

		for(u=1;u<deg_z;u++)	//u=1
		{
			for(m=0;m<2;m++)	//rs
			{
				for(z=0;z<4;z++)	//18>(max(deg_y) in encoding functions)*(rs-1)
				{
					for(v=0;v<4;v++)	//10>(max(deg_x) in encoding functions)*(rs-1)
					{
						p1[m+1][z][v]=expoly[m][z][v];	//z*poly
						p2[m][z+j][v+i]=mul(c, expoly[m][z][v]);	//a*xy*poly
					}
				}
			}
			
			//p1+p2
			for(m=0;m<2;m++)	//rs
				for(z=0;z<4;z++)	//18>(max(deg_y) in encoding functions)*(rs-1)
					for(v=0;v<4;v++)	//10>(max(deg_x) in encoding functions)*(rs-1)
						expoly[m][z][v]=add(p1[m][z][v], p2[m][z][v]);
		}
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
	int i,c[2],d[2],e[2],f;
	unsigned int mask=1;

	for(i=1;i>=0;i--)
		{
		if((fac1 & mask)>0)
			c[i]=1;
		else
			c[i]=0;

		mask=mask<<1; //shift 1 bit left
		}

	mask=1;
	for(i=1;i>=0;i--)
		{
		if((fac2 & mask)>0)
			d[i]=1;
		else
			d[i]=0;

		mask=mask<<1;
		}
	for(i=0;i<2;i++)
		e[i]=c[i]^d[i];


	f=e[0]*2+e[1]*1; //conver to decimal

	return f;
}

int power(int a, int b)
{
	int i,temp,pow_result;

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
