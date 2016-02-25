
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/********change element with different code************
n,k,q,p			m-->iter_num,tm,lm,mono_size,fa_size
a[n]
*****************************************************/

//***********var define*****************
#define n 15	//length of codeword
#define k 9 //length of message
#define q 16 //GF(q) q=n+1
#define p 4 //GF(2^p) = GF(q)
#define maxcost 300 //Designed maximal interpolation cost
#define m 1 // zero of multiplicity
#define iter_num 15  // C
#define tm 3
#define lm 1 //Desgined factorization output list size   (should be mended)
#define poly_Ysize 3  // Ysize==lm+1
#define poly_Xsize 20 // C+1 approximately 
#define poly_lod_maxsize 10000 //when k is very large, then maxsize should be modified
#define coeff_size 20 // the size of fac coeff
#define mono_size 20  // the size of inter_polynomial
#define fa_size 20	  // the number of layers of factorization'tree 
//********************************

//***************basic function statement**************
int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
void mono_table(void);
void encoder(void);
void modulation(void);
void channel(void);
void demodulation(void);
void interpolation(void);
void factorization(void);
//**********************************************


//**************my function statement***************	
int cal_delta(int,int,int);
int comb(int,int);
int cal_poly_lod(int);
void update_interpoly(int);
//*****
void DFS(int);
int cal_rootlist(int,int);
void cal_update_poly(int,int,int);
void output_message_poly(int);
void choose_correct_dec_codeword(void);
void re_encoder(int,int reEncode_codeword[]);
float cal_EucDis(int reEncode_codeword[],int bicodeword_temp[],int symbol_temp[][2]);
int cal_mindistan_seqnum(float distance_temp[]);



//*******
/*
void choose_correct_dec_codeword(void);
void re_encoder(int,int arr1[]);
float cal_EucDis(int arr1[],int arr2[],int arr3[][2]);
int cal_mindistan_seqnum(float arr1[]);
*/
//*********************************************

//**************basic var*****************
unsigned long int seq_num;	//number of input binary sequences
int SNR;
double BER, FER;
int mularray[]={1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};  //this array is used to the infinite field element mul
int a[15]={1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};      //caution: this array must change with different sequence of infinite element
int root[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};	//n+1		be used to the factorization

int mono_order[40][200];	//multiplexity matrix
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

//*********
int inter_point[2];  //interlolation points (pi,ri)
int inter_poly[(lm+1)][poly_Ysize][poly_Xsize];  //store the inter_poly of everygroup
int inter_poly_minlod[poly_Ysize][poly_Xsize];  // the poly with minimium used to factorization
int poly_lod[(lm+1)];			//store the lod of inter_poly of everygroup
int modify_flag[(lm+1)];      //store the flag if the poly should be continued to update 
int delta[lm+1];					   //store the group of delta, totalling


//int message_poly[fa_size][poly_Xsize];   //the array of poly used in factorization
int time;     // factorization time;
//factorization's element
int pai[fa_size],deg[fa_size];
int Q[fa_size][poly_Ysize][poly_Xsize+lm+1];	//poly_Ysize is fixed, but the Xsize of Q must be increase to (C+lm)
int Coeff[coeff_size];
int f[n][k];  //the total number of f is n is no sure,store the message polynomial after the fac, note that may be overflow
int f_num=0;
//int f_error_temp=0;

//*********


void main()
{
	int i, x, s, start, finish, value, count;
	unsigned long int j;
	unsigned int mask=1;
	long int error, ferror;
	double progress, b;

	srand(1977);

	mono_table();
/*
	//*********debug*********
	printf("mono_table:\n");
	for(j=0;j<mono_size;j++)
	{
		for(i=0;i<mono_size;i++)
			printf("\t%d",mono_order[j][i]);
		printf("\n");
	}
	//*******************
*/

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

		//printf("\n\n*************For the %dth frame*************:\n", j);
		//Generate binary message
			for(i=0;i<k*p;i++)
				bi_message[i]=rand()%2;

			//Convert to nonbinary
			for(i=0;i<k;i++)
				message[i]=bi_message[p*i]+2*bi_message[p*i+1]+4*bi_message[p*i+2]+8*bi_message[p*i+3];

/*			//******debug*******
			message[0]=6;	message[1]=9;	message[2]=2;	message[3]=6;
			message[4]=10;	message[5]=8;	message[6]=2;	message[7]=0;
			message[8]=0;
*/
			//*******************

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

			interpolation();

			factorization();
			
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
			if( d[j]<temp )
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

void interpolation(void)
{
	int i,j,h,temp,poly_lod_temp;
	int tk;   			// tk is the num of iteration, t is the number of inter_point
	int minlod_seq,minlod=0;	//minlod_seq is the j' with delta[j']!=0

	int alpha,beta;

	//********initialize****************
	//initialize group_inter_poly
    for(h=0;h<(lm+1);h++)
   	{
   		for(j=0;j<poly_Ysize;j++)
   		{
   			for(i=0;i<poly_Xsize;i++)
   			{
   				inter_poly[h][j][i]=0;
   			}
   		}
   	}

   	//initialize Group_0
	//initialize the flag of poly should be modified   1-->modified   0-->stop to modified
    for(j=0;j<(lm+1);j++)
    {
		inter_poly[j][j][0]=1;
		modify_flag[j]=1;
	}

	//cal the lod of G0
	for(h=0;h<(lm+1);h++)
		poly_lod[h]=cal_poly_lod(h);

	//*************************************

    //	statement:	complexity reduction have not been used, if it been done,
    //				the poly_lod must be calculate at first
    //				calculate the lod of poly and store it in poly_lod[i]
	//				for(h=0;h<(lm+1);h++)
  	//					poly_lod[h]=cal_poly_lod(h);

     //**********begin to interpolation**********************
	//iteration of modifying the inter_polynomial   -------> optimize to a function
    for(tk=0;tk<n;tk++)
    {
		inter_point[0]=a[tk]; 		// the first word is pi--->xi
		inter_point[1]=rx_codeword[tk];	// the second word is ri--->yi

		for(alpha=0;alpha<m;alpha++)
			for(beta=0;beta<(m-alpha);beta++)
			{
				//calculate delta_i and poly_lod
    			for(h=0;h<(lm+1);h++)
					if( modify_flag[h]==1 )
    					delta[h]=cal_delta(h,alpha,beta);

    			//choose the delta_i with delta != 0 and j' = minlod_seq
    			temp=poly_lod_maxsize;
    			for(h=0;h<(lm+1);h++)
    				if( (delta[h]!=0) && (modify_flag[h]==1) )
    					if( temp > poly_lod[h] )
  						{
 	 						temp = poly_lod[h] ;
  							minlod_seq=h;
  						}
  				
  				//calculate the new inter polynomial -------> update
  				update_interpoly(minlod_seq);

				//cal the lod of Group
				for(h=0;h<(lm+1);h++)
					if(modify_flag[h]==1)
					{	
						poly_lod[h]=cal_poly_lod(h);
						
						if(poly_lod[h]>iter_num)
							modify_flag[h]=0;
					}

			}
		}

	//find the interpolation polynomial of Group_C with minimium lod-----------> inter_poly[minlod]
 	temp=poly_lod_maxsize;
 	for(h=1;h<(lm+1);h++)
		if(modify_flag[h]==1)
		{
			poly_lod_temp=cal_poly_lod(h);

		 	if( poly_lod_temp<temp )
			{
		  		temp=poly_lod_temp;
				minlod_seq=h;
			}	
		}
   	
	// find out the inter_poly_minlod used to factorization
    for( j=0;j<poly_Ysize;j++)
		for( i=0;i<poly_Xsize;i++)
			inter_poly_minlod[j][i]=inter_poly[minlod_seq][j][i]; 

	//************debug***************
	//check the inter_poly_minlod must be over the every inter_point m=1 times
	int temp1[n];
	for(h=0;h<n;h++)
	{
		temp1[h]=0;
		temp=0;
		for( j=0;j<poly_Ysize;j++)
			for( i=0;i<poly_Xsize;i++)
				if(inter_poly_minlod[j][i]!=0)
				{
					temp=mul( power(a[h],i),power(rx_codeword[h],j) );
					temp=mul( temp,inter_poly_minlod[j][i] );
					temp1[h]= add(temp1[h],temp);
				}
	}

	temp=0;
	for(h=0;h<n;h++)
		if(temp1[h]!=0)
			temp++;

	if(temp!=0)
		printf("\n\ninterpolation is failed, error=%d\n\n",temp);

	//***********************************

	//*****interpolation debug***********
/*	//calculate the codeword score
	int flag2=0,flag3=0;
	for(i=0;i<n;i++)
		if( codeword[i] == rx_codeword[i] )
			flag2++;
	flag2=flag2*m;

	//calculate the deg of Q(x,y)
		for(j=0;j<poly_Ysize;j++)
			for(int i=0;i<poly_Xsize;i++)
				if( inter_poly_minlod[j][i] != 0 )
					if(	flag3 < ( i+j*(k-1) ) )
						flag3=( i+j*(k-1) );

	if( flag2>flag3 )
		printf("\n\nthe codeword score %d > degree Q(x,y) %d\n",flag2,flag3);
	else if(flag2<=flag3)
		printf("\n\nthe codeword score %d <= degree Q(x,y) %d\n",flag2,flag3);
*/	//****************************************

}		

void factorization(void)
{
	int i,j,h,flag,u;
	int error_1=0;	  //erro1 is the number of codeaword error after demodulation
	int error_2=0;    //error_2=1 means message_poly being worked out after fac
	int check[n];

	//*************initialize**************
	//initialize Q
	for(h=0;h<fa_size;h++)
		for(j=0;j<poly_Ysize;j++)
			for(i=0;i<(poly_Xsize+lm+1);i++)
				Q[h][j][i]=0;

	//initialize pai[] and deg[]
	for(i=0;i<fa_size;i++)
		pai[i]=deg[i]=-1;

	//initialize Coeff[]
	for(i=0;i<coeff_size;i++)
		Coeff[i]=-1;

	//intialize the result polynomial of factorization
	for(h=0;h<n;h++)
		for(i=0;i<k;i++)
			f[h][i]=-1;	  // 0 also is a valid codeword bit

	//initialize Q0=<<inter_poly_minlod>>   --> have been tested
	flag=0;
	h=0;
	for(i=0;i<poly_Xsize;i++)
	{
		for(j=0;j<poly_Ysize;j++)
			if( inter_poly_minlod[j][i] != 0)
			{
				flag=1;
				break;
			}

			if( flag == 1 )
			{
				h=i;
				break;
			}
	}

	for(j=0;j<poly_Ysize;j++)
		for(i=0;i<(poly_Xsize-h);i++)
			if( inter_poly_minlod[j][i+h] != 0 )
				Q[0][j][i]=inter_poly_minlod[j][i+h];

/*	//*********debug**************
	printf("\ninter_poly_minlod=\n");
	for(j=0;j<poly_Ysize;j++)
	{
		for(i=0;i<poly_Xsize;i++)
			printf("\t%d",inter_poly_minlod[j][i]);
		printf("\n");
	}

	printf("\nQ[0]=\n");
	for(j=0;j<poly_Ysize;j++)
	{
		for(i=0;i<(poly_Xsize+lm+1);i++)
			printf("\t%d",Q[0][j][i]);
		printf("\n");
	}
*/	//********************************

	//*********begin factorization*******************
	f_num=0;
	pai[0]=-1;
	deg[0]=-1;
	time=1;
	u=0;

	DFS(u);
	//****************************
	// choose the final dec_codeword
	if(f_num >= 1)
		choose_correct_dec_codeword();  
	else if(f_num == 0)
	{
		for(i=0;i<n*p;i++)
				dec_codeword[i]=rx_bicodeword[i];
	//	printf("\ndecoding failed!!");
	}
	//at last output the binary dec_codeword[i]   	

		printf("\nrx_codeword=\t");
			for(i=0;i<n;i++)
				printf("%d ",rx_codeword[i]);
			
			printf("\ncodeword=\t");
			for(i=0;i<n;i++)
				printf("%d ",codeword[i]);

			printf("\nmessage_poly=\t");
			for(i=0;i<k;i++)
				printf("%d ",message[i]);

			printf("\n\t\t\t\t\t\t\t\tf_num=\t%d",f_num);
			for(j=0;j<=f_num;j++)
			{
				printf("\nf(x)_%d=\t",j);
				for(i=0;i<k;i++)
					printf("%d ",f[j][i]);
			}
			printf("\n\n");




/*
	// choose the final dec_codeword-->cheating
	flag=0;
	for(j=0;j<f_num;j++)
	{
		for(i=0;i<k;i++)
			if(f[j][i]==message[i])
				flag++;
			else flag=0;

			//once the message_poly can be found in the fac selection, we break out the loop
			if(flag==k)	
			{
				for(i=0;i<n*p;i++)
					dec_codeword[i]=bi_codeword[i];
			//	printf("\ndecoding succeed!!!!");
				break;
			}
			else flag=0;

	}

	if(flag==0)
		{
			for(i=0;i<n*p;i++)
				dec_codeword[i]=rx_bicodeword[i];
			
			printf("\ndecoding failed!!");

	
		}
*/	


/*	//*****************debug*************
	//validate output
	error_1=0;
	for(h=0;h<n;h++)	
	if( codeword[h] != rx_codeword[h] )
		error_1++;
//	printf("\t\terror=%d\n",error_1);

	flag=0;
	if(error_1<=tm)
	{
		//check the fac selection f including the right message_poly or not 
		for(j=0;j<f_num;j++)
		{
			check[j]=0;
			for(i=0;i<k;i++)
				if(f[j][i]==message[i])
					check[j]++;

			if( check[j] == k )
				flag++;
		}



		if( flag != 1)
		{
			printf("\t\terror is happened, flag=%d\n",flag);
//			printf("\n\ncheck=");
//			for(i=0;i<n;i++)
//				printf("\t%d",check[i]);

//			for(h=0;h<f_num;h++)
//			{
//				printf("\nf[%d]=",h);
				
//				for(j=0;j<k;j++)
//					printf("\t%d",f[h][j]);
//			}

			printf("\n\n");
//			f_error_temp=1;
		}
		
	}
*/
	//******************************
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
		


int cal_delta(int h,int alpha,int beta)    //calculate the delta of inter_polynomial, c is inter_poly[i] and t is tk
{
	int delta=0,delta_temp,i,j,flag;
	int a,b,x,y;

	x=inter_point[0];			//xi=pi
	y=inter_point[1];			//yi=ri

	for(j=0;j<poly_Ysize; j++)
	{
		for(i=0;i<poly_Xsize;i++)
		{
			if( (inter_poly[h][j][i] != 0) && ( i >= alpha) && ( j >= beta) )
			{
				a=i;
				b=j;
				flag= comb(a,alpha)*comb(b,beta);
				// if flag is even number , delta_temp is equal to 0
				// if flag is odd number, delta_temp can be calculated to the delta
				if( (flag%2) !=0 )
				{
					delta_temp=mul(inter_poly[h][j][i],power(x,(a-alpha)));
					delta_temp=mul(delta_temp,power(y,(b-beta)));
					delta= add(delta,delta_temp);
				}
			}
		}	
	}

	return delta;
}


int comb(int a,int b)  //calculate combination(a_up,b_down)
{
	int i,temp1,temp2=1;

	for(i=0;i<b;i++)
	{
		temp1=(a-i)/(b-i);
		temp2*=temp1;
	}

	return temp2;
}


int cal_poly_lod(int h)   //input one inter_polynomial, output this poly's leading order
{
	int j,i;
	int max_ord=0;

	for(j=0;j<poly_Ysize;j++)
		{
			for(i=0;i<poly_Xsize;i++)
			{
				if( inter_poly[h][j][i] != 0 )
				{
					if(	max_ord < mono_order[j][i])
					{
						max_ord=mono_order[j][i];
					}
				}
			}
		}
	return max_ord;
}


void update_interpoly(int minlod_seq)
{
	int h,j,i,temp;
	int poly_temp[poly_Ysize][poly_Xsize];
	int poly_temp1[poly_Ysize][poly_Xsize+1],poly_temp2[poly_Ysize][poly_Xsize+1];

	for(j=0;j<poly_Ysize;j++)
		for(i=0;i<poly_Xsize;i++)
			poly_temp[j][i]=inter_poly[minlod_seq][j][i];


	//start to update
	for(h=0;h<(lm+1);h++)
		if( (delta[h]!=0) && (modify_flag[h]==1) )
		{
			//refresh poly_temp1 & poly_temp2
			for(j=0;j<poly_Ysize;j++)
   				for(i=0;i<(poly_Xsize+1);i++)
   					poly_temp1[j][i]=poly_temp2[j][i]=0;

			if( h!=minlod_seq)
			{
				// delta[minlod_seq] * inter_poly[h]
				for(j=0;j<poly_Ysize;j++)
					for(i=0;i<poly_Xsize;i++)
						if( inter_poly[h][j][i] != 0 )
							poly_temp1[j][i]=mul(inter_poly[h][j][i],delta[minlod_seq]);
				// delta[h] * inter_poly[minlod_seq]
				for(j=0;j<poly_Ysize;j++)
					for(i=0;i<poly_Xsize;i++)
						if( poly_temp[j][i] != 0 )
							poly_temp2[j][i]=mul(poly_temp[j][i],delta[h]);
				//  new inter_poly[i+1] = new inter_poly[i] + new inter_ploy[minlod_seq]
				for(j=0;j<poly_Ysize;j++)
					for(i=0;i<poly_Xsize;i++)
						inter_poly[h][j][i]=add(poly_temp1[j][i],poly_temp2[j][i]);	
			}
			else if( h==minlod_seq )
			{
				// delta[minlod_seq]*x*inter_poly[minlod]
				for(j=0;j<poly_Ysize;j++)
					for(i=0;i<poly_Xsize;i++)	
						if( poly_temp[j][i] != 0)
							poly_temp1[j][i+1]=mul(poly_temp[j][i],delta[h]);		
				
				// delta[minlod_seq]*xi*inter_poly[minlod_seq]
				temp=mul(delta[h],inter_point[0]);    //delta[minlod_seq]*xi	
				for(j=0;j<poly_Ysize;j++)
					for(i=0;i<poly_Xsize;i++)
       					if( poly_temp[j][i] != 0)
       						poly_temp2[j][i] = mul(poly_temp[j][i],temp);  

       			// inter_poly + inter_poly_temp
   				for(j=0;j<poly_Ysize;j++)
					for(i=0;i<poly_Xsize;i++)   
						inter_poly[h][j][i]=add(poly_temp1[j][i],poly_temp2[j][i]);
			}					 						     						

		}						
}



void DFS(int u)
{
	int v,i,h,temp,j,flag;
	int rootlist[lm];		//lm is the size of rootlist, because max deg_y(Q0)=lm

	//initialize rootlist
	for(h=0;h<lm;h++)
		rootlist[h]=-1;			//notification, rootlist=0 also is a effective codeword

	//judge the Qu(x,0)?=0
	flag=1;
	for(i=0;i<(poly_Xsize+lm+1);i++)
		if( Q[u][0][i] != 0 )		//if one of coeff !=0, flag=0;
		{
			flag=0;
			break;
		}

	//caution: serve the last degree of polynomial candidate, special state
/*	if( (flag==1) && (deg[u]==k-2) )
	{
		rootlist[0]=0;
		v=time;
		time=time+1;
		pai[v]=u;
		deg[v]=deg[u]+1;
		Coeff[v]=rootlist[0];
		u=v;
	}
*/	//************************

	//start the loop
	if( flag==1 && deg[u]%(k-1)==0 && deg[u]>0 )
		output_message_poly(u);
	else if(deg[u] < (k-1))		// the "D" of thesis represent k-1
	{
		//try n field element to calculate the y-root of Q(0,y)=0
		j=0;
		for(h=0;h<q;h++)		// q is the size of finite field
		{
			temp=cal_rootlist(u,h);
			if( temp == 0 )
			{	
				rootlist[j]=root[h];
				j++;
			}
		}


		for(i=0;i<lm;i++)
			if( rootlist[i]!=-1 )
			{	
				v=time;
				time=time+1;
				pai[v]=u;
				deg[v]=deg[u]+1;
				Coeff[v]=rootlist[i];
				cal_update_poly(v,u,rootlist[i]);  
				DFS(v);
			}

	}	

}


int cal_rootlist(int u,int h)
{
	int j,temp;

	//traverse all the field element root[i]
	temp=0;
	for(j=0;j<poly_Ysize;j++)
		if( Q[u][j][0] != 0)
		{
			temp=add( temp , mul( Q[u][j][0] , power(root[h],j) ) );
		}
	return temp;
}


void output_message_poly(int u)
{
	int h,temp,deg_temp;

	//output_message_poly
	temp=u;
	deg_temp=deg[u];
	while( deg_temp>(-1) && pai[temp]>(-1) && Coeff[temp]>(-1) )
	{
		f[f_num][deg_temp]=Coeff[temp];
		temp=pai[temp];
		deg_temp=deg[temp];
	}

	//check message_poly's validity
	temp=0;
	for(h=0;h<k;h++)
		if( f[f_num][h]>(-1) )
			temp++;
	//when poly f is k bits enough, f is a valid output message poly
	if( temp==k )
		f_num=f_num+1; 
}

void cal_update_poly(int v,int u,int a)
{
	int i,j,g,h,flag;
	int Q_update[lm+1][lm+1];  //the maxsize is (xy+a)^i
	int temp1[poly_Ysize][poly_Xsize*2],temp2[poly_Ysize][poly_Xsize*2];
	int update_temp[lm+1][poly_Xsize+lm+1];

	//clean Q_update[lm+1][lm+1], and temp1 & temp2
	for(j=0;j<(lm+1);j++)
		for(i=0;i<(lm+1);i++)	
			Q_update[j][i]=0;

	//clean update_temp[lm+1][poly_Xsize+lm+1]
	for(g=0;g<(lm+1);g++)
		for(h=0;h<(poly_Xsize+lm+1);h++)	
			update_temp[g][h]=0;

	//Q_update = a+xy;
	Q_update[0][0]=a;
	Q_update[1][1]=1;

	//**************begin factorization*********************
	//calculate the Q(x,xy+a)
	for (j=0; j < poly_Ysize ; j++) // when j=0, the poly does not change
	{
		//clean temp1 & temp2
		for(g=0;g<poly_Ysize;g++)
			for(h=0;h<poly_Xsize*2;h++)	
				temp1[g][h]=temp2[g][h]=0;

		if( j>=2 ) //calculate (xy+a)^j, and store it into the Q_update
		{
			//temp1 = Q_update * xy
			for(g=0;g<(lm+1);g++)
				for(h=0;h<(lm+1);h++)
					if( Q_update[g][h] != 0 )
						temp1[g+1][h+1]=Q_update[g][h];
			//temp2 = Q_update * a
			for(g=0;g<(lm+1);g++)
				for(h=0;h<(lm+1);h++)
					if( Q_update[g][h] != 0 )
						temp2[g][h] = mul(Q_update[g][h],a);
			// Q_update = temp1 + temp2 = (Q_update*a) + (Q_update*xy)
			for(g=0;g<(lm+1);g++)
				for(h=0;h<(lm+1);h++)
						Q_update[g][h] = add( temp1[g][h] , temp2[g][h] );
		}//caution Q_update[poly_Ysize][poly_Xsize] just enough, but it will overflow


		//begin to y=(xy+a)^j and calculate (bx^i)*y
		if(j==0)  //j=0, unchange
		{
			for(i=0;i<(poly_Xsize+lm);i++)
				if( Q[u][0][i] != 0)
					Q[v][0][i]=Q[u][0][i];
		}
		else if( j>=1 )
		{
			for(i=0;i<poly_Xsize;i++)	//traverse row j's every x^i
			{
				if( Q[u][j][i]!=0 )	//  coeff *x^i * (xy+a)^j
					for(g=0;g<(lm+1);g++)   //the maxsize of Q_update[g][h] is lm*lm  
							for(h=0;h<(lm+1);h++)
								if( Q_update[g][h]!=0 )
									update_temp[g][h+i]=mul( Q[u][j][i],Q_update[g][h] );
				
				for(g=0;g<(lm+1);g++)    //the maxsize of update_temp[g][h] is lm*(poly_Xsize+lm)
					for(h=0;h<(poly_Xsize+lm+1);h++)
							if( update_temp[g][h]!=0 )
								Q[v][g][h]=add( Q[v][g][h],update_temp[g][h] );

				//clean update_temp[lm+1][poly_Xsize+lm+1]
				for(g=0;g<(lm+1);g++)
					for(h=0;h<(poly_Xsize+lm+1);h++)	
						update_temp[g][h]=0;

			}
		}		

	}
	//*******************************

	//calculate <<Q(x,y)>> = Q(x,y)/xm
	flag=0;
	h=0;
	for(i=0;i<(poly_Xsize+lm+1);i++)  //the maxsize of Q[j][i] is poly_Ysize*(poly_Xsize+lm+1)
	{
			for(j=0;j<poly_Ysize;j++)
				if( Q[v][j][i] != 0)
				{
					flag=1;
					break;
				}

			if( flag == 1 )
			{
				h=i;
				break;
			}
		}
	
		if(h>0)
			for(j=0;j<poly_Ysize;j++)	//the maxsize of Q[j][i] is poly_Ysize*(poly_Xsize+lm+1)
				for(i=0;i<(poly_Xsize+lm+1-h);i++)
					{
						Q[v][j][i]=Q[v][j][i+h];
						Q[v][j][i+h]=0;
					}
	//**************************
}

//********find the correct f(x)*****************
void choose_correct_dec_codeword(void)
{
	int i,j,value,min_distance_seqnum;
	unsigned int mask=1;
	int reEncode_codeword[n];   //store the codeword of reencoding f(x)
	int bicodeword_temp[n*p];   // convert nonbinary reEncode_codeword into binary bicodeword_temp
	int symbol_temp[n*p][2];	//convert binary bicodeword_temp into coordinate symbol_temp
	float distance_temp[n];		//store the Euclidean distance between symbol_temp(int) and rx_symbol(float)

	//clean the distance_temp[n]
	for(i=0;i<n;i++)
		distance_temp[i]=0;    

	if(f_num>1)
	{
		for(i=0;i<f_num;i++)
		{
			re_encoder(i,reEncode_codeword);
			
			distance_temp[i]=cal_EucDis( reEncode_codeword , bicodeword_temp , symbol_temp );
		
		}
		//find out the min_distance_seqnum
		min_distance_seqnum=cal_mindistan_seqnum(distance_temp);   
	}
	else if(f_num==1)
		min_distance_seqnum=0;

	//after knowing the min_distance f(x)'i =  min_distance_seqnum, reencoder this fi(x)
	//get the final binary dec_codeword
	re_encoder(min_distance_seqnum,reEncode_codeword);

	//Convert the codeword into binary
	for(i=0;i<n;i++)
		{
			value=reEncode_codeword[i];
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


void re_encoder(int u,int reEncode_codeword[])
{
	int i,j,temp;
	
	//Encoding
	for(i=0;i<n;i++)
	{
		reEncode_codeword[i]=0;
		for(j=0;j<k;j++)
		{
			temp=mul( f[u][j],power(a[i],j) );
			reEncode_codeword[i]=add(reEncode_codeword[i],temp);
		}
	}
}

float cal_EucDis(int reEncode_codeword[],int bicodeword_temp[],int symbol_temp[][2])
{
	int i,j,value;
	unsigned int mask=1;
	float d1,d2;
	float temp=0.0;

	//Convert the nonbinary codeword into binary
	for(i=0;i<n;i++) 
	{
		value=reEncode_codeword[i];
		mask=1;
		for(j=0;j<p;j++)
		{
			if( (value & mask)>0 )
				bicodeword_temp[p*i+j]=1;
			else
				bicodeword_temp[p*i+j]=0;
			mask=mask<<1;
		}
	}

	//convert the bianry into coordinate
	for(i=0;i<(n*p)/2;i++)
	{
		//QPSK
		symbol_temp[i][0]=2*bicodeword_temp[i*2+0]-1;		//coordinate x 0-->-1
		symbol_temp[i][1]=2*bicodeword_temp[i*2+1]-1;		//coordinate y 1-->1
	}

	//calculate the distance
	for(i=0;i<(n*p)/2;i++)
	{	//caution: here must use compulsory convert, int-->float
		d1=( rx_symbol[i][0] - (float)symbol_temp[i][0] ) * ( rx_symbol[i][0] - (float)symbol_temp[i][0] );  // coordinate x 
		d2=( rx_symbol[i][1] - (float)symbol_temp[i][1] ) * ( rx_symbol[i][1] - (float)symbol_temp[i][1] );  // coordinate y
		temp=temp + d1 + d2 ;
	}

	return temp;
}

int cal_mindistan_seqnum(float distance_temp[])
{
	int i,j;
	float temp1;

	//find out the min distance
	temp1=distance_temp[0];
	j=0;
	for(i=0;i<f_num;i++)
		if(temp1>distance_temp[i])
		{
			temp1=distance_temp[i];
			j=i;
		}

	return j;
}

