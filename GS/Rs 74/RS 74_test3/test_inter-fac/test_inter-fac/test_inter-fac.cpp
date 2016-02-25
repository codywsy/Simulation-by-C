/****************************************
RS(7,4)code   GS algorithm  m=2   
*****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/********change element with different code************
n,k,q,p			m-->iter_num,tm,lm,mono_size,fa_size
a[n]
*****************************************************/
#define n 7	//length of codeword
#define k 4 //length of message
#define q 8 //GF(q)
#define p 3 //GF(2^p) = GF(q)
#define maxcost 300 //Designed maximal interpolation cost
#define m 1 // zero of multiplicity
#define iter_num 7
#define tm 1
#define lm 1   //Desgined factorization output list size   (should be mended)
#define mono_size 20  // the size of inter_polynomial
#define fa_size 30	  // the number of layers of factorization'tree 
//******************


int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
void mono_table(void);
void encoder(void);
void modulation(void);
void channel(void);
void demodulation(void);
//******
void decoding(void);
int cal_poly_lod(int);	
int cal_minlod(void);	
int cal_delta(int,int,int);
void cal_Newinter_poly(int,int);
void cal_Newpoly_minlod(int,int);
void cal_Newpoly(int,int);
void refresh_poly_temp(int array[][mono_size+1]);
int fac(int);
void DFS(int);
int cal_Qy_zero(int);
int cal_rootlist(int,int);
void cal_update_poly(int,int,int);
void output_message_poly(int);
void choose_correct_dec_codeword(void);
void re_encoder(int,int arr1[]);
float cal_EucDis(int arr1[],int arr2[],int arr3[][2]);
int cal_mindistan_seqnum(float arr1[]);

//******

unsigned long int seq_num;	//number of input binary sequences
int SNR;
double BER, FER;
int mularray[]={1,2,4,3,6,7,5};
int root[]={0,1,2,3,4,5,6,7};	//n+1	
int a[7]={1,2,3,4,5,6,7};      //caution: this array must change with different sequence of infinite element

int bi_message[k*p];	//binary input message
int message[k];			//decimal input message
int codeword[n];		//encoded code word
int bi_codeword[n*p];   //encoded binary codeword
int dec_codeword[n*p];  //decoded binary codeword
float pi=3.141593;
float tx_symbol[n*p][2];	//modulated tranmitted signal
float rx_symbol[n*p][2];	//received destorted signal

int mono_order[40][200], channel_type, kv_signal;	//multiplexity matrix
float N0, sgm;

//*********
int rx_bicodeword[n*p];   //demodulated binary codeword
int rx_codeword[n];   //demodulated codeword to be decoded
int inter_point[n][2];  //interlolation points (pi,ri)
int inter_poly[(lm+1)][mono_size][mono_size];  //store the inter_poly of everygroup
int inter_poly_minlod[mono_size][mono_size];  // the poly with minimium used to factorization
int poly_lod[(lm+1)];			//store the lod of inter_poly of everygroup
int poly_lod_temptfg[(lm+1)];		// lod temp
int alpha_beta[iter_num/n][2];     //store alpha and beta
int delta[lm+1];					   //store the group of delta, totalling

int message_poly[fa_size][mono_size];   //the array of poly used in factorization
int time;     // factorization time;

//factorization's element
int last[fa_size],deg[fa_size];
int Q[fa_size][mono_size][mono_size];
int Coeff[mono_size];
int f[n][k];		//store the message polynomial after the factorization, note that may be overflow
int f_i=0;

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

	printf("\nEnter start SNR: ");
	scanf("%d", &start);
	printf("\nEnter finish SNR: ");
	scanf("%d", &finish);
	printf("\nPlease input the number of k*p-bit sequence: ");
	scanf("%d", &seq_num);

		//*******debug*********
		printf("\nstart=%d	finish=%d	seq_num=%d	\n\n",start,finish,seq_num);
		//****************** 
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

void demodulation(void)
{
	int i;
	float d1,d2;
	for(i=0;i<n*p;i++)	
	{
		d1=(rx_symbol[i][0]-1)*(rx_symbol[i][0]-1)+rx_symbol[i][1]*rx_symbol[i][1];   // a*a= the square of a 
		d2=(rx_symbol[i][0]+1)*(rx_symbol[i][0]+1)+rx_symbol[i][1]*rx_symbol[i][1];
		if(d1<d2)
			rx_bicodeword[i]=0;
		else
			rx_bicodeword[i]=1;

	}

	//Convert rx_bicodeword to nonbinary's rx_codeword
	for(i=0;i<n;i++)
		rx_codeword[i] = rx_bicodeword[p*i]+2*rx_bicodeword[p*i+1]+4*rx_bicodeword[p*i+2];	

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

	
//*********
void decoding(void)
{
	int i,j,h,temp,u;
	int t,tk;   			// tk is the num of iteration, t is the number of inter_point
	int minlod_seq,minlod=0;	//minlod_seq is the j' with delta[j']!=0

	// interpolation points
	for(i=0;i<n;i++)					
	{
		inter_point[i][0]=a[i]; 		// the first word is pi--->xi
		inter_point[i][1]=rx_codeword[i];	// the second word is ri--->yi
	}

    //initialize group_inter_poly
    for(j=0;j<(lm+1);j++)
   	{
   		for(h=0;h<mono_size;h++)
   		{
   			for(i=0;i<mono_size;i++)
   			{
   				inter_poly[j][h][i]=0;
   			}
   		}
   	}


     //initialize alpha_beta points
    for(i=0,j=0;j<m;j++)     
    {
    	for(h=0;h<(m-j);h++)
    	{
    		if(i<(iter_num/n))
  			{
  				alpha_beta[i][0]=j;
   				alpha_beta[i][1]=h;
   				++i;
   			}
   		}
   	}

    //initialize Group_0
    for(i=0;i<(lm+1);i++)
    	inter_poly[i][i][0]=1;

	// work out the poly_lod of inter_poly[i] of Group[0]
	for(i=0;i<(lm+1);i++)      
		poly_lod[i]=cal_poly_lod(i);  
    
    //**********begin to interpolation**********************
	//iteration of modifying the inter_polynomial   -------> optimize to a function
    for(tk=0;tk<n;tk++)
    {

    	for(t=0;t<(iter_num/n);t++)
    	{
    		//calculate delta_i
    		for(i=0;i<(lm+1);i++)
    		{
    			delta[i]=cal_delta(i,tk,t);
    		}

       		//choose the delta_i with delta != 0 and j' = minlod_seq
  			minlod_seq=cal_minlod();

		 	//calculate the new inter polynomial -------> update
		 	cal_Newinter_poly(minlod_seq,tk);

		 	// work out the new poly_lod of inter_poly[i]
			for(i=0;i<(lm+1);i++)      
		 		poly_lod[i]=cal_poly_lod(i);  

 		 	//find the interpolation polynomial of Group_C with minimium lod-----------> inter_poly[minlod]
 		 	temp=poly_lod[0];
		 	for(i=0;i<(lm+1);i++)
		 		if( temp > poly_lod[i] )
		  		{
		  			temp=poly_lod[i];
		  			minlod=i;
		  		}
		}
    }

	// find out the inter_poly_minlod used to factorization
    for( j=0;j<mono_size;j++)
		for( h=0;h<mono_size;h++)
			inter_poly_minlod[j][h]=inter_poly[minlod][j][h];
	//********************************************************

	//***************begin factorization************************

	//initialize Q, and Q is used to 
	for(i=0;i<fa_size;i++)
		for(j=0;j<mono_size;j++)
			for(h=0;h<mono_size;h++)
				Q[i][h][j]=0;

	//initialize last[] and deg[]
	for(i=0;i<fa_size;i++)
	{
		last[i]=deg[i]=0;
	}

	//initialize Coeff[]
	for(i=0;i<mono_size;i++)
		Coeff[i]=0;

	//intialize the result polynomial of factorization
	for(i=0;i<n;i++)
		for(j=0;j<k;j++)
			f[i][j]=0;


	//initialize Q0
	for(h=0;h<mono_size;h++)
		for(j=0;j<mono_size;j++)
			Q[0][h][j]=inter_poly_minlod[h][j];

	f_i=0;
	last[0]=-1;
	deg[0]=-1;
	time=1;  u=0;

	DFS(u);

	// choose the final dec_codeword
	choose_correct_dec_codeword();  

	//at last output the binary dec_codeword[i]    

}	



int cal_poly_lod(int i)   //input one inter_polynomial, output this poly's leading order
{
	int j,h;
	int max_ord=0;

	if( poly_lod[i] <= iter_num )  //modification
		for(j=0;j<mono_size;j++)
			{
				for(h=0;h<mono_size;h++)
				{
					if( inter_poly[i][j][h] != 0 )
					{
						if(	max_ord < mono_order[j][h])
						{
							max_ord=mono_order[j][h];
						}
					}
				}
			}
	
	return max_ord;
}

int cal_minlod(void)				//calculate the minlod of this group
{
	int i,temp,flag=1,minlod_seqnum=0;

	//calculate the lod of poly and store it in poly_lod[i]
	for(i=0;i<(lm+1);i++)
  		poly_lod[i]=cal_poly_lod(i);

    //find out the minlod_seqnum of polynomial Group
  	for(i=0;i<(lm+1);i++)
  	  if( poly_lod[i] <= iter_num ) //modification
  		if( delta[i]!=0 )
  		{
  			if(flag==1)
  			{
  				temp = poly_lod[i];
				minlod_seqnum=i;
  				flag=0;
  			}

  			if( temp > poly_lod[i] )
  			{
  				temp = poly_lod[i];
  				minlod_seqnum=i;
  			}
  		}

  	return minlod_seqnum;

}



int cal_delta(int i,int tk,int t)    //calculate the delta of inter_polynomial, c is inter_poly[i] and t is tk
{
	int delta=0,delta_temp,h,j,flag;
	int a,b,alpha,beta,x,y;

	alpha = alpha_beta[t][0];    //alpha
	beta = alpha_beta[t][1];		//beta
	x=inter_point[tk][0];			//xi=pi
	y=inter_point[tk][1];			//yi=ri

	if( poly_lod[i] <= iter_num )  //modification
		for(j=0;j<mono_size; j++)
		{
			for(h=0;h<mono_size;h++)
			{
				if( (inter_poly[i][j][h] != 0) && ( h >= alpha) && ( j >= beta) )
				{
					a=h;
					b=j;
					flag= ( fac(a) / ( fac(a-alpha)*fac(alpha) ) ) * ( fac(b) / ( fac(b-beta)*fac(beta) ) );
					// if flag is even number , delta_temp is equal to 0
					// if flag is odd number, delta_temp can be calculated to the delta
					if( (flag%2) !=0 )
					{
						delta_temp=mul(inter_poly[i][j][h],power(x,(a-alpha)));
						delta_temp=mul(delta_temp,power(y,(b-beta)));
						delta= add(delta,delta_temp);
					}
				}
			}	
		}

	return delta;
}



void cal_Newinter_poly(int minlod_seq,int tk)
{
	int i;

	// start to update
	// (delta != 0)-----> update  
    for( i=0 ; i<(lm+1) ; i++ )
      if( poly_lod[i] <= iter_num )  //modification
        if( (delta[i]!=0) && (i!=minlod_seq) )
        	cal_Newpoly(i,minlod_seq);	

    // at last update the inter_poly[minlod_seq]
    if( poly_lod[i] <= iter_num )  //modification
		cal_Newpoly_minlod(minlod_seq,tk);
 
}



void cal_Newpoly_minlod(int minlod_seq,int tk)					
{
	int j,h,temp;
	int poly_temp1[mono_size+1][mono_size+1];
	int poly_temp2[mono_size+1][mono_size+1];

	refresh_poly_temp(poly_temp1);
	refresh_poly_temp(poly_temp2);

	
	// delta[i]*x*inter_poly[i]
	for( j=0 ; j<mono_size ; j++ )
   		for( h=0 ; h<mono_size ; h++ )
       		if( inter_poly[minlod_seq][j][h] != 0)
       			poly_temp1[j][h+1]=mul(inter_poly[minlod_seq][j][h],delta[minlod_seq]);

    // delta[i] * xi * inter_poly[i]
  	temp=mul(delta[minlod_seq],inter_point[tk][0]);    //delta[i] * xi
  
	for( j=0 ; j<mono_size ; j++ )
   		for( h=0 ; h<mono_size ; h++ )
       		if( inter_poly[minlod_seq][j][h] != 0)
       			poly_temp2[j][h] = mul(inter_poly[minlod_seq][j][h],temp);

    // inter_poly + inter_poly_temp
    for( j=0;j<mono_size;j++)
		for( h=0;h<mono_size;h++)
			if( (poly_temp1[j][h]==0) && (poly_temp2[j][h]==0) )  
				inter_poly[minlod_seq][j][h]=0;
			else
				inter_poly[minlod_seq][j][h]=add(poly_temp1[j][h],poly_temp2[j][h]);
}

void cal_Newpoly(int i,int minlod_seq)
{
	int h,j;
	int poly_temp1[mono_size+1][mono_size+1];
	int poly_temp2[mono_size+1][mono_size+1];

	refresh_poly_temp(poly_temp1);
	refresh_poly_temp(poly_temp2);

	// delta[j'] * inter_poly[i]
	for( j=0;j<mono_size;j++)
		for( h=0;h<mono_size;h++)  
			if( inter_poly[i][j][h] != 0 )
				poly_temp1[j][h]=mul(inter_poly[i][j][h],delta[minlod_seq]);
	// delta[i] * Q'(inter_poly_temp)
	for( j=0;j<mono_size;j++)
		for( h=0;h<mono_size;h++)  	
			if( inter_poly[minlod_seq][j][h] != 0 )
				poly_temp2[j][h]=mul(inter_poly[minlod_seq][j][h],delta[i]);
	//  new inter_poly[i+1] = new inter_poly[i] + new Q'
	for( j=0;j<mono_size;j++)
		for( h=0;h<mono_size;h++)  
			inter_poly[i][j][h]=add(poly_temp1[j][h],poly_temp2[j][h]);	


}


void refresh_poly_temp(int array[][mono_size+1])     //clean the poly_temp
{
	int h,j;

	for(j=0;j<(mono_size+1);j++)
   	{
   		for(h=0;h<(mono_size+1);h++)
   		{
   			array[j][h]=0;
   		}
   	}
}

int fac(int a)		//calculate the factorial of a
{
	int i,result=1;

	for(i=1;i<=a;i++)
		result=result*i;

	return result;
}

//************factorization part function***************
void DFS(int u)
{
	int v,i,h,temp,j=0;
	int rootlist[n];

	//initialize rootlist
	for(i=0;i<n;i++)
		rootlist[i]=0;

	if( cal_Qy_zero(u) )
			output_message_poly(u);     
	else if( deg[u] < (k-1) )		// the "D" of thesis represent k-1
	{
		//try n field element to calculate the y-root of Q(0,y)=0
		for(h=(n-1);h>=0;h--)		//field element from large to small
		{
			temp=cal_rootlist(u,h);
			if( temp == 0 )
			{	
				rootlist[j]=a[h];
				j++;
			}
		}

		for(i=0;i<n;i++)
			if( rootlist[i]!=0 )
			{
				v=time;
				time=time+1;
				last[v]=u;
				deg[v]=deg[u]+1;
				Coeff[v]=rootlist[i];
				cal_update_poly(v,u,rootlist[i]);  
				DFS(v);
			}

	}	

}

int cal_Qy_zero(int u)				//calculate the Q[u](x,0) ?= 0
{
	int i,flag=1;

	for(i=0;i<mono_size;i++)
		if( Q[u][0][i] != 0 )		//if one of coeff !=0, flag=0;
		{
			flag=0;
			break;
		} 

	return flag;
}

int cal_rootlist(int u,int h)
{
	int j,temp;

	//traverse all the field element a[i]
	temp=0;
	for(j=0;j<mono_size;j++)
		if( Q[u][j][0] != 0)
		{
			temp=add( temp , mul( Q[u][j][0] , power(a[h],j) ) );
		}
	return temp;
}

void output_message_poly(int u)
{
	int temp,deg_temp;		//use temp=last[temp] to iteration the last[u]

	temp=u;
	deg_temp=deg[u];
	while( deg_temp>(-1) )
	{
		f[f_i][deg_temp]=Coeff[temp];
		temp=last[temp];
		deg_temp=deg[temp];
	}
	f_i=f_i+1;       //f_i is the global variate
}


void cal_update_poly(int v,int u,int a)     // notification   update_temp[mono_size*2][mono_size*2]
{
	int i,j,g,h=0;
	int flag=0;
	int Q_update[mono_size][mono_size];
	int temp1[mono_size][mono_size],temp2[mono_size][mono_size];
	int update_temp1[mono_size*2][mono_size*2],update_temp2[mono_size*2][mono_size*2];

	//clean Q_update[mono_size][mono_size], and temp1¡¢temp2
	for(j=0;j<mono_size;j++)
		for(i=0;i<mono_size;i++)	
			Q_update[j][i]=temp1[j][i]=temp2[j][i]=0;
		
	//begin
	//Q_update = a+xy;
	Q_update[0][0]=a;
	Q_update[1][1]=1;

	//calculate the Q(x,xy+a)
	for (j=0; j < mono_size ; j++)    // when j=0, the poly does not change
	{	
		// clean up_date[mono_size*2][mono_size*2]
			for(h=0;h<(mono_size*2);h++)
				for(g=0;g<(mono_size*2);g++)	
					update_temp1[h][g]=update_temp2[h][g]=0;

		if( j>=2 )	//calculate (xy+a)^j, and store it into the Q_update 		
		{	
			//update_temp = Q_update * xy
			for(h=0;h<mono_size;h++)
				for(g=0;g<mono_size;g++)
					if( Q_update[h][g] != 0 )
						update_temp1[h+1][g+1]=Q_update[h][g];
			//temp1 = Q_update * a
			for(h=0;h<mono_size;h++)
				for(g=0;g<mono_size;g++)
						temp1[h][g] = mul(Q_update[h][g],a);	
			// Q_update = temp1 + update_temp = (Q_update*a) + (Q_update*xy)
			for(h=0;h<mono_size;h++)
				for(g=0;g<mono_size;g++)
					if( (temp1[h][g]!=0)||(update_temp1[h][g]!=0)||(Q_update[h][g]!=0) )
						Q_update[h][g] = add( temp1[h][g] , update_temp1[h][g] );

		}
		//caution Q_update[mono_size][mono_size] just enough, but it will overflow

		//begin to y=(xy+a)^j and calculate h(x)*y
		if(j==0)
		{
			for(i=0;i<mono_size;i++)
				if( Q[u][j][i] != 0)
					Q[v][j][i]=Q[u][j][i];
		}
		else if( j>=1 )
		{		
			for(i=0;i<mono_size;i++)	//  coeff *x^i * (xy+a)^j
				{
					if( Q[u][j][i]!=0 )    
					{
						for(h=0;h<mono_size;h++)     
							for(g=0;g<mono_size;g++)
								if( Q_update[h][g]!=0 )
									update_temp2[h][g+i]=mul( Q[u][j][i],Q_update[h][g] );
					}

					for(h=0;h<mono_size;h++)     
						for(g=0;g<mono_size;g++)
							if( update_temp2[h][g]!=0 )
								Q[v][h][g]=add( Q[v][h][g],update_temp2[h][g] );

					for(h=0;h<(mono_size*2);h++)
						for(g=0;g<(mono_size*2);g++)	
							if(update_temp2[h][g]!=0)
								update_temp2[h][g]=0;
				}
		}

	}
	//*******************************

	//calculate <<Q(x,y)>> = Q(x,y)/xm
		for(i=1;i<mono_size;i++)
		{
			for(j=0;j<mono_size;j++)
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
		{
			for(j=0;j<mono_size;j++)
				for(i=0;i<(mono_size-h);i++)
					if( Q[v][j][i+h] != 0 )
					{
						Q[v][j][i]=Q[v][j][i+h];
						Q[v][j][i+h]=0;
					}
		}
	//**************************

}

//********find the correct f(x)*****************
void choose_correct_dec_codeword(void)
{
	int i,min_distance_seqnum;
	int reEncode_codeword[n];   //store the codeword of reencoding f(x)
	int bicodeword_temp[n*p];   // convert nonbinary reEncode_codeword into binary bicodeword_temp
	int symbol_temp[n*p][2];	//convert binary bicodeword_temp into coordinate symbol_temp
	float distance_temp[n];		//store the Euclidean distance between symbol_temp(int) and rx_symbol(float)

	//clean the distance_temp[n]
	for(i=0;i<n;i++)
		distance_temp[i]=0;    

	for(i=0;i<f_i;i++)
	{
		re_encoder(i,reEncode_codeword);
		
		distance_temp[i]=cal_EucDis( reEncode_codeword , bicodeword_temp , symbol_temp );
		
	}
	//find out the min_distance_seqnum
	min_distance_seqnum=cal_mindistan_seqnum(distance_temp);   

	//after knowing the min_distance f(x)'i =  min_distance_seqnum, reencoder this fi(x)
	//get the final binary dec_codeword
	re_encoder(min_distance_seqnum,reEncode_codeword);
	for(i=0;i<n*p;i++)
	{
	
		dec_codeword[i]=bicodeword_temp[i];
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

	//Convert the codeword into binary
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
	for(i=0;i<n*p;i++)
	{
		symbol_temp[i][0]=-1*(2*bicodeword_temp[i]-1);       //coordinate x
		symbol_temp[i][1]=0; 								 //coordinate y
	}

	//calculate the distance
	for(i=0;i<n*p;i++)
	{	//caution: here must use compulsory convert, int-->float
		d1=( rx_symbol[i][0] - (float)symbol_temp[i][0] ) + ( rx_symbol[i][0] - (float)symbol_temp[i][0] );  // coordinate x 
		d2=( rx_symbol[i][1] - (float)symbol_temp[i][1] ) + ( rx_symbol[i][1] - (float)symbol_temp[i][1] );  // coordinate y
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
	for(i=0;i<f_i;i++)
		if(temp1>distance_temp[i])
		{
			temp1=distance_temp[i];
			j=i;
		}

	return j;
}

//*********

	
	
	
