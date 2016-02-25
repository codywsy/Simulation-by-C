	

	//*******demodulation debug**************
	int h,error_1=0;
	for(h=0;h<n*p;h++)
		if(bi_codeword[h]!=rx_bicodeword[h])
			error_1++;
	
	printf("\nbiterror=%d",error_1);

	//********************


	//************debug***************
	printf("\ncodeword=");
	for(h=0;h<n;h++)
		printf("\t%d",codeword[h]);

	printf("\nrx_codeword=");
	for(h=0;h<n;h++)
		printf("\t%d",rx_codeword[h]);
	//*************************


	//*****interpolation debug***********
	//predict the validate of factorization
	//calculate the codeword score
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
	//****************************************


	//debug: output message candidate poly
	printf("\n\n");
	for(h=0;h<f_num;h++)
	{
		printf("\nf[%d]=",h);
		for(i=0;i<k;i++)
			printf("\t%d",f[h][i]);
	}
