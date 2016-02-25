void New_zerobasis(int pointIndex, int LM_zb[][2])
{

	int i, j, u, v;
	int temp, flag, xi, yi, index_flag, temp2;
	int zb_temp1[zb_Xsize], zb_temp2[zb_Ysize][zb_Xsize+1];
	int zb_temp_a[w+1][zb_Xsize];
	int poly_temp1[zb_Yszie+w][zb_Xsize+w], poly_temp2[zb_Ysize+1][zb_Xsize];

	int num_deg_y, num_deg_x;

	//Initialisation
	for(i=0;i<zb_alpha;i++)
		for(u=0;u<zb_Ysize;u++)
			for(v=0;v<zb_Xsize;v++)
			{
				zb[i][u][v] = 0;
			}

	memset(zb_temp_a,0,sizeof(zb_temp_a));		

	xi = point[pointIndex][0];
	yi = point[pointIndex][1];
	flag = 0;

	num_deg_y = zb_alpha/(w+1);
	num_deg_x = zb_alpha%(w+1);

	//calculate Zb

	//1. pre-calculate zb_temp_a
	//zb_temp_a[0]
	zb_temp_a[0][0]=1;

	//zb_temp_a[1]
	zb_temp_a[1][0] = xi;
	zb_temp_a[1][1] = 1;

	//zb_temp_z[others]
	for(i=2;i<=w;i++)
	{
		//x*zb[i-1]
		for(v=0;v<zb_Xsize;v++)
			if(zb_temp_a[i-1][v]!=0)
			{
				zb_temp_a[i][v+1] = zb_temp_a[i-1][v];
			}

		//xi*zb[i-1]
		for(v=0;v<zb_Xsize;v++)
			if(zb_temp_a[i-1][v]!=0)
			{
				zb_temp_a[i][v] = add( zb_temp_a[i][v],zb_temp_a[i-1][v] );
			}
	}

	//2. calculate zb
	for(i=0;i<num_deg_y;i++)
	{
		if(i==0)
		{
			for(j=0;j<(w+1);j++)
			{
				for(v=0;v<zb_Xsize;v++)
					if(zb_temp_a[j][v]!=0)
					{
						zb[j][0][v] = zb_temp_a[j][v]; 
					}
			}
		}
		else if( i>0 && i<num_deg_y )
		{

			for(j=0;j<(w+1);j++)
			{
				temp_x = power( xi,w );	//xi^w, the coefficient of x
				temp_norm = add( y,mul(temp_x,xi) );	//xi^(w+1) + yi
				index_i_1 = (i-1)*(w+1);
				index_i = i*(w+1);

				//temp_norm*zb[(i-1)*(w+1)]
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb[index_i_1][u][v]!=0)
						{
							zb[index_i][u][v] = mul( temp_norm,zb[index_i_1][u][v] ); 
						}

				//y*zb[(i-1)*(w+1)]
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if()



			}


		}


	}

}