
void cal_update_poly(int v,int u,int a)     // notification   update_temp[mono_size*2][mono_size*2]
{
	int i,j,g,h;
	int flag=0;

	clean_Q_temp();    // clean Q_update and update_temp
	Q_update[0][0]=a;
	Q_updata[1][1]=1;

	for (j = 1; j < mono_size ; j++)    // when j=0, the poly does not change
	{	
		//update_temp = Q_update * xy
		for(h=0;h<mono_size;h++)
			for(g=0;g<mono_size;g++)
				if( Q_update[h][g] != 0 )
					update_temp[h+1][g+1]=Q_update[h][g];
		//Q_update = Q_update * a
		for(h=0;h<mono_size;h++)
			for(g=0;g<mono_size;g++)
				if( Q_update[h][g] != 0 )
					Q_update[h][g] = mul(Q_update,a);				
		// Q_update = Q_update + update_temp
		for(h=0;h<mono_size;h++)
			for(g=0;g<mono_size;g++)
				if( (Q_update[h][g]!=0)||(update_temp[h][g]!=0) )
				{
					Q_update[h][g] = add( Q_update[h][g] , update_temp[h][g] );
				//	update_temp[h][i]=Q_update[h][i];
				}

		for(i=0;i<mono_size;i++)
			if( Q[u][i][j]!=0 )
			{
				for(h=0;h<mono_size;h++)
					for(g=0;g<mono_size;g++)
					{
						update_temp[h+i][g]=Q_update[h][g];	
						update_temp[h+i][g]=mul( update_temp[h+i][g] , Q[u][i][j] );
					}

				for(h=0;h<mono_size;h++)
					for(g=0;g<mono_size;g++)
					{
						if( (h==i)&&(g==j) )
							Q[v][h][g]=update_temp[h][g];
						else
							Q[v][h][g]=add(Q[u][h][g],update_temp[h][g]);
					}			
			}		
	}

	//calculate <<Q(x,y)>> = Q(x,y)/xm
	for(i=0;i<mono_size;i++)
	{
		for(j=0;j<mono_size;j++)
			if( Q[v][i][j] != 0)
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

	for(i=0;i<(mono_size-h);i++)
		for(j=0;j<mono_size;j++)
			if( (Q[v][i+h][j] != 0) && (h>0) )
				Q[v][i][j]=Q[v][i+h][j];

}