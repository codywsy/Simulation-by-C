#ifndef _H_main_
#define _H_main_

//main.h
#define w 4
#define weiz 39
#define iterNum 64	//when m=1, C is equal to n
#define lm 1
#define able_correct 12
#define pointNum 2
#define interval 0.25

//test_vec_construction
#define choose_num 2	//fixed

//tgorder()
#define tg_size 275	//more than k, tg_size represent the probably used pole basis num, (w+1)*w/2 + (w+1) * ( (w+interpoly_Ysize)-w+1 ) 
//mono_table()
#define weight_Zsize 6
#define weight_XYsize tg_size	//equals to the tg_size
#define mono_ordinSize 700	//choose the value with debugging, equals to max weigtdegree in weight[monoTable_Zsize][weight_XYsize]
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 67	//large than interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize tg_size //>index of term[max(deg_y)][w] in the pole basis + 1, nearly equals to tg_size
//interpolation()
#define init_polyNum (w*(lm+1))	//change with diff lm, the poly num of the init polyGroup
#define interpoly_Zsize (lm+1)	//maxValue of z is lm=1, so the Zsize should add 1 more.
#define interpoly_Ysize 52	//maxdeg of y is (w-1) + w*(n/(w+1)), so the Ysize should add 1 more.
#define interpoly_Xsize	(w+1)	//maxdeg of x is w, so the Xsize should add 1 more.
//factorization()
#define facpoly_Zsize interpoly_Zsize// same as the interpoly_Zsize
#define facpoly_Ysize interpoly_Ysize		// same as the interpoly_Ysize
#define facpoly_Xsize interpoly_Xsize	// same as the interpoly_Xsize
//rcs()
#define rcspoly_Ysize 67	//degY+1, degY = interpoly_Ysize + max(j_1) + w*( (max(i_1)+w)/(w+1) )
#define rcspoly_Xsize 11	//degX+1, degX = max(i_1) + w
#define faiMax_Ysize 8	//change with diff w and k, the max degY of probably used polebasis
#define faiMax_Xsize 4	//change with diff w and k, the max degX of probably used polebasis
//expoly()-->expanded polynomial
#define expoly_Ysize (lm*faiMax_Ysize+1)	
#define expoly_Xsize (lm*faiMax_Xsize+1)


//conditional compile
//#define _GS_Normal_

#ifndef _GS_Normal_
#define eta_min 1
#define eta_max 4
#define test_vec_num_max 16	//the num of test_vec is 2^eta
#define r_thread 0.5

#else
#define eta 0
#define test_vec_num 1 //the num of test_vec is 2^eta
#endif

//#define _Complexity_
#define _NoReductionCom_
#define _NoReductionUncom_
//#define _PolyCoeffNumUncom_
//#define _PolyCoeffNumFac_
//#define PrintfMonoTable
//#define _PrintRThread_
#define cheatingEncoding
//#define cheatingFac

#define _ResultComplexity_

#ifdef _ResultComplexity_
	#define FrameError 5000
#else
	#define cheatingDec
	#define FrameError 309
#endif


#define OpenFile fp=fopen("[1,4]New-LCC_Herm(64£¬34)thread=0.5.txt","a")
//#define OpenFile fp=fopen("LCC_Herm(64£¬34)_¦Ç=3.txt","a")


//#define _PrintGroupLod_
//#define _PrintLCCLod_
#define _PrintGammaValue_


#endif