#ifndef _H_NewInter_
#define _H_NewInter_

//NewInter()
#define New_interpoly_Zsize (lm+1)
#define New_interpoly_Ysize w
#define New_interpoly_Xsize (iterNum+1)

//Function
int MonoOrderConvert(int degz, int degy, int degx);
//void NewInter(int g[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int interpoint[][n - eta]);
//void NewInter();
void NewInter_1();
void NewInter_2();
void NewInter_3();
void InterOnce(int poly[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int interpoint[3]);
void InterOnce_2(int poly[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int interpoint[3], int poly_min_order[eta], int location, int poly_min_order_all[n], int location_all);
void ConvertX2Y(int polySour[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int polyDes[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int PolyNum, int omga);
int ChooseTheMinPoly(int polySour[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int PolyNum);

#endif