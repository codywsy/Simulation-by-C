#ifndef _H_NewInter_
#define _H_NewInter_

//NewInter()
#define New_interpoly_Zsize (lm+1)
#define New_interpoly_Ysize w
#define New_interpoly_Xsize (iterNum+1)

//Function
int MonoOrderConvert(int degz, int degy, int degx);
void NewInter(int g[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int interpoint[][n - eta]);

#endif