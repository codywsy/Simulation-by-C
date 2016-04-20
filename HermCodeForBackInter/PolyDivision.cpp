/* 
 * POJ_2527.cpp 
 * 
 *  Created on: 2013��10��26�� 
 *      Author: Administrator 
 */ 
   
#include <iostream>  
#include <cstdio>  
   
using namespace std;  
   
const int maxn = 10010;  
int main(){  
    int n,k;  
   
    int val[maxn];  
    while(scanf("%d%d",&n,&k)!=EOF,n!=-1 || k !=-1){  
        int i;  
        for(i = 0 ; i <= n ; ++i){  
            scanf("%d",&val[i]);  
        }  
   
        //���г�������  
        for(i = n ; i >= k ; --i){  
            if(val[i] == 0){  
                continue;  
            }  
   
            val[i-k] = val[i-k] - val[i];  
            val[i] = 0;  
        }  
   
        //�������鳤��,����λ��0�������  
        int t = n;  
        while(val[t] == 0 && t > 0){  
            --t;  
        }  
   
        for(i = 0 ; i < t ; ++i){  
            printf("%d ",val[i]);  
        }  
        printf("%d\n",val[t]);  
   
   
    }  
   
    return 0;  
}  
