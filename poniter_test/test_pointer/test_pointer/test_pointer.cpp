/* commandLine-2.c */
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>


int main(int argc, char *argv[])
{
	if(argc == 2)
	{
		switch(*(argv[1]+1))
			{
				case 'd':d
					system("dir/w c:\\");
					break;
				case 't':
					system("type c:\\test_pointer.cpp");
					break;
				default:
					printf("Using commandLine -d or -t");
			}
	}
 	else
 		printf("Using commandLine -d or -t");
 	printf("\n");
 	getch();
 	return 0;
 }
 