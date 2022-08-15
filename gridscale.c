#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gridtools.c"

int main(int argc, char *argv[])
{
    void CheckInput(int argc, char *argv[]);
    void multgrid2(int numfiles, char *infilenames[]);
    
    CheckInput(argc,argv);
    multgrid2(argc-1,argv);
}

void multgrid2(int numfiles, char *infilenames[])
{
    void readgrid(gridtype *grid, char *filename);
    void IndexForm2BlockForm(gridtype *grid);
    void writegrid(gridtype *grid, char *filename);
    void printstats(gridtype *grid);
    void gridstats(gridtype *grid);
    gridtype grid1;
    int ii,i,j,k;
    float mval;
    
    mval = atof(infilenames[1]);                                    //takes first input and assigns it to the float multiplier
    
    for(ii=2; ii<=numfiles; ii++)                                   //runs program as many times as there are files
    {
        if(mval>0)                                                  //prevents multiplication by 0 which would destroy our file
        {                                                           //and prevents multiplication by negative values
            readgrid(&grid1,infilenames[ii]);
            IndexForm2BlockForm(&grid1);
            fprintf(stdout,"\nOld version of %s.\n",infilenames[ii]);
            gridstats(&grid1);
            printstats(&grid1);
            for (i=1; i<=grid1.imax; i++)                           //runs through every point in file and mulitplies by the multiplier
            {
                for (j=1; j<=grid1.jmax; j++)
                {
                    for (k=1; k<=grid1.kmax; k++)
                    {
                        grid1.data[i][j][k] = grid1.data[i][j][k] * mval;
                    }
                }
            }
            sprintf(grid1.filename,"%s.%gxscale",grid1.filename,mval);  //determines new file name
            writegrid(&grid1,grid1.filename);                           //creates the new file
            fprintf(stdout,"\nNew version of %s multiplied by %f.\n",infilenames[ii],mval);
            gridstats(&grid1);
            printstats(&grid1);
        }
        else                                                        //terminates program if scalar factor is 0
        {
            fprintf(stdout,"\n\nCareful! Multiplying your grid by 0 will\n");
            fprintf(stdout,"permanently set all cells 0, rendering the\n");
            fprintf(stdout,"file useless. Terminating program.\n\n");
            exit(0);
        }
    }
}

void CheckInput(int argc, char *argv[])
{
    void CheckInputGrids(char *argv[], int first, int last, int *error, char *finalerrmsg);
    void usage(char *argv[],char *errmsg);
    void beginerror(char *arg, char *errmsg);
    char errmsg[100];
    char finalerrmsg[500];
    int error, ii;

    error = 0;                                                         // every error ecountered increases this value by 1
    beginerror(argv[0],finalerrmsg);
    
    if (argc < 3)                                                      // minimum of 2 argument
    {
        sprintf(errmsg,"    - %i input parameters passed in: 2 or more expected!\n",argc-1);
        strcat(finalerrmsg,errmsg);
        usage(argv,finalerrmsg);
    }
    
    CheckInputGrids(argv,2,argc-1,&error,finalerrmsg);
    
    for(ii=2; ii<=argc; ii++)                                           //if multiple files are put in, this loops through all
    {                                                                   //possible file names to ensure file doesn't already exist
        sprintf(argv[argc+ii-1],"%s.%sxscale",argv[ii],argv[1]);
        CheckOutputFile(argv[argc+ii-1],&error,finalerrmsg);
    }
    
   
    if (error >= 1)                                                     //if any errors occur, displays errors and terminates
    {
        usage(argv,finalerrmsg);
    }
}

void usage(char *argv[], char *finalerrmsg)
{
    fprintf(stdout,"%s",finalerrmsg);
    fprintf(stdout,"Program %s multiplies a file that is \n",argv[0]);
    fprintf(stdout,"a binary 3D grids of densities by a factor\n");
    fprintf(stdout,"input by user.\n\n");
    fprintf(stdout,"    Example 1:   %s 5 Ell30\n\n",argv[0]);
    fprintf(stdout,"This example would mulitply all values in\n");
    fprintf(stdout,"the file Ell30 by 5.\n\n");
    fprintf(stdout,"    Example 2:   %s 2.5 Ell20 Ell30\n\n",argv[0]);
    fprintf(stdout,"This example would mulitply all values in\n");
    fprintf(stdout,"the files Ell20 & Ell30 by 2.5.\n\n");
    fflush(stdout);
    exit(1);
}
