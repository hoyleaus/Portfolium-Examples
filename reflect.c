#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gridtools.c"
#include <unistd.h>

int main(int argc, char *argv[])
{
    void CheckInput(int argc, char *argv[], int nn);
    void reflectgrid(int numfiles, char *infilenames[]);
    
    int nn=0;
    
    if(argc > 1)
    {
        float VALUE = atof(argv[1]);
        
        for(nn=2;nn<=argc;nn++)                                 // quick for loop to look if the .ref(n) file already exists
        {
            if(VALUE == 1) { sprintf(argv[argc+nn-1],"%s.refx",argv[nn]); }
            if(VALUE == 2) { sprintf(argv[argc+nn-1],"%s.refy",argv[nn]); }
            if(VALUE == 3) { sprintf(argv[argc+nn-1],"%s.refz",argv[nn]); }
            CheckInput(argc,argv, nn);                          // checks for proper input, which includes look for an already existing file
        }
    }
    else
    {
        CheckInput(argc,argv, nn);
    }
    reflectgrid(argc-1,argv);                                   //we are going to convert agrc into "numfiles". We subtract one to get rid
}                                                               //of the executable name

void reflectgrid(int numfiles, char *infilenames[])             // takes user input - infilenames[1] - and decides which branch of code to run
{
    void readgrid(gridtype *grid, char *filename);
    void IndexForm2BlockForm(gridtype *grid);
    void NewBlockGrid(gridtype *grid, short int imax, short int jmax, short int kmax, double xcell, double ycell, double zcell);
    void reflectX(gridtype *grid1,gridtype *grid2);
    void reflectY(gridtype *grid1,gridtype *grid2);
    void reflectZ(gridtype *grid1,gridtype *grid2);
    gridtype grid1,grid2;
    void writegrid(gridtype *grid, char *filename);
    void printstats(gridtype *grid);
    void gridstats(gridtype *grid);
    
    int ii;
    
    for (ii=2;ii<=numfiles;ii++)                                //infilenames[1] is the reflection parameter so we start at 2
    {
        readgrid(&grid1,infilenames[ii]);                       //reads file to give us grid 1 parameters (imax jmax etc..)
        IndexForm2BlockForm(&grid1);                            //converts forms if needed
        NewBlockGrid(&grid2, grid1.imax, grid1.jmax, grid1.kmax, grid1.xcell,grid1.ycell,grid1.zcell);
                                                                // ^^^^^uses info from readgrid to create new grid with the same dimensions
        float VAL = atof(infilenames[1]);                       // takes reflection parameter and converts it to a float variable "VAL"
        
        if(VAL == 1) { sprintf(infilenames[numfiles+2],"%s.refx",infilenames[ii]); } //logic to create new file name
        if(VAL == 2) { sprintf(infilenames[numfiles+2],"%s.refy",infilenames[ii]); }
        if(VAL == 3) { sprintf(infilenames[numfiles+2],"%s.refz",infilenames[ii]); }
       
        sprintf(grid2.filename,"%s",infilenames[numfiles+2]);
        
        if(VAL == 1) { reflectX(&grid1,&grid2); }               //logic to determine which reflectN branch to run
        if(VAL == 2) { reflectY(&grid1,&grid2); }
        if(VAL == 3) { reflectZ(&grid1,&grid2); }

        writegrid(&grid2,grid2.filename);                       //Writes a new file with .refN extension
    }
}

void reflectZ(gridtype *grid1,gridtype *grid2)                  //reflects along z axis (option 3)
{                                                               //reflectZ, reflectY and reflectX have the exact same logic
    void printstats(gridtype *grid);                            //just different variables
    void gridstats(gridtype *grid);
    
    int i,j,k,n;
    
    n = grid1->kmax;
    
    for (k=1; k<=grid1->kmax; k++, n--)
    {
        for (j=1; j<=grid1->jmax; j++)
        {
            for (i=1; i<=grid1->imax; i++)
            {
                grid2->data[i][j][n] = grid1->data[i][j][k];    //loops to replace on cell's value with its "opposite" cell value
            }
        }
    }
    gridstats(grid2);
    printstats(grid2);
}

void reflectY(gridtype *grid1,gridtype *grid2)                  //reflects along y axis (option 2)
{
    void printstats(gridtype *grid);
    void gridstats(gridtype *grid);
    
    int i,j,k,n;
    
    n = grid1->jmax;
    
    for (j=1; j<=grid1->jmax; j++, n--)
    {
        for (i=1; i<=grid1->imax; i++)
        {
            for (k=1; k<=grid1->kmax; k++)
            {
                grid2->data[i][n][k] = grid1->data[i][j][k];
            }
        }
    }
    gridstats(grid2);
    printstats(grid2);
}

void reflectX(gridtype *grid1,gridtype *grid2)                  //reflects along x axis (option 1)
{
    void printstats(gridtype *grid);
    void gridstats(gridtype *grid);
    
    int i,j,k,n;
    
    n = grid1->imax;
    
    for (i=1; i<=grid1->imax; i++, n--)
    {
        for (j=1; j<=grid1->jmax; j++)
        {
            for (k=1; k<=grid1->kmax; k++)
            {
                grid2->data[n][j][k] = grid1->data[i][j][k];
            }
        }
    }
    gridstats(grid2);
    printstats(grid2);
}

void CheckInput(int argc, char *argv[], int nn)                 //logic to determine error count
{
    void CheckInputGrids(char *argv[], int first, int last, int *error, char *finalerrmsg);
    void CheckOutputFile(char *filename, int *error,char *finalerrmsg);
    void CheckOutputFile(char *filename, int *error,char *finalerrmsg);
    void usage(char *argv[],char *errmsg);
    void beginerror(char *arg, char *errmsg);
    char errmsg[100];
    char finalerrmsg[500];
    int error;

    error = 0;                                                  //for every error encountered, this value increases by one
    beginerror(argv[0],finalerrmsg);
    
    if (argc < 3) // minimum of 2 argument
    {
        sprintf(errmsg,"    - %i input parameters passed in: 2 or more expected!\n",argc-1);
        strcat(finalerrmsg,errmsg);
        usage(argv,finalerrmsg);
    }
    
    CheckInputGrids(argv,2,argc-1,&error,finalerrmsg);
    CheckOutputFile(argv[argc+nn-1],&error,finalerrmsg);
    
    if (error >= 1)                                             //if error is greater than one, we will get an error message and terminate
    {
        usage(argv,finalerrmsg);
    }
}
void usage(char *argv[], char *finalerrmsg)                     //Brief explanation of how to operate the program. 
{
    fprintf(stdout,"%s",finalerrmsg);
    fprintf(stdout,"Program %s takes reflection parameter first:\n",argv[0]);
    fprintf(stdout,"3 = reflection along z axis\n");
    fprintf(stdout,"2 = reflection along y axis\n");
    fprintf(stdout,"1 = reflection along x axis\n");
    fprintf(stdout,"Program %s then also takes at least one file\n",argv[0]);
    fprintf(stdout,"that is/are binary 3D grids of densities.\n");
    fprintf(stdout,"Each grid is reflected along chosen axis and\n");
    fprintf(stdout,"written to a new file with a .ref(n) extension.  \n");
    fprintf(stdout,"(n) in .ref(n) will be the axis along which the\n");
    fprintf(stdout,"reflection took place. \n\n");
    fprintf(stdout,"    Example:   %s 2 Ell30 Ell20 \n\n",argv[0]);
    fprintf(stdout,"This example would reflect the files Ell30 & Ell20\n");
    fprintf(stdout,"along the y axis, and write them to new files called\n");
    fprintf(stdout,"Ell30.refy and Ell20.refy\n\n");
    fflush(stdout);
    exit(1);
}
