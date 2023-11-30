/* Atomic Program by Sverker Edvardsson & Daniel Aberg
   July 1998 */

/* Edited in 2023 by Juan David Lizarazo*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <zlib.h>

#define PI 3.141592654
#define IARITHIF(x)  (x > 0 ? 1 : x < 0 ? -1 : 0)

/* OBS! only equivalent electrons allowed ! */
/*                       N                 */
/* Program general for nl  configuration except for the CI */
/* parts of Trees parameters i.e. functions t2T2-t8T8 */
/* Those are valid for f-electrons only. */


/* atom/material dependecies are in files:
   param.inp, pauli.inp, cfp.inp and vk1k2k3.gz */


int STATES, electrons;

typedef double doublereal;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int integer;

double clebsg(double a, double b, double c, double xx, double yy, double zz, 
double gg, double hh, double pp,int mode);
double sign(double a,double b);

struct t_fgercm {
        int ierr, ierct;
	}       fgercm;


struct state {
double l;
double ml;
double ms;
};

// const char *lanthanides[] = {"Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm"};

struct state **states, **states_tmp;

int signe, *un_paired;

double lx(double la, double mla, double lb, double mlb);
double ly(double la, double mla, double lb, double mlb);

double sx(double msa, double msb);
double sy(double msa, double msb);

double Spin_Orbit(double la,double mla, double msa, double lb,double mlb, double msb);

double Spin_Orbit_Xi(void);

double Matrix_element_SO(int method,int i, int j);

void Matrix_element_Hcf(int method,int i, int j);

void Crystal_Field(double la,double mla, double msa, double lb,double mlb, double msb);

double Atp[20][20], Btp[20][20];

void readatp(void);

double RE_Atp(int t,int p);

double IM_Atp(int t,int p);

double Rt(int t);

double Sigmat(int t);

double RE_Matrix_element_Hcf=0.0, RE_CFP=0.0, IM_Matrix_element_Hcf=0.0, IM_CFP=0.0;

int diag1(integer n);
void diag2(integer n);

double ck(double k,double la,double mla, double lb,double mlb);

double General(double la,double mla, double msa, double lb,double mlb, double msb, double lc,double mlc, double msc, double ld,double mld, double msd);

double Matrix_element_Hm(int method, int i, int j);

double H_m(double la,double mla, double msa, double lb,double mlb, double msb);

double RE_Matrix_element_Hm=0.0, IM_Matrix_element_Hm=0.0;

double RE_H=0.0, IM_H=0.0;

double Matrix_element_alpha(int method, int i, int j);

double Matrix_element_beta(int method, int i, int j);

double Matrix_element_gamma1(int method, int i, int j);

double Matrix_element_gamma2(double power, int method, int i, int j);

double General_gamma(double power, double la,double mla, double msa, double lb,double mlb, double msb, double lc,double mlc, double msc, double ld,double mld, double msd);

double General_Trees_small(double power, double q, double la,double mla, double msa, double lb,double mlb, double msb);

double Reduce_General_Trees_small(double power1, double q1, double l1, double ml1, double ms1, double l2, double ml2, double ms2, double power2, double q2, double l3, double ml3, double ms3, double l4, double ml4, double ms4, double power3, double q3, double l5, double ml5, double ms5, double l6, double ml6, double ms6);

double General_Trees(double power1, double power2, double power3, double q1, double q2, double q3, double la,double mla, double msa, double lb,double mlb, double msb, double lc,double mlc, double msc, double ld,double mld, double msd,double le,double mle, double mse, double lf,double mlf, double msf );

double Matrix_element_Trees(double power1, double power2, double power3, double q1, double q2, double q3, int method, int i, int j);

double Vkkpkb(double power1, double power2, double power3, int method, int i, int j);

double Vkkpkbq1q2q3(double power1, double power2, double power3, double q1, double q2, double q3, int method, int i, int j);

double V222(int method, int i, int j);
double V444(int method, int i, int j);
double V666(int method, int i, int j);
double V224(int method, int i, int j);
double V244(int method, int i, int j);
double V446(int method, int i, int j);
double V266(int method, int i, int j);
double V466(int method, int i, int j);

void calc_Vkkpkb(int method, int i, int j);

double t2T2(void);
double t3T3(void);
double t4T4(void);
double t6T6(void);
double t7T7(void);
double t8T8(void);

double Matrix_element(int method,int i, int j);

double Delta(double a, double b);

double Rk(double k);

void load_pauli(void);

void save_Vkkpkb(int i, int j);

void load_Vkkpkb(int i, int j);

void load_parameters(void);

double R2, R4, R6, S2, S4, S6;

double F2, F4, F6, Xi;

/* CI parameters */

double Alpha, Beta, Gamma;

double T2,T3,T4,T6,T7,T8;

double Bx, By, Bz;

double Vkkpkb222,Vkkpkb224,Vkkpkb244,Vkkpkb246,Vkkpkb444;
double Vkkpkb446,Vkkpkb266,Vkkpkb466,Vkkpkb666;

gzFile pek_save;

int main(int argc, char **argv)
{

    int i,j,method,elements=0;
    int check, load, index;
    double RE_matrix, IM_matrix;
    gzFile pek;

    // Check for proper argument count
    if (argc != 4) {
        printf("Input format: load STATES electrons\n");
        exit(EXIT_FAILURE);
    }

    // Read command line arguments
    // load=1 read vkkpkb and diag matrix
    // load=0 calc. vkkpkb
    load = atoi(argv[1]);
    STATES = atoi(argv[2]);
    electrons = atoi(argv[3]);

    printf("\n\nload=%d #STATES=%d #electrons=%d\n\n", load, STATES, electrons);

    // Memory allocation for unpaired electrons
    if(!(un_paired = malloc(electrons * sizeof(int)))) {
      printf("No mem_un_paired!"); 
      exit(0);
      }

    if(!(states = malloc(STATES * sizeof(double)))) {
      printf("No mem_states!"); 
      exit(0);
      }

    for(index=0;index < STATES; index++)
      if(!( states[index]=malloc(electrons * sizeof(struct state)))) {
        printf("No mem_states!");
        exit(0);
        }

    // Memory allocation for temporary states
    if(!(states_tmp = malloc(STATES * sizeof(double)))) {
      printf("No mem_states_tmp!"); 
      exit(0);
      }

    for(index=0;index < STATES; index++)
    if(!( states_tmp[index]=malloc(electrons * sizeof(struct state)))) {
      printf("No mem_states_tmp!"); 
      exit(0);
      }

    // If load is 0 and the number of three-particle elements is large, notify the user
    if (load == 0 && 0.5 * STATES * (STATES + 1.0) > 50000.0) {
        printf("\n\nCalculating %ld three-particle elements...\n", (long)(0.5 * STATES * (STATES + 1.0)));
        printf("Go get coffee!\n\n\n");
    }

    // Initialize required data and parameters
    printf("Reading free-ion params ...\n");
    readatp();
    printf("Loading pauli.inp ...\n");
    load_pauli();
    printf("Reading cf params ...\n");
    load_parameters();

    // Open matrix file for writing or reading
    if (load == 1) {
        pek = gzopen("matrix.gz", "wb");
        if (!pek) {
            puts("Error in matrix.gz!");
            exit(EXIT_FAILURE);
        }
    }

    // Open vk1k2k3 file for writing or reading based on 'load'
    if (load == 1) {
        pek_save = gzopen("vk1k2k3.gz", "rb");
        if (!pek_save) {
            puts("Error in vk1k2k3.gz!");
            exit(EXIT_FAILURE);
        }
    } else {
        pek_save = gzopen("vk1k2k3.gz", "wb");
        if (!pek_save) {
            puts("Error in vk1k2k3.gz!");
            exit(EXIT_FAILURE);
        }
    }

    // Matrix calculation loop
    for (j = 0; j < STATES; j++) 
    {
        for (i = 0; i <= j; i++)
        {
            // Calculate upper triangle of the matrix
            if (load == 0)
            {
                elements++;
                if (elements % 500 == 0)
                {
                    printf("#elements = %ld\n", elements);
                }
                method = method_judd(i, j);
                calc_Vkkpkb(method, i, j);
                save_Vkkpkb(i, j); // Enormously CPU intensive, so save it
            } else
            {
                elements++;
                method = method_judd(i, j);
                Matrix_element_Hcf(method, i, j);
                Matrix_element_Hm(method, i, j);
                load_Vkkpkb(i, j); // If it's already saved, get it from here!
                // Calculate real and imaginary parts of the matrix element
                RE_matrix = Matrix_element(method, i, j) 
                            + Matrix_element_SO(method, i, j) 
                            + Matrix_element_gamma1(method, i, j) 
                            + Matrix_element_beta(method, i, j) 
                            + Matrix_element_alpha(method, i, j)
                            + RE_Matrix_element_Hcf 
                            + RE_Matrix_element_Hm 
                            + t2T2() 
                            + t3T3() 
                            + t4T4() 
                            + t6T6() 
                            + t7T7() 
                            + t8T8();

                IM_matrix = IM_Matrix_element_Hcf + IM_Matrix_element_Hm;

                // Zero out matrix elements below a threshold
                if (fabs(RE_matrix) < 1e-10) RE_matrix = 0.0;
                if (fabs(IM_matrix) < 1e-10) IM_matrix = 0.0;

                // Write matrix elements to file
                gzprintf(pek, "%12.10lf %12.10lf\n", RE_matrix, IM_matrix);

                // Progress logging
                if (elements % 1000 == 0) {
                    printf("#elements = %ld\n", elements);
                }
            }
        }
    }

    // Cleanup and close files
    if (load == 1)
    {
        gzclose(pek);
    }
    gzclose(pek_save);

    // Free allocated memory
    free(un_paired);
    free(states);
    free(states_tmp);

    // Attempt to diagonalize the matrix if 'load' is 1
    if (load == 1) {
        printf("Using diag1 ... \n");
        check = diag1(STATES); // If enough memory try divide and conquer algorithm

        if (check < 0) {
            printf("Using diag2 ... \n");
            diag2(STATES); // Otherwise use the standard routine
        }
    }

    return EXIT_SUCCESS;
}


void load_parameters(void)
{
    char comment[81];
    FILE *pek_param;

    if((pek_param=fopen("param.inp","r"))==NULL) {
            puts("Error in param.inp!");
            exit(0);
            }

    fscanf(pek_param,"%lf%lf%lf%lf%lf%lf", &R2, &R4, &R6, &S2, &S4, &S6);
    printf("%lf %lf %lf %lf %lf %lf\n", R2, R4, R6, S2, S4, S6);

    fscanf(pek_param,"%lf%lf%lf%lf", &F2, &F4, &F6, &Xi);
    printf("%lf %lf %lf %lf\n", F2, F4, F6, Xi);

    fscanf(pek_param,"%lf%lf%lf", &Alpha, &Beta, &Gamma);
    printf("%lf %lf %lf\n", Alpha, Beta, Gamma);

    fscanf(pek_param,"%lf%lf%lf%lf%lf%lf", &T2, &T3, &T4, &T6, &T7, &T8);
    printf("%lf %lf %lf %lf %lf %lf\n", T2, T3, T4, T6, T7, T8);

    fscanf(pek_param,"%lf%lf%lf", &Bx, &By, &Bz);
    printf("%lf %lf %lf\n", Bx, By, Bz);

    fclose(pek_param);
    } 


void load_Vkkpkb(int i, int j)
{
    int gn,a,b,gi,gj,gk;
    char buf[1], line[200], flyt[20];
    double CI[11];

    for(gj=0; gj<200; gj++)
        line[gj]=' ';
        gi=-1;
        while((buf[0]=gzgetc(pek_save))!='\n')
        {
            gi++;
            if(gi>199)
            {printf("Index problem in load_Vkkpkb!\n"); exit(0);}
            line[gi]=buf[0];
        }

    gj=0;
    for(gn=0; gn<11; gn++)
    {
    for(gk=0; gk<20; gk++)
        flyt[gk]=' ';
        gi=-1;
        while( ((int)line[gj]) != 32 )
        {
            gi++;
            if( (gi>19) || (gj> 199) )
                {printf("Index problem in load_Vkkpkb!\n"); exit(0);}
            flyt[gi]=line[gj];
            gj++;
        }
    if(gn==0)
        a=atoi(flyt);
    else if(gn==1)
        b=atoi(flyt);
    else
        CI[gn]=atof(flyt);
    gj++;
    }

    Vkkpkb222=CI[2];
    Vkkpkb224=CI[3];
    Vkkpkb244=CI[4];
    Vkkpkb246=CI[5];
    Vkkpkb444=CI[6];
    Vkkpkb446=CI[7];
    Vkkpkb266=CI[8];
    Vkkpkb466=CI[9];
    Vkkpkb666=CI[10];

    if( (a==i)&&(b==j) ) return;
    else {printf("Problems in load_Vkkpkb!\n"); exit(0);}
}



void save_Vkkpkb(int i, int j)
{
    gzprintf(pek_save,"%d %d ",i,j);

    if(Vkkpkb222==0.0) 
        gzprintf(pek_save,"0 ");
    else
        gzprintf(pek_save,"%le ",Vkkpkb222);
    
    if(Vkkpkb224==0.0)
        gzprintf(pek_save,"0 ");
    else 
        gzprintf(pek_save,"%le ",Vkkpkb224);
    
    if(Vkkpkb244==0.0)
        gzprintf(pek_save,"0 ");
    else
        gzprintf(pek_save,"%le ",Vkkpkb244);

    if(Vkkpkb246==0.0)
        gzprintf(pek_save,"0 ");
    else
        gzprintf(pek_save,"%le ",Vkkpkb246);
    
    if(Vkkpkb444==0.0)
        gzprintf(pek_save,"0 ");
    else
        gzprintf(pek_save,"%le ",Vkkpkb444);

    if(Vkkpkb446==0.0)
        gzprintf(pek_save,"0 ");
    else 
        gzprintf(pek_save,"%le ",Vkkpkb446);
    
    if(Vkkpkb266==0.0)
        gzprintf(pek_save,"0 ");
    else 
        gzprintf(pek_save,"%le ",Vkkpkb266);
    
    if(Vkkpkb466==0.0)
        gzprintf(pek_save,"0 ");
    else 
        gzprintf(pek_save,"%le ",Vkkpkb466);
    
    if(Vkkpkb666==0.0)
        gzprintf(pek_save,"0\n");
    else 
        gzprintf(pek_save,"%le\n",Vkkpkb666);

    return;

}


double Matrix_element_gamma1(int method, int i, int j)
{
    double sum=0.0, power;
    int t;

    for(t=1; t<=states[0][0].l; t++)
        {
        power = 2.0*(double)t -1.0;
        sum = sum + (2.0*power+1.0)*Matrix_element_gamma2(power, method, i, j);
        }

    return(Gamma*(2.0/5.0)*sum);

}


double Matrix_element_beta(int method, int i, int j)
{

 return(Beta * (1.5*Matrix_element_gamma2(1.0, method, i, j) + 
                5.5*Matrix_element_gamma2(5.0, method, i, j)));

}


double Matrix_element_alpha(int method, int i, int j)
{
    double l,C;
    l = states[0][0].l;
    C = 2.0*l*(l+1.0)*(2.0*l+1.0);
    return(Alpha * C * Matrix_element_gamma2(1.0, method, i, j));

}





double Matrix_element_gamma2(double power, int method, int i, int j)
{

    int k,l;
    double sum=0.0;


    if(method==0)
    {
    for(k=0; k<electrons-1; k++)
    {
    for(l=k+1; l<electrons; l++)
        {
        sum=sum + General_gamma(power, states[i][l].l, states[i][l].ml,states[i][l].ms, states[i][k].l, states[i][k].ml,states[i][k].ms,states[i][l].l, states[i][l].ml,states[i][l].ms,states[i][k].l, states[i][k].ml,states[i][k].ms) 
                    - 

    General_gamma(power, states[i][l].l, states[i][l].ml,states[i][l].ms, states[i][k].l, states[i][k].ml,states[i][k].ms,states[i][k].l, states[i][k].ml,states[i][k].ms,states[i][l].l, states[i][l].ml,states[i][l].ms);


        }
    }


    return(sum*pow(-1.0,(double)signe));


    }

    else if(method==1)
    {

    for(k=0; k<electrons; k++)
    sum=sum + General_gamma(power, states[i][k].l, states[i][k].ml,states[i][k].ms, states[i][un_paired[1]].l, states[i][un_paired[1]].ml,states[i][un_paired[1]].ms,states[i][k].l, states[i][k].ml,states[i][k].ms,states_tmp[j][un_paired[1]].l, states_tmp[j][un_paired[1]].ml,states_tmp[j][un_paired[1]].ms)                  

                - 

    General_gamma(power, states[i][k].l, states[i][k].ml,states[i][k].ms, states[i][un_paired[1]].l, states[i][un_paired[1]].ml,states[i][un_paired[1]].ms,states_tmp[j][un_paired[1]].l, states_tmp[j][un_paired[1]].ml,states_tmp[j][un_paired[1]].ms,states[i][k].l, states[i][k].ml,states[i][k].ms);

    return(sum*pow(-1.0,(double)signe));

    }

    else if(method==2)
    {

    return((General_gamma(power, states[i][un_paired[1]].l, states[i][un_paired[1]].ml,states[i][un_paired[1]].ms, states[i][un_paired[2]].l, states[i][un_paired[2]].ml,states[i][un_paired[2]].ms,states_tmp[j][un_paired[1]].l, states_tmp[j][un_paired[1]].ml,states_tmp[j][un_paired[1]].ms,states_tmp[j][un_paired[2]].l, states_tmp[j][un_paired[2]].ml,states_tmp[j][un_paired[2]].ms)
                - 

    General_gamma(power, states[i][un_paired[1]].l, states[i][un_paired[1]].ml,states[i][un_paired[1]].ms, states[i][un_paired[2]].l, states[i][un_paired[2]].ml,states[i][un_paired[2]].ms,states_tmp[j][un_paired[2]].l, states_tmp[j][un_paired[2]].ml,states_tmp[j][un_paired[2]].ms,states_tmp[j][un_paired[1]].l, states_tmp[j][un_paired[1]].ml,states_tmp[j][un_paired[1]].ms)) * pow(-1.0,(double)signe));

    }

    else return(0.0);

}





double General_gamma(double power, double la,double mla, double msa, double lb,double mlb, double msb, double lc,double mlc, double msc, double ld,double mld, double msd)
{
    double q;
    double trej1,trej2;

    if(mla-mlc != mld-mlb) return(0.0);

    if((msa !=msc) || (msb != msd)) return(0.0);

    q=mla-mlc;

    trej1 = clebsg( la, power, lb, -mla, q, mlc,0, 0, 0, 1);
    trej2 = clebsg( la, power, lb, -mlb, -q, mld, 0, 0, 0, 1);

    return(pow(-1,mla+mld)*trej1*trej2);
}


void calc_Vkkpkb(int method, int i, int j)
{

    if(electrons<3) return;

    Vkkpkb222=V222(method, i, j);
    Vkkpkb224=V224(method, i, j);
    Vkkpkb244=V244(method, i, j);

    Vkkpkb246=Vkkpkb( 2.0, 4.0, 6.0, method, i, j);

    Vkkpkb444=V444(method, i, j);
    Vkkpkb446=V446(method, i, j);
    Vkkpkb266=V266(method, i, j);
    Vkkpkb466=V466(method, i, j);
    Vkkpkb666=V666(method, i, j);

    return;

}



double t2T2(void)
{
    if(electrons<3) return(0.0);

    return(T2*
    (sqrt(605.0/5292.0)       * Vkkpkb222 -
    sqrt(6760.0/43659.0)      * Vkkpkb224 -
    sqrt(1805.0/391314.0)     * Vkkpkb244 -
    sqrt(4160.0/754677.0)     * Vkkpkb246 +
    sqrt(55016.0/717409.0)    * Vkkpkb444 -
    sqrt(195.0/204974.0)      * Vkkpkb446 +
    sqrt(1625.0/143748.0)     * Vkkpkb266 +
    sqrt(88400.0/1185921.0)   * Vkkpkb466 +
    sqrt(29393.0/790614.0)    * Vkkpkb666));
}

double t3T3(void)
{
    if(electrons<3) return(0.0);

    return(T3*
    (sqrt(32761.0/889056.0)    * Vkkpkb222 +
    sqrt(33.0/1372.0)         * Vkkpkb224 -
    sqrt(4.0/33957.0)         * Vkkpkb244 -
    sqrt(13.0/264.0)          * Vkkpkb246 +
    sqrt(49972.0/622545.0)    * Vkkpkb444 +
    sqrt(52.0/1089.0)         * Vkkpkb446 +
    sqrt(325.0/199584.0)      * Vkkpkb266 -
    sqrt(442.0/12705.0)       * Vkkpkb466 +
    sqrt(205751.0/784080.0)   * Vkkpkb666));
}

double t4T4(void)
{
    if(electrons<3) return(0.0);

    return(T4*
    (sqrt(3575.0/889056.0)     * Vkkpkb222 -
    sqrt(325.0/37044.0)       * Vkkpkb224 -
    sqrt(54925.0/373527.0)    * Vkkpkb244 +
    sqrt(625.0/26136.0)       * Vkkpkb246 +
    sqrt(92480.0/1369599.0)   * Vkkpkb444 +
    sqrt(529.0/11979.0)       * Vkkpkb446 +
    sqrt(6889.0/2195424.0)    * Vkkpkb266 -
    sqrt(10880.0/251559.0)    * Vkkpkb466 -
    sqrt(79135.0/1724976.0)   * Vkkpkb666));
}

double t6T6(void)
{
    if(electrons<3) return(0.0);


    return(T6*
    (-sqrt(1573.0/8232.0)      * Vkkpkb222 -
    sqrt(15028.0/305613.0)    * Vkkpkb224 +
    sqrt(4693.0/12326391.0)   * Vkkpkb244 +
    sqrt(1568.0/107811.0)     * Vkkpkb246 -
    sqrt(297680.0/5021863.0)  * Vkkpkb444 -
    sqrt(49.0/395307.0)       * Vkkpkb446 -
    sqrt(1.0/223608.0)        * Vkkpkb266 -
    sqrt(174080.0/8301447.0)  * Vkkpkb466 +
    sqrt(79135.0/175692.0)    * Vkkpkb666));

}

double t7T7(void)
{
    if(electrons<3) return(0.0);

    return(T7*
    (sqrt(264407.0/823200.0)   * Vkkpkb222 +
    sqrt(28717.0/2778300.0)   * Vkkpkb224 -
    sqrt(1273597.0/28014525.0)* Vkkpkb244 +
    sqrt(841.0/1960200.0)     * Vkkpkb246 -
    sqrt(719104.0/2282665.0)  * Vkkpkb444 -
    sqrt(1369.0/35937.0)      * Vkkpkb446 +
    sqrt(625.0/81312.0)       * Vkkpkb266 -
    sqrt(8704.0/3773385.0)    * Vkkpkb466 +
    sqrt(15827.0/319440.0)    * Vkkpkb666));

}

double t8T8(void)
{
    if(electrons<3) return(0.0);

    return(T8*
    (sqrt(21879.0/274400.0)    * Vkkpkb222 -
    sqrt(37349.0/926100.0)    * Vkkpkb224 +
    sqrt(849524.0/9338175.0)  * Vkkpkb244 -
    sqrt(17.0/653400.0)       * Vkkpkb246 -
    sqrt(73644.0/2282665.0)   * Vkkpkb444 +
    sqrt(68.0/11979.0)        * Vkkpkb446 +
    sqrt(1377.0/27104.0)      * Vkkpkb266 -
    sqrt(103058.0/1257795.0)  * Vkkpkb466 -
    sqrt(8379.0/106480.0)     * Vkkpkb666));

}


/* The below stuff speed up the Vkkpkb calc's (no approx.) */

double V222(int method, int i, int j)
{
    return(
    6.0 * Vkkpkbq1q2q3(2.0, 2.0, 2.0, -2.0, 0.0, 2.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(2.0, 2.0, 2.0, -2.0, 1.0, 1.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(2.0, 2.0, 2.0, -1.0, -1.0, 2.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(2.0, 2.0, 2.0, -1.0, 0.0, 1.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 2.0, 2.0, 0.0, 0.0, 0.0, method, i, j));

}


double V444(int method, int i, int j)
{
    return(
    6.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -4.0, 0.0, 4.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -4.0, 1.0, 3.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -4.0, 2.0, 2.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -3.0, -1.0, 4.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -3.0, 0.0, 3.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -3.0, 1.0, 2.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -2.0, -2.0, 4.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -2.0, -1.0, 3.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -2.0, 0.0, 2.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -2.0, 1.0, 1.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -1.0, -1.0, 2.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, -1.0, 0.0, 1.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 4.0, 4.0, 0.0, 0.0, 0.0, method, i, j));

}


double V666(int method, int i, int j)
{

    return(
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -6.0, 0.0, 6.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -6.0, 1.0, 5.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -6.0, 2.0, 4.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -6.0, 3.0, 3.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -5.0, -1.0, 6.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -5.0, 0.0, 5.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -5.0, 1.0, 4.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -5.0, 2.0, 3.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -4.0, -2.0, 6.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -4.0, -1.0, 5.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -4.0, 0.0, 4.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -4.0, 1.0, 3.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -4.0, 2.0, 2.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -3.0, -3.0, 6.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -3.0, -2.0, 5.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -3.0, -1.0, 4.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -3.0, 0.0, 3.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -3.0, 1.0, 2.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -2.0, -2.0, 4.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -2.0, -1.0, 3.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -2.0, 0.0, 2.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -2.0, 1.0, 1.0, method, i, j)  +
    3.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -1.0, -1.0, 2.0, method, i, j)  +
    6.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, -1.0, 0.0, 1.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(6.0, 6.0, 6.0, 0.0, 0.0, 0.0, method, i, j));

}

double V224(int method, int i, int j)
{

    return(
    1.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, -2.0, -2.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, -2.0, -1.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, -2.0, 0.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, -2.0, 1.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, -2.0, 2.0, 0.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, -1.0, -1.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, -1.0, 0.0, 1.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, 0.0, 0.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, 0.0, 1.0, -1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, 0.0, 2.0, -2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, 1.0, -1.0, 0.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, 1.0, 1.0, -2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, 1.0, 2.0, -3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, 2.0, -1.0, -1.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 2.0, 4.0, 2.0, 2.0, -4.0, method, i, j));


}

double V244(int method, int i, int j)
{

    return(
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 2.0, -4.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 1.0, -4.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 0.0, -4.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 2.0, -3.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 1.0, -3.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 0.0, -3.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, -1.0, -3.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 2.0, -2.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 1.0, -2.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 0.0, -2.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, -1.0, -2.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, -2.0, -2.0, 4.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 2.0, -1.0, -1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 1.0, -1.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 0.0, -1.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, -1.0, -1.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, -2.0, -1.0, 3.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, 0.0, 0.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, -1.0, 0.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, -2.0, 0.0, 2.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 4.0, 4.0, -2.0, 1.0, 1.0, method, i, j));


}


double V446(int method, int i, int j)
{

    return(
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -4.0, -2.0, 6.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -4.0, -1.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -4.0, 0.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -4.0, 1.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -4.0, 2.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -4.0, 3.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -4.0, 4.0, 0.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -3.0, -3.0, 6.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -3.0, -2.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -3.0, -1.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -3.0, 0.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -3.0, 1.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -3.0, 2.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -3.0, 3.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -3.0, 4.0, -1.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -2.0, -2.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -2.0, -1.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -2.0, 0.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -2.0, 1.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -2.0, 2.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -2.0, 3.0, -1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -2.0, 4.0, -2.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -1.0, -1.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -1.0, 0.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -1.0, 1.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -1.0, 2.0, -1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -1.0, 3.0, -2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, -1.0, 4.0, -3.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 0.0, 0.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 0.0, 1.0, -1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 0.0, 2.0, -2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 0.0, 3.0, -3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 0.0, 4.0, -4.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 1.0, 1.0, -2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 1.0, 2.0, -3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 1.0, 3.0, -4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 1.0, 4.0, -5.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 2.0, 2.0, -4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 2.0, 3.0, -5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 2.0, 4.0, -6.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 4.0, 6.0, 3.0, 3.0, -6.0, method, i, j));

}



double V266(int method, int i, int j)
{

    return(
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 2.0, -6.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 1.0, -6.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 0.0, -6.0, 6.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 2.0, -5.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 1.0, -5.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 0.0, -5.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -1.0, -5.0, 6.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 2.0, -4.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 1.0, -4.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 0.0, -4.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -1.0, -4.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -2.0, -4.0, 6.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 2.0, -3.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 1.0, -3.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 0.0, -3.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -1.0, -3.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -2.0, -3.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 2.0, -2.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 1.0, -2.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 0.0, -2.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -1.0, -2.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -2.0, -2.0, 4.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 2.0, -1.0, -1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 1.0, -1.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 0.0, -1.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -1.0, -1.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -2.0, -1.0, 3.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, 0.0, 0.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -1.0, 0.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -2.0, 0.0, 2.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(2.0, 6.0, 6.0, -2.0, 1.0, 1.0, method, i, j));

}



double V466(int method, int i, int j)
{

    return(
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 4.0, -6.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 3.0, -6.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 2.0, -6.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 1.0, -6.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 0.0, -6.0, 6.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 4.0, -5.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 3.0, -5.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 2.0, -5.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 1.0, -5.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 0.0, -5.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -1.0, -5.0, 6.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 4.0, -4.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 3.0, -4.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 2.0, -4.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 1.0, -4.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 0.0, -4.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -1.0, -4.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -2.0, -4.0, 6.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 4.0, -3.0, -1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 3.0, -3.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 2.0, -3.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 1.0, -3.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 0.0, -3.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -1.0, -3.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -2.0, -3.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -3.0, -3.0, 6.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 4.0, -2.0, -2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 3.0, -2.0, -1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 2.0, -2.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 1.0, -2.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 0.0, -2.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -1.0, -2.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -2.0, -2.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -3.0, -2.0, 5.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -4.0, -2.0, 6.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 2.0, -1.0, -1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 1.0, -1.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 0.0, -1.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -1.0, -1.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -2.0, -1.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -3.0, -1.0, 4.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -4.0, -1.0, 5.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, 0.0, 0.0, 0.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -1.0, 0.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -2.0, 0.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -3.0, 0.0, 3.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -4.0, 0.0, 4.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -2.0, 1.0, 1.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -3.0, 1.0, 2.0, method, i, j)  +
    2.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -4.0, 1.0, 3.0, method, i, j)  +
    1.0 * Vkkpkbq1q2q3(4.0, 6.0, 6.0, -4.0, 2.0, 2.0, method, i, j));

}



double Vkkpkbq1q2q3(double power1, double power2, double power3, double q1, double q2, double q3, int method, int i, int j)
{

    double trej;

    trej = clebsg( power1, power2, power3, q1, q2, q3,0, 0, 0, 1);

    if(fabs(trej)>0.0)
    return( 6.0*sqrt( (2.0*power1+1.0)*(2.0*power2+1.0)*(2.0*power3+1.0) ) * trej * Matrix_element_Trees( power1, power2, power3, q1, q2, q3, method, i, j));
    else
    return(0.0);

}

double Vkkpkb(double power1, double power2, double power3, int method, int i, int j)
{

    int q1, q2, q3;
    double sum=0.0,trej;



    for(q1=-(int)power1; q1<=(int)power1; q1++)
    {
    for(q2=-(int)power2; q2<=(int)power2; q2++)
        {
        for(q3=-(int)power3; q3<=(int)power3; q3++)
        {
        if(q1+q2+q3==0)
        {
        /* if( (power1==2.0)&&(power2==2.0)&&(power3==4.0) )
    {printf("%d %d %d\n",q1,q2,q3); getchar();} */

        trej = clebsg( power1, power2, power3, (double)q1, (double)q2, (double)q3,0, 0, 0, 1);

            if(fabs(trej)>0.0)
        sum = sum + trej *

            Matrix_element_Trees( power1, power2, power3, (double)q1, (double)q2, (double)q3, method, i, j);
            
        }

        }
        }
    }

    return( 6.0*sqrt( (2.0*power1+1.0)*(2.0*power2+1.0)*(2.0*power3+1.0) )*sum );

    /* factor 6 because a tripple has 3! permutations */

}





double General_Trees_small(double power, double q, double la,double mla, double msa, double lb,double mlb, double msb)
{

    double trej;


    if(msa != msb) return(0.0);

    if(-mla+q+mlb != 0.0) return(0.0);


    trej = clebsg( la, power, la, -mla, q, mlb,0, 0, 0, 1);


    return(pow(-1.0,la-mla)*trej);

    }

    double Reduce_General_Trees_small(double power1, double q1, double l1, double ml1, double ms1, double l2, double ml2, double ms2, double power2, double q2, double l3, double ml3, double ms3, double l4, double ml4, double ms4, double power3, double q3, double l5, double ml5, double ms5, double l6, double ml6, double ms6)

    {
    double a,b,c;

    if(-ml1+q1+ml2!=0.0 )
    return(0.0);
    if(-ml3+q2+ml4!=0.0 )
    return(0.0);
    if(-ml5+q3+ml6!=0.0 )
    return(0.0);


    a=General_Trees_small(power1, q1, l1, ml1, ms1, l2, ml2,  ms2);


    b=General_Trees_small(power2, q2, l3, ml3, ms3, l4, ml4,  ms4);


    c=General_Trees_small(power3, q3, l5, ml5, ms5, l6, ml6,  ms6);


    return(a*b*c);
    }


    double General_Trees(double power1, double power2, double power3, double q1, double q2, double q3, double la,double mla, double msa, double lb,double mlb, double msb, double lc,double mlc, double msc, double ld,double mld, double msd,double le,double mle, double mse, double lf,double mlf, double msf )
    {

    double t1, t2, t3, t4, t5, t6;


    t1 = Reduce_General_Trees_small(power1, q1, lc, mlc, msc, lf, mlf,  msf, power2, q2, lb, mlb, msb, le, mle,  mse, power3, q3, la, mla, msa, ld, mld,  msd) 

    -    
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, le, mle,  mse, power2, q2, lb, mlb, msb, lf, mlf,  msf, power3, q3, la, mla, msa, ld, mld,  msd)
    
    -    
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, lf, mlf,  msf, power2, q2, lb, mlb, msb, ld, mld,  msd, power3, q3, la, mla, msa, le, mle,  mse) 

    +    
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, ld, mld,  msd, power2, q2, lb, mlb, msb, lf, mlf,  msf, power3, q3, la, mla, msa, le, mle,  mse) 

    +    
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, le, mle,  mse, power2, q2, lb, mlb, msb, ld, mld,  msd, power3, q3, la, mla, msa, lf, mlf,  msf) 

    -    
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, ld, mld,  msd, power2, q2, lb, mlb, msb, le, mle,  mse, power3, q3, la, mla, msa, lf, mlf,  msf);



    

    t2 = -
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, lf, mlf,  msf, power2, q2, lc, mlc, msc, le, mle,  mse, power3, q3, la, mla, msa, ld, mld,  msd) 

    +    
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, le, mle,  mse, power2, q2, lc, mlc, msc, lf, mlf,  msf, power3, q3, la, mla, msa, ld, mld,  msd)  
    
    +    
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, lf, mlf,  msf, power2, q2, lc, mlc, msc, ld, mld,  msd, power3, q3, la, mla, msa, le, mle,  mse) 

    -
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, ld, mld,  msd, power2, q2, lc, mlc, msc, lf, mlf,  msf, power3, q3, la, mla, msa, le, mle,  mse) 

    -    
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, le, mle,  mse, power2, q2, lc, mlc, msc, ld, mld,  msd, power3, q3, la, mla, msa, lf, mlf,  msf) 

    +    
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, ld, mld,  msd,power2, q2, lc, mlc, msc, le, mle,  mse, power3, q3, la, mla, msa, lf, mlf,  msf);



    t3 = -
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, lf, mlf,  msf, power2, q2, la, mla, msa, le, mle,  mse, power3, q3, lb, mlb, msb, ld, mld,  msd) 

    +    
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, le, mle,  mse, power2, q2, la, mla, msa, lf, mlf,  msf, power3, q3, lb, mlb, msb, ld, mld,  msd)  
    
    +    
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, lf, mlf,  msf, power2, q2, la, mla, msa, ld, mld,  msd, power3, q3, lb, mlb, msb, le, mle,  mse) 

    -    
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, ld, mld,  msd, power2, q2, la, mla, msa, lf, mlf,  msf, power3, q3, lb, mlb, msb, le, mle,  mse) 

    -    
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, le, mle,  mse, power2, q2, la, mla, msa, ld, mld,  msd, power3, q3, lb, mlb, msb, lf, mlf,  msf) 

    +    
    Reduce_General_Trees_small(power1, q1, lc, mlc, msc, ld, mld,  msd, power2, q2, la, mla, msa, le, mle,  mse, power3, q3, lb, mlb, msb, lf, mlf,  msf);



    t4 = 
    Reduce_General_Trees_small(power1, q1, la, mla, msa, lf, mlf,  msf, power2, q2, lc, mlc, msc, le, mle,  mse, power3, q3, lb, mlb, msb, ld, mld,  msd) 

    -    
    Reduce_General_Trees_small(power1, q1, la, mla, msa, le, mle,  mse, power2, q2, lc, mlc, msc, lf, mlf,  msf, power3, q3, lb, mlb, msb, ld, mld,  msd)  
    
    -    
    Reduce_General_Trees_small(power1, q1, la, mla, msa, lf, mlf,  msf, power2, q2, lc, mlc, msc, ld, mld,  msd, power3, q3, lb, mlb, msb, le, mle,  mse) 

    +    
    Reduce_General_Trees_small(power1, q1, la, mla, msa, ld, mld,  msd, power2, q2, lc, mlc, msc, lf, mlf,  msf, power3, q3, lb, mlb, msb, le, mle,  mse) 

    +    
    Reduce_General_Trees_small(power1, q1, la, mla, msa, le, mle,  mse, power2, q2, lc, mlc, msc, ld, mld,  msd, power3, q3, lb, mlb, msb, lf, mlf,  msf) 

    -    
    Reduce_General_Trees_small(power1, q1, la, mla, msa, ld, mld,  msd, power2, q2, lc, mlc, msc, le, mle,  mse, power3, q3, lb, mlb, msb, lf, mlf,  msf);


    t5 = 
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, lf, mlf,  msf, power2, q2, la, mla, msa, le, mle,  mse, power3, q3, lc, mlc, msc, ld, mld,  msd) 

    -    
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, le, mle,  mse, power2, q2, la, mla, msa, lf, mlf,  msf, power3, q3, lc, mlc, msc, ld, mld,  msd)  
    
    -    
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, lf, mlf,  msf, power2, q2, la, mla, msa, ld, mld,  msd, power3, q3, lc, mlc, msc, le, mle,  mse) 

    +    
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, ld, mld,  msd, power2, q2, la, mla, msa, lf, mlf,  msf, power3, q3, lc, mlc, msc, le, mle,  mse) 

    +    
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, le, mle,  mse, power2, q2, la, mla, msa, ld, mld,  msd, power3, q3, lc, mlc, msc, lf, mlf,  msf) 

    -    
    Reduce_General_Trees_small(power1, q1, lb, mlb, msb, ld, mld,  msd, power2, q2, la, mla, msa, le, mle,  mse, power3, q3, lc, mlc, msc, lf, mlf,  msf);


    t6 =- 
    Reduce_General_Trees_small(power1, q1, la, mla, msa, lf, mlf,  msf, power2, q2, lb, mlb, msb, le, mle,  mse, power3, q3, lc, mlc, msc, ld, mld,  msd) 

    +    
    Reduce_General_Trees_small(power1, q1, la, mla, msa, le, mle,  mse, power2, q2, lb, mlb, msb, lf, mlf,  msf, power3, q3, lc, mlc, msc, ld, mld,  msd)  
    
    +    
    Reduce_General_Trees_small(power1, q1, la, mla, msa, lf, mlf,  msf, power2, q2, lb, mlb, msb, ld, mld,  msd, power3, q3, lc, mlc, msc, le, mle,  mse) 

    -    
    Reduce_General_Trees_small(power1, q1, la, mla, msa, ld, mld,  msd, power2, q2, lb, mlb, msb, lf, mlf,  msf, power3, q3, lc, mlc, msc, le, mle,  mse) 

    -    
    Reduce_General_Trees_small(power1, q1, la, mla, msa, le, mle,  mse, power2, q2, lb, mlb, msb, ld, mld,  msd, power3, q3, lc, mlc, msc, lf, mlf,  msf) 

    +    
    Reduce_General_Trees_small(power1, q1, la, mla, msa, ld, mld,  msd, power2, q2, lb, mlb, msb, le, mle,  mse, power3, q3, lc, mlc, msc, lf, mlf,  msf);


    return( (1.0/6.0) * ( t1+t2+t3+t4+t5+t6 ) );

}



/* Rules for three particle operators */

double Matrix_element_Trees(double power1, double power2, double power3, double q1, double q2, double q3, int method, int i, int j)
{

    int k,l,m;
    double sum=0.0;


    if(method==0)
    {

    for(k=0; k<electrons-2; k++)
    {
        for(l=k+1; l<electrons-1; l++)
        {
            for(m=l+1; m<electrons; m++)
            {
            sum = (sum 
                + General_Trees(power1,
                                power2,
                                power3,
                                q1,
                                q2,
                                q3,
                                states[i][k].l,
                                states[i][k].ml,
                                states[i][k].ms,
                                states[i][l].l,
                                states[i][l].ml,
                                states[i][l].ms,
                                states[i][m].l,
                                states[i][m].ml,
                                states[i][m].ms, 
                                states[i][k].l,
                                states[i][k].ml,
                                states[i][k].ms,
                                states[i][l].l,
                                states[i][l].ml,
                                states[i][l].ms,
                                states[i][m].l,
                                states[i][m].ml,
                                states[i][m].ms));
            }
        }
    }

    return(sum*pow(-1.0,(double)signe));


    }

    else if(method==1)
    {
        for(k=0; k<electrons-1; k++)
        {
            for(l=k+1; l<electrons; l++)
            {

            if(l!=un_paired[1])
                sum = (sum 
                + General_Trees(power1,
                                power2,
                                power3,
                                q1,
                                q2,
                                q3,
                                states[i][k].l,
                                states[i][k].ml,
                                states[i][k].ms,
                                states[i][l].l, 
                                states[i][l].ml,
                                states[i][l].ms,
                                states[i][un_paired[1]].l,
                                states[i][un_paired[1]].ml,
                                states[i][un_paired[1]].ms,
                                states[i][k].l,
                                states[i][k].ml,
                                states[i][k].ms,
                                states[i][l].l,
                                states[i][l].ml,
                                states[i][l].ms,
                                states_tmp[j][un_paired[1]].l,
                                states_tmp[j][un_paired[1]].ml,
                                states_tmp[j][un_paired[1]].ms));

            }
        }


        return(sum*pow(-1.0,(double)signe));

    }

    else if(method==2)
    {

    for(k=0; k<electrons; k++)
    {
        if((k!=un_paired[1]) && (k!=un_paired[2]))
            sum = (sum 
                + General_Trees(power1,
                                power2, 
                                power3,
                                q1,
                                q2,
                                q3,
                                states[i][k].l,
                                states[i][k].ml,
                                states[i][k].ms,
                                states[i][un_paired[1]].l,
                                states[i][un_paired[1]].ml,
                                states[i][un_paired[1]].ms,
                                states[i][un_paired[2]].l,
                                states[i][un_paired[2]].ml,
                                states[i][un_paired[2]].ms, 
                                states[i][k].l,
                                states[i][k].ml,
                                states[i][k].ms,
                                states_tmp[j][un_paired[1]].l,
                                states_tmp[j][un_paired[1]].ml,
                                states_tmp[j][un_paired[1]].ms,
                                states_tmp[j][un_paired[2]].l,
                                states_tmp[j][un_paired[2]].ml,
                                states_tmp[j][un_paired[2]].ms));

    }
    return(sum*pow(-1.0,(double)signe));

    }

    else if(method==3)
    {


    return(General_Trees(power1,
                         power2,
                         power3,
                         q1,
                         q2,
                         q3,
                         states[i][un_paired[1]].l,
                         states[i][un_paired[1]].ml,
                         states[i][un_paired[1]].ms,
                         states[i][un_paired[2]].l,
                         states[i][un_paired[2]].ml,
                         states[i][un_paired[2]].ms,
                         states[i][un_paired[3]].l,
                         states[i][un_paired[3]].ml,
                         states[i][un_paired[3]].ms,
                         states_tmp[j][un_paired[1]].l,
                         states_tmp[j][un_paired[1]].ml,
                         states_tmp[j][un_paired[1]].ms,
                         states_tmp[j][un_paired[2]].l,
                         states_tmp[j][un_paired[2]].ml,
                         states_tmp[j][un_paired[2]].ms,
                         states_tmp[j][un_paired[3]].l,
                         states_tmp[j][un_paired[3]].ml,
                         states_tmp[j][un_paired[3]].ms) * pow(-1.0,(double)signe));


    }

    else return(0.0);

}

/* ordering of Determinantal Product States */

int method_judd(int i, int j)
{

    int k,l,m=0;
    double ml,ms;

    for (k=0; k<electrons; k++)
    {
        states_tmp[j][k].l=states[j][k].l;
        states_tmp[j][k].ml=states[j][k].ml;
        states_tmp[j][k].ms=states[j][k].ms;
    }

    signe=0;

    for (k=0; k<electrons; k++)
    {
        for (l=0;l<electrons; l++)
        {
            if((states[i][k].ml==states_tmp[j][l].ml)
                &&(states[i][k].ms==states_tmp[j][l].ms
            ))
            {
        
            if(l!=k)
        { 
        ml=states_tmp[j][l].ml;
        ms=states_tmp[j][l].ms;

        states_tmp[j][l].ml=states_tmp[j][k].ml;
        states_tmp[j][l].ms=states_tmp[j][k].ms;

        states_tmp[j][k].ml=ml;
        states_tmp[j][k].ms=ms;
        signe++;
        }
        break;
        }
    
    } 
    
    }

    for(k=0; k<electrons; k++)
    {
        if((states[i][k].ml!=states_tmp[j][k].ml)||(states[i][k].ms!=states_tmp[j][k].ms))
        {
            m++;  
            un_paired[m]=k;
        }
    }

    return(m);
}


double Matrix_element_Hm(int method, int i, int j)
{

    int k;
    double sumr=0.0,sumi=0.0;

    if(method==0)
    {
        for(k=0; k<electrons; k++)
        {

            H_m( states[i][k].l,states[i][k].ml,states[i][k].ms, states_tmp[j][k].l,states_tmp[j][k].ml, states_tmp[j][k].ms);
            
            sumr=sumr+RE_H;
            sumi=sumi+IM_H;

        }

        RE_Matrix_element_Hm=sumr*pow(-1.0,(double)signe);
        IM_Matrix_element_Hm=sumi*pow(-1.0,(double)signe);
        return;

    }

    else if(method==1)
    {

        H_m(states[i][un_paired[1]].l,states[i][un_paired[1]].ml,states[i][un_paired[1]].ms, states_tmp[j][un_paired[1]].l,states_tmp[j][un_paired[1]].ml, states_tmp[j][un_paired[1]].ms);

        RE_Matrix_element_Hm=RE_H*pow(-1.0,(double)signe);
        IM_Matrix_element_Hm=IM_H*pow(-1.0,(double)signe);

        return;

    }
    
    RE_Matrix_element_Hm=0.0;
    IM_Matrix_element_Hm=0.0;
    
    return;

}

double H_m(double la,double mla, double msa, double lb,double mlb, double msb)
{
    double myb=0.5, gs=2*1.00115962, tesla=0.0000042551*219475.644;
    double sumr=0.0, sumi=0.0, sum=0.0;

    /*  Bz  */

    sumr=(Delta(mla,mlb)*Delta(msa,msb)*myb*Bz*tesla*(mlb+gs*msb));

    /*  Bx  */

    sum=Delta(msa,msb) * ( Delta(mla,mlb+1) * sqrt( lb*(lb+1)-mlb*(mlb+1) ) +
                            Delta(mla,mlb-1) * sqrt( lb*(lb+1)-mlb*(mlb-1) ) );
    sum=sum+ (gs * Delta(mla,mlb) * ( Delta(msb,-0.5)*Delta(msa,0.5) +
                                        Delta(msb,0.5)*Delta(msa,-0.5) ));
    sum=sum*myb*Bx*tesla*0.5;
    sumr=sumr+sum;
    
    /* By */
    
    sumi=Delta(msa,msb) * ( Delta(mla,mlb-1) * sqrt( lb*(lb+1)-mlb*(mlb-1) ) -
                            Delta(mla,mlb+1) * sqrt( lb*(lb+1)-mlb*(mlb+1) ) );
    sumi=sumi+( gs * Delta(mla,mlb) * ( Delta(msb,0.5)*Delta(msa,-0.5) -
                                        Delta(msb,-0.5)*Delta(msb,0.5) ));
    sumi=sumi*myb*By*tesla*0.5;
    
    RE_H=sumr;
    IM_H=sumi;

    return;

}

void Matrix_element_Hcf(int method, int i, int j)
{

    int k;
    double sumr=0.0, sumi=0.0;


    if(method==0)
    {
        for(k=0; k<electrons; k++)
        {
            Crystal_Field(states[i][k].l,
                          states[i][k].ml,
                          states[i][k].ms,
                          states_tmp[j][k].l,
                          states_tmp[j][k].ml,
                          states_tmp[j][k].ms);
            sumr=sumr+RE_CFP;
            sumi=sumi+IM_CFP;
        }

        RE_Matrix_element_Hcf=sumr*pow(-1.0,(double)signe);
        IM_Matrix_element_Hcf=sumi*pow(-1.0,(double)signe);

        return;

    }

    else if(method==1)
    {


    Crystal_Field(states[i][un_paired[1]].l,
                  states[i][un_paired[1]].ml,
                  states[i][un_paired[1]].ms,
                  states_tmp[j][un_paired[1]].l,
                  states_tmp[j][un_paired[1]].ml,
                  states_tmp[j][un_paired[1]].ms);

    RE_Matrix_element_Hcf=RE_CFP*pow(-1.0,(double)signe);
    IM_Matrix_element_Hcf=IM_CFP*pow(-1.0,(double)signe);

    return;

    }
    
    else
    {
        RE_Matrix_element_Hcf=0.0;
        IM_Matrix_element_Hcf=0.0; 
        return;
    }

}



/* obs! Complex */
void Crystal_Field(double la,double mla, double msa, double lb,double mlb, double msb)
{

    int t,p;
    double trej1, trej2;
    double sumr=0.0, sumi=0.0;

    if(la!=lb) 
    {
        printf("Error in Crystal_Field!\n");
        exit(0);
    }


    if(msa!=msb) 
    {
        RE_CFP=0.0; 
        IM_CFP=0.0;
        return;
    }

    for (t=2; t<=2*(int)la; t=t+2)
    {
    
    for (p=-t;p<=t; p++)
        {
            trej1 = clebsg( la, (double)t, la, 0, 0, 0, 0, 0, 0, 1);
            trej2 = clebsg( la, (double)t, la, -mla, (double)p, mlb, 0, 0, 0, 1);

            sumr = sumr + RE_Atp(t,p) * Rt(t) * pow(-1.0, mla) * (2.0*la+1.0) 
                        * trej1
                        * trej2;

            sumi = sumi + IM_Atp(t,p) * Rt(t) * pow(-1.0, mla) * (2.0*la+1.0)
                        * trej1
                        * trej2;

        }
    }


    RE_CFP=sumr;
    IM_CFP=sumi;

    return;

}

double RE_Atp(int t,int p)
{

    return(Atp[t][p+7]);

}



double IM_Atp(int t,int p)
{

    return(Btp[t][p+7]);

}


void readatp(void)
{
    FILE *pekatp;
    int dummy1,dummy2;
    int t,p;
    double real,imag;
    extern double Atp[20][20], Btp[20][20];

    if((pekatp=fopen("cfp.inp","r"))==NULL) {puts("Error in cfp.inp!"); exit(0);}


    for(t=0;t<=7; t++)
    {
    for(p=-t;p<=t;p++)
        {
        if(fscanf( pekatp, "%d%d%lf%lf",&dummy1,&dummy2,&real,&imag))
        {

            /* pseudo cm-1, will be upon mult. with <r**t> */
            Atp[t][p+7]=(1.0-Sigmat(t))*real;
            Btp[t][p+7]=(1.0-Sigmat(t))*imag;

        }
        else 
        {
            puts("End of file.");
            exit(0);}
        }
    }
    fclose(pekatp);

    if((dummy1!=7)||(dummy2!=7))
    {
        printf("End of file or error in cfp.inp!!!\n");
        exit(0);
    }

}




double Rt(int t)
{
if(t==2) return(R2);
else
if(t==4) return(R4);
else
if(t==6) return(R6);
else {printf("Error in Rt!\n"); exit(0);}

}


double Sigmat(int t)
{

if(t==2) return(S2);
else
if(t==4) return(S4);
else
if(t==6) return(S6);
else return(0.0);

}


double Matrix_element_SO(int method, int i, int j)
{

    int k;
    double sum=0.0;


    if(method==0)
    {
        for(k=0; k<electrons; k++)
            sum=(sum 
                + Spin_Orbit(states[i][k].l,
                            states[i][k].ml,
                            states[i][k].ms,
                            states_tmp[j][k].l,
                            states_tmp[j][k].ml, 
                            states_tmp[j][k].ms));

        return(sum*pow(-1.0,(double)signe));

    }

    else if(method==1)
    {
        return(pow(-1.0,(double)signe) 
            * Spin_Orbit(states[i][un_paired[1]].l,
                         states[i][un_paired[1]].ml,
                         states[i][un_paired[1]].ms,
                         states_tmp[j][un_paired[1]].l,
                         states_tmp[j][un_paired[1]].ml,
                         states_tmp[j][un_paired[1]].ms));
    }
    
    else return(0.0);

}



double Spin_Orbit_Xi(void)
{

    return(Xi);
}



double Spin_Orbit(double la,double mla, double msa, double lb,double mlb, double msb)
{
    double offdiag, diag;

    offdiag = (lx(la, mla, lb, mlb) * sx(msa, msb) 
             - ly(la, mla, lb, mlb)*sy(msa, msb));

    /* negative sign since i**2 = -1 */

    diag = mlb * msb * Delta(mla,mlb) * Delta(msa,msb);

    return( Spin_Orbit_Xi() * ( offdiag + diag) );

}


double lx(double la, double mla, double lb, double mlb)
{

    if(la!=lb)
    {
        printf("Error in lx!\n");
        exit(0);
    }

    if(mla==mlb+1.0)
        return(0.5*sqrt(la*(la+1.0) - mlb*(mlb+1.0)));
    else if(mla==mlb-1.0)
        return(0.5*sqrt(la*(la+1.0) - mlb*(mlb-1.0)));
    else return(0.0);

}

double ly(double la, double mla, double lb, double mlb)
{

  /* obs! this is the real ly without i */

    if(la!=lb)
    {
        printf("Error in ly!\n");
        exit(0);
    }

    if(mla==mlb+1.0)
        return(-0.5*sqrt(la*(la+1.0) - mlb*(mlb+1.0)));
    else if(mla==mlb-1.0)
        return(0.5*sqrt(la*(la+1.0) - mlb*(mlb-1.0)));
    else return(0.0);

}

double sx(double msa, double msb)
{

    if(msa==msb+1.0)
        return(0.5*sqrt(0.75 - msb*(msb+1.0)));
    else if(msa==msb-1.0)
        return(0.5*sqrt(0.75 - msb*(msb-1.0)));
    else return(0.0);

}

double sy(double msa, double msb)
{

    if(msa==msb+1.0)
        return(-0.5*sqrt(0.75 - msb*(msb+1.0)));
    else if(msa==msb-1.0)
        return(0.5*sqrt(0.75 - msb*(msb-1.0)));
    else return(0.0);

}

int diag1(integer n)
{
    char jobz = 'V';
    char uplo = 'U';
    integer lwork=2*n, lrwork, liwork;
    integer *iwork;
    integer ldz=n, m, info, size;
    doublecomplex *AP, *z, *work;
    doublereal *w, *rwork;
    integer i;
    gzFile pek, pek3;
    FILE *pek2;
    char buf[1], line[60], flyt[30];
    int gn,a,b,gi,gj,gk;
    double dummy;

    dummy = (((double)n*((double)n+1.0))/2.0);
    size = (integer)dummy;

    dummy=log((double)(n))/log(2.0) +1.0;

    lrwork=1+4*n+2*n*(integer)(dummy)+3*n*n;

    liwork=3+5*n;

    if(!(AP = malloc( size*sizeof(doublecomplex))))    
    {
        printf("No mem_ap!");
        free(AP);
        return(-1);
    }

    if(!(z = malloc( ldz*n*sizeof(doublecomplex))))    
    {
        printf("No mem_z!");
        free(AP);
        free(z);
        return(-1);
    }

    if(!(work = malloc(lwork*sizeof(doublecomplex)))) 
    {
        printf("No mem_work!");
        free(AP);
        free(z);
        free(work);
        return(-1);
    }

    if(!(w = malloc(n*sizeof(doublereal))))           
    {
        printf("No mem_w!");
        free(AP);
        free(z);
        free(work);
        free(w);
        return(-1);
    }

    if(!(rwork = malloc(lrwork*sizeof(doublereal)))) 
    {
        printf("No mem_rwork!");
        free(AP);
        free(z);
        free(work);
        free(w);
        free(rwork);
        return(-1);
    }

    if(!(iwork = malloc(liwork*sizeof(integer)))) 
    {
        printf("No mem_iwork!");
        free(AP);
        free(z);
        free(work);
        free(w);
        free(rwork);
        free(iwork);
        return(-1);
    }

    if((pek=gzopen("matrix.gz","rb"))==NULL) 
    {
        puts("Error in matrix.gz!");
        exit(0);
    }

    if((pek3=gzopen("eigenfunctions.gz","wb"))==NULL)
    {
        puts("Error in eigenfunctions.gz!");
        exit(0);
    }

    if((pek2=fopen("energies","w"))==NULL)
    {
        puts("Error in energies!");
        exit(0);
    }


    /* ------------ read the gzipped matrix -------------- */

    for(i=0; i < (int)((n*n+n)*0.5); i++)
    {
        for(gj=0; gj<60; gj++)
            line[gj]=' ';

        gi=-1;
        while((buf[0]=gzgetc(pek))!='\n')
        {
            gi++;
            if(gi>59)
            {
                printf("Index problem in diag!\n");
                exit(0);
            }
            line[gi]=buf[0];
        }

        gj=0;
        for(gn=0; gn<2; gn++)
        {
            for(gk=0; gk<30; gk++)
                flyt[gk]=' ';

            gi=-1;
            while( ((int)line[gj]) != 32 )
            {
                gi++;
                if( (gi>29) || (gj> 59) )
                    {
                        printf("Index problem in diag!\n");
                        exit(0);
                    }
                flyt[gi]=line[gj];  
                gj++;
            }
            if(gn==0)
                AP[i].r=atof(flyt);
            else if(gn==1)
                AP[i].i=atof(flyt);
            
            gj++;
        }
    }

    /* ------------ end read of gzipped matrix -------------- */

    zhpevd_(&jobz, &uplo, &n, AP, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    if(info!=0)
    {
        printf("problems with diagonalize! info = %ld\n",info);
        exit(0);
    }


    for(i=0; i < n*n; i++)
        gzprintf(pek3,"%le %le\n",z[i].r,z[i].i);

    printf("Saving eigenvalues ...\n");

    for(i=0; i < n; i++)
        fprintf(pek2,"%lf\n",w[i]);

    gzclose(pek);
    fclose(pek2);
    gzclose(pek3);
    free(AP);
    free(z);
    free(work);
    free(w);
    free(rwork);
    free(iwork);

    return(0);
}

void diag2(integer n)
{

    char jobz = 'V';
    char uplo = 'U';
    integer lwork=2*n-1;
    integer ldz=n, m, info, size;
    doublecomplex *AP, *z, *work;
    doublereal *w, *rwork;
    integer i;
    gzFile pek, pek3;
    FILE *pek2;
    char buf[1], line[60], flyt[30];
    int gn,a,b,gi,gj,gk;
    double dummy;

    dummy = (((double)n*((double)n+1.0))/2.0);
    size = (integer)dummy;

    if(!(AP = malloc( size*sizeof(doublecomplex))))
    {
        printf("No mem_ap!");
        exit(0);
    }

    if(!(z = malloc( ldz*n*sizeof(doublecomplex))))
    {
        printf("No mem_z!");
        exit(0);
    }

    if(!(work = malloc(lwork*sizeof(doublecomplex))))
    {
        printf("No mem_work!");
        exit(0);
    }
    if(!(w = malloc(n*sizeof(doublereal))))
    {
        printf("No mem_w!");
        exit(0);
    }
    if(!(rwork = malloc((3*n-2)*sizeof(doublereal))))
    {
        printf("No mem_rwork!");
        exit(0);
    }

    if((pek=gzopen("matrix.gz","rb"))==NULL)
    {
        puts("Error in matrix.gz!");
        exit(0);
    }

    if((pek3=gzopen("eigenfunctions.gz","wb"))==NULL)
    {
        puts("Error in eigenfunctions.gz!");
        exit(0);
    }

    if((pek2=fopen("energies","w"))==NULL)
    {
        puts("Error in energies!");
        exit(0);
    }


    /* ------------ read the gzipped matrix -------------- */

    for(i=0; i < (int)((n*n+n)*0.5); i++)
    {
        for(gj=0; gj<60; gj++)
            line[gj]=' ';

        gi=-1;
        while((buf[0]=gzgetc(pek))!='\n')
        {
            gi++;

            if(gi>59)
            {
                printf("Index problem in diag!\n");
                exit(0);
            }

            line[gi]=buf[0];
        }

        gj=0;
        for(gn=0; gn<2; gn++)
        {

            for(gk=0; gk<30; gk++)
                flyt[gk]=' ';

            gi=-1;
            while( ((int)line[gj]) != 32 )
            {
                gi++;

                if( (gi>29) || (gj> 59) )
                {
                    printf("Index problem in diag!\n");
                    exit(0);
                }

                flyt[gi]=line[gj];
                gj++;
            }

            if(gn==0)
                AP[i].r=atof(flyt);
            else if(gn==1)
                AP[i].i=atof(flyt);
            
            gj++;

        }

    }

    /* ------------ end read of gzipped matrix -------------- */


    zhpev_(&jobz, &uplo, &n, AP, w, z, &ldz, work, rwork, &info);
    if(info!=0)
    {
        printf("problems with diagonalize! info = %ld\n",info);
        exit(0);
    }


    for(i=0; i < n*n; i++)
        gzprintf(pek3,"%le %le\n",z[i].r,z[i].i);

    printf("Eigenvalues:\n");

    for(i=0; i < n; i++)
        fprintf(pek2,"%lf\n",w[i]);

    gzclose(pek);
    fclose(pek2);
    gzclose(pek3);
    free(AP);
    free(z);
    free(work);
    free(w);
    free(rwork);

}

void load_pauli(void)
{

    int i,j;
    FILE *pek;

    if((pek=fopen("pauli.inp","r"))==NULL)
    {
        puts("Error in pauli.inp!");
        exit(0);
    }

    for(i=0; i<STATES; i++)
    { 
        for(j=0; j<electrons; j++)
        {
        fscanf(pek,
            "%lf%lf%lf",
            &states[i][j].l,
            &states[i][j].ml,
            &states[i][j].ms);
        }
    }

    fclose(pek);
}


double Matrix_element(int method, int i, int j)
{

    int k,l;
    double sum=0.0;

    if(method==0)
    {
        for(k=0; k<electrons-1; k++)
        {
            for(l=k+1; l<electrons; l++)
            {
                sum=sum + General(states[i][l].l,
                                states[i][l].ml,
                                states[i][l].ms,
                                states[i][k].l,
                                states[i][k].ml,
                                states[i][k].ms,
                                states[i][l].l,
                                states[i][l].ml,
                                states[i][l].ms,
                                states[i][k].l,
                                states[i][k].ml,
                                states[i][k].ms) 
                        - General(states[i][l].l,
                                states[i][l].ml,
                                states[i][l].ms,
                                states[i][k].l,
                                states[i][k].ml,
                                states[i][k].ms,
                                states[i][k].l,
                                states[i][k].ml,
                                states[i][k].ms,
                                states[i][l].l,
                                states[i][l].ml,
                                states[i][l].ms);

            }
        }

        return(sum*pow(-1.0,(double)signe));

    }

    else if(method==1)
    {

    for(k=0; k<electrons; k++)
        sum = sum + General(states[i][k].l,
                            states[i][k].ml,
                            states[i][k].ms,
                            states[i][un_paired[1]].l,
                            states[i][un_paired[1]].ml,
                            states[i][un_paired[1]].ms,
                            states[i][k].l,
                            states[i][k].ml,
                            states[i][k].ms,
                            states_tmp[j][un_paired[1]].l,
                            states_tmp[j][un_paired[1]].ml,
                            states_tmp[j][un_paired[1]].ms)                  
                  - General(states[i][k].l,
                            states[i][k].ml,
                            states[i][k].ms, 
                            states[i][un_paired[1]].l,
                            states[i][un_paired[1]].ml,
                            states[i][un_paired[1]].ms,
                            states_tmp[j][un_paired[1]].l,
                            states_tmp[j][un_paired[1]].ml,
                            states_tmp[j][un_paired[1]].ms,
                            states[i][k].l,
                            states[i][k].ml,
                            states[i][k].ms);

    return(sum*pow(-1.0,(double)signe));

    }

    else if(method==2)
    {
        return((General(states[i][un_paired[1]].l,
                        states[i][un_paired[1]].ml,
                        states[i][un_paired[1]].ms,
                        states[i][un_paired[2]].l,
                        states[i][un_paired[2]].ml,
                        states[i][un_paired[2]].ms,
                        states_tmp[j][un_paired[1]].l,
                        states_tmp[j][un_paired[1]].ml,
                        states_tmp[j][un_paired[1]].ms,
                        states_tmp[j][un_paired[2]].l,
                        states_tmp[j][un_paired[2]].ml,
                        states_tmp[j][un_paired[2]].ms)
                - General(states[i][un_paired[1]].l,
                          states[i][un_paired[1]].ml,
                          states[i][un_paired[1]].ms,
                          states[i][un_paired[2]].l,
                          states[i][un_paired[2]].ml,
                          states[i][un_paired[2]].ms,
                          states_tmp[j][un_paired[2]].l,
                          states_tmp[j][un_paired[2]].ml,
                          states_tmp[j][un_paired[2]].ms,
                          states_tmp[j][un_paired[1]].l,
                          states_tmp[j][un_paired[1]].ml,
                          states_tmp[j][un_paired[1]].ms)
                ) * pow(-1.0,(double)signe));

    }

    else return(0.0);

}


double General(double la,double mla, double msa, double lb,double mlb, double 
msb, double lc,double mlc, double msc, double ld,double mld, double msd)
{
    double sum=0;
    int k;

    if(Delta(msa,msc)==0) return(0.0);
    if(Delta(msb,msd)==0) return(0.0);
    if(Delta(mla+mlb, mlc+mld)==0) return(0.0);

    /* k=0 is excluded since it's a central contribution equal for all levels */

    for(k=2; k<=(int)(la+lc); k=k+2)
        sum=sum+ck((double)k,la,mla,lc,mlc)
                  * ck((double)k,ld,mld,lb,mlb)
                  * Rk((double)k);

    return(sum);
}



double Rk(double k)
{
    if(k==2)
        return(F2);
    else if(k==4)
        return(F4);
    else if(k==6)
        return(F6);
    else
    {
        printf("Problem in Rk!\n");
        exit(0);
    }
}


double ck(double k,double la,double mla, double lb,double mlb)
{
    double trej1,trej2;

    trej1 = clebsg( la, k, lb, 0, 0, 0, 0, 0, 0, 1);
    trej2 = clebsg( la, k, lb, -mla, mla-mlb, mlb, 0, 0, 0, 1);

    return(pow(-1,mla)
          *sqrt((2*la+1)*(2*lb+1))
          *trej1
          *trej2);

}



double Delta(double a, double b)
{
    if(a==b) return(1.0);
    else return(0.0);
}


/* Calculation of 3j- and 6j-symbols */

double clebsg(double a, double b, double c, double xx, double yy, double zz, double gg, double hh, double pp,int mode)
{

    int i, i_,iay[4],ij, ijpar, iyy, iz, jq, js, jspar,
        k, k1, k2, k3, k4, k5, k6, ka, kb, kc, key, keyrac, keytri, kk1,
        kk2, kk3, kk4, kk5, kk6, kk7, kk8, kk9, kup, m1, m10, m2, m3,
        m4, m5, m6, m7, m8, m9, ma, mb, mc, md, mm1, mm2, mm3, mm4, mm5,
        n4, n5, n5par;
    double anine,ay[4],clebsg_v, clebsh, factor, ra,
        racah, rb, x, y, z;

    static int ix,iy,j[101];
    static double h[101];
    static int jjj = 0;
    static double signz = 1.;


        /*
        *     WIGN3J- WIGNER 3-J SYMBOL
        *     CLEBSG- CLEBSCH-GORDAN COEFFICIENT
        *     WIGN6J- WIGNER 6-J SYMBOL
        *     RACAHC- RACAH COEFFICIENT
        *     JAHNUF- U-FUNCTION (JAHN)
        *     WIGN9J- WIGNER 9-J SYMBOL
        */
    #define INTPTF(q)       (int)((q) + (q) + sign( .10e0, (q) ))
    #define IPARF(i)        (int)(4*((i)/4) - (i) + 1)

        if(mode==1) goto wign3j;
        else
        if(mode==2) goto wign6j;
        else
        {
        printf("ERROR!!!");
        exit(1);
        }

    wign3j:
            if(xx + yy +zz != 0.0) return(0.0);  /* speed up of 3j symbols */
        fgercm.ierr = 0;
        /*     CALL NOARG(NARG)
        *     IF(NARG.EQ.3) GOTO 2 */
        signz = -1.;
        key = 1;
        goto L_1;
        /*   2 KEY=3
        *     GOTO 1
        */


    wign6j:
        key = 11;
        fgercm.ierr = 0;
        goto L_100;




    L_1:
        k1 = INTPTF( a );
        k2 = INTPTF( b );
        k3 = INTPTF( c );
        if( key >= 3 )
            goto L_100;
        k4 = INTPTF( xx );
        k5 = INTPTF( yy );
        k6 = INTPTF( zz*signz );

    L_100:
        if( jjj != 0 )
            goto L_500;
        jjj = 1;
        fgercm.ierct = 0;
        h[0] = 1.0;
        j[0] = 0;
        x = 0.;
        for( i = 2; i <= 101; i++ ){
            i_ = i - 1;
            x = x + 1.0;
            h[i_] = h[i_ - 1]*x;
            j[i_] = j[i_ - 1];
    L_200:
            if( h[i_] < 10.0 )
                goto L_400;
            h[i_] = 0.01*h[i_];
            j[i_] = j[i_] + 2;
            goto L_200;
    L_400:
            ;
            }

    L_500:
        if( key < -5 )
            goto L_750;
            if( key >= 3 )
                    goto L_320;
            if( (k4 + k5 - k6) != 0 )
                    goto L_710;
            m1 = k1 + k2 - k3;
        m2 = k2 + k3 - k1;
            m3 = k3 + k1 - k2;
            m4 = k1 + k4;
        m5 = k1 - k4;
            m6 = k2 + k5;
            m7 = k2 - k5;
            m8 = k3 + k6;
            m9 = k3 - k6;
            m10 = k1 + k2 + k3 + 2;

            if( m1 < 0 )
                    goto L_710;
            if( m2 < 0 )
            goto L_710;
        if( m3 < 0 )
                    goto L_710;
            if( m4 < 0 )
                    goto L_710;
            if( m5 < 0 )
                    goto L_710;
        if( m6 < 0 )
                    goto L_710;
            if( m7 < 0 )
            goto L_710;
            if( m8 < 0 )
                    goto L_710;
            if( m9 < 0 )
                    goto L_710;
            if( (m4 - (m4/2) - (m4/2)) != 0 )
                    goto L_710;
            if( (m6 - (m6/2) - (m6/2)) != 0 )
                    goto L_710;
            if( (m8 - (m8/2) - (m8/2)) != 0 )
            goto L_710;
        if( (m10 - (m10/2) - (m10/2)) != 0 )
                    goto L_710;

            y = k3 + 1;
            m1 = m1/2 + 1;
            m2 = m2/2 + 1;
        m3 = m3/2 + 1;
            m4 = m4/2 + 1;
            m5 = m5/2 + 1;
        m6 = m6/2 + 1;
            m7 = m7/2 + 1;
            m8 = m8/2 + 1;
            m9 = m9/2 + 1;
            m10 = m10/2 + 1;

            y = sqrt( y*h[m1 - 1]*h[m2 - 1]*h[m3 - 1]*h[m4 - 1]*h[m5 - 1]*
            h[m6 - 1]*h[m7 - 1]*h[m8 - 1]*h[m9 - 1]/h[m10 - 1] );
            iy = j[m1 - 1] + j[m2 - 1] + j[m3 - 1] + j[m4 - 1] + j[m5 - 1] +
            j[m6 - 1] + j[m7 - 1] + j[m8 - 1] + j[m9 - 1] - j[m10 - 1];
        iy=(int)((float)(iy)/2.0);

            n4 = m1;
            if( n4 > m5 )
                    n4 = m5;
            if( n4 > m6 )
                    n4 = m6;
        n4 = n4 - 1;
            m2 = k2 - k3 - k4;
            m3 = k1 + k5 - k3;
        n5 = 0;
            if( n5 < m2 )
                    n5 = m2;
            if( n5 < m3 )
                    n5 = m3;
            n5par = IPARF( n5 );
            n5 = n5/2;
            z = 0.0;
            goto L_610;

    L_700:
        mm1 = m1 - n5;
            mm2 = m5 - n5;
            mm3 = m6 - n5;
            mm4 = n5 - (m2/2) + 1;
            mm5 = n5 - (m3/2) + 1;

        x = 1./(h[mm1 - 1]*h[mm2 - 1]*h[mm3 - 1]*h[mm4 - 1]*h[mm5 - 1]*
                h[n5]);

        ix = -j[mm1 - 1] - j[mm2 - 1] - j[mm3 - 1] - j[mm4 - 1] - j[mm5 - 1] -
            j[n5];

    L_800:
            switch( IARITHIF(ix + iy) ){
                    case -1: goto L_900;
                    case  0: goto L_210;
                    case  1: goto L_110;
                    }
    L_900:
        x = 0.1*x;
        ix = ix + 1;
            goto L_800;
    L_110:
            x = 10.0*x;
            ix = ix - 1;
            goto L_800;

    L_210:
            if( n5par < 0 )
            x = -x;
            z = z + x;
            n5par = -n5par;
            n5 = n5 + 1;

    L_610:
            switch( IARITHIF(n5 - n4) ){
                    case -1: goto L_700;
                    case  0: goto L_700;
                    case  1: goto L_810;
            }

    L_710:
            clebsh = 0.0;
            fgercm.ierr = 1;
            fgercm.ierct = fgercm.ierct + 1;
            goto L_220;

    L_810:
            clebsh = z*y;
        switch( key ){
                    case 1: goto L_120;
                    case 2: goto L_220;
                    }

    L_220:
            clebsg_v = clebsh;
            return( clebsg_v );

    L_120:
        js = k1 - k2 + k6;
        if( js < 0 )
                    js = -js;
            jspar = IPARF( js );
        clebsg_v = jspar*clebsh/sqrt( k3 + 1.0 );
        signz = 1.;
            return( clebsg_v );

    L_320:
            if( key >= 10 )
            goto L_130;
            key = key - 2;
            if( (k1 - (k1/2) - (k1/2)) != 0 )
                    goto L_420;
            if( (k2 - (k2/2) - (k2/2)) != 0 )
                    goto L_420;
            if( (k3 - (k3/2) - (k3/2)) != 0 )
                    goto L_420;
            ij = k1 + k2 + k3;
            ijpar = IPARF( ij );
        if( ijpar <= 0 )
            goto L_420;
            m1 = ij - k1 - k1;
            m2 = ij - k2 - k2;
            m3 = ij - k3 - k3;
            m4 = ij + 2;
            if( m1 < 0 )
            goto L_420;
            if( m2 < 0 )
                    goto L_420;
        if( m3 < 0 )
                    goto L_420;
            m1 = m1/2 + 1;
            m2 = m2/2 + 1;
            m3 = m3/2 + 1;
            m4 = ij/2 + 2;
            y = sqrt( h[m1 - 1]*h[m2 - 1]*h[m3 - 1]/h[m4 - 1] );
            iy = (j[m1 - 1] + j[m2 - 1] + j[m3 - 1] - j[m4 - 1])/2;
            ij = ij/2;
            ijpar = IPARF( ij );
        ij = ij/2 + 1;
        m1 = m1/2 + 1;
            m2 = m2/2 + 1;
            m3 = m3/2 + 1;
            z = h[ij - 1]/(h[m1 - 1]*h[m2 - 1]*h[m3 - 1]);
            iz = j[ij - 1] - j[m1 - 1] - j[m2 - 1] - j[m3 - 1];
            iz = iz + iy;
        clebsh = ijpar*y*z*pow(10.0, iz);
            switch( key ){
                    case 1: goto L_220;
            case 2: goto L_720;
                    }

    L_720:
            jq = k2 - k1;
            if( jq < 0 )
                    jq = -jq;
            ijpar = IPARF( jq );
            clebsg_v = clebsh*ijpar*sqrt( k3 + 1.0 );
            return( clebsg_v );

    L_420:
            clebsh = 0.0;
            fgercm.ierr = 1;
            fgercm.ierct = fgercm.ierct + 1;
            switch( key ){
                    case 1: goto L_220;
            case 2: goto L_720;
                    }

    L_130:
            if( key == 11 )
                    goto L_450;
            if( key > 19 )
                    goto L_750;
        k1 = INTPTF( a );
        k2 = INTPTF( b );
        k3 = INTPTF( yy );
        k4 = INTPTF( xx );
        k5 = INTPTF( c );
        k6 = INTPTF( zz );

    L_750:
            ka = k1;
            kb = k2;
            kc = k3;
            keytri = 1;
        goto L_630;

    L_230:
        ka = k4;
            kb = k5;
            keytri = 2;
            goto L_630;

    L_330:
            kb = k2;
            kc = k6;
            keytri = 3;
            goto L_630;

    L_430:
            ka = k1;
            kb = k5;
            keytri = 4;
            goto L_630;

    L_530:
        y = ay[0]*ay[1]*ay[2]*ay[3];
            iyy = iay[0] + iay[1] + iay[2] + iay[3];
        m1 = (k1 + k2 + k4 + k5)/2 + 2;
            m2 = (k1 + k2 - k3)/2 + 1;
            m3 = (k4 + k5 - k3)/2 + 1;
            m4 = (k1 + k5 - k6)/2 + 1;
            m5 = (k2 + k4 - k6)/2 + 1;
            m6 = k1 + k4 - k3 - k6;
            m7 = k2 + k5 - k3 - k6;

            n4 = m1;
            if( n4 > m2 )
            n4 = m2;
        if( n4 > m3 )
                    n4 = m3;
            if( n4 > m4 )
                    n4 = m4;
            if( n4 > m5 )
                    n4 = m5;
        n4 = n4 - 1;
            n5 = 0;
            if( n5 < m6 )
            n5 = m6;
            if( n5 < m7 )
                    n5 = m7;
            n5par = IPARF( n5 );
            n5 = n5/2;
            m6 = m6/2 - 1;
            m7 = m7/2 - 1;
            z = 0.0;
            goto L_730;

    L_140:
        x = h[m1 - n5 - 1]/(h[n5]*h[m2 - n5 - 1]*h[m3 - n5 - 1]*h[m4 - n5 - 1]*
            h[m5 - n5 - 1]*h[n5 - m6 - 1]*h[n5 - m7 - 1]);
            ix = j[m1 - n5 - 1] - j[n5] - j[m2 - n5 - 1] - j[m3 - n5 - 1] -
            j[m4 - n5 - 1] - j[m5 - n5 - 1] - j[n5 - m6 - 1] - j[n5 - m7 - 1];
    L_240:
            switch( IARITHIF(ix + iyy) ){
            case -1: goto L_340;
                    case  0: goto L_440;
                    case  1: goto L_540;
            }
    L_340:
            x = 0.1*x;
            ix = ix + 1;
            goto L_240;
    L_540:
            x = 10.0*x;
            ix = ix - 1;
            goto L_240;
    L_440:
        if( n5par < 0 )
            x = -x;
            z = z + x;
            n5par = -n5par;
            n5 = n5 + 1;

    L_730:
        if( n5 <= n4 )
                    goto L_140;

        racah = z*y;
    L_840:
            if( key < -5 )
                    goto L_160;
            key = key - 10;
            switch( key ){
                    case 1: goto L_150;
                    case 2: goto L_250;
                    case 3: goto L_350;
                    }

    L_830:
            racah = 0.0;
            fgercm.ierr = 1;
            fgercm.ierct = fgercm.ierct + 1;
            goto L_840;

    L_150:
            ijpar = IPARF( k1 + k2 + k4 + k5 );
            if( ijpar < 0 )
            racah = -racah;
    L_250:
            clebsg_v = racah;
            return( clebsg_v );

    L_350:
            factor = sqrt( (k3 + 1.0)*(k6 + 1) );
            clebsg_v = factor*racah;
            return( clebsg_v );
    L_450:
        k1 = INTPTF( a );
        k2 = INTPTF( b );
        k3 = INTPTF( c );
        k4 = INTPTF( xx );
        k5 = INTPTF( yy );
        k6 = INTPTF( zz );
            goto L_750;

            /*     TRIANGLE FUNCTION
            */
    L_630:
            ma = ka + kb - kc;
            mb = ka - kb + kc;
            mc = -ka + kb + kc;
            md = ka + kb + kc + 2;
            if( ma < 0 )
                    goto L_830;
            if( mb < 0 )
                    goto L_830;
            if( mc < 0 )
            goto L_830;
        if( (md - (md/2) - (md/2)) != 0 )
                    goto L_830;
            ma = ma/2 + 1;
            mb = mb/2 + 1;
            mc = mc/2 + 1;
            md = md/2 + 1;
        ay[keytri - 1] = sqrt( h[ma - 1]*h[mb - 1]*h[mc - 1]/h[md - 1] );


        iay[keytri - 1] = j[ma - 1] + j[mb - 1] + j[mc - 1] - j[md - 1];

            iay[keytri - 1]=(int)((float)(iay[keytri - 1])/2.0);


            switch( keytri ){
                    case 1: goto L_230;
                    case 2: goto L_330;
                    case 3: goto L_430;
                    case 4: goto L_530;
            }


            key = -10;
            fgercm.ierr = 0;

        kk1 = INTPTF( a );
        kk2 = INTPTF( b );
        kk3 = INTPTF( c );
        kk4 = INTPTF( xx );
        kk5 = INTPTF( yy );
        kk6 = INTPTF( zz );
        kk7 = INTPTF( gg );
        kk8 = INTPTF( hh );
        kk9 = INTPTF( pp );

        kup = kk1 + kk9;
        m1 = kk4 + kk8;
        m2 = kk2 + kk6;
        if( kup > m1 )
            kup = m1;
        if( kup > m2 )
            kup = m2;

        k = kk1 - kk9;
        if( k < 0 )
                    k = -k;
            m1 = kk4 - kk8;
            if( m1 < 0 )
                    m1 = -m1;
        m2 = kk2 - kk6;
            if( m2 < 0 )
                    m2 = -m2;
            if( k < m1 )
                    k = m1;
            if( k < m2 )
                    k = m2;

            anine = 0.0;

    L_660:
        if( k > kup )
                    goto L_260;
            k1 = kk1;
            k2 = kk4;
        k3 = kk7;
            k4 = kk8;
            k5 = kk9;
            k6 = k;
            keyrac = 1;
        goto L_100;

    L_160:
            switch( keyrac ){
                    case 1: goto L_360;
                    case 2: goto L_460;
                    case 3: goto L_560;
                    }

    L_360:
        ra = racah;
        k1 = kk2;
            k2 = kk8;
            k3 = kk5;
            k4 = kk4;
        k5 = kk6;
            keyrac = 2;
            goto L_750;

    L_460:
        rb = racah;
            k1 = kk9;
            k2 = kk6;
            k3 = kk3;
            k4 = kk2;
            k5 = kk1;
            keyrac = 3;
            goto L_750;

    L_560:
        anine = anine + ra*rb*racah*(k + 1);
        k = k + 2;
            goto L_660;

    L_260:
        clebsg_v = anine;
            return( clebsg_v );
    #undef  IPARF
    #undef  INTPTF
} 


double sign(double a,double b)
{
    if(b>=0.0) return(fabs(a));
    else return(-fabs(a));
}



