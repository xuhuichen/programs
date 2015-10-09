/* post processing porgram for the compton2d Monte Carlo code to read the leaked
photon files and collect output the Flux v.s. Frequency information.
Edited by Xuhui Chen on 2/13/09 to allow more than one input files. */
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#define ch_max 200
#define st_max 1024
#define t_max 90

#define qu(x) ((x)*(x))

/* ------------------------------------------------------------ */
double max(x, y)
double x, y;
{
if (x > y) return (x);
else return (y);
}

/* ------------------------------------------------------------ */
int digit(c)
char c;
{
if (c == '0' || c == '1' || c == '2' || c == '3' || c == '4'
 || c == '5' || c == '6' || c == '7' || c == '8' || c == '9')
    return (1);
else return (0);
}

/* ------------------------------------------------------------ */
void add_dat(str)
char *str;
{
int len;

len = strlen(str);
if (str[len-4] != '.')
  { str[len] = '.';
    str[len+1] = 'd';
    str[len+2] = 'a';
    str[len+3] = 't';
    str[len+4] = '\0';
    }
return;
}
/*--------------------------------------------------------------*/
void file_append(str,substr)
char *str, *substr;
{
int len, sublen, i;
len= strlen(str); sublen=strlen(substr);
for(i=0;i<=sublen+4;i++){
  if(i<5){  str[len+sublen-i]=str[len-i];}
  else{str[len+sublen-i]=substr[sublen+4-i];}
}
return;
}
/* ------------------------------------------------------------ */
void dir_append(str,substr)
char *str, *substr;
{
int len, sublen, i;
len=strlen(str);sublen=strlen(substr);
for(i=0;i<sublen;i++){
str[len+i]=substr[i];}
str[len+sublen]='\0';
return;
}
/* ------------------------------------------------------------ */
double dinput(x)
double x;
{
char a[30];

gets(a);
if (strlen(a)) return (atof(a));
else return (x);
}

/* ------------------------------------------------------------ */
int iinput(x)
int x;
{
char a[30];

gets(a);
if (strlen(a)) return (atoi(a));
else return (x);
}

/* ------------------------------------------------------------ */
void file_error(fname)
char *fname;
{
printf ("\n\n File error: Can not open file '%s'! - Exit.\n\n", fname);
exit(0);
}

/* ------------------------------------------------------------ */
main()
{
double E0[ch_max], E1[ch_max], ph[ch_max];
double time, E, ew, dt, F[t_max][ch_max];
double F_ec[t_max][ch_max],F_ssc[t_max][ch_max],F_ssc2[t_max][ch_max];
double E_lower=1.e-7, E_upper=1.e10, dE;
double t0[t_max], t1[t_max];
double rmax=1.e16, t_bound, r, z, mu, mu_obs, mu0, mu1, mup0, mup1;
double phi, cdt, gam_bulk=33., doppler, betta;
int emitype;
int k, l, m, n, regions=1, reg, n_bins=0, n_r=100;
int n_t=30, n_file, factor=0, max_file=5120 ;
int i, dir;
char infile[30]="p001_evb.dat", outfile[30]="seds_30.dat", a[30], b[30];
char ecfile[30], sscfile[30], ssc2file[30], directory[30]="photons/", filepath[30];
FILE *fpin, *fpout, *fpec, *fpssc, *fpssc2;

/* particle_sum is the sum of the particles in a sertain energy range. 
   n_sum k1_sum and k2_sum are the n (angular bin) and k (energy bins) of the bins we
   are trying to find the number of particles for.
   Added by JF July 2003  */
int particle_sum[t_max][ch_max];
int n_sum, k1_sum, k2_sum;
char test;

printf ("\n\n Input from event file: [%s]  ", infile);
gets(a);
if (strlen(a)) {
   if (digit(a[0])) {
      for (k=0; k<strlen(a); k++) infile[2+k] = a[k];
      infile[2+k] = '\0';
   } else {
      strcpy(infile, a);
   } 
   add_dat(infile);
}

//if (!(fpin=fopen(infile, "rb"))) file_error(infile);

l = strlen(infile);

printf ("\n gam_bulk = [%e]  ", gam_bulk);
gam_bulk = dinput(gam_bulk);
betta   = sqrt(1.-1./(gam_bulk*gam_bulk));
printf ("\n rmax = [%e]  ", rmax);
rmax = dinput(rmax);

printf ("\n Output file: [%s]  ", outfile);
gets(a);
if (strlen(a)) {
   strcpy(outfile, a);
   add_dat(outfile);
}
   strcpy(ecfile,outfile);
   file_append(ecfile,"_ec");
   strcpy(sscfile,outfile);
   file_append(sscfile,"_ssc");
   strcpy(ssc2file,sscfile);
   file_append(ssc2file,"2"); 

printf ("\n Number of time bins: [%1d]  ", n_t);
n_t = iinput(n_t);

if (n_t > t_max)  {
   printf ("\n Number of angular bins set to maximum number = %d.\n", t_max);
   n_t = t_max;
}
t0[0]=1.6e4; t1[n_t]=6e4;
printf ("\n start from time: [%e]  ", t0[0]);
t0[0] = dinput(t0[0]);
printf ("\n end at time: [%e]  ", t1[n_t]);
t1[n_t] = dinput(t1[n_t]);
dt = (t1[n_t]-t0[0])/n_t;

for(n=0; n<n_t; n++) {
   t0[n] = t0[0] + n*dt;
   t1[n] = t0[0] + n*dt+dt;
}

mu0 = -0.12;
mu1 = 0.12;
printf ("\n mu_min = [%e]  ", mu0 );
mu0 = dinput(mu0);
printf (" mu_max = [%e]  ", mu1);
mu1 = dinput(mu1);
// get the angle in the observer's frame, to be used in the conversion to flux
mup0=(mu0-betta)/(1.-mu0*betta);
mup1=(mu1-betta)/(1.-mu1*betta);

m1:
printf ("\n Energy grid: \n");
printf ("\n Number of energy regions: [%d]  ", regions);
regions = iinput(regions);

for (reg=0; reg<regions; reg++) {
    printf ("   Region no. %d:", reg+1);

    printf ("\n   Lower boundary: [%e] keV  ", E_lower);
    E_lower = dinput(E_lower);
    printf ("   Upper boundary: [%e] keV  ", E_upper);
    E_upper = dinput(E_upper);

    printf ("   Number of energy bins in this region: [%d]  ", n_r);
    n_r = iinput(n_r);
    if( (n_r + n_bins) > ch_max ) {
       printf ("\n Not more than %d energy channels!", ch_max);
       printf ("\n Re-specify energy grid! \n");
       goto m1;
    }

    printf ("   Logarithmic (0) or linear (1) bins? [0]  ");
    gets(a);
    if (a[0] == '1') {
        dE = (E_upper - E_lower)/((double)(n_r));
        E0[n_bins] = E_lower;
        for (k=n_bins; k<(n_bins + n_r - 1); k++) {
           E1[k] = E0[k+1] = E0[k] + dE;
	}
        E1[n_bins + n_r - 1] = E0[n_bins + n_r - 1] + dE;
    } else {
        dE = exp(log(E_upper/E_lower)/((double)(n_r)));
        E0[n_bins] = E_lower;
        for (k=n_bins; k<(n_bins + n_r - 1); k++) {
           E1[k] = E0[k+1] = E0[k]*dE;
	}
        E1[n_bins + n_r - 1] = E0[n_bins + n_r - 1]*dE;
    }
    n_bins+=n_r;
    E_lower = E_upper; //GF :?
}

printf ("\n Show energy binning? [n]  ");
gets(a);
test = a[0];
// printf ("\n number of input files: [%1d]", max_file);
// max_file = iinput(max_file);
/*
if (a[0] == 'y')
  { printf ("\n Energy binning:\n");
    for (k=0; k<n_bins; k++)
       printf ("\n %d. E0 = %e keV, E1 = %e keV", k, E0[k], E1[k]);
    printf ("\n");
    }
printf ("\n");
*/

for (k=0; k<ch_max; k++) {
   for (n=0; n<t_max; n++) {
      F[n][k] = 0.0;
      F_ec[n][k] = 0.0;
      F_ssc[n][k] = 0.0;
      F_ssc2[n][k] = 0.0;
      particle_sum[n][k] = 0.0;
   }
}

/* Added by JF July 2003. */
/*particle_sum=0; */
n_sum  = 1;
k1_sum = 22;
k2_sum = 23;

m = 0;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
for(dir = 0;dir<10;dir++){
 if(dir==0) {
  directory[6]='s';
 }else{
  directory[6]='0'+dir;
 }

 n_file = 1;
 while(n_file <= max_file) {
   strcpy(filepath,directory);
   dir_append(filepath,infile);
   if( !(fpin=fopen(filepath, "rb")) ) {
	do {
	 infile[7]++;
         infile[1]='0';
         infile[2]='0';
         infile[3]='1';
	 strcpy(filepath,directory);
	 dir_append(filepath,infile);
	}while((!(fpin=fopen(filepath,"rb"))) && infile[7]<'l');

	if( infile[7]>='l' ){
	 infile[7]='b';
	 break;
	}
   }
   if( infile[1]=='0' && infile[2]=='0' && infile[3]=='1' && directory[6]=='s')factor++;
   //file_error(infile);
   printf ("\n Input from another event file [%s%s] \n",directory,infile);
  
//   while ((a[0]=fgetc(fpin)) != EOF) {
  
      /* it reads one photon at a time */
//      fseek(fpin, -1, 1);
//      fscanf(fpin, "%lf %lf %lf %lf %lf %lf %lf %i\n", &t_bound, &E, &ew, &r, &z, &mu, &phi, &emitype);
      while(!feof(fpin)){
            fread(&t_bound,sizeof(double),1,fpin);
            fread(&E,sizeof(double),1,fpin);
            fread(&ew,sizeof(double),1,fpin);
            fread(&r,sizeof(double),1,fpin);
            fread(&z,sizeof(double),1,fpin);
            fread(&mu,sizeof(double),1,fpin);
            fread(&phi,sizeof(double),1,fpin);
            fread(&emitype,sizeof(int),1,fpin);  
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // relativistic beaming
      //    if (z > zmax*(1.-1.e-7)) continue;          

      mu      = -mu;
      if (mu < mu0 || mu >= mu1) continue;
      mu = (mu0+mu1)/2.; // force the photons within the angle bin to have the same angle
      mu_obs  = (mu+betta)/(1.+mu*betta);
      doppler = 1./gam_bulk/(1.-mu_obs*betta);      //doppler = gam_bulk*(1.+mu*betta);
      t_bound = (t_bound-betta*z*3.33333333e-11)/doppler;
      E       = E*doppler;
      ew      = ew*doppler;
      cdt     = z*mu_obs/gam_bulk + sqrt(1.-mu_obs*mu_obs)*(rmax-r*cos(phi));
      time = t_bound + 3.33333333e-11*cdt;
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
      if (m++ < 10 ) {
          printf (" t = %e, E = %e, ew = %e\n", t_bound, E, ew);
          printf (" r = %e, z = %e, mu = %e, phi = %e\n", r, z, mu, phi);
      }
  
      /* cdt = (rmax - r*cos(phi))*sqrt(1. - qu(mu)) + (zmax - z)*mu; */
  

      

  
      /* if (m > 1000000 && m < 1000100) printf ("time = %e s\n", time); */
  
      for (k=0; k<n_bins; k++) {
         if (E >= E0[k] && E < E1[k]) break;
      }
  
      for (n=0; n<n_t; n++) {
         if (time >= t0[n] && time < t1[n]) break;
      }
  
      /* if (m > 1000000 && m < 1000100) 
          printf ("E-bin no. %d;  mu-bin no. %d.\n", k, n);  */
      
      if (k < n_bins && n < n_t) {
         F[n][k]+=ew;
	 if(emitype==11){F_ec[n][k]+=ew;}
         if(emitype==1){F_ssc[n][k]+=ew;}
         if(emitype==2){F_ssc2[n][k]+=ew;}
         /* if( (n==n_sum) && ( (k==k1_sum)||(k==k2_sum) )  ) */
         particle_sum[n][k]++;
      }
   }
   /* end of while on input photon file */
  
   fclose(fpin);
  
   if(infile[3]=='9') {
      if(infile[2]=='9') {
         infile[1]++;
         infile[2]='0';
      } else {
         infile[2]++;
      }
  
      infile[3]='0';
   } else {
      infile[3]++;
   }
  
   n_file++;
 } //different infiles
} //different directories
/* end of loop over input files */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Xuhui */

for (k=0; k<n_bins; k++) {
   for (n=0; n<n_t; n++) {
      F[n][k]/=(dt*(E1[k] - E0[k])*(mup1-mup0)/2.);
      F_ec[n][k]/=(dt*(E1[k] - E0[k])*(mup1-mup0)/2.);
      F_ssc[n][k]/=(dt*(E1[k] - E0[k])*(mup1-mup0)/2.);
      F_ssc2[n][k]/=(dt*(E1[k] - E0[k])*(mup1-mup0)/2.);
      /* F[n][k]/=(dt*(E1[k] - E0[k])); */
   }
}

if( !(fpout=fopen(outfile, "w")) ) file_error(outfile);
if( !(fpec=fopen(ecfile, "w")) ) file_error(ecfile);
if( !(fpssc=fopen(sscfile, "w")) ) file_error(sscfile);
if( !(fpssc2=fopen(ssc2file, "w")) ) file_error(ssc2file);

if (test == 'y') printf ("\n Energy binning:\n");

fprintf (fpout, "#time(s): %e %e dt(s): %e\n", t0[0], t1[n_t], dt);
fprintf (fpout, "#angle: %f %f\n", mu0, mu1);
fprintf (fpout, "#factor: %i\n", factor);
fprintf (fpout, "#Energy(keV)   Luminosity(erg/s/keV)\n");
fprintf (fpec, "#time(s): %e %e dt(s): %e\n", t0[0], t1[n_t], dt);
fprintf (fpec, "#angle: %f %f\n", mu0, mu1);
fprintf (fpec, "#factor: %i\n", factor);
fprintf (fpec, "#Energy(keV)   Luminosity(erg/s/keV)\n");
fprintf (fpssc, "#time(s): %e %e dt(s): %e\n", t0[0], t1[n_t], dt);
fprintf (fpssc, "#angle: %f %f\n", mu0, mu1);
fprintf (fpssc, "#factor: %i\n", factor);
fprintf (fpssc, "#Energy(keV)   Luminosity(erg/s/keV)\n");
fprintf (fpssc2, "#time(s): %e %e dt(s): %e\n", t0[0], t1[n_t], dt);
fprintf (fpssc2, "#angle: %f %f\n", mu0, mu1);
fprintf (fpssc2, "#factor: %i\n", factor);
fprintf (fpssc2, "#Energy(keV)   Luminosity(erg/s/keV)\n");


for (k=0; k<n_bins; k++) {
   if (test == 'y') {
      printf ("\n %d. E0 = %e keV, E1 = %e keV N(E) = %d", k, E0[k], E1[k], particle_sum[n_sum][k]);
   }

   fprintf (fpout, "%f ", log10(max(1.e-20, sqrt(E0[k]*E1[k]))));
   fprintf (fpec, "%f ", log10(max(1.e-20, sqrt(E0[k]*E1[k]))));
   fprintf (fpssc, "%f ", log10(max(1.e-20, sqrt(E0[k]*E1[k]))));
   fprintf (fpssc2, "%f ", log10(max(1.e-20, sqrt(E0[k]*E1[k]))));
   for (n=0; n<n_t-1; n++) {
      fprintf (fpout, "%f ", log10(max(1.e-20, F[n][k])));
      fprintf (fpec, "%f ", log10(max(1.e-20, F_ec[n][k])));
      fprintf (fpssc, "%f ", log10(max(1.e-20, F_ssc[n][k])));
      fprintf (fpssc2, "%f ", log10(max(1.e-20, F_ssc2[n][k])));
   }

   fprintf (fpout, "%f %i\n", log10(max(1.e-20, F[n][k])),particle_sum[n][k]);
   fprintf (fpec, "%f\n", log10(max(1.e-20, F_ec[n][k])));
   fprintf (fpssc, "%f\n", log10(max(1.e-20, F_ssc[n][k])));
   fprintf (fpssc2, "%f\n", log10(max(1.e-20, F_ssc2[n][k])));
}
fclose(fpout);
fclose(fpec);
fclose(fpssc);
fclose(fpssc2);
printf("\n");

/*printf("Number of particle bundles:  %i", particle_sum); */

}


