/*------------------------------------------------------------------------------
 * antenna.c : ANTEX 1.4 receiver antenna phase centre model
 *----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "nrtk/antenna.h"

static void trim(char *s)
{
    char *p=s+strlen(s);
    while (p>s&&isspace((unsigned char)p[-1])) *--p='\0';
    for (p=s;*p&&isspace((unsigned char)*p);p++) ;
    if (p!=s) memmove(s,p,strlen(p)+1);
}
static int is_satid(const char *s)
{
    return strlen(s)==3&&strchr("GRECJIS",s[0])&&isdigit((unsigned char)s[1])&&
           isdigit((unsigned char)s[2]);
}
static void free_ant(nrtk_ant_t *a)
{
    int i;
    for (i=0;i<a->nf;i++) {free(a->f[i].noazi); free(a->f[i].azi);}
    a->nf=0;
}
static int add_ant(nrtk_atx_t *atx, const nrtk_ant_t *a)
{
    nrtk_ant_t *p;
    if (atx->n>=atx->nmax) {
        atx->nmax=atx->nmax<=0?256:atx->nmax*2;
        if (!(p=(nrtk_ant_t *)realloc(atx->ant,sizeof(nrtk_ant_t)*atx->nmax))) return 0;
        atx->ant=p;
    }
    atx->ant[atx->n++]=*a;
    return 1;
}
static int read_vals(const char *buff, int start, double *v, int n)
{
    int i,len=(int)strlen(buff);
    char s[16];
    for (i=0;i<n;i++) {
        if (start+8*i>=len) return i;
        memcpy(s,buff+start+8*i,8); s[8]='\0';
        v[i]=atof(s)*1E-3;  /* mm -> m */
    }
    return n;
}
extern int nrtk_atx_read(const char *file, nrtk_atx_t *atx)
{
    FILE *fp;
    char buff[1024],serial[32],code[4];
    nrtk_ant_t a;
    nrtk_antfrq_t *f=NULL;
    double zen1=0,zen2=90,dzen=5,dazi=0,neu[3];
    int state=0,skip=0,i,nrow;

    if (!(fp=fopen(file,"r"))) return 0;
    memset(&a,0,sizeof(a));

    while (fgets(buff,sizeof(buff),fp)) {
        const char *label=strlen(buff)>60?buff+60:"";

        if (strstr(label,"START OF ANTENNA")) {
            memset(&a,0,sizeof(a)); state=1; skip=0; f=NULL;
            continue;
        }
        if (strstr(label,"END OF ANTENNA")) {
            if (state&&!skip&&a.nf>0) {if (!add_ant(atx,&a)) free_ant(&a);}
            else free_ant(&a);
            state=0; continue;
        }
        if (!state||skip) continue;

        if (strstr(label,"TYPE / SERIAL NO")) {
            memcpy(a.type,buff,16); a.type[16]='\0'; trim(a.type);
            memcpy(a.radome,buff+16,4); a.radome[4]='\0'; trim(a.radome);
            if (!*a.radome) strcpy(a.radome,"NONE");
            memcpy(serial,buff+20,20); serial[20]='\0'; trim(serial);
            if (is_satid(serial)) skip=1;    /* satellite antenna */
        }
        else if (strstr(label,"ZEN1 / ZEN2 / DZEN")) {
            sscanf(buff,"%lf %lf %lf",&zen1,&zen2,&dzen);
        }
        else if (strstr(label,"DAZI")) {
            sscanf(buff,"%lf",&dazi);
        }
        else if (strstr(label,"START OF FREQUENCY")&&!strstr(label,"RMS")) {
            if (a.nf>=NRTK_ANT_MAXFRQ) {f=NULL; continue;}
            f=&a.f[a.nf++];
            memset(f,0,sizeof(*f));
            memcpy(code,buff+3,3); code[3]='\0';
            if (code[1]==' ') code[1]='0';
            strcpy(f->code,code);
            f->zen1=zen1; f->zen2=zen2; f->dzen=dzen>0.0?dzen:5.0; f->dazi=dazi;
            f->nzen=(int)((zen2-zen1)/f->dzen+0.5)+1;
            f->nazi=dazi>0.0?(int)(360.0/dazi+0.5)+1:0;
            f->noazi=(double *)calloc(f->nzen,sizeof(double));
            f->azi=f->nazi>0?(double *)calloc((size_t)f->nazi*f->nzen,sizeof(double)):NULL;
        }
        else if (strstr(label,"END OF FREQUENCY")&&!strstr(label,"RMS")) {
            f=NULL;
        }
        else if (strstr(label,"START OF FREQ RMS")) {
            f=NULL; state=2;
        }
        else if (strstr(label,"END OF FREQ RMS")) {
            state=1;
        }
        else if (state==1&&f&&strstr(label,"NORTH / EAST / UP")) {
            if (sscanf(buff,"%lf %lf %lf",neu,neu+1,neu+2)==3) {
                for (i=0;i<3;i++) f->neu[i]=neu[i]*1E-3;
            }
        }
        else if (state==1&&f&&!strncmp(buff+3,"NOAZI",5)) {
            read_vals(buff,8,f->noazi,f->nzen);
        }
        else if (state==1&&f&&f->nazi>0&&strlen(buff)>8&&
                 (isdigit((unsigned char)buff[5])||buff[5]=='.')) {
            double az=atof(buff);
            nrow=(int)(az/f->dazi+0.5);
            if (nrow>=0&&nrow<f->nazi) read_vals(buff,8,f->azi+(size_t)nrow*f->nzen,f->nzen);
        }
    }
    fclose(fp);
    return 1;
}
extern void nrtk_atx_free(nrtk_atx_t *atx)
{
    int i;
    for (i=0;i<atx->n;i++) free_ant(atx->ant+i);
    free(atx->ant);
    atx->ant=NULL; atx->n=atx->nmax=0;
}
extern const nrtk_ant_t *nrtk_atx_find(const nrtk_atx_t *atx, const char *name)
{
    char type[64]="",radome[16]="NONE";
    int i,pass;

    if (!atx||!name) return NULL;
    if (sscanf(name,"%63s %15s",type,radome)<1) return NULL;
    for (pass=0;pass<2;pass++) {
        for (i=0;i<atx->n;i++) {
            if (strcmp(atx->ant[i].type,type)) continue;
            if (!strcmp(atx->ant[i].radome,pass?"NONE":radome)) return atx->ant+i;
        }
    }
    return NULL;
}
/* ANTEX frequency code for system/frequency ---------------------------------*/
static void freq_code(int sys, double freq, char *code)
{
    const double MHz=freq/1E6;
    char c='G';
    int n=1;

    switch (sys) {
        case SYS_GAL: c='E'; break;
        case SYS_CMP: c='C'; break;
        case SYS_QZS: c='J'; break;
        case SYS_GLO: c='R'; break;
        default: c='G';
    }
    if      (fabs(MHz-1575.42 )<1.0) n=1;
    else if (fabs(MHz-1227.60 )<1.0) n=2;
    else if (fabs(MHz-1176.45 )<1.0) n=5;
    else if (fabs(MHz-1278.75 )<1.0) n=6;
    else if (fabs(MHz-1207.14 )<1.0) n=7;
    else if (fabs(MHz-1191.795)<1.0) n=8;
    else if (fabs(MHz-1561.098)<1.0) n=2;  /* BDS B1I */
    else if (fabs(MHz-1268.52 )<1.0) n=6;  /* BDS B3I */
    sprintf(code,"%c%02d",c,n);
}
static const nrtk_antfrq_t *find_frq(const nrtk_ant_t *ant, int sys, double freq)
{
    char code[8];
    int i;
    const double MHz=freq/1E6;

    freq_code(sys,freq,code);
    for (i=0;i<ant->nf;i++) if (!strcmp(ant->f[i].code,code)) return ant->f+i;

    /* fall back to the GPS frequency closest in band */
    strcpy(code,MHz>1500.0?"G01":(MHz>1200.0?"G02":"G05"));
    for (i=0;i<ant->nf;i++) if (!strcmp(ant->f[i].code,code)) return ant->f+i;
    for (i=0;i<ant->nf;i++) if (!strcmp(ant->f[i].code,"G01")) return ant->f+i;
    return ant->nf>0?ant->f:NULL;
}
static double interp_pcv(const nrtk_antfrq_t *f, double az, double el)
{
    double zen=90.0-el*R2D,x,y,a,v0,v1;
    int iz,ia;

    if (f->nzen<=0) return 0.0;
    x=(zen-f->zen1)/f->dzen;
    if (x<0.0) x=0.0;
    if (x>f->nzen-1) x=f->nzen-1;
    iz=(int)x; if (iz>=f->nzen-1) iz=f->nzen-2;
    if (iz<0) return f->noazi[0];
    a=x-iz;

    if (f->nazi<=1||!f->azi) {
        return f->noazi[iz]*(1.0-a)+f->noazi[iz+1]*a;
    }
    az=fmod(az*R2D,360.0); if (az<0.0) az+=360.0;
    y=az/f->dazi;
    ia=(int)y; if (ia>=f->nazi-1) ia=f->nazi-2;
    y-=ia;
    v0=f->azi[(size_t)ia*f->nzen+iz]*(1.0-a)+f->azi[(size_t)ia*f->nzen+iz+1]*a;
    v1=f->azi[(size_t)(ia+1)*f->nzen+iz]*(1.0-a)+f->azi[(size_t)(ia+1)*f->nzen+iz+1]*a;
    return v0*(1.0-y)+v1*y;
}
extern double nrtk_ant_dr(const nrtk_ant_t *ant, int sys, double freq,
                          const double *azel, const double *del, int pcv)
{
    const nrtk_antfrq_t *f=NULL;
    double e[3],off[3]={0},cosel=cos(azel[1]);
    int i;

    e[0]=sin(azel[0])*cosel;
    e[1]=cos(azel[0])*cosel;
    e[2]=sin(azel[1]);

    if (ant&&freq>0.0) f=find_frq(ant,sys,freq);
    if (f) {
        off[0]=f->neu[1]; off[1]=f->neu[0]; off[2]=f->neu[2]; /* NEU -> ENU */
    }
    if (del) for (i=0;i<3;i++) off[i]+=del[i];

    return -dot3(off,e)+(f&&pcv?interp_pcv(f,azel[0],azel[1]):0.0);
}
