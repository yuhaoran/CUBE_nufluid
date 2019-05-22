#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef	struct 
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  double Hz;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8 -8];  /* fills to 256 Bytes */
} IO_HEADER;

#define HUBBLE0 0.1

int read_gadget_header(FILE *fp,IO_HEADER *h);

int main(int argc, char **argv)
{	
	float Hubble=0.1;
	int i,dummy,dummy2,n,ns;
	long int Ntotal=0;

	FILE *fp;
	IO_HEADER h;
	
	char snapfile[1024];

	if(argc!=2)
	{
	printf("usage: %s [filename]\n ",argv[0]);
	exit(1);
	}
	strcpy(snapfile,argv[1]);

	/*=====reading snapshot========*/
	if((fp=fopen(snapfile,"r"))==NULL)
	{
		fprintf(stderr,"error: file open failed for %s!\n",snapfile);
		fflush(stderr);
		exit(1);
	}
	
	read_gadget_header(fp,&h);	
		
	printf("partnum\t\tpartnum_total\t\tpartmass\n");
	for(i=0;i<6;i++)
	{
	printf("%d\t\t%u\t\t%g\n",h.npart[i],h.npartTotal[i],h.mass[i]);
	Ntotal+=h.npartTotal[i];
	}
	printf("---------------------------------------\n");
	printf("Ntotal:\t\t%ld\n",Ntotal);
	printf("NFiles:\t\t%d\n",h.num_files);
	printf("SFR: %d\tFeedback: %d\tCooling: %d\n",h.flag_sfr,h.flag_feedback,h.flag_cooling);
	printf("Boxsize:%g\nOmega0:%f\t OmegaL0:%f\t h: %f\t Hz: %f\n",h.BoxSize,h.Omega0,h.OmegaLambda,h.HubbleParam,h.Hz);
	printf("Redshift: %3.2f\tscaleFac: %g\n",h.redshift,h.time);
	return 0;
}

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | (((x) << 8) & 0x00ff0000) | (((x) >> 8) & 0x0000ff00) | ((unsigned) (x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_LONG(x) (*(unsigned *)&(x) = SWAP_4(*(unsigned *)&(x)))

void swap_Nbyte(void *data2swap,size_t nel,size_t mbyte)
/*This function is used to switch endian, for data2swap[nel] with element size mbyte*/
{
  size_t i,j;
  char *data, *old_data;//by definition, sizeof(char)=1, one byte
  
  data=(char *)data2swap;
  
  switch(mbyte)
  {
	  case 1 :break;
	  case 2 :
	  		for(j=0;j<nel;j++)
			FIX_SHORT(data[j*2]);
			break;
	  case 4 :
	  		for(j=0;j<nel;j++)
			FIX_LONG(data[j*4]);
			break;
	  default :
			old_data=malloc(mbyte);
			for(j=0;j<nel;j++)
			{
			  memcpy(&old_data[0],&data[j*mbyte],mbyte);
			  for(i=0;i<mbyte;i++)
				{
				  data[j*mbyte+i]=old_data[mbyte-i-1];
				}
			}
			free(old_data);
	}
}
size_t fread_swap(void *buf,size_t Nsize,size_t Nbuf,FILE *fp, int FlagByteSwap)
{
	size_t Nread;
	Nread=fread(buf,Nsize,Nbuf,fp);
	if(FlagByteSwap)
	swap_Nbyte(buf,Nbuf,Nsize);
	return Nread;
}

int read_gadget_header(FILE *fp,IO_HEADER * h)
{//read the header part, assign header extensions, and do several consistency check
//return ByteOrder
	int dummy,dummy2,n,ns,ByteOrder;
	
    #define myfread(a,b,c,d) fread_swap(a,b,c,d,ByteOrder)

	n=sizeof(IO_HEADER);
	ns=n;
	swap_Nbyte((char *)&ns,1,sizeof(ns));
			
	fread(&dummy,sizeof(dummy),1,fp);
	
	/*=====determine byteorder======*/
	if(dummy==n)
	 ByteOrder=0;
	else if(dummy==ns)
	{
		ByteOrder=1;
		printf(" --------------------------------------------------------\n");
		printf("|WARNING: Different Endianness detected! Doing ByteSwap..|\n");
		printf(" --------------------------------------------------------\n");
	}
	else
	{
		printf("header size not expected:%d;%d,%d\n",dummy,n,ns);
		exit(1);
	}
	
	dummy=n;

	myfread(h->npart,sizeof(int),6,fp);
	myfread(h->mass,sizeof(double),6,fp);
	myfread(&h->time,sizeof(double),1,fp);
	myfread(&h->redshift,sizeof(double),1,fp);
	myfread(&h->flag_sfr,sizeof(int),1,fp);
	myfread(&h->flag_feedback,sizeof(int),1,fp);
	myfread(h->npartTotal,sizeof(int),6,fp);
	myfread(&h->flag_cooling,sizeof(int),1,fp);
	myfread(&h->num_files,sizeof(int),1,fp);
	myfread(&h->BoxSize,sizeof(double),1,fp);
	myfread(&h->Omega0,sizeof(double),1,fp);
	myfread(&h->OmegaLambda,sizeof(double),1,fp);
	myfread(&h->HubbleParam,sizeof(double),1,fp);
	fseek(fp,n+sizeof(int),SEEK_SET);
	myfread(&dummy2,sizeof(dummy2),1,fp);
    if(dummy!=dummy2)
	{
		fprintf(stderr,"error!record brackets not match for header!\t%d,%d\n",dummy,dummy2);
		exit(1);
	} 
	
  
	h->Hz=HUBBLE0 * sqrt(h->Omega0 / (h->time * h->time * h->time) 
			+ (1 - h->Omega0 - h->OmegaLambda) / (h->time * h->time)
			+ h->OmegaLambda);//Hubble param for the current catalogue;
		
	return ByteOrder;
	
}
