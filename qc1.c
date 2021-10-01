#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <string.h>
#define MIN(a, b) (a < b ? a : b)


KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	FILE *resFileWriter;
	gzFile fp;
	kseq_t *seq;
	int l;
	
	int i,j;
	
	int readLength = 0, Q, Qi;
	double thisReadQuality=0;
	double sumQualityReads=0;
	int thisReadQualitySum=0;
	int maxreadLen=0;
	unsigned long long QD[4096][4];
	unsigned long long BD[4096][5];
	unsigned long long numberOfReads=0;
	unsigned long long rQD[4];

	/*unsigned long long Qlt10;
	unsigned long long Q1020;
	unsigned long long Q2030;
	unsigned long long Qgt30;*/

	unsigned long long totalGC=0, totalBases=0;
	unsigned long long GC[10];
	int thisGC=0, thisGCpercentIndex=0;
	double thisGCpercent=0;
	
	for (i=0;i<4;i++)
	{
		rQD[i]=0;
	}

	for (i=0;i<10;i++)
	{
		GC[i]=0;
	}

	for (i=0;i<4096;i++)
	{
		for (j=0;j<4;j++)
		{
			QD[i][j]=0;
		}
	}
	
	//printf("%d\n",argc);
	if (argc != 6)
	{
		fprintf(stderr, "Usage: %s <in.seq> <qualityDistribution.out> <baseDistribution.out> <gcDistribution.out> <stats.out>\n", argv[0]);
		return 1;
	}
	
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	
	while ((l = kseq_read(seq)) >= 0)
	{
		/*printf("name: %s\n", seq->name.s);
		if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		printf("seq: %s\n", seq->seq.s);
		if (seq->qual.l) printf("qual: %s\n", seq->qual.s);*/
		
		//printf("%s\n",seq->qual.s);
		
		readLength=(int)strlen(seq->qual.s);
		totalBases+=readLength;
		//printf("%d\n",readLength);
	    if(readLength>maxreadLen){
			maxreadLen = readLength;
		}	
		numberOfReads++;
		thisReadQualitySum=0;
		thisGC=0;
		for (i=0;i<readLength;i++)
		{
			if (seq->seq.s[i]=='A')
			{
				BD[i][0]++;
			}


			if (seq->seq.s[i]=='C')
			{
				BD[i][1]++;
				totalGC++;
				thisGC++;
			}

			if (seq->seq.s[i]=='G')
			{
				BD[i][2]++;
				totalGC++;
				thisGC++;
			}

			if (seq->seq.s[i]=='T')
			{
				BD[i][3]++;
			}

			if (seq->seq.s[i]=='N')
			{
				BD[i][4]++;
			}

			Q=(int)seq->qual.s[i]-33;

			if (Q<10)          {rQD[0]++;}
			if (Q>=10 && Q<20) {rQD[1]++;}
			if (Q>=20 && Q<30) {rQD[2]++;}
			if(Q>=30)          {rQD[3]++;}

			thisReadQualitySum+=Q;
			Qi=MIN (3, Q/10);
			QD[i][Qi]++;
		}
		thisReadQuality=(double)thisReadQualitySum/(double)readLength;
		sumQualityReads+=thisReadQuality;
		
		/*rQD[MIN(3,(int)thisReadQuality/10)]++;*/
		
		thisGCpercent = ((double)thisGC/(double)readLength)*100;
		thisGCpercentIndex=MIN(9,(int)(thisGCpercent/10));
		GC[thisGCpercentIndex]++;
		//return 0;
	}
	
	
	resFileWriter = fopen(argv[2], "w");
	fprintf(resFileWriter,"Position\tQ<10\t10<Q<20\t20<Q<30\tQ>30\n");
//	printf("%d\n",maxreadLen);
	for (i=0;i<maxreadLen;i++)
	{
		fprintf(resFileWriter,"%d\t",i+1);
		for (j=0;j<4;j++)
		{
			fprintf(resFileWriter,"%f\t",(double)QD[i][j]/numberOfReads*100);
		}fprintf(resFileWriter,"\n");
	}
	fclose(resFileWriter);

	
	resFileWriter = fopen(argv[3], "w");
	fprintf(resFileWriter,"Position\tA\tC\tG\tT\tN\n");
	for (i=0;i<maxreadLen;i++)
	{
		fprintf(resFileWriter,"%d\t",i+1);
		for (j=0;j<5;j++)
		{
			fprintf(resFileWriter,"%f\t",(double)BD[i][j]/numberOfReads*100);
		}fprintf(resFileWriter,"\n");
	}
	fclose(resFileWriter);
	
	resFileWriter = fopen(argv[4], "w");
	for (i=0;i<10;i++)
	{
		fprintf(resFileWriter,"\"%d-%d\", %f\n",i*10,(i+1)*10,(double)GC[i]/(double)numberOfReads*100);
	}
	fclose(resFileWriter);

	resFileWriter = fopen(argv[5], "w");
        fprintf(resFileWriter,"Sample Name                 : %s\n", argv[1]);
	fprintf(resFileWriter,"Number of reads             : %Lu\n",numberOfReads);
        fprintf(resFileWriter,"Number of reads in Million  : %f\n", ((double)numberOfReads/1000000));
	fprintf(resFileWriter,"Number of bases             : %Lu\n",totalBases);
        fprintf(resFileWriter,"Total Data in MB            : %.2f\n",((double)totalBases/1000000.0));
	fprintf(resFileWriter,"Total Data in GB            : %f\n",((double)totalBases/1000000000.0));
        fprintf(resFileWriter,"Percent GC                  : %.2f\n",((double)totalGC/(double)totalBases)*100);
	fprintf(resFileWriter,"Average quality             : %.2f\n",sumQualityReads/numberOfReads);
	fprintf(resFileWriter,"%% bases qual < 10           : %.2f\n",((double)rQD[0]/(double)totalBases)*100);
	fprintf(resFileWriter,"%% bases qual >= 10 and < 20 : %.2f\n",((double)rQD[1]/(double)totalBases)*100);
	fprintf(resFileWriter,"%% bases qual >= 20 and < 30 : %.2f\n",((double)rQD[2]/(double)totalBases)*100);
	fprintf(resFileWriter,"%% bases qual >= 30          : %.2f\n",((double)rQD[3]/(double)totalBases)*100);
	fprintf(resFileWriter,"Average Read length          : %.2f\n",((double)totalBases/(double)numberOfReads));
	fclose(resFileWriter);

	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
