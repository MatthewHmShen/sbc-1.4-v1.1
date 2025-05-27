#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../sbc/sbc.h" // 假设sbc.h是SBC编解码库的头文件

#define MAXCHANNELS 2
#define DEFACCURACY 7
#define BUFFER_SIZE 1024

static double sampletobits(short sample16, int verbose)
{
    double bits = 0;
    unsigned short bit;
    int i;

    if (verbose)
        printf("===> sampletobits(%hd, %04hX)\n", sample16, sample16);

    /* Bit 0 is MSB */
    if (sample16 < 0)
        bits = -1;

    if (verbose)
        printf("%d", (sample16 < 0) ? 1 : 0);

    /* Bit 15 is LSB */
    for (i = 1; i < 16; i++) {
        bit = (unsigned short) sample16;
        bit >>= 15 - i;
        bit %= 2;

        if (verbose)
            printf("%d", bit);

        if (bit)
            bits += (1.0 / pow(2.0, i));
    }

    if (verbose)
        printf("\n");

    return bits;
}

static int calculate_rms_level(FILE *ref, FILE *tst, int channels, int frames, int accuracy, char *csvname)
{
    short refsample[MAXCHANNELS], tstsample[MAXCHANNELS];
    double refbits, tstbits;
    double rms_accu[MAXCHANNELS];
    double rms_level[MAXCHANNELS];
    double rms_limit = 1.0 / (pow(2.0, accuracy - 1) * pow(12.0, 0.5));
    FILE *csv = NULL;
    int i, j, verdict;

    if (csvname)
        csv = fopen(csvname, "wt");

    if (csv) {
        fprintf(csv, "num;");
        for (j = 0; j < channels; j++)
            fprintf(csv, "ref channel %d;tst channel %d;", j, j);
        fprintf(csv, "\r\n");
    }

    memset(rms_accu, 0, sizeof(rms_accu));
    memset(rms_level, 0, sizeof(rms_level));

    for (i = 0; i < frames; i++) {
        if (csv)
            fprintf(csv, "%d;", i);

        for (j = 0; j < channels; j++) {
            if (fread(&refsample[j], sizeof(short), 1, ref) != 1) {
                printf("Failed to read reference data\n");
                if (csv)
                    fclose(csv);
                return -1;
            }

            if (fread(&tstsample[j], sizeof(short), 1, tst) != 1) {
                printf("Failed to read test data\n");
                if (csv)
                    fclose(csv);
                return -1;
            }

            if (csv)
                fprintf(csv, "%d;%d;", refsample[j], tstsample[j]);

            refbits = sampletobits(refsample[j], 0);
            tstbits = sampletobits(tstsample[j], 0);

            rms_accu[j] += pow(tstbits - refbits, 2.0);
        }

        if (csv)
            fprintf(csv, "\r\n");
    }

    printf("Limit: %f\n", rms_limit);

    for (j = 0; j < channels; j++) {
        printf("Channel %d\n", j);
        printf("Accumulated %f\n", rms_accu[j]);
        rms_accu[j] /= (double) frames;
        printf("Accumulated / %f = %f\n", (double) frames, rms_accu[j]);
        rms_level[j] = sqrt(rms_accu[j]);
        printf("Level = %f (%f x %f = %f)\n", rms_level[j], rms_level[j], rms_level[j], rms_level[j] * rms_level[j]);
    }

    verdict = 1;

    for (j = 0; j < channels; j++) {
        printf("Channel %d: %f\n", j, rms_level[j]);

        if (rms_level[j] > rms_limit)
            verdict = 0;
    }

    printf("%s return %d\n", __FUNCTION__, verdict);

    return verdict;
}

static int check_absolute_diff(FILE *ref, FILE *tst, int channels, int frames, int accuracy)
{
    short refsample[MAXCHANNELS], tstsample[MAXCHANNELS];
    double refbits, tstbits;
    double rms_absolute = 1.0 / (pow(2, accuracy - 2));
    int i, j, verdict;

    verdict = 1;

    printf("Absolute max: %f\n", rms_absolute);
    for (i = 0; i < frames; i++) {
        for (j = 0; j < channels; j++) {
            if (fread(&refsample[j], sizeof(short), 1, ref) != 1) {
                printf("Failed to read reference data\n");
                return -1;
            }

            if (fread(&tstsample[j], sizeof(short), 1, tst) != 1) {
                printf("Failed to read test data\n");
                return -1;
            }

            refbits = sampletobits(refsample[j], 0);
            tstbits = sampletobits(tstsample[j], 0);

            double cur_diff = fabs(tstbits - refbits);

            if (cur_diff > rms_absolute) {
                printf("Channel %d exceeded : fabs(%f - %f) = %f > %f\n", j, tstbits, refbits, cur_diff, rms_absolute);
                verdict = 0;
            }
        }
    }

    printf("%s return %d\n", __FUNCTION__, verdict);

    return verdict;
}

static void usage(void)
{
    // printf("SBC conformance test ver %s\n", VERSION);
    printf("Copyright (c) 2007-2010  Marcel Holtmann\n");
    printf("Copyright (c) 2007-2008  Frederic Dalleau\n\n");

    printf("Usage:\n"
           "\tsbctester reference.pcm checkfile.pcm\n"
           "\tsbctester integer\n"
           "\n");

    printf("To test the encoder:\n");
    printf("\tUse a reference codec to encode original.pcm to reference.sbc\n");
    printf("\tUse sbcenc to encode original.pcm to checkfile.sbc\n");
    printf("\tDecode both file using the reference decoder\n");
    printf("\tRun sbctester with these two pcm files to get the result\n\n");

    printf("\tA file called out.csv is generated to use the data in a\n");
    printf("\tspreadsheet application or database.\n\n");
}

int main(int argc, char *argv[])
{
    FILE *ref = NULL;
    FILE *tst = NULL;
    char *refname, *tstname;
    int pass_rms, pass_absolute, pass, accuracy;
    int frames, channels;

    if (argc == 2) {
        double db;

        printf("Test sampletobits\n");
        db = sampletobits((short) atoi(argv[1]), 1);
        printf("db = %f\n", db);
        exit(0);
    }

    if (argc < 3) {
        usage();
        exit(1);
    }

    refname = argv[1];
    tstname = argv[2];

    printf("opening reference %s\n", refname);

    ref = fopen(refname, "rb");
    if (!ref) {
        printf("Failed to open reference file\n");
        exit(1);
    }

    printf("opening testfile %s\n", tstname);
    tst = fopen(tstname, "rb");
    if (!tst) {
        printf("Failed to open test file\n");
        fclose(ref);
        exit(1);
    }

    // 假设输入文件是16位线性PCM格式，单声道
    channels = 1;
    frames = 44100; // 假设采样率为44100 Hz，1秒的数据

    accuracy = DEFACCURACY;
    printf("Accuracy: %d\n", accuracy);

    pass_rms = calculate_rms_level(ref, tst, channels, frames, accuracy, "out.csv");
    if (pass_rms < 0)
        goto error;

    rewind(ref);
    rewind(tst);

    pass_absolute = check_absolute_diff(ref, tst, channels, frames, accuracy);
    if (pass_absolute < 0)
        goto error;

    pass = pass_rms && pass_absolute;
    printf("Verdict: %s\n", pass ? "pass" : "fail");

    fclose(ref);
    fclose(tst);
    return 0;

error:
    fclose(ref);
    fclose(tst);
    exit(1);
}

