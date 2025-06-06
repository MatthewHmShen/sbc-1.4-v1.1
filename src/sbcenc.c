/*
 *
 *  Bluetooth low-complexity, subband codec (SBC) encoder
 *
 *  Copyright (C) 2008-2010  Nokia Corporation
 *  Copyright (C) 2004-2010  Marcel Holtmann <marcel@holtmann.org>
 *  Copyright (C) 2012-2013  Intel Corporation
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>

#include "../sbc/sbc.h"

static int verbose = 0;

#define BUF_SIZE 32768
static unsigned char input[BUF_SIZE], output[BUF_SIZE + BUF_SIZE / 4];

static void encode(char *filename, int subbands, int bitpool, int joint,
                   int dualchannel, int snr, int blocks, bool msbc)
{
    sbc_t sbc;
    int fd, size, srate, codesize, nframes;
    ssize_t encoded;
    ssize_t len;
    unsigned char input[1024 * 1024]; // 假设输入缓冲区大小为1MB
    unsigned char output[1024 * 1024]; // 假设输出缓冲区大小为1MB

    if (strcmp(filename, "-")) {
        fd = open(filename, O_RDONLY);
        if (fd < 0) {
            fprintf(stderr, "Can't open file %s: %s\n", filename, strerror(errno));
            return;
        }
    } else {
        fd = fileno(stdin);
    }

    if (!msbc) {
        sbc_init(&sbc, 0L);

        // 假设采样率为44100 Hz
        srate = 44100;
        sbc.frequency = SBC_FREQ_44100;

        sbc.subbands = subbands == 4 ? SBC_SB_4 : SBC_SB_8;

        if (joint && !dualchannel)
            sbc.mode = SBC_MODE_JOINT_STEREO;
        else if (!joint && dualchannel)
            sbc.mode = SBC_MODE_DUAL_CHANNEL;
        else if (!joint && !dualchannel)
            sbc.mode = SBC_MODE_STEREO;
        else {
            fprintf(stderr, "Both joint and dualchannel mode have been specified\n");
            goto done;
        }

        sbc.endian = SBC_BE;
        sbc.bitpool = bitpool;
        sbc.allocation = snr ? SBC_AM_SNR : SBC_AM_LOUDNESS;
        sbc.blocks = blocks == 4 ? SBC_BLK_4 :
                   blocks == 8 ? SBC_BLK_8 :
                   blocks == 12 ? SBC_BLK_12 : SBC_BLK_16;
    } else {
        // 假设mSBC要求16位，16kHz，单声道输入
        srate = 16000;
        sbc_init_msbc(&sbc, 0);
        sbc.endian = SBC_BE;
    }

    codesize = sbc_get_codesize(&sbc);
    nframes = sizeof(input) / codesize;

    while (1) {
        unsigned char *inp, *outp;
        size = read(fd, input, codesize * nframes);
        if (size < 0) {
            perror("Can't read audio data");
            break;
        }
        if (size < codesize) {
            break;
        }

        inp = input;
        outp = output;
        while (size >= codesize) {
            len = sbc_encode(&sbc, inp, codesize, outp, sizeof(output) - (outp - output), &encoded);
            if (len != codesize || encoded <= 0) {
                fprintf(stderr, "sbc_encode fail, len=%zd, encoded=%lu\n", len, (unsigned long) encoded);
                break;
            }
            size -= len;
            inp += len;
            outp += encoded;
        }

        len = write(fileno(stdout), output, outp - output);
        if (len != outp - output) {
            perror("Can't write SBC output");
            break;
        }
        if (size != 0) {
            break;
        }
    }

    sbc_finish(&sbc);

done:
    if (fd > fileno(stderr)) {
        close(fd);
    }
}

static void usage(void)
{
	// printf("SBC encoder utility ver %s\n", VERSION);
	printf("Copyright (c) 2004-2010  Marcel Holtmann\n\n");

	printf("Usage:\n"
		"\tsbcenc [options] file(s)\n"
		"\n");

	printf("Options:\n"
		"\t-h, --help           Display help\n"
		"\t-v, --verbose        Verbose mode\n"
		"\t-m, --msbc           mSBC codec\n"
		"\t-s, --subbands       Number of subbands to use (4 or 8)\n"
		"\t-b, --bitpool        Bitpool value (default is 32)\n"
		"\t-j, --joint          Joint stereo\n"
		"\t-d, --dualchannel    Dual channel\n"
		"\t-S, --snr            Use SNR mode (default is loudness)\n"
		"\t-B, --blocks         Number of blocks (4, 8, 12 or 16)\n"
		"\n");
}

static struct option main_options[] = {
	{ "help",	0, 0, 'h' },
	{ "verbose",	0, 0, 'v' },
	{ "msbc",	0, 0, 'm' },
	{ "subbands",	1, 0, 's' },
	{ "bitpool",	1, 0, 'b' },
	{ "joint",	0, 0, 'j' },
	{ "dualchannel",0, 0, 'd' },
	{ "snr",	0, 0, 'S' },
	{ "blocks",	1, 0, 'B' },
	{ 0, 0, 0, 0 }
};

int main(int argc, char *argv[])
{
	int i, opt, subbands = 8, bitpool = 32, joint = 0, dualchannel = 0;
	int snr = 0, blocks = 16;
	bool msbc = false;

	while ((opt = getopt_long(argc, argv, "+hmvs:b:jdSB:",
						main_options, NULL)) != -1) {
		switch(opt) {
		case 'h':
			usage();
			exit(0);

		case 'v':
			verbose = 1;
			break;

		case 'm':
			msbc = true;
			break;

		case 's':
			subbands = atoi(optarg);
			if (subbands != 8 && subbands != 4) {
				fprintf(stderr, "Invalid subbands\n");
				exit(1);
			}
			break;

		case 'b':
			bitpool = atoi(optarg);
			break;

		case 'j':
			joint = 1;
			break;

		case 'd':
			dualchannel = 1;
			break;

		case 'S':
			snr = 1;
			break;

		case 'B':
			blocks = atoi(optarg);
			if (blocks != 16 && blocks != 12 &&
						blocks != 8 && blocks != 4) {
				fprintf(stderr, "Invalid blocks\n");
				exit(1);
			}
			break;

		default:
			usage();
			exit(1);
		}
	}

	argc -= optind;
	argv += optind;
	optind = 0;

	if (argc < 1) {
		usage();
		exit(1);
	}

	for (i = 0; i < argc; i++)
		encode(argv[i], subbands, bitpool, joint, dualchannel,
							snr, blocks, msbc);

	return 0;
}
