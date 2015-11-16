#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <libgen.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>

#define PI 3.14159265358979323846

void cube_value(double *p, double *v)
{
	double abs_x = fabs(p[0]), abs_y = fabs(p[1]), abs_z = fabs(p[2]);
	if (abs_x < 1.0 && abs_y < 1.0 && abs_z < 1.0)
		*v = 65535.0 * (fmax(abs_x, fmax(abs_y, abs_z)));
	else
		*v = 0.0;
}

void cube_gradient(double *p, double *g)
{
	(void)p;

	g[0] = 0.0;
	g[1] = 0.0;
	g[2] = 0.0;
}

void octahedron_value(double *p, double *v)
{
	double abs_x = fabs(p[0]), abs_y = fabs(p[1]), abs_z = fabs(p[2]);
	if ((abs_x + abs_y + abs_z) < 1.0)
		*v = 65535.0 * (1.0 - (abs_x + abs_y + abs_z));
	else
		*v = 0.0;
}

void octahedron_gradient(double *p, double *g)
{
	(void)p;

	g[0] = 0.0;
	g[1] = 0.0;
	g[2] = 0.0;
}

void sphere_value(double *p, double *v)
{
	double len = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

	if (len < 1.0)
		*v = (1.0 - len) * 65535.0;
	else
		*v = 0.0;
}

void sphere_gradient(double *p, double *g)
{
	double len2 = 2.0 * sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

	g[0] = 65535.0 * (p[0] / len2 + 0.5);
	g[1] = 65535.0 * (p[1] / len2 + 0.5);
	g[2] = 65535.0 * (p[2] / len2 + 0.5);
}

#define ML_fm 6
#define ML_alpha 0.25

void marschnerlobb_value(double *p, double *v)
{
	double rhor = cos(2.0 * M_PI * ML_fm * cos(M_PI * sqrt(p[0] * p[0] + p[1] * p[1]) / 2.0));
	*v = 65535.0*(1-sin(M_PI * p[2] / 2.0) + ML_alpha *(1.0 + rhor)) / (2.0 * (1.0 + ML_alpha));
}

void marschnerlobb_gradient(double *p, double *g)
{
    g[0] = ML_alpha*sin(2.0*M_PI*ML_fm*cos(M_PI*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/2.0))*M_PI*M_PI*ML_fm*sin(M_PI*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/2.0)/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])*p[0]/(1.0+ML_alpha)/2.0;
    g[1] = ML_alpha*sin(2.0*M_PI*ML_fm*cos(M_PI*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/2.0))*M_PI*M_PI*ML_fm*sin(M_PI*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/2.0)/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])*p[1]/(1.0+ML_alpha)/2.0;
    g[2] = M_PI*(-cos(M_PI*p[2]/2.0)*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])+2.0*ML_alpha*sin(2.0*M_PI*ML_fm*cos(M_PI*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/2.0))*M_PI*ML_fm*sin(M_PI*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/2.0)*p[2])/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/(1.0+ML_alpha)/4.0;

	double len2 = 2.0 * sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);

    g[0] = 65535.0 * ((-g[0] / len2) + 0.5);
    g[1] = 65535.0 * ((-g[1] / len2) + 0.5);
    g[2] = 65535.0 * ((-g[2] / len2) + 0.5);
}

#define LS_c 2

void linnerstrand_value(double *p, double *v)
{
	*v = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) < 1.0 ? 65535.0 * cos(LS_c * LS_c * (p[0]*p[0] + p[1]*p[1] + p[2]*p[2])) + 1 : 0.0;
}

void linnerstrand_gradient(double *p, double *g)
{
	(void)p;

	g[0] = 0.0;
	g[1] = 0.0;
	g[2] = 0.0;
}

FILE *open_raw(const char *base, int index, const char *ext)
{
	FILE *f;
	char path[strlen(base) + 7];
	sprintf(path, "%s_%d.%3s", base, index, ext);

	printf("opening %s...\n", path);

	if ((f = fopen(path, "w")) == NULL)
	{
		perror("fopen");
		exit(-1);
	}

	return f;
}

void usage()
{
	fputs(
		"usage: mkvol -g GRID -f FUNCTION -v VOXELS -o BASENAME [-i]\n"
		"where:\n"
		"  -g GRID     Grid type. GRID = cc | bcc | fcc\n"
		"  -f FUNCTION Function. FUNCTION = sphere | cube | octahedron | marschnerlobb | linnerstrand\n"
		"  -v VOXELS   Create a volume of this many voxels. VOXELS > 0\n"
		"  -o BASENAME The basename of the output files. BASENAME = valid file path\n"	
		"  -i          Use interleaved format.\n",
		stderr);
}

enum 
{
	GRID_NONE = 0,
	GRID_CC,
	GRID_BCC,
	GRID_FCC
};

int main(int argc, char **argv)
{
	void (*getValue)(double *, double *) = NULL;
	void (*getGradient)(double *, double *) = NULL;

	char *option_basename = NULL;
	long  option_voxels = 0;
	int   option_interleaved = 0;
	int   option_grid = GRID_NONE;
	int   option_gradient = 0;
	int c;

	while ((c = getopt (argc, argv, "g:f:v:o:iG")) != -1) {
		switch (c) {
			case 'g':
				if (strcmp("cc", optarg) == 0)
					option_grid = GRID_CC;
				else if (strcmp("bcc", optarg) == 0)
					option_grid = GRID_BCC;
				else if (strcmp("fcc", optarg) == 0)
					option_grid = GRID_FCC;
				else {
					usage();
					exit(1);
				}
				break;

			case 'f':
				if (strcmp("sphere", optarg) == 0) {
					getValue = &sphere_value;
					getGradient = &sphere_gradient;
				}
				else if (strcmp("cube", optarg) == 0) {
					getValue = &cube_value;
					getGradient = &cube_gradient;
				}
				else if (strcmp("octahedron", optarg) == 0) {
					getValue = &octahedron_value;
					getGradient = &octahedron_gradient;
				}
				else if (strcmp("marschnerlobb", optarg) == 0) {
					getValue = &marschnerlobb_value;
					getGradient = &marschnerlobb_gradient;
				}
				else if (strcmp("linnerstrand", optarg) == 0) {
					getValue = &linnerstrand_value;
					getGradient = &linnerstrand_gradient;
				}
				else {
					usage();
					exit(1);
				}
				break;

			case 'v':
				option_voxels = atol(optarg);
				break;

			case 'i':
				option_interleaved = 1;
				break;

			case 'G':
				option_gradient = 1;
				break;

			case 'o':
				option_basename = optarg;
				break;

			default:
				fprintf(stderr, "error: Unknown option: -%c\n", c);

				usage();
				exit(0);
				break;
		}
	}

	/* Sanity checks. */
	if (option_interleaved && option_gradient) {
		fputs("error: Can't use interleaved format with gradients enabled.\n", stderr);
		usage();
		exit(1);
	}
	if (option_grid == GRID_NONE) {
		fputs("error: No grid type selected.\n", stderr);
		usage();
		exit(1);
	}
	if (option_voxels <= 0) {
		fputs("error: Please specify a positive number of voxels.\n", stderr);
		usage();
		exit(1);
	}

	int sx = 0, sy = 0, sz = 0;
	double x, y, z;
	double step_size = 1.0;

	FILE *f1 = NULL, *f2 = NULL, *f3 = NULL, *f4 = NULL;

	f1 = open_raw(option_basename, 1, "raw");

	if (option_grid == GRID_CC)
	{
		step_size = 2.0 / pow(option_voxels, 1.0 / 3.0);
	}
	else if (option_grid == GRID_BCC)
	{
		if (option_interleaved)
			f2 = f1;
		else
			f2 = open_raw(option_basename, 2, "raw");

		step_size = 2.0 / pow(option_voxels / 2.0, 1.0 / 3.0);
	}
	else if (option_grid == GRID_FCC)
	{
		if (option_interleaved) {
			f2 = f1;
			f3 = f1;
			f4 = f1;
		}
		else {
			f2 = open_raw(option_basename, 2, "raw");
			f3 = open_raw(option_basename, 3, "raw");
			f4 = open_raw(option_basename, 4, "raw");
		}

		step_size = 2.0 / pow(option_voxels / 4.0, 1.0 / 3.0);		
	}


	for (z = -1.0; z < 1.0; z += step_size)
	{
		for (sy = 0, y = -1.0; y < 1.0; y += step_size)
		{
			for (sx = 0, x = -1.0; x < 1.0; x += step_size)
			{
				double v[3];
				double p[3];
				double g[3];
				uint16_t data[16];

				p[0] = x;
				p[1] = y;
				p[2] = z;

				getValue(p, v);
				data[option_gradient ? 3 : 0] = (uint16_t)v[0];

				if (option_gradient) {
					getGradient(p, v);
					data[0] = (uint16_t)v[0];
					data[1] = (uint16_t)v[1];
					data[2] = (uint16_t)v[2];
				}

				fwrite(data, sizeof(uint16_t), option_gradient ? 4 : 1, f1);

				if (option_grid == GRID_BCC) {
					p[0] = x + (step_size / 2.0);
					p[1] = y + (step_size / 2.0);
					p[2] = z + (step_size / 2.0);
					getValue(p, v);
					data[option_gradient ? 3 : 0] = (uint16_t)v[0];

					if (option_gradient) {
						getGradient(p, v);
						data[0] = (uint16_t)v[0];
						data[1] = (uint16_t)v[1];
						data[2] = (uint16_t)v[2];
					}

					fwrite(data, sizeof(uint16_t), option_gradient ? 4 : 1, f2);
				}

				if (option_grid == GRID_FCC) {
					p[0] = x;
					p[1] = y + (step_size / 2.0);
					p[2] = z + (step_size / 2.0);
					getValue(p, v);
					data[option_gradient ? 3 : 0] = (uint16_t)v[0];
					if (option_gradient) {
						getGradient(p, g);
						data[0] = (uint16_t)g[0];
						data[1] = (uint16_t)g[1];
						data[2] = (uint16_t)g[2];
					}					
					fwrite(data, sizeof(uint16_t), option_gradient ? 4 : 1, f2);

					p[0] = x + (step_size / 2.0);
					p[1] = y;
					p[2] = z;
					getValue(p, v);
					data[option_gradient ? 3 : 0] = (uint16_t)v[0];
					if (option_gradient) {
						getGradient(p, g);
						data[0] = (uint16_t)g[0];
						data[1] = (uint16_t)g[1];
						data[2] = (uint16_t)g[2];
					}					
					fwrite(data, sizeof(uint16_t), option_gradient ? 4 : 1, f3);
					
					p[0] = x + (step_size / 2.0);
					p[1] = y + (step_size / 2.0);
					p[2] = z;
					getValue(p, v);
					data[option_gradient ? 3 : 0] = (uint16_t)v[0];
					if (option_gradient) {
						getGradient(p, g);
						data[0] = (uint16_t)g[0];
						data[1] = (uint16_t)g[1];
						data[2] = (uint16_t)g[2];
					}					
					fwrite(data, sizeof(uint16_t), option_gradient ? 4 : 1, f4);
				}

				++sx;
			}
			++sy;
		}
		++sz;
	}

	if (option_grid == GRID_CC)
		printf("%d×%d×%d = %d voxels generated\n", sx, sy, sz, sx*sy*sz);
	if (option_grid == GRID_BCC)
		printf("2×%d×%d×%d = %d voxels generated\n", sx, sy, sz, 2*sx*sy*sz);
	if (option_grid == GRID_FCC)
		printf("4×%d×%d×%d = %d voxels generated\n", sx, sy, sz, 4*sx*sy*sz);

	char *object_model = NULL;

	if (option_interleaved) {
		if (option_grid == GRID_CC)
			object_model = "I";
		else if (option_grid == GRID_BCC)
			object_model = "II";
		else if (option_grid == GRID_FCC)
			object_model = "IIII";
	}
	else if (option_gradient) {
		object_model = "RGBA";
	}
	else {
		object_model = "I";
	}

	const char *objfile = basename(option_basename);

	FILE *datf1 = open_raw(option_basename, 1, "dat");
	fprintf(datf1,
		"ObjectFileName: %s_%d.raw\n"
		"Resolution:     %d %d %d\n"
		"SliceThickness: 1 1 1\n"
		"Format:         USHORT\n"
		"ObjectModel:    %s\n",
		objfile, 1, sx, sy, sz, object_model);
	fclose(datf1);

	if (!option_interleaved && (option_grid == GRID_BCC || option_grid == GRID_FCC)) {
		FILE *datf2 = open_raw(option_basename, 2, "dat");
		fprintf(datf2,
			"ObjectFileName: %s_%d.raw\n"
			"Resolution:     %d %d %d\n"
			"SliceThickness: 1 1 1\n"
			"Format:         USHORT\n"
			"ObjectModel:    %s\n",
			objfile, 2, sx, sy, sz, object_model);
		fclose(datf2);
	}

	if (!option_interleaved && option_grid == GRID_FCC) {
		FILE *datf3 = open_raw(option_basename, 3, "dat");
		fprintf(datf3,
			"ObjectFileName: %s_%d.raw\n"
			"Resolution:     %d %d %d\n"
			"SliceThickness: 1 1 1\n"
			"Format:         USHORT\n"
			"ObjectModel:    %s\n",
			objfile, 3, sx, sy, sz, object_model);
		fclose(datf3);

		FILE *datf4 = open_raw(option_basename, 4, "dat");
		fprintf(datf4,
			"ObjectFileName: %s_%d.raw\n"
			"Resolution:     %d %d %d\n"
			"SliceThickness: 1 1 1\n"
			"Format:         USHORT\n"
			"ObjectModel:    %s\n",
			objfile, 4, sx, sy, sz, object_model);
		fclose(datf4);
	}

	return 0;
}

