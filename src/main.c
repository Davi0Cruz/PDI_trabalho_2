#include <stdio.h>
#include <stdlib.h>
#include "tiffio.h"
#include <string.h>
#include <math.h>

int discovery_dimensions(const char* filename, uint32_t* w, uint32_t* h) {
	char* filenametiff = (char*) malloc(strlen(filename) + 5);
	strcpy(filenametiff, filename);
	strcat(filenametiff, ".tif");
	TIFF* tiff = TIFFOpen(filenametiff, "r");
	TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, w);
	TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, h);
	TIFFClose(tiff);
	free(filenametiff);
	return 1;
}

int readTiff(const char* filename, uint8_t** pixels) {
	char* filenametiff = (char*) malloc(strlen(filename) + 5);
	strcpy(filenametiff, filename);
	uint32_t w, h;
	TIFF* tiff = TIFFOpen(strcat(filenametiff, ".tif"), "r");
	TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
	TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
	uint32_t npixels = w * h;
	*pixels = (uint8_t*) malloc(npixels * sizeof(uint8_t));
	uint32_t* raster = (uint32_t*) _TIFFmalloc(npixels * sizeof(uint32_t));
	TIFFReadRGBAImage(tiff, w, h, raster, 0);
	for (uint32_t row = 0; row < h; row++) {
		for (uint32_t col = 0; col < w; col++) {
			uint32_t pixel = raster[row * w + col];
			uint8_t r = TIFFGetR(pixel);
			uint8_t g = TIFFGetG(pixel);
			uint8_t b = TIFFGetB(pixel);
			(*pixels)[row * w + col] = (uint8_t)(0.299 * r + 0.587 * g + 0.114 * b);
		}
	}
	_TIFFfree(raster);
	TIFFClose(tiff);
	free(filenametiff);
	return 1;
}

int writeTiff(char* filename, uint32_t w, uint32_t h, uint8_t* pixels) {
	char* filenametiff = (char*) malloc(strlen(filename) + 5);
	strcpy(filenametiff, filename);

	TIFF* tiff = TIFFOpen(strcat(filenametiff, ".tif"), "w");
	TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
	TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
	TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(tiff, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
	TIFFSetField(tiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
	for (uint32_t row = 0; row < h; row++) {
		TIFFWriteScanline(tiff, &pixels[row * w], row, 0);
	}
	TIFFClose(tiff);
	free(filenametiff);
	return 1;
}

void histogram(const char* filename, uint8_t* pixels, uint32_t w, uint32_t h) {
	uint32_t hist[256];
	for (uint32_t i = 0; i < 256; i++) hist[i] = 0;
	for (uint32_t row = 0; row < h; row++) {
		for (uint32_t col = 0; col < w; col++) {
			hist[pixels[row * w + col]]++;
		}
	}
	char* filenamedat = (char*) malloc(strlen(filename) + 5);
	strcpy(filenamedat, filename);
	FILE *f = fopen(strcat(filenamedat, ".dat"), "w");
	for (uint32_t i = 0; i < 256; i++) {
		if (hist[i] > 0)
		fprintf(f, "%d %d\n", i, hist[i]);
	}
	fclose(f);	
	free(filenamedat);

	char *filenamemean = (char*) malloc(strlen(filename) + 9);
	strcpy(filenamemean, filename);
	FILE *fmean = fopen(strcat(filenamemean, "_mean.dat"), "w");
	uint32_t sum = 0;
	uint32_t count = 0;
	for (uint32_t i = 0; i < 256; i++) {
		sum += i * hist[i];
		count += hist[i];
	}
	fprintf(fmean, "%f\n", (float)sum / count);
	fclose(fmean);
	free(filenamemean);

	char *filenamecount = (char*) malloc(strlen(filename) + 11);
	strcpy(filenamecount, filename);
	FILE *fcount = fopen(strcat(filenamecount, "_count.dat"), "w");
	fprintf(fcount, "%d\n", count);
	fclose(fcount);
	free(filenamecount);
}

void start_histogram(char* filename) {
	uint32_t w, h;
	uint8_t* pixels;
	discovery_dimensions(filename, &w, &h);
	pixels = (uint8_t*) malloc(w * h * sizeof(uint8_t));
	readTiff(filename, &pixels);
	histogram(filename, pixels, w, h);
	free(pixels);
}

void laplacian(uint8_t * pixels, float k, uint32_t w, uint32_t h) {
	int mask[3][3] = {
		{0, 1, 0},
		{1, -4, 1},
		{0, 1, 0}
	};
	uint8_t* pixels2 = (uint8_t*) malloc(w * h * sizeof(uint8_t));
	for (uint32_t row = 1; row < h - 1; row++) {
		for (uint32_t col = 1; col < w - 1; col++) {
			int sum = 0;
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					sum += pixels[(row + i) * w + (col + j)] * mask[i + 1][j + 1];
				}
			}
			sum = pixels[row * w + col]+k*sum;
			if (sum < 0) sum = 0;
			if (sum > 255) sum = 255;
			pixels2[row * w + col] = (uint8_t) sum;
		}
	}
	for (uint32_t row = 1; row < h - 1; row++) {
		for (uint32_t col = 1; col < w - 1; col++) {
			pixels[row * w + col] = pixels2[row * w + col];
		}
	}
	free(pixels2);
}

void start_laplacian(char* filename, float k) {
	uint32_t w, h;
	uint8_t* pixels;
	discovery_dimensions(filename, &w, &h);
	pixels = (uint8_t*) malloc(w * h * sizeof(uint8_t));
	readTiff(filename, &pixels);
	laplacian(pixels, k, w, h);
	char* filenameout = (char*) malloc(strlen(filename) + 5);
	strcpy(filenameout, filename);
	strcat(filenameout, "_laplacian");
	writeTiff(filenameout, w, h, pixels);
	free(pixels);
}

void high_boost(uint8_t * pixels, float k, uint32_t w, uint32_t h) {
	// gaussian mask
	float sigma = 5.0;
	float K = 1.0;
	float mask[3][3] = {
		{K*exp(-1/sigma), K*exp(-1/sigma/2), K*exp(-1/sigma)},
		{K*exp(-1/sigma/2), K, K*exp(-1/sigma/2)},
		{K*exp(-1/sigma), K*exp(-1/sigma/2), K*exp(-1/sigma)}
	};
	float sum_mask = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			sum_mask += mask[i][j];
		}
	}
	uint8_t* pixels2 = (uint8_t*) malloc(w * h * sizeof(uint8_t));
	for (uint32_t row = 1; row < h - 1; row++) {
		for (uint32_t col = 1; col < w - 1; col++) {
			double sum = 0;
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					sum +=  pixels[(row + i) * w + (col + j)] * mask[i + 1][j + 1] / sum_mask;
				}
			}
			sum = pixels[row*w+col] - sum;
			sum = pixels[row*w+col] + k*sum;
			if (sum < 0) sum = 0;
			if (sum > 255) sum = 255;
			pixels2[row * w + col] = (uint8_t)sum;
		}
	}
	for (uint32_t row = 1; row < h - 1; row++) {
		for (uint32_t col = 1; col < w - 1; col++) {
			pixels[row * w + col] = pixels2[row*w+col];
		}
	}
	free(pixels2);
}

void start_high_boost(char* filename, float k ) {
	uint32_t w, h;
	uint8_t* pixels;
	discovery_dimensions(filename, &w, &h);
	pixels = (uint8_t*) malloc(w * h * sizeof(uint8_t));
	readTiff(filename, &pixels);
	high_boost(pixels, k, w, h);
	char* filenameout = (char*) malloc(strlen(filename) + 5);
	strcpy(filenameout, filename);
	strcat(filenameout, "_high-boost");
	writeTiff(filenameout, w, h, pixels);
	free(pixels);
}

int main(int argc, char* argv[]) {
	
	char* command = argv[1];
	char* filename = argv[2];
	if (strcmp("hist",command)==0){
		start_histogram(filename);
	} else if (strcmp("laplacian",command)==0){
		char* karg = argv[3];
		start_laplacian(filename, atof(karg));
	} else if (strcmp("high-boost",command)==0){
		char* karg = argv[3];
		start_high_boost(filename, atof(karg));
	} else {
		printf("Command not found\n");
	}

	return 0;
}
