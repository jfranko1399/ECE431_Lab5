#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#define ITERATIONS 30

void add_contour(unsigned char *image, int cols, int rows, int *col_loc, int *row_loc, int length, int name);
float *normalize_energy(float *energy);

#define SQR(x) ((x)*(x))

int main(int argc, char *argv[])

{
  FILE *fpt;
  FILE *fpt2;
  FILE *sobel_ppm;
  FILE *normal_ppm;
  unsigned char *image;
  unsigned char *initial;
  unsigned char *sobel;
  unsigned char *normal_image;
  int *init_cols, *init_rows, *final_cols, *final_rows;
  char header[20];
  int ROWS, COLS, BYTES;
  int contours, r, c, i, j, k, l, sumx, sumy, index;
  int max, min, oldRange, max2, min2, oldRange2, max3, min3, oldRange3;
  float minf;
  float averageDis, *firstInternalEnergy, *secondInternalEnergy, *externalEnergy;
  float *normalFirst, *normalSecond, *normalExternal, *sumEnergy;
  char ch;

  if (argc != 3)
  {
    printf("Usage: ./contour [image.ppm] [contour_points.txt]\n");
    exit(0);
  }

  /* read image */
  if ((fpt = fopen(argv[1], "rb")) == NULL)
  {
    printf("Unable to open hawk.ppm for reading\n");
    exit(0);
  }

  fscanf(fpt, "%s %d %d %d", header, &COLS, &ROWS, &BYTES);

  /* error checking PPM format */
  if (strcmp(header, "P5") != 0  ||  BYTES != 255)
  {
    printf("Not a greyscale 8-bit PPM image\n");
    exit(0);
  }

  /* allocate memory for stored image */
  image = (unsigned char *)calloc(ROWS * COLS, sizeof(unsigned char));
  header[0] = fgetc(fpt);	/* read white-space character that separates header */
  fread(image, 1, COLS * ROWS, fpt);
  fclose(fpt);

  if ((fpt2 = fopen(argv[2], "rb")) == NULL)
  {
    printf("Unable to open hawk_init.txt for reading\n");
    exit(0);
  }

  contours = 0;
  while (!feof(fpt2))
  {
    ch = fgetc(fpt2);
    if (ch == '\n')
    {
      contours++;
    }
  }

  rewind(fpt2);

  init_cols = (int *)calloc(contours, sizeof(int));
  init_rows = (int *)calloc(contours, sizeof(int));

  i = 0;
  while ((fscanf(fpt2, "%d %d\n", &init_cols[i], &init_rows[i])) != EOF)
  {
    i++;
  }

  fclose(fpt2);

  add_contour(image, COLS, ROWS, init_cols, init_rows, contours, 0);

  sobel = (unsigned char *)calloc(ROWS * COLS, sizeof(unsigned char));

  int sobelx[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
  int sobely[] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};

  for (r = 1; r < ROWS - 1; r++)
  {
    for (c = 1; c < COLS - 1; c++)
    {
      sumx = 0;
      sumy = 0;
      for (i = -1; i < 2; i++)
      {
        for (j = -1; j < 2; j++)
        {
          sumx += sobelx[(i + 1) * 3 + (j + 1)] * image[(r + i) * COLS + (c + j)];
          sumy += sobely[(i + 1) * 3 + (j + 1)] * image[(r + i) * COLS + (c + j)];
        }
      }

      sobel[r * COLS + c] = sqrt(SQR(sumx) + SQR(sumy));
    }
  }

  // sobel_ppm = fopen("sobel.ppm", "w");
  // fprintf(sobel_ppm, "P5 %d %d 255\n", COLS, ROWS);
  // fwrite(sobel, COLS * ROWS, 1, sobel_ppm);
  // fclose(sobel_ppm);

  max = 0;
  min = INT_MAX;

  /* find the old maximum and minimum pixel values */
  for (r = 0; r < ROWS; r++)
  {
    for (c = 0; c < COLS; c++)
    {
      if (sobel[r * COLS + c] > max)
      {
        max = sobel[r * COLS + c];
      }
      if (sobel[r * COLS + c] < min)
      {
        min = sobel[r * COLS + c];
      }
    }
  }

  /* determine the old range of pixel values */
  oldRange = max - min;

  normal_image = (unsigned char *)calloc(ROWS * COLS, sizeof(unsigned char));

  /* map all old pixels to new range of 0-255 */
  for (i = 0; i < ROWS * COLS; i++)
  {
    normal_image[i] = (((sobel[i] - min) * 255) / oldRange);
  }

  /* write out normalized image to see result */
  // normal_ppm = fopen("normal.ppm", "w");
  // fprintf(normal_ppm, "P5 %d %d 255\n", COLS, ROWS);
  // fwrite(normal_image, COLS * ROWS, 1, normal_ppm);
  // fclose(normal_ppm);

  firstInternalEnergy = (float *)calloc(49, sizeof(float));
  secondInternalEnergy = (float *)calloc(49, sizeof(float));
  externalEnergy = (float *)calloc(49, sizeof(float));

  normalFirst = (float *)calloc(49, sizeof(float));
  normalSecond = (float *)calloc(49, sizeof(float));
  normalExternal = (float *)calloc(49, sizeof(float));

  sumEnergy = (float *)calloc(49, sizeof(float));

  final_cols = (int *)calloc(contours, sizeof(int));
  final_rows = (int *)calloc(contours, sizeof(int));

  for (i = 0; i < ITERATIONS; i++)
  {
    averageDis = 0.0;
    for (j = 0; j < contours - 1; j++)
    {
      averageDis += sqrt(SQR(init_rows[j] - init_rows[j + 1]) + SQR(init_cols[j] - init_cols[j + 1]));
    }

    averageDis += sqrt(SQR(init_rows[j] - init_rows[0]) + SQR(init_cols[j] - init_cols[0]));
    averageDis /= contours;
    //printf("%f\n", averageDis);

    for (j = 0; j < contours - 1; j++)
    {
      for (r = -3; r < 4; r++)
      {
        for (c = -3; c < 4; c++)
        {
          firstInternalEnergy[(r + 3) * 7 + (c + 3)] = SQR((init_rows[j] + r) - init_rows[j + 1]) + SQR((init_cols[j] + c) - init_cols[j + 1]);
          secondInternalEnergy[(r + 3) * 7 + (c + 3)] = SQR(averageDis - sqrt(firstInternalEnergy[(r + 3) * 7 + (c + 3)]));
          externalEnergy[(r + 3) * 7 + (c + 3)] = SQR(255 - sobel[(init_rows[j] + r) * COLS + (init_cols[j] + c)]);
        }
      }
    }

    firstInternalEnergy[(r + 3) * 7 + (c + 3)] = SQR((init_rows[j] + r) - init_rows[0]) + SQR((init_cols[j] + c) - init_cols[0]);
    secondInternalEnergy[(r + 3) * 7 + (c + 3)] = SQR(averageDis - sqrt(firstInternalEnergy[(r + 3) * 7 + (c + 3)]));
    externalEnergy[(r + 3) * 7 + (c + 3)] = SQR(255 - sobel[(init_rows[j] + r) * COLS + (init_cols[j] + c)]);

    normalFirst = normalize_energy(firstInternalEnergy);
    normalSecond = normalize_energy(secondInternalEnergy);
    normalExternal = normalize_energy(externalEnergy);

    minf = FLT_MAX;

    for (j = 0; j < 7; j++)
    {
      for (k = 0; k < 7; k++)
      {
        printf("%f %f %f %f\n", normalFirst[j * 7 + k], normalSecond[j * 7 + k], normalExternal[j * 7 + k], minf);
        sumEnergy[j * 7 + k] = normalFirst[j * 7 + k] + normalSecond[j * 7 + k] + normalExternal[j * 7 + k];
        if (sumEnergy[j * 7 + k] < minf)
        {
          minf = sumEnergy[j * 7 + k];
          r = j;
          c = k;
        }
      }
    }

    for (j = 0; j < contours; j++)
    {
      //printf("%d %d %d %d\n", init_cols[j], init_rows[j], c, r);
      init_cols[j] = init_cols[j] + c;
      init_rows[j] = init_rows[j] + r;
    }
  }

  add_contour(image, COLS, ROWS, init_cols, init_rows, contours, 1);
}


void add_contour(unsigned char *image, int cols, int rows, int *col_loc, int *row_loc, int length, int name)
{
  FILE *contour_ppm;
  unsigned char *output;
  int i, j;

  output = (unsigned char *)calloc(cols * rows, sizeof(unsigned char));

  for (i = 0; i < cols * rows; i++)
  {
    output[i] = image[i];
  }

  for (i = 0; i < length - 1; i++)
  {
    for (j = -3; j < 4; j++)
    {
      output[(row_loc[i] + j) * cols + col_loc[i]] = 0;
      output[row_loc[i] * cols + (col_loc[i] + j)] = 0;
    }
  }

  if (name == 0)
  {
    contour_ppm = fopen("init_contour.ppm", "w");
  }
  else
  {
    contour_ppm = fopen("final_contour.ppm", "w");
  }

  fprintf(contour_ppm, "P5 %d %d 255\n", cols, rows);
  fwrite(output, cols * rows, 1, contour_ppm);
  fclose(contour_ppm);
}


float *normalize_energy(float *energy)
{
  int r, c, i, max = 0, min = INT_MAX, oldRange;
  float *normal;

  /* find the old maximum and minimum pixel values */
  for (r = 0; r < 7; r++)
  {
    for (c = 0; c < 7; c++)
    {
      if (energy[r * 7 + c] > max)
      {
        max = energy[r * 7 + c];
      }
      if (energy[r * 7 + c] < min)
      {
        min = energy[r * 7 + c];
      }
    }
  }

  /* determine the old range of pixel values */
  oldRange = max - min;

  normal = (float *)calloc(49, sizeof(float));

  /* map all old pixels to new range of 0-255 */
  for (i = 0; i < 49; i++)
  {
    normal[i] = (((energy[i] - min) * 1) / oldRange);
  }

  return normal;
}
