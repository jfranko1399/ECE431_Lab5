#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

void add_contour(unsigned char *image, int cols, int rows, int *col_loc, int *row_loc, int name);
float *normalize_energy(float *energy);

#define SQR(x) ((x)*(x))

int main(int argc, char *argv[])

{
  FILE *fpt;
  FILE *fpt2;
  FILE *sobel_ppm;
  FILE *normal_ppm;
  FILE *final;
  unsigned char *image;
  unsigned char *initial;
  int *sobel;
  unsigned char *normal_image;
  int *init_cols, *init_rows, *add_cols, *add_rows;
  char header[20];
  int ROWS, COLS, BYTES;
  int r, c, i, j, k, l, sumx, sumy, index;
  int max, min, oldRange;
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

  // read initial contours
  if ((fpt2 = fopen(argv[2], "rb")) == NULL)
  {
    printf("Unable to open hawk_init.txt for reading\n");
    exit(0);
  }

  init_cols = (int *)calloc(42, sizeof(int));
  init_rows = (int *)calloc(42, sizeof(int));

  // store contour points
  i = 0;
  while ((fscanf(fpt2, "%d %d\n", &init_cols[i], &init_rows[i])) != EOF)
  {
    i++;
  }

  fclose(fpt2);

  // draw initial contour
  add_contour(image, COLS, ROWS, init_cols, init_rows, 0);

  sobel = (int *)calloc(ROWS * COLS, sizeof(int));

  int sobelx[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
  int sobely[] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};

  // calculate sobel filter
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

  firstInternalEnergy = (float *)calloc(49, sizeof(float));
  secondInternalEnergy = (float *)calloc(49, sizeof(float));
  externalEnergy = (float *)calloc(49, sizeof(float));

  normalFirst = (float *)calloc(49, sizeof(float));
  normalSecond = (float *)calloc(49, sizeof(float));
  normalExternal = (float *)calloc(49, sizeof(float));

  sumEnergy = (float *)calloc(49, sizeof(float));

  add_cols = (int *)calloc(42, sizeof(int));
  add_rows = (int *)calloc(42, sizeof(int));

  // 30 iterations of contour
  i = 0;
  while (i < 30)
  {
    // average distance between points
    averageDis = 0.0;
    for (j = 0; j < 41; j++)
    {
      averageDis += sqrt(SQR(init_rows[j] - init_rows[j + 1]) + SQR(init_cols[j] - init_cols[j + 1]));
    }

    averageDis += sqrt(SQR(init_rows[j] - init_rows[0]) + SQR(init_cols[j] - init_cols[0]));
    averageDis /= 42;

    // calculate energies for each contour point
    for (j = 0; j < 42; j++)
    {
      // check in 7x7 box
      for (r = -3; r < 4; r++)
      {
        for (c = -3; c < 4; c++)
        {
          if ((j + 1) == 42)
          {
            firstInternalEnergy[(r + 3) * 7 + (c + 3)] = SQR((init_rows[j] + r) - init_rows[0]) + SQR((init_cols[j] + c) - init_cols[0]);
            secondInternalEnergy[(r + 3) * 7 + (c + 3)] = SQR(averageDis - sqrt(firstInternalEnergy[(r + 3) * 7 + (c + 3)]));
            externalEnergy[(r + 3) * 7 + (c + 3)] = SQR(255 - sobel[(init_rows[j] + r) * COLS + (init_cols[j] + c)]);
          }
          else
          {
            firstInternalEnergy[(r + 3) * 7 + (c + 3)] = SQR((init_rows[j] + r) - init_rows[j + 1]) + SQR((init_cols[j] + c) - init_cols[j + 1]);
            secondInternalEnergy[(r + 3) * 7 + (c + 3)] = SQR(averageDis - sqrt(firstInternalEnergy[(r + 3) * 7 + (c + 3)]));
            externalEnergy[(r + 3) * 7 + (c + 3)] = SQR(255 - sobel[(init_rows[j] + r) * COLS + (init_cols[j] + c)]);
          }
        }
      }

      // normalize energies
      normalFirst = normalize_energy(firstInternalEnergy);
      normalSecond = normalize_energy(secondInternalEnergy);
      normalExternal = normalize_energy(externalEnergy);

      minf = FLT_MAX;
      r = c = 0;

      // sum energies, store minimum point in addition to row and col location
      for (l = -3; l < 4; l++)
      {
        for (k = -3; k < 4; k++)
        {
          sumEnergy[(l + 3) * 7 + (k + 3)] = 3 * (normalFirst[(l + 3) * 7 + (k + 3)] + normalSecond[(l + 3) * 7 + (k + 3)]) + (2 * normalExternal[(l + 3) * 7 + (k + 3)]);
          if (sumEnergy[(l + 3) * 7 + (k + 3)] < minf)
          {
            minf = sumEnergy[(l + 3) * 7 + (k + 3)];
            add_rows[j] = l;
            add_cols[j] = k;
          }
        }
      }
    }

    // update contour points with new pixel locations
    for (j = 0; j < 42; j++)
    {
      init_cols[j] += add_cols[j];
      init_rows[j] += add_rows[j];
    }

    i++;
  }

  // print final contour points in csv file
  i = 0;
  final = fopen("hawk_final.csv", "w");
  while (i < 42)
  {
    fprintf(final, "%d,%d\n", init_cols[i], init_rows[i]);
    i++;
  }
  fclose(final);

  // draw final contour
  add_contour(image, COLS, ROWS, init_cols, init_rows, 1);
}

/* Draw contour points on hawk image */
void add_contour(unsigned char *image, int cols, int rows, int *col_loc, int *row_loc, int name)
{
  FILE *contour_ppm;
  unsigned char *output;
  int i, j;

  output = (unsigned char *)calloc(cols * rows, sizeof(unsigned char));

  // copy image
  for (i = 0; i < cols * rows; i++)
  {
    output[i] = image[i];
  }

  // draw plus at each contour point
  for (i = 0; i < 42; i++)
  {
    for (j = -3; j < 4; j++)
    {
      output[(row_loc[i] + j) * cols + col_loc[i]] = 0;
      output[row_loc[i] * cols + (col_loc[i] + j)] = 0;
    }
  }

  // change file name depending on time function called
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

/* Normalize energy input to a new range of 0 - 1 */
float *normalize_energy(float *energy)
{
  int r, c, i, max = 0, min = INT_MAX, oldRange;
  float *normal;

  /* find the old maximum and minimum energy values */
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

  /* determine the old range of energy values */
  oldRange = max - min;

  normal = (float *)calloc(49, sizeof(float));

  /* map all old pixels to new range of 0-1 */
  for (i = 0; i < 49; i++)
  {
    normal[i] = (((energy[i] - min) * 1) / oldRange);
  }

  return normal;
}
