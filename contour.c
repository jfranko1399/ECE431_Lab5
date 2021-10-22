#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

void add_contour(unsigned char *image, int cols, int rows, int *col_loc, int *row_loc, int length, int name);

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
  int *init_cols, *init_rows;
  char header[20];
  int ROWS, COLS, BYTES;
  int contours, r, c, i, j, k, l, sumx, sumy, max, min, oldRange;
  float distance;
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

  for (i = 0; i < 31; i++)
  {
    distance = 0.0;
    for (j = 0; j < contours - 1; j++)
    {
      distance += sqrt(SQR(init_rows[j] - init_rows[j + 1]) + SQR(init_cols[i] - init_rows[i + 1]));
    }

    distance += sqrt(SQR(init_rows[j] - init_rows[0]) + SQR(init_cols[i] - init_rows[0]));
    distance /= contours;

    for (j = 0; j < contours; j++)
    {

    }
  }
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
