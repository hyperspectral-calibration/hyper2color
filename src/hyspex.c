/*
    Hyspex header routines

    Copyright (c) 2015-2021, Ruven <ruven@users.sourceforge.net>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hyspex.h"

#define HYSPEX_MAGIC "HYSPEX\0\0"
#define HYSPEX_SIZE 8
#define HYSPEX_BANDS 1961
#define HYSPEX_WIDTH 1965
#define HYSPEX_SCANLINES 2073
#define HYSPEX_WAVELENGTHS 2181



int is_hyspex( FILE* s, hyspex_header *header )
{
  /* Rewind our file if necessary
   */
  rewind( s );

  /* Hyspex Magic
   */
  unsigned char magic[8];
  if( fread( magic, 1, 8, s ) != 8 ){
    printf("Unable to read header\n");
    return 1;
  }

  if( memcmp( magic, HYSPEX_MAGIC, 8 ) == 0 ) return 1;
  else return 0;

}



int parse_envi_header( FILE *s, hyspex_header *header )
{
  ssize_t read;
  char *p;
  size_t len = 0;
  char* line = NULL;

  rewind( s );

/*   while( (read = getline(&line, &len, s) ) != -1 ){ */

/*     char *name, *value; */
/*     fscanf( s, "%s = %s",  ); */

/*     p = NULL; */
/*     p = strtok( line, "\t" ); */
/*     while( p ){ */
/*       p = strtok(NULL,"\t"); */
/*       if(p){ */
/* 	radiometric[wavelengths][n] = atof(p); */
/* 	n++; */
/*       } */
/*     } */
/*     n = 0; */
/*     wavelengths++; */
/*   } */



}


int parse_hyspex_header( FILE *s, hyspex_header *header )
{
  long int hh;

  /* Rewind our file if necessary
   */
  rewind( s );

  /* Hyspex Magic
   */
  unsigned char magic[8];
  if( fread( magic, 1, 8, s ) != 8 ){
    printf("Unable to read header\n");
    return 1;
  }
  if( memcmp( magic, HYSPEX_MAGIC, 8 ) != 0 ){
    printf("%s: Not a Hyspex file\n", magic );
  }
  
  /* Header size
   */
  fseek( s, HYSPEX_SIZE, SEEK_SET );
  if( fread( &hh, 4, 1, s ) != 1 ){
    printf("Unable to read header\n");
    return 1;
  }
  header->size = hh;

  /* Number of bands
   */
  fseek( s, HYSPEX_BANDS, SEEK_SET );
  if( fread( &hh, 4, 1, s ) != 1 ){
    printf("Unable to read header\n");
    return 1;
  }
  header->bands = hh;

  /* Number of samples
   */
  if( fread( &hh, 4, 1, s ) != 1 ){
    printf("Unable to read header\n");
    return 1;
  }
  header->samples = hh;

  /* Number of scanlines
   */
  fseek( s, HYSPEX_SCANLINES, SEEK_SET );
  if( fread( &hh, 4, 1, s ) != 1 ){
    printf("Unable to read header\n");
    return 1;
  }
  header->scanlines = hh;

  /* Number of bits per pixel
   */
  if( fread( &hh, 4, 1, s ) != 1 ){
    printf("Unable to read header\n");
    return 1;
  }
  header->bpp = 2;


  /* Extract list of wavelengths
   */
  header->wavelengths = malloc( sizeof(double) * header->bands );
  fseek( s, HYSPEX_WAVELENGTHS, SEEK_SET );
  int n;
  double w;
  for( n=0; n<header->bands; n++ ){
    if( fread( &w, 8, 1, s ) != 1 ){
      printf("Unable to read header\n");
      return 1;
    }
    header->wavelengths[n] = w;
  }


  /* Extract responsivity per pixel element
   */
  /*  header->responsivities = malloc( sizeof(double) * header->bands * header->samples );*/
  for( n=0; n<header->bands*header->samples; n++ ){
    if( fread( &w, 8, 1, s ) != 1 ){
      printf("Unable to read header\n");
      return 1;
    }
/*     header->responsivities[n] = w; */
  }

  /* Extract quantum efficiency per band
   */
  header->QE = malloc( sizeof(double) * header->bands );
  for( n=0; n<header->bands; n++ ){
    if( fread( &w, 8, 1, s ) != 1 ){
      printf("Unable to read header\n");
      return 1;
    }
    header->QE[n] = w;
  }

  /* Extract background
   */
  /* header->background = malloc( sizeof(double) * header->bands * header->samples ); */
  for( n=0; n<header->bands*header->samples; n++ ){
    if( fread( &w, 8, 1, s ) != 1 ){
      printf("Unable to read header\n");
      return 1;
    }
/*     header->background[n] = w; */
  }


  /* Rewind back to the beginning
   */
  rewind( s );


  return 0;
}



/* Load a spectral curve. Assume BIL
 */
int load_hyspex_pixel( FILE* s, hyspex_header *header, double *spectrum, int x, int y )
{
  unsigned short pixel;
  size_t index = header->size + ((size_t)y * (size_t)header->bpp * (size_t)header->samples * (size_t)header->bands) + (x * (size_t)header->bpp);

  unsigned int n = 0;

  for( n=0; n<header->bands; n++ ){
    fseek( s, index, SEEK_SET );
    if( fread( &pixel, (size_t)header->bpp, 1, s ) != 1 ){
      printf("Unable to read pixel data for band %d\n", n);
      return 1;
    }
    if( header->bpp == 2 ) spectrum[n] = (double)pixel / 65535.0;
    else spectrum[n] = (double)pixel;

    index += (size_t)header->samples * (size_t)header->bpp;
  }

  return 0;
}



/* Load a spectral curve. Assume BIL. Return number of pixels loaded
 */
int load_hyspex_bil( FILE* s, hyspex_header *header, void *buffer, int y )
{
  int size;

  size_t index = header->size + ((size_t)y * (size_t)header->bpp * (size_t)header->samples * (size_t)header->bands);

  fseek( s, index, SEEK_SET );
  size = fread( buffer, (size_t)header->bpp, header->samples*header->bands, s );

  return size;
}


void update_width( hyspex_header *header, int new_width )
{
  unsigned char *ptr = (unsigned char*) header;
  ptr[HYSPEX_WIDTH] = (int) new_width;
}


/* Free memory
 */
void free_hyspex( hyspex_header *header ){

  /* Free our memory
   */
  free( header->wavelengths );
  free( header->QE );
}
