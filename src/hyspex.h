/*
    Hyspex header structure

    Copyright (c) 2015-2018, Ruven <ruven@users.sourceforge.net>

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


#include <stdio.h>


/* Hyspex header structure
 */
typedef struct {
  int size;
  unsigned int bands;
  unsigned int samples;
  unsigned int scanlines;
  unsigned int bpp;
  double *wavelengths;

  double *responsivities;
  double *QE;
  double *background;

} hyspex_header;


int parse_hyspex_header( FILE*, hyspex_header* );
int load_hyspex_pixel( FILE*, hyspex_header*, double*, int, int );
int load_hyspex_bil( FILE*, hyspex_header*, void*, int );
void free_hyspex( hyspex_header* );
