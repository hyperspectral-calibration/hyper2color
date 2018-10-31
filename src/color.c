/*
    Copyright (C) 2015-2018 Ruven Pillay <ruven@users.sourceforge.net>

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

#include <math.h>
#include "color.h"



/* CIE XYZ->LAB conversion
 */
void XYZ2LAB( float X, float Y, float Z, float *L, float *a, float *b ){

  float var_X = X / X_D65;
  float var_Y = Y / Y_D65;
  float var_Z = Z / Z_D65;

  if ( var_X > 0.008856 ) var_X = cbrtf(var_X);
  else                    var_X = ( 7.787 * var_X ) + ( 16.0 / 116.0 );

  if ( var_Y > 0.008856 ) var_Y = cbrtf(var_Y);
  else                    var_Y = ( 7.787 * var_Y ) + ( 16.0 / 116.0 );

  if ( var_Z > 0.008856 ) var_Z = cbrtf(var_Z);
  else                    var_Z = ( 7.787 * var_Z ) + ( 16.0 / 116.0 );

  *L = ( 116.0 * var_Y ) - 16.0;
  *a = 500.0 * ( var_X - var_Y );
  *b = 200.0 * ( var_Y - var_Z );

  // Clip our L values to 0->100
  *L = (*L <= 0.0 ? 0.0 : *L >= 100.0 ? 100.0 : *L);
}



/* Apply sRGB gamma of 2.4
 */
static float _sRGB_Gamma( float c ){
  float alpha = 0.055;
  if( c <= 0.0031308 ) return 12.92 * c;
  else return (1 + alpha) * powf(c, 1/2.4) - alpha;
}
void sRGB_Gamma( float *R, float *G, float *B ){
  *R = _sRGB_Gamma(*R);
  *G = _sRGB_Gamma(*G);
  *B = _sRGB_Gamma(*B);
}


/* Apply AdobeRGB gamma of 563/256 (2.199)
 */
void AdobeRGB_Gamma( float *R, float *G, float *B ){
  float gamma = 256.0/563.0;
  *R = powf(*R,gamma);
  *G = powf(*G,gamma);
  *B = powf(*B,gamma);
}



/* Convert XYZ to RGB using a 3x3 Matrix
 */
void XYZ2RGB( float matrix[][3], float X, float Y, float Z, float *R, float *G, float *B ){

  /* XYZ -> linear RGB matrix
   */
  X = X / 100.0;    /* X from 0 to 95.047 for D65 */
  Y = Y / 100.0;    /* Y from 0 to 100.000 for D65 */
  Z = Z / 100.0;    /* Z from 0 to 108.883 for D65 */
 
  /* Convert to normalized linear RGB
   */
  float r = X * matrix[0][0] + Y * matrix[0][1] + Z * matrix[0][2];
  float g = X * matrix[1][0] + Y * matrix[1][1] + Z * matrix[1][2];
  float b = X * matrix[2][0] + Y * matrix[2][1] + Z * matrix[2][2];

  /* Clip to range 0.0 -> 1.0
   */
  *R = ( r <= 0.0 ? 0.0 : r >= 1.0 ? 1.0 : r );
  *G = ( g <= 0.0 ? 0.0 : g >= 1.0 ? 1.0 : g );
  *B = ( b <= 0.0 ? 0.0 : b >= 1.0 ? 1.0 : b );
}

