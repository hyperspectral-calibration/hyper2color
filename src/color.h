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



/* Define CIE XYZ tristimulus values for reference white
   for various color temperatures for a standard 2 degree observer
 */

#define X_D50 96.422
#define Y_D50 100.0
#define Z_D50 82.521

#define X_D65 95.047
#define Y_D65 100.0
#define Z_D65 108.88

#define X_D75 94.972
#define Y_D75 100.0
#define Z_D75 122.638

/* Incandescent lamp
 */
#define X_CIE_A 109.850
#define Y_CIE_A 100.0
#define Z_CIE_A 35.585

#define X_CIE_C 98.074
#define Y_CIE_C 100.0
#define Z_CIE_C 118.232


static float XYZ_sRGB_matrix_D65[3][3] = {
  { 3.240479, -1.537150, -0.498535},
  {-0.969256,  1.875992,  0.041556},
  { 0.055648, -0.204043,  1.057311}
};

static float XYZ_sRGB_matrix_D50[3][3] = {
  { 3.1338561, -1.6168667, -0.4906146},
  {-0.9787684,  1.9161415,  0.0334540},
  { 0.0719453, -0.2289914,  1.4052427}
};

static float XYZ_AdobeRGB_matrix_D65[3][3] = {
  { 2.0413690, -0.5649464, -0.3446944},
  {-0.9692660,  1.8760108,  0.0415560},
  { 0.0134474, -0.1183897,  1.0154096}
};

static float XYZ_AdobeRGB_matrix_D50[3][3] = {
  { 1.9624274, -0.6105343, -0.3413404},
  {-0.9787684,  1.9161415,  0.0334540},
  { 0.0286869, -0.1406752,  1.3487655}
};

/*static float XYZ_DCIP3_matrix_D65[3][3] = {
  	  x 	  y 	  z
R 	0.68 	0.32 	0.00
G 	0.265 	0.69 	0.045
B 	0.15 	0.06 	0.79
*/

/* Color conversion functions
 */
void XYZ2LAB( float, float, float, float*, float*, float* );
void XYZ2RGB( float m[][3], float, float, float, float*, float*, float* );
void sRGB_Gamma( float*, float*, float* );
void AdobeRGB_Gamma( float*, float*, float* );
