/*
    Generate a true color CIE L*a*b*, sRGB or AdobeRGB image from a hyperspectral cube
    using a given illuminant.


    Copyright (C) 2015-2019 Ruven Pillay <ruven@users.sourceforge.net>

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



#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

#include <math.h>
#include <sys/mman.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "tiffio.h"

/* Load our colorimetric matrices and functions
 */
#include "color.h"
#include "CIE-tristimulus.c"
#include "CIE-D65.c"
#include "CIE-A.c"
#include "icc.c"


/* Load our hyspex header library
 */
#include "hyspex.h"



/* Print help message
 */
void help( void ){
  printf( "\n \
 Generate a true color CIE L*a*b*, sRGB or AdobeRGB image from a hyperspectral cube\n \
 using a given illuminant.\n\n \
 Copyright (c) 2015-2018, Ruven <ruven@users.sourceforge.net>\n \
 Licence: GNU GENERAL PUBLIC LICENSE (GPL)\n\n \
 Generate colorimetric rendering of hyperspectral reflectance data with \n \
 requested illuminant or color temperature in sRGB, AdobeRGB or CIE L*a*b* color \n \
 output space. Output is in TIFF format using requested compression.\n\n \
 Requested color temperature can be a standard temperature (D65, D50, D75, D93, A)\n \
 or a temperature in degrees Kelvin (eg: 5000 for a temperature of 5000K)\n \
 Alternatively, an artibrary illuminant power spectrum can be provided \n \
 containing a list of wavelengths and power values\n\n \
 Output bits per channel can be 8, 16 or 32 bits, where 8 and 16 are encoded\n \
 as unsigned integer and 32 is encoded as floating point\n\n \
 eg: hyper2color -i data.img -o calibrated_color.tif -t D65 \n\n \
 Options:\n\n \
  --input,       -i:  input hyperspectral cube\n \
  --output,      -o:  output TIFF image\n \
  --temperature, -t:  output color temperature: D65 (default), D50 or temperature in K)\n \
  --colorspace,  -s:  output color space: CIELAB, sRGB (default) or AdobeRGB\n \
  --power,       -p:  illuminant power spectrum file (optional)\n \
  --bits,        -b:  output bits per channel: 8 (default), 16 or 32 bits\n \
  --width,       -x:  hyperspectral image width\n \
  --height,      -y:  hyperspectral image height\n \
  --channels,    -c:  number of bands in hyperspectral cube\n \
  --compression  -m:  TIFF output compression: none (default), deflate, lzw or jpeg\n \
  --help,        -h:  this help message\n \
  --verbose,     -v:  verbose output\n\n\n" );
}




/* Function to calculate the blackbody spectrum for a particular temperature
   for the desired wavelength, which should be given in nm
*/
double calculate_power_spectrum( int bbTemp, int wavelength ){

  double power;
  double wlm = wavelength * 1e-9;  /* Convert from nm to meters */
  power = 3.74177152e-16  /
    ( pow( wlm, 5 ) * ( exp( 0.0143877696 / (wlm * bbTemp) ) - 1.0 ) );

  return power;
}



void calculate_XYZ( hyspex_header *header, double *spectrum, double **power_spectrum, double **cie_color_match )
{

  /* Interpolate using GSL
   */
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc( gsl_interp_cspline, header->bands );
  gsl_spline_init( spline, header->wavelengths, spectrum, header->bands );


  int firstwav = 410;

  int tr = firstwav - cie_color_match[0][0];
  int te = firstwav - power_spectrum[0][0];


  double X, Y, Z, norm;
  X = Y = Z = norm = 0.0;

  int k;
  for( k=0; k<=830-firstwav; k+=5 ){

    double val = gsl_spline_eval(spline, k, acc);
    if( val < 0 ) val = 0.0;

    X += val * cie_color_match[tr+k][1] * power_spectrum[te+k][1];
    Y += val * cie_color_match[tr+k][2] * power_spectrum[te+k][1];
    Z += val * cie_color_match[tr+k][3] * power_spectrum[te+k][1];
    norm += cie_color_match[tr+k][2] * power_spectrum[te+k][1];
  }

  /* Free our interpolator
   */
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  float XX = (float)( (X * 100.0) / (norm) );
  float YY = (float)( (Y * 100.0) / (norm) );
  float ZZ = (float)( (Z * 100.0) / (norm) );

}

typedef struct { gsl_spline* s; gsl_interp_accel *a; double *cie; double *power; } my_f_params;


double f( double x, void *p ){
   my_f_params *params = (my_f_params*) p;
   gsl_spline* spline = (gsl_spline*) params->s;
   gsl_interp_accel* acc = (gsl_interp_accel*) params->a;
   double* match = (double*) params->cie;
   double* power = (double*) params->power;

   double result = gsl_spline_eval(spline, x, acc);
   result = result * match[(int)x-360] * power[(int)x-360];

   return result;
}


int main( int argc, char *argv[] )
{
  int c;
  int digit_optind = 0;
  int i, n, j, k;

  int width = 0;
  int height = 0;
  int bands = 0;

  int verbose = 0;
  FILE *in = NULL;
  TIFF *out = NULL;

  /* Output color space (default: sRGB)
   */
  uint16 colorspace = PHOTOMETRIC_RGB;

  /* Output color profile - 0: None, 1: sRGB, 2: AdobeRGB
   */
  unsigned int icc_profile = 1;

  /* Rendered color temperature (default: D65, 6504K)
   */
  int temperature = 6504;

  /* Output bits per channel (default: 8)
   */
  int bpc = 0;

  /* Output compression
   */
  short compression = COMPRESSION_NONE;

  /* Parse our options
   */
  while( 1 ) {

    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"input", 1, 0, 'i'},
      {"output", 1, 0, 'o'},
      {"temperature", 1, 0, 't'},
      {"colorspace", 1, 0, 's'},
      {"bits", 1, 0, 'b'},
      {"width", 1, 0, 'x'},
      {"height", 1, 0, 'y'},
      {"channels", 1, 0, 'c'},
      {"compression", 1, 0, 'm'},
      {"help", 0, 0, 'h'},
      {"verbose", 0, 0, 'v'},
      {0, 0, 0, 0}
    };

    c = getopt_long( argc, argv, "i:o:t:s:b:x:y:c:m:vh", long_options, &option_index );

    if( c == -1 ){
      break;
    }

    switch( c ){

    case 'i':
      /* Our input image
       */
      if( ! ( in = fopen( optarg, "rb" ) ) ){
	help();
	printf( "Unable to open input image file: '%s'\n\n", optarg );
	exit( 1 );
      }
      break;

    case 'o':
      /* Out output image
       */
      if( ! ( out = TIFFOpen( optarg, "w" ) ) ){
	help();
	printf( "Unable to open output image file: '%s'\n\n", optarg );
	exit( 1 );
      }
      break;

    case 's':
      /* Output color space - sRGB by default
       */
      if( strncasecmp( optarg, "CIELAB", 64 ) == 0 ){
	colorspace = PHOTOMETRIC_CIELAB;
	icc_profile = 0;
      }
      else if( strncasecmp( optarg, "AdobeRGB", 64 ) == 0 ){
	colorspace = PHOTOMETRIC_RGB;
	icc_profile = 2;
      }
      else{
	colorspace = PHOTOMETRIC_RGB;
	icc_profile = 1;
      }

      break;

    case 't':
      /* Requested color temperature
       */
      if( strcmp( optarg, "D93" ) == 0 ) temperature = 9300;
      else if( strcmp( optarg, "D75" ) == 0 ) temperature = 7500;
      else if( strcmp( optarg, "D65" ) == 0 ) temperature = 6504;
      else if( strcmp( optarg, "D50" ) == 0 ) temperature = 5000;
      else if( strcmp( optarg, "A" ) == 0 ) temperature = -1;
      else temperature = atoi( optarg );
      break;

    case 'b':
      /* Output bits per channel
       */
      n = atoi( optarg );
      if( n == 8 || n == 16 || n == 32 ){
	bpc = n;
      }
      else printf( "Unsupported bit depth '%s': defaulting to 8\n", optarg );
      break;

    case 'x':
      /* Width
       */
      width = atoi( optarg );
      break;

    case 'y':
      /* Height
       */
      height = atoi( optarg );
      break;

    case 'c':
      /* Channels
       */
      bands = atoi( optarg );
      break;

    case 'm':
      /* Output compression
       */
      if( strcasecmp( optarg, "deflate" ) == 0 ) compression = COMPRESSION_DEFLATE;
      if( strcasecmp( optarg, "lzw" ) == 0 ) compression = COMPRESSION_LZW;
      if( strcasecmp( optarg, "jpeg" ) == 0 ) compression = COMPRESSION_JPEG;
      break;

    case 'h':
      help();
      exit( 0 );

    case 'v':
      verbose = 1;
      break;

    default:
      help();
      printf ("?? getopt returned character code 0%o ??\n", c);
      exit( 1 );
    }

  }


  /* Make sure we have properly intialized some stuff
   */
  if( !in || !out ){
    help();
    if( !in ) printf( "No input image specified\n" );
    if( !out ) printf( "No output image specified\n" );
    printf( "\n" );
    exit( 1 );
  }



  /* Extract info from hyspex header
   */
  hyspex_header header;

  if( width==0 && height==0 && bands==0 ){
    parse_hyspex_header( in, &header );
  }


  const double wavelengths40[40] = {412.880826, 427.454604, 442.028383, 456.602161, 471.175940, 485.749718,
				    500.323497, 514.897275, 529.471054, 544.044832, 558.618611, 573.192389,
				    587.766168, 602.339946, 616.913725, 631.487503, 646.061282, 660.635060,
				    675.208839, 689.782617, 704.356396, 718.930174, 733.503953, 748.077731,
				    762.651510, 777.225288, 791.799067, 806.372845, 820.946624, 835.520402,
				    850.094181, 864.667959, 879.241737, 893.815516, 908.389294, 922.963073,
				    937.536851, 952.110630, 966.684408, 981.258187 };

  const double wavelengths80[80] = { 416.524261, 423.811157, 431.098053, 438.384949, 445.671814, 452.958710,
				     460.245605, 467.532501, 474.819397, 482.106262, 489.393158, 496.680054,
				     503.966949, 511.253845, 518.540710, 525.827637, 533.114502, 540.401367,
				     547.688293, 554.975159, 562.262085, 569.548950, 576.835815, 584.122742,
				     591.409607, 598.696472, 605.983398, 613.270264, 620.557190, 627.844055,
				     635.130920, 642.417847, 649.704712, 656.991638, 664.278503, 671.565369,
				     678.852295, 686.139160, 693.426086, 700.712952, 707.999817, 715.286743,
				     722.573608, 729.860535, 737.147400, 744.434265, 751.721191, 759.008057,
				     766.294983, 773.581848, 780.868713, 788.155640, 795.442505, 802.729370,
				     810.016296, 817.303162, 824.590088, 831.876953, 839.163818, 846.450745,
				     853.737610, 861.024536, 868.311401, 875.598267, 882.885193, 890.172058,
				     897.458984, 904.745850, 912.032715, 919.319641, 926.606506, 933.893433,
				     941.180298, 948.467163, 955.754089, 963.040955, 970.327881, 977.614746,
				     984.901611, 992.188538 };

  const double wavelengths160[160] = { 414.702548, 418.345993, 421.989437, 425.632882, 429.276326, 432.919771,
				       436.563216, 440.206660, 443.850105, 447.493550, 451.136994, 454.780439,
				       458.423883, 462.067328, 465.710773, 469.354217, 472.997662, 476.641106,
				       480.284551, 483.927996, 487.571440, 491.214885, 494.858330, 498.501774,
				       502.145219, 505.788663, 509.432108, 513.075553, 516.718997, 520.362442,
				       524.005887, 527.649331, 531.292776, 534.936220, 538.579665, 542.223110,
				       545.866554, 549.509999, 553.153444, 556.796888, 560.440333, 564.083777,
				       567.727222, 571.370667, 575.014111, 578.657556, 582.301001, 585.944445,
				       589.587890, 593.231334, 596.874779, 600.518224, 604.161668, 607.805113,
				       611.448558, 615.092002, 618.735447, 622.378891, 626.022336, 629.665781,
				       633.309225, 636.952670, 640.596115, 644.239559, 647.883004, 651.526448,
				       655.169893, 658.813338, 662.456782, 666.100227, 669.743672, 673.387116,
				       677.030561, 680.674005, 684.317450, 687.960895, 691.604339, 695.247784,
				       698.891229, 702.534673, 706.178118, 709.821562, 713.465007, 717.108452,
				       720.751896, 724.395341, 728.038786, 731.682230, 735.325675, 738.969119,
				       742.612564, 746.256009, 749.899453, 753.542898, 757.186343, 760.829787,
				       764.473232, 768.116676, 771.760121, 775.403566, 779.047010, 782.690455,
				       786.333900, 789.977344, 793.620789, 797.264233, 800.907678, 804.551123,
				       808.194567, 811.838012, 815.481457, 819.124901, 822.768346, 826.411790,
				       830.055235, 833.698680, 837.342124, 840.985569, 844.629014, 848.272458,
				       851.915903, 855.559347, 859.202792, 862.846237, 866.489681, 870.133126,
				       873.776571, 877.420015, 881.063460, 884.706904, 888.350349, 891.993794,
				       895.637238, 899.280683, 902.924128, 906.567572, 910.211017, 913.854461,
				       917.497906, 921.141351, 924.784795, 928.428240, 932.071685, 935.715129,
				       939.358574, 943.002018, 946.645463, 950.288908, 953.932352, 957.575797,
				       961.219242, 964.862686, 968.506131, 972.149575, 975.793020, 979.436465,
				       983.079909, 986.723354, 990.366799, 994.010243 };


  if( width != 0 ) header.samples = width;
  if( height != 0 ) header.scanlines = height;
  if( bands != 0 ){
    header.bands = bands;
    //header.wavelengths = malloc( sizeof(float) * bands );
    if( bands == 40 ) header.wavelengths = (double*) wavelengths40;
    else if( bands == 80 ) header.wavelengths = (double*) wavelengths80;
    else if( bands == 160 ) header.wavelengths = (double*) wavelengths160;
    header.bpp = 2;
  }


  if( verbose ){
    printf( "Hyspex header size %d bytes\n", header.size );
    printf( "Hyperspectral data cube: %dx%d pixels, %d bands\n", header.samples, header.scanlines, header.bands );
    unsigned char* space = "CIE L*a*b*";
    if( icc_profile == 1 ) space = "sRGB";
    else if( icc_profile == 2 ) space = "AdobeRGB";
    printf( "Output color space: %s\n", space );
    printf( "Output color temperature: %d Kelvin\n", temperature );
    //    printf( "Output bits per pixel: %d\n", bpc );
  }


  unsigned short *scanline_spectrum;
  scanline_spectrum = malloc( header.samples * sizeof(unsigned short) * header.bands );

  double spectrum[320];

  /* Load up our illuminant power spectrum
   */
  double power_spectrum[531][2];
  for( k=300; k<=830; k++ ){

    power_spectrum[k-300][0] = k;

    if( temperature == -1 ) power_spectrum[k-300][1] = CIE_A[k-300][1];
    else if( temperature == 6504 ) power_spectrum[k-300][1] = D65[k-300][1];
    else power_spectrum[k-300][1] = calculate_power_spectrum( temperature, k );
  }


  int firstwav = ceil( header.wavelengths[0] );

  int tr = firstwav - cie_color_match[0][0];
  int te = firstwav - power_spectrum[0][0];


  /* Calculate CIE normalization
   */
  double norm = 0.0;
  for( k=0; k<=830-firstwav; k++ ){
    norm += cie_color_match[tr+k][2] * power_spectrum[te+k][1];
  }
  double power[471];
  for( k=0; k<471; k++ ){
    power[k] = power_spectrum[k+60][1];
  }


  double X_match[471];
  for( k=0; k<471; k++ ){
    X_match[k] = cie_color_match[k][1];
  }
  double Y_match[471];
  for( k=0; k<471; k++ ){
    Y_match[k] = cie_color_match[k][2];
  }
  double Z_match[471];
  for( k=0; k<471; k++ ){
    Y_match[k] = cie_color_match[k][3];
  }



  /* First define our bits per sample and format
   */
  int bits_per_sample = 8;
  short sample_format = SAMPLEFORMAT_UINT;

  if( bpc == 32 ){
    bits_per_sample = 32;
    sample_format = SAMPLEFORMAT_IEEEFP;
  }
  else if( bpc == 16 ){
    bits_per_sample = 16;
    sample_format = SAMPLEFORMAT_UINT;
  }
  else{
    bits_per_sample = 8;
    sample_format = SAMPLEFORMAT_UINT;
  }


  /* Allocate memory for the color values of a single scan line
   */
  void *calculated_color = NULL;

  if( bits_per_sample == 32 ) calculated_color = malloc( sizeof(float)*header.samples*3 );
  else if( bits_per_sample == 16 ) calculated_color = malloc( sizeof(unsigned short)*header.samples*3 );
  else calculated_color = malloc( sizeof(unsigned char)*header.samples*3 );


  /* Set basic TIFF metadata tags
   */
  TIFFSetField( out, TIFFTAG_IMAGEWIDTH, header.samples );          // set the width of the image
  TIFFSetField( out, TIFFTAG_IMAGELENGTH, header.scanlines );       // set the height of the image
  TIFFSetField( out, TIFFTAG_SAMPLESPERPIXEL, 3 );                  // set number of channels per pixel
  TIFFSetField( out, TIFFTAG_BITSPERSAMPLE, bits_per_sample );      // set the size of the channels
  TIFFSetField( out, TIFFTAG_SAMPLEFORMAT, sample_format );         // Floating point precision
  TIFFSetField( out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );    // set the origin of the image.
  TIFFSetField( out, TIFFTAG_RESOLUTIONUNIT, RESUNIT_CENTIMETER );  // set resolution to cm
  TIFFSetField( out, TIFFTAG_XRESOLUTION, 150.0 );                  // 150 pixels per cm
  TIFFSetField( out, TIFFTAG_YRESOLUTION, 150.0 );                  // 150 pixels per cm
  //  TIFFSetField( out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out,header.scanlines*header.samples) );
  TIFFSetField( out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize( out, header.samples*3 ) );
  TIFFSetField( out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG );
  TIFFSetField( out, TIFFTAG_PHOTOMETRIC, colorspace );
  TIFFSetField( out, TIFFTAG_COMPRESSION, compression );
  TIFFSetField( out, TIFFTAG_SOFTWARE, "hyper2color" );
  TIFFSetField( out, TIFFTAG_IMAGEDESCRIPTION, "Color rendering of hyperspectral image cube" );

  // TIFFTAG_COLORSPACE is an EXIF field - leave for now
  //  if( colorspace == PHOTOMETRIC_RGB ) TIFFSetField( out, TIFFTAG_COLORSPACE, 1 );
  //  else TIFFSetField( out, TIFFTAG_COLORSPACE, 65535 );

  // Set ICC profile if necessary
  if( icc_profile > 0 ){
    unsigned char* buffer = NULL;
    size_t len;
    if( icc_profile == 1 ){
      buffer = sRGB_ICC;
      len = sRGB_ICC_size;
    }
    else if( icc_profile == 2 ){
      buffer = AdobeRGB_ICC;
      len = AdobeRGB_ICC_size;
    }
    TIFFSetField( out, TIFFTAG_ICCPROFILE, len, buffer );
  }


  // We set the strip size of the file to be size of one row of pixels
  //  TIFFSetField( out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize( out, header.samples*3 ) );


  /* Set up our integration function
   */
  gsl_function F;
  F.function = &f;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc( 2000 );


  /* Loop through our pixels and calculate the CIE XYZ
   */
  for( j=0; j<header.scanlines; j++ ){

    /* Load an entire line in BIL (Band Interleaved Line) format
     */
    load_hyspex_bil( in, &header, scanline_spectrum, j );

    for( i=0; i<header.samples; i++ ){

      /* Extract the spectral values for pixel i, j
       */
      for( k=0; k<header.bands; k++ ){
	n = i + header.samples*k;
	if( header.bpp == 2 ) spectrum[k] = (double)scanline_spectrum[n] / 65535.0;
	else spectrum[k] = (double)scanline_spectrum[n];
      }

      /* Interpolate using GSL
       */
      gsl_interp_accel *acc = gsl_interp_accel_alloc();
      gsl_spline *spline = gsl_spline_alloc( gsl_interp_linear, header.bands );
      gsl_spline_init( spline, header.wavelengths, spectrum, header.bands );


      double X, Y, Z;
      X = Y = Z = 0.0;


      //#pragma omp parallel for
      for( k=0; k<=830-firstwav; k++ ){

	double val = gsl_spline_eval(spline, k+firstwav, acc);
	if( val < 0 ) val = 0.0;

	X += val * cie_color_match[tr+k][1] * power_spectrum[te+k][1];
	Y += val * cie_color_match[tr+k][2] * power_spectrum[te+k][1];
	Z += val * cie_color_match[tr+k][3] * power_spectrum[te+k][1];
      }


      /* Free our interpolator
       */
      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);


      float XX = (float)( (X * 100.0) / (norm) );
      float YY = (float)( (Y * 100.0) / (norm) );
      float ZZ = (float)( (Z * 100.0) / (norm) );


      /* Convert to RGB (sRGB or AdobeRGB)
       */
      if( colorspace == PHOTOMETRIC_RGB ){

	float R, G, B;

	// Use appropriate color conversion matrix depending on color space and color temperature
	if( icc_profile == 2 ){
	  if( temperature == 5000 ) XYZ2RGB(XYZ_AdobeRGB_matrix_D50,XX,YY,ZZ,&R,&G,&B);
	  else                      XYZ2RGB(XYZ_AdobeRGB_matrix_D65,XX,YY,ZZ,&R,&G,&B);

	  AdobeRGB_Gamma( &R, &G, &B );
	}
	else{
	  if( temperature == 5000 ) XYZ2RGB(XYZ_sRGB_matrix_D50,XX,YY,ZZ,&R,&G,&B);
	  else                      XYZ2RGB(XYZ_sRGB_matrix_D65,XX,YY,ZZ,&R,&G,&B);

	  sRGB_Gamma( &R, &G, &B );
	}

	if( bits_per_sample == 32 ){
	  ((float*)calculated_color)[i*3]     = R;
	  ((float*)calculated_color)[i*3 + 1] = G;
	  ((float*)calculated_color)[i*3 + 2] = B;
	}
	else if( bits_per_sample == 16 ){
	  ((unsigned short*)calculated_color)[i*3]     = (unsigned short)( R * 65535.0 );
	  ((unsigned short*)calculated_color)[i*3 + 1] = (unsigned short)( G * 65535.0 );
	  ((unsigned short*)calculated_color)[i*3 + 2] = (unsigned short)( B * 65535.0 );
	}
	else{
	  ((unsigned char*)calculated_color)[i*3]     = (unsigned char)( R * 255.0 );
	  ((unsigned char*)calculated_color)[i*3 + 1] = (unsigned char)( G * 255.0 );
	  ((unsigned char*)calculated_color)[i*3 + 2] = (unsigned char)( B * 255.0 );
	}
      }

      // CIE L*a*b* color space
      else{
	float L, a, b;
	XYZ2LAB(XX,YY,ZZ,&L,&a,&b);

	if( bits_per_sample == 8 ){
	  ((unsigned char*)calculated_color)[i*3]     = (unsigned char)( L * 2.55 );
	  ((unsigned char*)calculated_color)[i*3 + 1] = (signed char) (a);
	  ((unsigned char*)calculated_color)[i*3 + 2] = (signed char) (b);
	}
	else if( bits_per_sample == 16 ){
	  ((unsigned short*)calculated_color)[i*3]     = (unsigned short)( L * 655.35 );
	  ((unsigned short*)calculated_color)[i*3 + 1] = (signed short)( a * 255.0 );
	  ((unsigned short*)calculated_color)[i*3 + 2] = (signed short)( b * 255.0 );
	}
	else{
          ((float*)calculated_color)[i*3] = L;
          ((float*)calculated_color)[i*3 + 1] = a;
	  ((float*)calculated_color)[i*3 + 2] = b;	  
	}
      }
    }


    /* Write out a whole scanline
     */
    if( TIFFWriteScanline(out, (void*) calculated_color, j, 0) == -1 ){
      printf( "TIFF write error at scanline %d \n", j );
      break;
    }


    /* Report progress
     */
    if( verbose ){
      printf( "Processing: %3d\%\r", (int)(j*100.0/header.scanlines) );
      fflush( stdout );
    }

  }

  /* Free our line of color output values
   */
  free( calculated_color );
  free( scanline_spectrum );


  /* Free our integration workspace
   */
  gsl_integration_workspace_free( w );


  /* Close our files
   */
  if( in ) fclose( in );
  if( out ) TIFFClose( out );

  return 0;
}
