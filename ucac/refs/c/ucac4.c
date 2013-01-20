#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ucac4.h"

/* Basic access functions for UCAC-4.  Public domain.  Please contact
pluto (at) projectpluto.com with comments/bug fixes.  */

/* History: */

/*   2012 Mar 18:  (BJG) Created using the UCAC3 code as a template. */

/* UCAC4 consists of 900 zones,  each .2 degrees high in declination. */

#ifdef CUMULATIVE_DESIGNATORS_USED
      /* ...unless we end up with a designation system that is a */
      /* cumulative count across the entire catalog,  instead of */
      /* a zone/number-within-zone designation scheme.  I don't  */
      /* know which will be used yet,  so I've commented this    */
      /* out but am not deleting it yet:                         */

static const int32_t ucac4_offsets[901] = {
        0,       206,       866,      2009,      3622,      5571,
     8119,     11083,     14581,     18489,     22819,     27668,
    32894,     38573,     44607,     51083,     57959,     65122,
    73128,     81479,     90380,     99902,    109805,    120449,
   131069,    141935,    153543,    165254,    177944,    190397,
   203265,    217164,    230956,    245739,    260436,    275441,
   291131,    306759,    323263,    339786,    356856,    374663,
   392636,    411953,    430914,    450375,    470619,    490371,
   511931,    533948,    556980,    580928,    605167,    630973,
   656223,    682372,    710516,    738355,    768228,    797955,
   828196,    859161,    889896,    922517,    954434,    987871,
  1022311,   1057120,   1094390,   1131881,   1171406,   1212583,
  1254514,   1299068,   1343566,   1389869,   1437370,   1485272,
  1537752,   1590554,   1645964,   1703267,   1761235,   1824574,
  1888291,   1954427,   2020545,   2087100,   2158181,   2229282,
  2305819,   2382310,   2458399,   2539228,   2619455,   2703320,
  2791647,   2880146,   2978698,   3078455,   3183901,   3291044,
  3399522,   3517174,   3631529,   3749099,   3869671,   3989408,
  4116811,   4242528,   4372612,   4501990,   4630724,   4766601,
  4900430,   5043092,   5190015,   5334648,   5490087,   5644280,
  5806690,   5977882,   6149759,   6333387,   6511715,   6693107,
  6879918,   7058706,   7239823,   7414877,   7596441,   7776823,
  7955439,   8142451,   8321547,   8501071,   8682641,   8862674,
  9057390,   9249537,   9452538,   9662243,   9869267,  10088246,
 10301374,  10528856,  10764546,  10997480,  11233989,  11455728,
 11681153,  11909327,  12124701,  12350926,  12572014,  12795436,
 13012238,  13222762,  13443183,  13658499,  13880831,  14102818,
 14314906,  14532699,  14741997,  14956716,  15172778,  15380595,
 15598750,  15806866,  16018075,  16229926,  16437454,  16655825,
 16866560,  17086363,  17310404,  17533912,  17766260,  17989304,
 18219078,  18440914,  18659549,  18885704,  19102986,  19323150,
 19540437,  19748650,  19961298,  20165181,  20373302,  20582550,
 20783195,  20989743,  21186775,  21385982,  21584871,  21776272,
 21972232,  22157647,  22350854,  22546489,  22736746,  22931113,
 23116091,  23304438,  23492534,  23677999,  23871500,  24055348,
 24245864,  24438253,  24627079,  24831459,  25023554,  25219448,
 25417851,  25610344,  25808561,  25996908,  26189924,  26379076,
 26563662,  26757932,  26947542,  27140183,  27341349,  27539669,
 27742006,  27934620,  28127345,  28314380,  28494206,  28683526,
 28867915,  29057414,  29242460,  29419129,  29601907,  29778955,
 29963139,  30151977,  30336478,  30530769,  30718133,  30909666,
 31100339,  31285931,  31479450,  31664777,  31856760,  32053192,
 32245001,  32442855,  32630361,  32825001,  33021120,  33214516,
 33417171,  33611502,  33812537,  34012212,  34208558,  34410672,
 34601527,  34798895,  34995096,  35187473,  35388276,  35580630,
 35780449,  35979256,  36177626,  36387369,  36591643,  36807134,
 37025815,  37243104,  37475827,  37698299,  37925218,  38150060,
 38370268,  38597638,  38811749,  39035564,  39259069,  39475312,
 39702650,  39920680,  40146857,  40375027,  40598240,  40834389,
 41064908,  41306434,  41544953,  41779844,  42027725,  42276859,
 42544638,  42820639,  43097301,  43376972,  43649874,  43936633,
 44223548,  44498981,  44777602,  45044307,  45311896,  45577796,
 45837837,  46102751,  46354839,  46604248,  46850161,  47089350,
 47334866,  47567152,  47800465,  48034689,  48258171,  48484058,
 48703272,  48924962,  49134462,  49336043,  49544964,  49749124,
 49958415,  50172955,  50387447,  50613268,  50833890,  51055687,
 51264896,  51470509,  51681290,  51882668,  52088557,  52295264,
 52496319,  52705052,  52907451,  53112636,  53311236,  53504657,
 53702628,  53895625,  54092142,  54288125,  54478742,  54674913,
 54869875,  55073788,  55277754,  55479681,  55690209,  55897312,
 56112255,  56319448,  56522797,  56724822,  56911977,  57104100,
 57292101,  57475865,  57667867,  57855066,  58045405,  58232617,
 58414752,  58603745,  58784456,  58968537,  59148928,  59323159,
 59499108,  59666285,  59835396,  60006720,  60171236,  60337699,
 60496827,  60658849,  60812631,  60960518,  61109621,  61249081,
 61388880,  61524571,  61655527,  61791729,  61921423,  62053876,
 62190195,  62323601,  62462614,  62597982,  62737390,  62870650,
 63002358,  63142211,  63276593,  63413539,  63555908,  63695591,
 63841731,  63982652,  64127964,  64265998,  64402435,  64551415,
 64695997,  64843498,  64996236,  65146242,  65301430,  65450368,
 65600245,  65743693,  65882206,  66024395,  66159006,  66292509,
 66427050,  66553501,  66680021,  66799980,  66922092,  67042185,
 67160865,  67283938,  67403459,  67525135,  67645929,  67764731,
 67889823,  68014244,  68140479,  68269657,  68399630,  68534896,
 68667463,  68804968,  68946515,  69084129,  69225206,  69363053,
 69501421,  69634836,  69759606,  69886083,  70006664,  70129161,
 70253717,  70376410,  70505720,  70632203,  70763001,  70895156,
 71023425,  71158961,  71291445,  71428721,  71566984,  71702674,
 71844316,  71980829,  72118027,  72251613,  72381864,  72516495,
 72645490,  72777581,  72912081,  73043773,  73182303,  73315870,
 73453777,  73590973,  73725822,  73866983,  74002648,  74140226,
 74280406,  74419211,  74564959,  74707161,  74852948,  75005061,
 75153213,  75305697,  75451463,  75598671,  75747614,  75892780,
 76043759,  76190156,  76339798,  76488674,  76634382,  76784765,
 76931496,  77080036,  77228311,  77375178,  77527606,  77675768,
 77825495,  77970062,  78109440,  78254528,  78396584,  78540380,
 78688178,  78830813,  78978369,  79122902,  79269264,  79416466,
 79561951,  79714331,  79863992,  80014458,  80167197,  80317763,
 80472563,  80622769,  80774939,  80924495,  81070768,  81220355,
 81365073,  81513125,  81664670,  81813156,  81966524,  82117939,
 82271689,  82422213,  82570771,  82722687,  82870182,  83019300,
 83168091,  83312922,  83462845,  83610177,  83759610,  83903733,
 84043779,  84184355,  84316706,  84453813,  84595155,  84733909,
 84879114,  85022743,  85166973,  85310552,  85450975,  85591772,
 85726443,  85865377,  85997980,  86128168,  86263466,  86394862,
 86526209,  86659574,  86789742,  86926075,  87059553,  87196375,
 87330165,  87462684,  87599249,  87726329,  87855697,  87989149,
 88120449,  88255327,  88387574,  88521530,  88654434,  88785614,
 88922537,  89057337,  89197471,  89337989,  89477554,  89622144,
 89760316,  89900524,  90043815,  90184212,  90326902,  90466149,
 90609411,  90752343,  90892033,  91034943,  91172644,  91311310,
 91449841,  91585099,  91726461,  91861689,  91998957,  92132668,
 92263230,  92400963,  92533468,  92669353,  92807937,  92945608,
 93088072,  93222770,  93359805,  93498210,  93631255,  93766785,
 93896843,  94028739,  94160038,  94290246,  94424231,  94554567,
 94688997,  94826105,  94960640,  95096976,  95226447,  95357687,
 95486349,  95612360,  95738734,  95859895,  95984034,  96112786,
 96238538,  96365556,  96490151,  96616133,  96737014,  96854193,
 96976089,  97093620,  97211356,  97331700,  97450285,  97573204,
 97692810,  97811776,  97928630,  98044212,  98162393,  98277717,
 98394471,  98513752,  98630482,  98751581,  98871789,  98994428,
 99117358,  99236264,  99359070,  99478692,  99600615,  99723109,
 99843973,  99970116, 100095278, 100225651, 100355481, 100485073,
100620812, 100753408, 100887796, 101022844, 101157240, 101295313,
101427331, 101560501, 101694327, 101825340, 101958059, 102088564,
102220275, 102350267, 102478211, 102610306, 102739290, 102869786,
102999872, 103127790, 103258649, 103385394, 103514541, 103641316,
103765125, 103890749, 104011157, 104133400, 104253451, 104372695,
104495855, 104616365, 104739017, 104860983, 104979547, 105101052,
105217572, 105334669, 105452060, 105568395, 105687016, 105801577,
105916899, 106032335, 106146860, 106264074, 106376355, 106488910,
106599663, 106707327, 106817932, 106923848, 107030933, 107135924,
107239646, 107348770, 107454737, 107561296, 107669088, 107774161,
107878965, 107977751, 108076372, 108174958, 108270708, 108364656,
108452733, 108541106, 108630830, 108717179, 108806291, 108892147,
108978812, 109064547, 109149065, 109236525, 109320479, 109406010,
109489445, 109570239, 109651768, 109729456, 109807640, 109884948,
109960688, 110037996, 110114562, 110192074, 110268905, 110343407,
110416551, 110484908, 110552200, 110617712, 110681497, 110744680,
110803097, 110860027, 110915955, 110968945, 111021532, 111070721,
111119000, 111165227, 111209791, 111254569, 111297610, 111340866,
111383787, 111425666, 111469402, 111512046, 111554826, 111600005,
111644195, 111688980, 111731472, 111773517, 111814410, 111854126,
111895147, 111934334, 111974133, 112014255, 112053072, 112093290,
112131676, 112170831, 112208531, 112245060, 112282544, 112318561,
112354797, 112390557, 112425011, 112459768, 112492430, 112524745,
112556391, 112586924, 112618156, 112649072, 112680236, 112710897,
112740453, 112770354, 112797835, 112825364, 112851786, 112877701,
112904019, 112929142, 112954530, 112979406, 113003202, 113027028,
113049299, 113071829, 113093972, 113115571, 113138157, 113159919,
113181445, 113203017, 113224071, 113244866, 113264853, 113284879,
113304606, 113323850, 113343510, 113361896, 113380253, 113398208,
113415740, 113432675, 113448765, 113464551, 113479492, 113494194,
113508899, 113522894, 113536655, 113549834, 113562609, 113575291,
113587163, 113599160, 113610752, 113621762, 113632588, 113642685,
113652637, 113662010, 113671156, 113680187, 113688547, 113697000,
113704961, 113712710, 113720022, 113726470, 113732697, 113738459,
113743893, 113749088, 113753913, 113758382, 113762417, 113766129,
113769622, 113772648, 113775600, 113778330, 113780735, 113782892,
113784770, 113786420, 113787918, 113789245, 113790079, 113790591,
113790762 };
#endif     // CUMULATIVE_DESIGNATORS_USED

#ifdef __BYTE_ORDER
#if __BYTE_ORDER == __BIG_ENDIAN
static void swap_32( int32_t *ival)
{
   int8_t temp, *zval = (int8_t *)ival;

   temp = zval[0];
   zval[0] = zval[3];
   zval[3] = temp;
   temp = zval[1];
   zval[1] = zval[2];
   zval[2] = temp;
}

static void swap_16( int16_t *ival)
{
   int8_t temp, *zval = (int8_t *)ival;

   temp = zval[0];
   zval[0] = zval[1];
   zval[1] = temp;
}

void flip_ucac4_star( UCAC4_STAR *star)
{
   int i;

   swap_32( &star->ra);
   swap_32( &star->spd);
   swap_16( &star->mag1);
   swap_16( &star->mag2);
   swap_16( &star->epoch_ra);
   swap_16( &star->epoch_dec);
   swap_32( &star->pm_ra);
   swap_32( &star->pm_dec);
   swap_32( &star->twomass_id);
   swap_16( &star->mag_j);
   swap_16( &star->mag_h);
   swap_16( &star->mag_k);
   for( i = 0; i < 5; i++)
      swap_16( &star->apass_mag[i]);
   swap_32( &star->catalog_flags);
   swap_32( &star->id_number);
   swap_16( &star->ucac2_zone);
   swap_32( &star->ucac2_number);
}
#endif                   // #if __BYTE_ORDER == __BIG_ENDIAN
#endif                   // #ifdef __BYTE_ORDER

/* The following function writes out a UCAC4 star in the same ASCII */
/* format as the FORTRAN code.                                      */

int write_ucac4_star_fortran_style( char *obuff, const UCAC4_STAR *star)
{
/*
          WRITE (line,'(2i10,2i6,i3,i2,i3,2i4,3i3,2i6,2i7,2i4
     .                ,i11,3i6,6i3,5i6,5i4,i2,9i2,2i3,i10,i4,i7)')

 451313731,   544918,16176,16187,12,0, 0, 58, 76, 4, 4, 2, 9762, 9709,   -47,   121,139,140, 863323727,13612,12869,12680, 5, 5, 5, 3, 3, 3,20000,20000,20000,20000,20000,  0,  0,  0,  0,  0,4,0,0,0,0,0,0,0,0,1, 0, 0,  1000284,  0,     0

*/
   int i;
   char *tptr;

   sprintf( obuff, "%10d%10d%6d%6d%3d%2d%3d%4d%4d%3d%3d%3d%6d%6d%7d%7d%4d%4d",
        star->ra, star->spd, star->mag1, star->mag2,
        star->mag_sigma, star->obj_type, star->double_star_flag,
        star->ra_sigma + 128, star->dec_sigma + 128, star->n_ucac_total,
        star->n_ucac_used, star->n_cats_used,
        star->epoch_ra, star->epoch_dec,
        star->pm_ra, star->pm_dec,
        star->pm_ra_sigma + 128, star->pm_dec_sigma + 128);
   sprintf( obuff + strlen( obuff), "%11d%6d%6d%6d%3d%3d%3d%3d%3d%3d",
        star->twomass_id, star->mag_j, star->mag_h, star->mag_k,
        star->icq_flag[0], star->icq_flag[1], star->icq_flag[2],
        star->e2mpho[0], star->e2mpho[1], star->e2mpho[2]);
   for( i = 0; i < 5; i++)
      sprintf( obuff + strlen( obuff), "%6d", star->apass_mag[i]);
   for( i = 0; i < 5; i++)
      sprintf( obuff + strlen( obuff), "%4d", star->apass_mag_sigma[i]);
   sprintf( obuff + strlen( obuff), "%2d ", star->yale_gc_flags);
               /* Show catalog flags as separate digits: */
   tptr = obuff + strlen( obuff);
   sprintf( tptr, "%09d", star->catalog_flags);
   for( i = 8; i >= 0; i--)
      {
      tptr[i + i] = tptr[i];
      tptr[i + i + 1] = ' ';
      }
   tptr[17] = '\0';
   sprintf( obuff + strlen( obuff), "%3d%3d%10d%4d%7d\n",
       star->leda_flag, star->twomass_ext_flag,
       star->id_number, star->ucac2_zone, star->ucac2_number);
   return( 0);
}

/* This function,  like the above,  writes out the data for a UCAC-4 star */
/* into an ASCII buffer,  but it tries for a more human-readable format.  */

int write_ucac4_star( const int zone, const long offset, char *obuff,
                     const UCAC4_STAR *star, const int output_format)
{
   const long epoch_ra  = 190000 + star->epoch_ra;
   const long epoch_dec = 190000 + star->epoch_dec;
   int i;

   if( output_format & UCAC4_FORTRAN_STYLE)
      return( write_ucac4_star_fortran_style( obuff, star));

   sprintf( obuff, "%03d-%06ld %12.8lf %12.8lf %2d.%03d %2d.%03d %3d ",
               zone, offset,
               (double)star->ra / 3600000., (double)star->spd / 3600000. - 90.,
               star->mag1 / 1000, abs( star->mag1 % 1000),
               star->mag2 / 1000, abs( star->mag2 % 1000),
               star->mag_sigma);

   sprintf( obuff + strlen( obuff), "%2d %2d ",
               star->obj_type, star->double_star_flag);

   sprintf( obuff + strlen( obuff), "%4d.%02d %4d.%02d ",
               (int)epoch_ra / 100, (int)epoch_ra % 100,
               (int)epoch_dec / 100, (int)epoch_dec % 100);

   sprintf( obuff + strlen( obuff),
            "%3d %3d %2d %2d %2d ",
            star->ra_sigma + 128, star->dec_sigma + 128,
            (int)star->n_ucac_total, (int)star->n_ucac_used,
            (int)star->n_cats_used);

   if( star->pm_ra || star->pm_dec || !(output_format & UCAC4_WRITE_SPACES))
      sprintf( obuff + strlen( obuff), "%6d %6d %3d %3d ",
            star->pm_ra, star->pm_dec,
            star->pm_ra_sigma + 128, star->pm_dec_sigma + 128);
   else        /* no proper motion given,  keep these fields blank */
      strcat( obuff, "                      ");

   if( star->twomass_id || !(output_format & UCAC4_WRITE_SPACES))
      {
      sprintf( obuff + strlen( obuff),
            "%10ld %2d.%03d %2d.%03d %2d.%03d ",
            (long)star->twomass_id,
            star->mag_j / 1000, abs( star->mag_j % 1000),
            star->mag_h / 1000, abs( star->mag_h % 1000),
            star->mag_k / 1000, abs( star->mag_k % 1000));

      sprintf( obuff + strlen( obuff), "%03d %03d %03d %02d %02d %02d ",
            star->e2mpho[0], star->e2mpho[1], star->e2mpho[2],
            star->icq_flag[0], star->icq_flag[1], star->icq_flag[2]);
      }
   else        /* no 2MASS data given;  keep these fields blank */
      {
      memset( obuff + 116, ' ', 53);
      obuff[169] = '\0';
      }

   for( i = 0; i < 5; i++)
      if( (star->apass_mag[i] && star->apass_mag[i] != 20000)
                   || !(output_format & UCAC4_WRITE_SPACES))
         sprintf( obuff + strlen( obuff), "%2d.%03d ",
                star->apass_mag[i] / 1000, star->apass_mag[i] % 1000);
         else
            strcat( obuff, "       ");
   for( i = 0; i < 5; i++)
      if( star->apass_mag_sigma[i] || !(output_format & UCAC4_WRITE_SPACES))
         sprintf( obuff + strlen( obuff), "0.%03d ",
                    star->apass_mag_sigma[i]);
      else
         strcat( obuff, "      ");

   sprintf( obuff + strlen( obuff), "%09d", star->catalog_flags);
   sprintf( obuff + strlen( obuff),
            " %2d %03d %03d %9d", star->yale_gc_flags,
            star->leda_flag, star->twomass_ext_flag, star->id_number);
   if( star->ucac2_zone || !(output_format & UCAC4_WRITE_SPACES))
      sprintf( obuff + strlen( obuff), " %03d-%06d\n",
               star->ucac2_zone, star->ucac2_number);
   else
      strcat( obuff, "           \n");
   return( 0);
}

#if defined( linux) || defined( unix)
   static const char *path_separator = "/", *read_only_permits = "r";
#else
   static const char *path_separator = "\\", *read_only_permits = "rb";
#endif

/* The layout of UCAC-4 is such that files for the north (zones 380-900,
corresponding to declinations -14.2 to the north celestial pole) are in
the 'u4n' folder of one DVD.  Those for the south (zones 1-379,
declinations -14.2 and south) are in the 'u4s' folder of the other DVD.
People may copy these retaining the path structure,  or maybe they'll
put all 900 files in one folder.  So if you ask this function for,  say,
zone_number = 314 and files in the folder /data/ucac4,  the function will
look for the data under the following four names:

z314         (i.e.,  all data copied to the current folder)
u4s/z314     (i.e.,  you've copied everything to two subfolders of the current)
/data/ucac4/z314
/data/ucac4/u4s/z314

   ...stopping when it finds a file.  This will,  I hope,  cover all
likely situations.  If you make things any more complicated,  you've
only yourself to blame.   */

static FILE *get_ucac4_zone_file( const int zone_number, const char *path)
{
   FILE *ifile;
   char filename[80];

   sprintf( filename, "u4%c%sz%03d", (zone_number >= 380 ? 'n' : 's'),
                     path_separator, zone_number);
            /* First,  look for file in current path: */
   ifile = fopen( filename + 4, read_only_permits);
   if( !ifile)
      ifile = fopen( filename, read_only_permits);
         /* If file isn't there,  use the 'path' passed in as an argument: */
   if( !ifile && *path)
      {
      char filename2[80], *endptr;
      int i;

      strcpy( filename2, path);
      endptr = filename2 + strlen( filename2);
      if( endptr[-1] != *path_separator)
         *endptr++ = *path_separator;
      for( i = 0; !ifile && i < 2; i++)
         {
         strcpy( endptr, filename + 4 * (1 - i));
         ifile = fopen( filename2, read_only_permits);
         }
      }
   return( ifile);
}

int extract_ucac4_info( const int zone, const long offset, UCAC4_STAR *star,
                     const char *path)
{
   int rval;

   if( zone < 1 || zone > 900)     /* not a valid sequential number */
      rval = -1;
   else
      {
      FILE *ifile = get_ucac4_zone_file( zone, path);

      if( ifile)
         {
         if( fseek( ifile, (offset - 1) * sizeof( UCAC4_STAR), SEEK_SET))
            rval = -2;
         else if( !fread( star, sizeof( UCAC4_STAR), 1, ifile))
            rval = -3;
         else           /* success! */
            {
            rval = 0;
#ifdef __BYTE_ORDER
#if __BYTE_ORDER == __BIG_ENDIAN
            flip_ucac4_star( star);
#endif
#endif
            }
         fclose( ifile);
         }
      else
         rval = -4;
      }
   return( rval);
}

static FILE *get_ucac4_index_file( const char *path)
{
   FILE *index_file;
   const char *idx_filename = "u4index.asc";

                     /* Look for the index file in the local directory... */
   index_file = fopen( idx_filename, read_only_permits);
                     /* ...and if it's not there,  look for it in the same */
                     /* directory as the data: */
   if( !index_file)
      {
      char filename[100];

      strcpy( filename, path);
      if( filename[strlen( filename) - 1] != path_separator[0])
         strcat( filename, path_separator);
      strcat( filename, idx_filename);
      index_file = fopen( filename, read_only_permits);
      }
   return( index_file);
}

/* The layout of the ASCII index is a bit peculiar.  There are 1440
lines per dec zone (of which there are,  of course,  900). Each line
contains 21 bytes,  except for the first,  which includes the dec
and is therefore six bytes longer. */

static long get_index_file_offset( const int zone, const int ra_start)
{
   int rval = (zone - 1) * (1440 * 21 + 6) + ra_start * 21;

   if( ra_start)
      rval += 6;
   return( rval);
}

/* RA, dec, width, height are in degrees */

/* A note on indexing:  within each zone,  we want to locate the stars
within a particular range in RA.  If an index is unavailable,  then
we have things narrowed down to somewhere between the first and
last records.  If an index is available,  our search can take
place within a narrower range.  But in either case,  the range is
refined by doing a secant search which narrows down the starting
point to within 'acceptable_limit' records,  currently set to
40;  i.e.,  it's possible that we will read in forty records that
are before the low end of the desired RA range.  The secant search
is slightly modified to ensure that each iteration knocks off at
least 1/8 of the current range.

   Records are then read in 'buffsize' stars at a time and,  if
they're in the desired RA/dec rectangle,  written out to 'ofile'. */

#include <time.h>

clock_t time_searching = 0;

int extract_ucac4_stars( FILE *ofile, const double ra, const double dec,
                  const double width, const double height, const char *path,
                  const int output_format)
{
   const double dec1 = dec - height / 2., dec2 = dec + height / 2.;
   const double ra1 = ra - width / 2., ra2 = ra + width / 2.;
   const double zone_height = .2;    /* zones are .2 degrees each */
   int zone = (int)( (dec1  + 90.) / zone_height) + 1;
   const int end_zone = (int)( (dec2 + 90.) / zone_height) + 1;
   const int index_ra_resolution = 1440;  /* = .25 degrees */
   int ra_start = (int)( ra1 * (double)index_ra_resolution / 360.);
   int rval = 0;
   const int buffsize = 400;     /* read this many stars at a try */
   FILE *index_file = NULL;
   UCAC4_STAR *stars = (UCAC4_STAR *)calloc( buffsize, sizeof( UCAC4_STAR));

   if( !stars)
      rval = -1;
   if( zone < 1)
      zone = 1;
   if( ra_start < 0)
      ra_start = 0;
   while( rval >= 0 && zone <= end_zone)
      {
      FILE *ifile = get_ucac4_zone_file( zone, path);

      if( ifile)
         {
         int keep_going = 1;
         int i, n_read;
         const int32_t max_ra  = (int32_t)( ra2 * 3600. * 1000.);
         const int32_t min_ra  = (int32_t)( ra1 * 3600. * 1000.);
         const int32_t min_spd = (int32_t)( (dec1 + 90.) * 3600. * 1000.);
         const int32_t max_spd = (int32_t)( (dec2 + 90.) * 3600. * 1000.);
         uint32_t offset, end_offset;
         const uint32_t acceptable_limit = 40;
         long index_file_offset = get_index_file_offset( zone, ra_start);
         clock_t t0 = clock( );
         static long cached_index_data[5] = {-1L, 0L, 0L, 0L, 0L};
         const uint32_t ra_range = (uint32_t)( 360 * 3600 * 1000);
         uint32_t ra_lo = (uint32_t)( ra_start * (ra_range / index_ra_resolution));
         uint32_t ra_hi = ra_lo + ra_range / index_ra_resolution;

         if( index_file_offset == cached_index_data[0])
            {
            offset = cached_index_data[1];
            end_offset = cached_index_data[2];
            }
         else
            {
            if( !index_file)
               index_file = get_ucac4_index_file( path);
            if( index_file)
               {
               char ibuff[50];

               fseek( index_file, index_file_offset, SEEK_SET);
               fgets( ibuff, sizeof( ibuff), index_file);
               sscanf( ibuff, "%d%d", &offset, &end_offset);
               end_offset += offset;
               cached_index_data[0] = index_file_offset;
               cached_index_data[1] = offset;
               cached_index_data[2] = end_offset;
               }
            else     /* no index:  binary-search within entire zone: */
               {
               offset = 0;
               fseek( ifile, 0L, SEEK_END);
               end_offset = ftell( ifile) / sizeof( UCAC4_STAR);
//             end_offset = ucac4_offsets[zone] - ucac4_offsets[zone - 1];
               ra_lo = 0;
               ra_hi = ra_range;
               }
            }
//       printf( "Seeking RA=%u between offsets %u to %u\n", min_ra, offset, end_offset);
//                   /* Secant-search within the known limits: */
         while( end_offset - offset > acceptable_limit)
            {
            UCAC4_STAR star;
            uint32_t delta = end_offset - offset, toffset;
            uint32_t minimum_bite = delta / 8 + 1;
            uint64_t tval = (uint64_t)delta *
                        (uint64_t)( min_ra - ra_lo) / (uint64_t)( ra_hi - ra_lo);

            if( tval < minimum_bite)
               tval = minimum_bite;
            else if( tval > delta - minimum_bite)
               tval = delta - minimum_bite;
            toffset = offset + (uint32_t)tval;
            fseek( ifile, toffset * sizeof( UCAC4_STAR), SEEK_SET);
            fread( &star, 1, sizeof( UCAC4_STAR), ifile);
//          printf( "At offset %u: RA=%u; range %u to %u\n",
//                   toffset, star.ra, offset, end_offset);
            if( star.ra < min_ra)
               {
               offset = toffset;
               ra_lo = star.ra;
               }
            else
               {
               end_offset = toffset;
               ra_hi = star.ra;
               }
            }
         time_searching += clock( ) - t0;
         fseek( ifile, offset * sizeof( UCAC4_STAR), SEEK_SET);

         while( (n_read = fread( stars, sizeof( UCAC4_STAR), buffsize, ifile)) > 0
                                                   && keep_going)
            for( i = 0; i < n_read && keep_going; i++)
               {
               UCAC4_STAR star = stars[i];

#ifdef __BYTE_ORDER
#if __BYTE_ORDER == __BIG_ENDIAN
               flip_ucac4_star( &star);
#endif
#endif
               if( star.ra > max_ra)
                  keep_going = 0;
               else if( star.ra > min_ra && star.spd > min_spd
                                           && star.spd < max_spd)
//                if( !(output_format & UCAC4_OMIT_TYCHO_STARS) ||
//                         !star.catalog_flags[UCAC4_CATFLAG_TYCHO])
//                   if( star.twomass_id ||
//                         (output_format & UCAC4_INCLUDE_DOUBTFULS))
                        {
                        rval++;
                        if( ofile)
                           {
                           if( output_format & UCAC4_RAW_BINARY)
                              fwrite( &star, 1, sizeof( UCAC4_STAR), ofile);
                           else
                              {
                              char buff[UCAC4_ASCII_SIZE];

                              write_ucac4_star( zone, offset + 1, buff, &star,
                                                            output_format);
                              fwrite( buff, 1, strlen( buff), ofile);
                              }
                           }
                        }
               offset++;
               }
         fclose( ifile);
         }
      zone++;
      }
   if( index_file)
      fclose( index_file);
   free( stars);

            /* We need some special handling for cases where the area
               to be extracted crosses RA=0 or RA=24: */
   if( rval >= 0 && ra > 0. && ra < 360.)
      {
      if( ra1 < 0.)      /* left side crosses over RA=0h */
         rval += extract_ucac4_stars( ofile, ra+360., dec, width, height,
                                          path, output_format);
      if( ra2 > 360.)    /* right side crosses over RA=24h */
         rval += extract_ucac4_stars( ofile, ra-360., dec, width, height,
                                          path, output_format);
      }
   return( rval);
}
