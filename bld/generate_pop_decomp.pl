#!/usr/bin/env perl
#=======================================================================
#
#  This is a script to return the decomposition information for
#  either CICE or POP.
#
# Usage:
#
# generate_cice_pop_decomp [options]
#
# To get help on options and usage:
#
# generate_cice_pop_decomp -help
#
#=======================================================================

use Cwd;
use strict;
#use diagnostics;
use Getopt::Long;
use English;

#-----------------------------------------------------------------------------------------------

#Figure out where configure directory is and where can use the XML/Lite module from
my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!; # name of program
my $ProgDir = $1;                         # name of directory where program lives

my $cwd = getcwd();  # current working directory
my $cfgdir;
my $utilroot = $ENV{'UTILROOT'};

if ($ProgDir) { $cfgdir = $ProgDir; }
else { $cfgdir = $cwd; }

#-----------------------------------------------------------------------------------------------
# Add $cfgdir to the list of paths that Perl searches for modules
my @dirs = ( $cfgdir, "$cfgdir/perl5lib", "$cfgdir/../../../../scripts/ccsm_utils/Tools/perl5lib", "$utilroot/Tools/perl5lib" );
unshift @INC, @dirs;
my $result = eval "require Decomp::Config";
if ( ! defined($result) ) {
   die <<"EOF";
** Cannot find perl module \"Decomp::Config\" from directories: @dirs **
EOF
}
require Decomp::Config;

my $model    = "pop";
my $platform = "XT";
my $res      = "gx1v6";
my $output   = "all";

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]
OPTIONS
     -res <resolution>    (or -r)   Horizontal resolution (gx1v6 etc.). (default $res))
     -model <model>       (or -m)   Model type (pop or cice).	        (default $model)
     -platform <platform> (or -p)   Platform type (XT etc.).		(defualt $platform)
     -nproc <number>      (or -n)   Number of processors to use.	(required)
     -thrds <number>      (or -t)   Number of threads per processor	(default 1)
     -output <type>	  (or -o)   Either output: all, maxblocks, bsize-x, bsize-y, or decomptype
			  (default: $output)
EXAMPLES

   $ProgName -res gx1v6 -model cice -platform XT -nproc 80 -output maxblocks

   will return a single value -- the optimum max number of blocks to use.

EOF
}

#------------------------------------------------------------------------------------------------

  my %opts = (
                res      => $res,
                model    => $model,
                platform => $platform,
                nproc    => undef,
                thrds    => 1,
                output   => $output,
                printing => 1,
                help     => 0,
                file     => "$cfgdir/pop_decomp.xml",
           );

  my $cmdline = @ARGV;
  GetOptions( 
              "r|res=s"      => \$opts{'res'},
              "m|model=s"    => \$opts{'model'},
              "p|platform=s" => \$opts{'platform'},
              "n|nproc=i"    => \$opts{'nproc'},
              "t|thrds=i"    => \$opts{'thrds'},
              "o|output=s"   => \$opts{'output'},
              "h|elp"        => \$opts{'help'},
          ) or usage();

  # Check for unparsed arguments
  if (@ARGV) {
      print "ERROR: unrecognized arguments: @ARGV\n";
      usage();
  }
  if ( $opts{'help'} ) {
      usage();
  }

  foreach my $key ( keys( %opts ) ) {
     if ( $key ne "help" && ! defined($opts{$key}) ) {
        print "ERROR: required input $key was not set\n";
        usage();
     }
  }

  $opts{'ProgName'} = $ProgName;
  $opts{'ProgDir'}  = $cfgdir;
  $opts{'cmdline'}  = $cmdline;

# redefine nproc to be total procs, nproc*thrds
  $opts{'nproc'} = $opts{'nproc'} * $opts{'thrds'};

# try to read from the xml file
  my $dcmp = Decomp::Config->new( \%opts );
  my %decomp = ( maxblocks=>0, bsize_x=>0, bsize_y=>0, decomptype=>"",
                 nlats=>0, nlons=>0 );
  my $matches = $dcmp->ReadXML( $opts{'file'}, \%decomp );

# if no xml entry, try to generate something
  if ( $decomp{'maxblocks'} == 0) {
     %decomp = CalcDecompInfo( $decomp{'nlats'}, $decomp{'nlons'}, \%opts);
  }
 
# adjust maxblocks to take into account threading
  $decomp{'maxblocks'} = $decomp{'maxblocks'} * $opts{'thrds'};

  if ( $decomp{'maxblocks'} == 0 ) {
     printf "%d %s",-1, "ERROR:($ProgName) No Decomp Created \n";
  } else {
     if (      $opts{'output'} eq "all"       ) {
       printf "%d %d %d %d %d %s", $decomp{'nlons'}, $decomp{'nlats'}, 
          $decomp{'bsize_x'}, $decomp{'bsize_y'}, $decomp{'maxblocks'}, $decomp{'decomptype'};
      } elsif ( $opts{'output'} eq "maxblocks" ) {
        print $decomp{'maxblocks'};
      } elsif ( $opts{'output'} eq "bsize_x"   ) {
        print $decomp{'bsize_x'};
      } elsif ( $opts{'output'} eq "bsize_y"   ) {
        print $decomp{'bsize_y'};
      } elsif ( $opts{'output'} eq "decomptype") {
        print $decomp{'decomptype'};
      } else {
        print "ERROR:($ProgName) bad argument to output option $opts{'output'}\n";
        usage();
      }
      print "\n";
  }


sub CalcDecompInfo {
#
# Calculate decomposition information
# Tries to first find an even cartesian decomposition (set = 1)
# If can't find an even decomp, tries a space filling curve (set = 2)
#
  my $nlats    = shift;
  my $nlons    = shift;
  my $opts_ref = shift;

  my %opts   = %$opts_ref;
  my $nprocs = $opts{'nproc'};
  my $model  = $opts{'model'};

  my ($maxblocks,$bsize_x,$bsize_y,$decomptype);
  my %decomp;
  my $found = 0;
  my $done = 0;
  my $nprocsx = 0;
  my $nprocsy = 0;
  my $nx = 0;
  my $ny = 0;
  my $nn = 0;
  my $tmp = 0;
  my $nscore = 0.0 ;
  my $bscore = $nlons * $nlats * $nprocs ;

# find an even 2d decomposition that has the most square blocks
# for a given total number of processors, nproc, such that
#   nx can be 1 to nproc
#   nx*ny = nproc
#   mod(nlats,ny) = mod(nlons,nx) = 0
#   nscore is the minimum value where
#     nscore = tmp * tmp where
#       tmp is (bsize_x/bsize_y - 1)
#   we want bsize_x/bsize_y to be closest to 1 so subtract
#   1 and "square it" to create a function that maximizes
#   squareness when it's minimum.
# found indicates a decomp has been found, but need to continue
# to search to find the best decomp.

  $nn = 0;
  do {
     $nn = $nn + 1;
     $ny = $nn;
     $nx = int($nprocs/$ny);
     if ($ny * $nx == $nprocs &&
         $nlats % $ny == 0 &&
         $nlons % $nx == 0) {

        $tmp = ($nlons/$nx * $ny/$nlats) - 1.0 ;
        $nscore = $tmp * $tmp;
        if ($nscore < $bscore) {
          $bscore = $nscore;
          $nprocsx = $nx;
          $nprocsy = $ny;
          $found = 1;
       }
     }

#     print "debug $nn $nx $ny $nprocsx $nprocsy $nscore $bscore $found \n";

  } until ($nn == $nprocs);

# space filling curves
  if ($found == 0) {
#    what do we do here?     
#    $found = 2
  }


  if ($found == 1) {
    $decomp{'nlats'}      = $nlats;
    $decomp{'nlons'}      = $nlons;
    $decomp{'maxblocks'}  = 1;
    $decomp{'decomptype'} = "cartesian";
    $decomp{'bsize_x'}    = int( $nlons / $nprocsx );
    $decomp{'bsize_y'}    = int( $nlats / $nprocsy );
  }

  return(%decomp);
}
