#!/usr/bin/perl
#
# SYNOPSIS:
#   make_html.pl [OPTS] 
#
# DESCRIPTION:
#   Create an html page that provides a link to the latest generated
#   BioNetGen distribution.
#
# OPTIONS:
#   --platform   PLAT  : Choices are Linux, MacOSX, Windows


use strict;
use warnings;

# Perl Core Modules
use FindBin;
use File::Spec;
use Getopt::Long;
use Cwd ("getcwd");
use Config;
use File::Path qw(remove_tree);


# distribution version (default undefined)
my $version = '';
# platform choices are: (MacOSX or Linux or Windows)
my $platform = '';
# path to version file
my $path_to_version_file = '.';


GetOptions( 'help|h'        => sub { display_help(); exit(0); },
            'platform=s'    => \$platform);


print "Platform: $platform\n";


my $zip_type  = '';
my $travis_os = '';
if ($platform eq "linux") {
  $zip_type = ".tar.gz";  $travis_os = "Linux";
  my $archive_file = "./dist/NFsim-source-".$platform.$zip_type;

  print " pwd = \n";
  system("pwd");
  
  print "\nCreating distribution archive:\n";
#  print "tar cvzf ${archive_file} ${dist_dir}\n";
  system("tar cvzf ${archive_file}  bin doc models src test tools validate CMakeLists.txt LICENSE.txt README.txt makefile NFsim-manual_v1.12.pdf ");
  system("ls -lt dist");

} else {
  if ($platform eq "osx") {
  $zip_type = ".tar.gz";  $travis_os = "MacOSX";
  } else {
    if ($platform eq "windows") {
      $zip_type = ".zip";  $travis_os = "Windows";
    } else {
      print "Invalid platform: ".$platform."\n";
      exit;
    }
  }
}


exit;


# ########################################################################
#   HELP 
# ########################################################################
# display help menu
sub display_help
{
print <<END_HELP
make_html.pl
SYNOPSIS:
   make_html.pl [OPTS] 
DESCRIPTION:
   Create an html page that provides a link to the latest generated
   BioNetGen distribution.
OPTIONS:
   --platform   PLAT  : Choices are Linux, MacOSX, Windows
END_HELP

}
