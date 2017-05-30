#!/usr/bin/perl
#
# SYNOPSIS:
#   make_html.pl [OPTS] 
#
# DESCRIPTION:
#   Create an html page that provides a link to the latest generated
#   NFsim distribution.
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

    print "platform: $platform\n";
    &write_html($platform,$version);

exit;

# ########################################################################
#      WRITE HTML
# ########################################################################
sub write_html {
  my $travis_os   = $_[0];


  my $zip_type  = '';
  if ($travis_os eq "linux") {
    $zip_type = ".tar.gz";  $platform = "Linux";
  } else {
    if ($travis_os eq "osx") {
    $zip_type = ".tar.gz";  $platform = "MacOSX";
    } else {
      if ($travis_os eq "windows") {
        $zip_type = ".zip";  $platform = "Windows";
      } else {
        print "Invalid platform: ".$travis_os."\n";
        exit;
      }
    }
  }
  
  my $ofile_name = "./dist/NFsim-".$travis_os.".html";


my $date_stamp; my $sec; my $min; my $hour; my $mday; my $mon; my $year; my $wday; my $yday; my $isdst;
open(FNEW,">$ofile_name");
  print FNEW "<html>\n";
  print FNEW "<head>\n";
  print FNEW "<title>".$platform." Beta Site for NFsim</title>\n";
  print FNEW "<META NAME=\"ROBOTS\" CONTENT=\"NOINDEX,NOFOLLOW\">\n";


  print FNEW "<style>\n";
  print FNEW "body {\n";
  print FNEW "   margin-top: 100px;\n";
  print FNEW "   margin-bottom: 100px;\n";
  print FNEW "   margin-right: 150px;\n";
  print FNEW "   margin-left: 80px;\n";
  print FNEW "}\n";
  print FNEW "</style>\n";


  print FNEW "</head>\n";
  print FNEW "<body bgcolor=\"LightSkyBlue\">\n";

  print FNEW "<font size=\"4\">\n"; 
  print FNEW "<center><h1>".$platform." Beta Site for NFsim</h1></center>\n";
  print FNEW "<br>\n";
  print FNEW "<font color=\"red\">WARNING:</font> This is not the download \n";
  print FNEW "site for BioNetGen.  If you wish to download the latest \n";
  print FNEW "version of BioNetGen, please visit: \n";
  print FNEW "<center><h1><a href=\"http://bionetgen.org\">bionetgen.org</a></h1></center>\n";
  print FNEW "<br>\n";
  print FNEW "<br>\n";
  print FNEW "<br>\n";

  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
  $date_stamp          = sprintf "%4d-%02d-%02d %02d:%02d\n",$year+1900,$mon+1,$mday,$hour-4,$min;

  print FNEW "If you are an NFsim developer, and you wish to get \n";
  print FNEW "access to the latest build for ".$platform.", please click here: <br>\n";
  print FNEW "<center><h1>\n";
  print FNEW "<a href=\"NFsim-".$travis_os."\">\n";
  print FNEW           "NFsim-".$travis_os."</a>\n";
  print FNEW "</h1></center>\n";
  print FNEW "<center>\n";
  print FNEW "<a href=\"NFsim-".$travis_os."\">\n";
  print FNEW "Time Stamp: ".$date_stamp."</a>\n";
  print FNEW "</center>\n";
  print FNEW "<br>\n";
  print FNEW "<center> Dont forget that after doing the download, you will need to do: \n";
  print FNEW "<font color=\"red\"> chmod +x NFsim-".$travis_os."</center>";
  print FNEW "<br>\n";
  print FNEW "<br>\n";

  print FNEW "</font></body></html>\n";
  
  close(FNEW);

  return;
}


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
   NFsim distribution.
OPTIONS:
   --platform   PLAT  : Choices are Linux, MacOSX, Windows
END_HELP

}
