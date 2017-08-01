

curl -v --ftp-ssl  -u $FTP_USER:$FTP_PASSWORD -T ./dist/NFsim-$TRAVIS_OS_NAME  ftp://ftp.box.com/BioNetGen_Beta/d_data/d_travis/  



#    The FTP_USER and FTP_PASSWORD strings are actually defined in the env: global: -secure: strings in travis.yml.  Those strings
#  Were generated using a program that may be installed on Linux.  Documentation can be found at: https://docs.travis-ci.com/user/encryption-keys/
#  To make a long story short, the program may be installed on Linux with:      gem install travis
#  the command to generate the FTP_PASSWORD string is
#
#              travis encrypt  FTP_USER="password" --skip-version-check -r RuleWorld/nfsim
#
#  The password should be surroundedd by double quotes, and the -r parameter indicates the repository for which the password
#  will be used.
#


