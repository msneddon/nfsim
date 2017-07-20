if [ "${TRAVIS_OS_NAME}" = "linux" ]; then

curl -T ./dist/NFsim-$TRAVIS_OS_NAME               -u $FTP_USER_0:$FTP_PASSWORD ftp://ftp.drivehq.com/d_data/d_travis/  
curl -T ./dist/NFsim-source-$TRAVIS_OS_NAME.tar.gz -u $FTP_USER:$FTP_PASSWORD ftp://ftp.drivehq.com/d_data/d_travis/

else

curl -T ./dist/NFsim-$TRAVIS_OS_NAME               -u $FTP_USER_1:$FTP_PASSWORD ftp://ftp.drivehq.com/d_data/d_travis/  
curl -T ./dist/NFsim-source-$TRAVIS_OS_NAME.tar.gz -u $FTP_USER:$FTP_PASSWORD ftp://ftp.drivehq.com/d_data/d_travis/

fi




#    The FTP_USER and FTP_PASSWORD strings are actually defined in the env: global: -secure: strings above.  Those strings
#  Were generated using a program that may be installed on Linux.  Documentation can be found at: https://docs.travis-ci.com/user/encryption-keys/
#  To make a long story short, the program may be installed on Linux with:      gem install travis
#  the command to generate the FTP_PASSWORD string is
#
#              travis encrypt  FTP_USER="password" --skip-version-check -r RuleWorld/nfsim
#
#  The password should be surroundedd by double quotes, and the -r parameter indicates the repository for which the password
#  will be used.
#


