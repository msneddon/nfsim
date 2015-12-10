import shutil
import platform
import os

destdir = '{0}-{1}'.format(platform.system(), platform.architecture()[0])
os.makedirs(destdir)
if platform.system() != 'Windows':
    shutil.move('NFsim', os.path.join(destdir, 'NFsim'))
else:
    shutil.move('NFsim.exe', os.path.join(destdir, 'NFsim.exe'))
