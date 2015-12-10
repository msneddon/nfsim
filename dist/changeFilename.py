import shutil
import platform

if platform.system() != 'Windows':
    shutil.move('NFsim', 'NFsim-{0}-{1}'.format(platform.system(), platform.architecture()[0]))
else:
    shutil.move('NFsim.exe', 'NFsim-{0}-{1}.exe'.format(platform.system(), platform.architecture()[0]))
