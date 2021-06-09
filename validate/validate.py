import unittest
import os
import numpy as np
import subprocess
import re
import fnmatch
import sys
import bionetgen

# bngPath = os.path.join('.', 'BioNetGen-2.2.6-stable', 'BNG2.pl')
# bngPath = '/home/boltzmann/apps/BioNetGen-2.3.1/BNG2.pl'  # <<< SET YOUR BIONETGEN PATH HERE <<<
# bngPath = '/home/boltzmann/apps/BioNetGen-2.6.0/BNG2.pl'  # <<< SET YOUR BIONETGEN PATH HERE <<<
# nfsimPath = '/home/boltzmann/apps/BioNetGen-2.6.0/bin/NFsim'  # <<< SET YOUR BIONETGEN PATH HERE <<<
bngPath = os.path.join(bionetgen.defaults.bng_path, "BNG2.pl")
nfsimPath = os.path.join('..', 'build', 'NFsim')


class ParametrizedTestCase(unittest.TestCase):

    """ TestCase classes that want to be parametrized should
        inherit from this class.
    """

    def __init__(self, methodName='runTest', param=None):
        super(ParametrizedTestCase, self).__init__(methodName)
        self.param = param

    @staticmethod
    def parametrize(testcase_klass, param=None):
        """ Create a suite containing all tests taken from the given
            subclass, passing them the parameter 'param'.
        """
        testloader = unittest.TestLoader()
        testnames = testloader.getTestCaseNames(testcase_klass)
        suite = unittest.TestSuite()
        for name in testnames:
            suite.addTest(testcase_klass(name, param=param))
        return suite


def loadResults(fileName, split):
    try:
        with open(fileName) as dataInput:
            timeCourse = []
            # remove spaces
            line = dataInput.readline().strip()
            headers = re.sub('\s+', ' ', line).split(split)

            for line in dataInput:
                nline = re.sub('\s+', ' ', line.strip()).split(' ')
                try:
                    timeCourse.append([float(x) for x in nline])
                except:
                    print('++++', nline)
        return headers, np.array(timeCourse)
    except IOError:
        print('no file')
        return [], []


class TestNFSimFile(ParametrizedTestCase):
    

    # XXX:ideally this should be done through the console but I'm doing the quick and dirty version right now
    def BNGtrajectoryGeneration(self, outputDirectory, fileNumber):
        bngFileName = os.path.join(outputDirectory, 'v{0}.bngl'.format(fileNumber))
        with open(os.devnull, "w") as fnull:
            subprocess.check_call(['perl', bngPath, '-outdir', outputDirectory, bngFileName], stdout=fnull)

    def NFsimtrajectoryGeneration(self, outputDirectory, fileNumber, runOptions):
        runOptions = [x.strip() for x in runOptions.split(' ')]
        with open(os.devnull, "w") as fnull:

            subprocess.check_call([nfsimPath, '-xml', os.path.join(outputDirectory, 'v{0}.xml'.format(fileNumber)),
                                   '-o', os.path.join(outputDirectory, 'v{0}_nf.gdat'.format(fileNumber))] + runOptions,
                                  stdout=fnull)

    def loadConfigurationFile(self, outputDirectory, fileNumber):
        with open(os.path.join(outputDirectory, 'r{0}.txt').format(fileNumber), 'r') as f:
            return f.readlines()

    def test_nfsim(self):
        tol = 0.35 # this is the error tolerance when comparing nfsim's run to the ssa where 0.35 = 35%
        (modelName, runOptions) = self.loadConfigurationFile(self.param['odir'], self.param['num'])
        print('Processing model "{0}"'.format(modelName.strip()))
        ssaDiff = nfDiff = 0
        for index in range(self.param['iterations']):
            self.BNGtrajectoryGeneration(self.param['odir'], self.param['num'])
            self.NFsimtrajectoryGeneration(self.param['odir'], self.param['num'], runOptions)
            odeh, ode = loadResults(os.path.join(self.param['odir'], 'v{0}_ode.gdat'.format(self.param['num'])), ' ')
            ssah, ssa = loadResults(os.path.join(self.param['odir'], 'v{0}_ssa.gdat'.format(self.param['num'])), ' ')
            nfh, nf = loadResults(os.path.join(self.param['odir'], 'v{0}_nf.gdat'.format(self.param['num'])), ' ')

            #square root difference
            ssaDiff += pow(sum(pow(ode[:, 1:] - ssa[:, 1:], 2)), 0.5)
            nfDiff += pow(sum(pow(ode[:, 1:] - nf[:, 1:], 2)), 0.5)
            #nfDiff += sum(abs(ode[:, 1:] - nf[:, 1:]))

        ssaDiff = ssaDiff / self.param['iterations']
        nfDiff = nfDiff / self.param['iterations']

        #rdiff = (nfDiff - ssaDiff) / nfDiff
        rdiff = nfDiff - ssaDiff - (tol * ssaDiff)
        # relative difference should be less than 'tol'
        for element in rdiff:
            # self.assertTrue(element < tol)
            self.assertTrue(element <= 0)


def getTests(directory):
    """
    Gets a list of bngl files that could be correctly translated in a given 'directory'
    """
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*txt'):
            matches.append(''.join(filename.split('.')[0][1:]))
    return sorted(matches)

if __name__ == "__main__":
    suite = unittest.TestSuite()
    if len(sys.argv) > 1:
        os.chdir(sys.argv[1])
    testFolder = './basicModels'
    tests = getTests(testFolder)
    for index in tests:
        suite.addTest(ParametrizedTestCase.parametrize(TestNFSimFile, param={'num': index,
                      'odir': 'basicModels', 'iterations': 30}))
    result = unittest.TextTestRunner(verbosity=2).run(suite)

    ret = (list(result.failures) == [] and list(result.errors) == [])
    ret = 0 if ret else 1
    if ret > 0:
        sys.exit("Validation return an error code")
    else:
        sys.exit()
