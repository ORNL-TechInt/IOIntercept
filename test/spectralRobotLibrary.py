#    spectralRobotLibrary.py

import os
import subprocess
import sys
import json
import time
from robot.libraries.BuiltIn import BuiltIn

TIMESTAMP = time.ctime

class spectralRobotLibrary(object):

    def __init__(self, workdir):
        self._work_path= workdir
        self._result = ''

    def status_should_be(self, expected_status):
        if expected_status.strip() != self._result.strip():
            raise AssertionError("Expected result to be '%s' but was '%s'."
                                 % (expected_status, self._result))

    def status_should_not_contain(self, error_indicator):
	if error_indicator.lower() in self._result.lower():
	    raise AssertionError("Result contains '%s' indicating failure." 
				% (error_indicator))

    def build_system(self, working_dir):
	#Get Current Location
	cwd = os.getcwd()

	#Directory probably already exists
	try:
    	    os.mkdir(working_dir)
	except OSError:
	    pass
    	
	#Move to working directory
	os.chdir(working_dir)	

	#Execute Cmake from base directory
	command = ['cmake3',cwd + "/../"]
	process = subprocess.Popen(command, universal_newlines=True, stdout=subprocess.PIPE,
				   stderr = subprocess.STDOUT)

	self._result = process.communicate()[0].strip()

	if 'Generating done' not in self._result:
		self._result = 'failure'
		return

	command = ['make']
	process = subprocess.Popen(command, universal_newlines=True, stdout=subprocess.PIPE,
				   stderr = subprocess.STDOUT)

		
	self._result = process.communicate()[0]



    def build_gtc(self, working_dir):
	#Get Current Location
	cwd = os.getcwd()

	#Directory probably already exists
	try:
    	    os.mkdir(working_dir)
	except OSError:
	    pass
    	
	#Move to working directory
	os.chdir(working_dir)	

	command = ['make','clean']
	process = subprocess.Popen(command, universal_newlines=True, stdout=subprocess.PIPE,
				   stderr = subprocess.STDOUT)

		
	command = ['make']
	process = subprocess.Popen(command, universal_newlines=True, stdout=subprocess.PIPE,
				   stderr = subprocess.STDOUT)

	self._result = process.communicate()[0]

    def run_simple(self, ver, working_dir, run_cnt):
	os.chdir(working_dir)
	pfs_dir = os.environ['PFS_DIR']
	if ver == 'Static':
	    command = ['./test_static',run_cnt]
	else:
	    command = ['env','LD_PRELOAD=../src/libspectral.so','./test_dynamic',run_cnt]

	process = subprocess.Popen(command, universal_newlines=True, stdout=subprocess.PIPE,
	    			    stderr = subprocess.STDOUT)

    	self._result = process.communicate()[0]
	cnt = str(len([name for name in os.listdir(pfs_dir) if os.path.isfile(os.path.join(pfs_dir, name))]))
	self._result = cnt

    def run_mpi(self, working_dir, mpi_bin):
	mpi_run = mpi_bin + "/mpirun"
	os.chdir(working_dir)
	pfs_dir = os.environ['PFS_DIR']
	command = [mpi_run.strip(),"--allow-run-as-root","-np", "1","./test_mpi_static"]

	sys.stderr.write(command)
	sys.stderr.write('\n')
	process = subprocess.Popen(command, universal_newlines=True, stdout=subprocess.PIPE,
	    			    stderr = subprocess.STDOUT)

	self._result = process.communicate()[0]
	sys.stderr.write(self._result)
	cnt = str(len([name for name in os.listdir(pfs_dir) if os.path.isfile(os.path.join(pfs_dir, name))]))
	sys.stderr.write("%s" % (cnt))
	self._result = cnt
