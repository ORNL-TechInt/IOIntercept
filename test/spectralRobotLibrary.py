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
        if expected_status != self._result:
            raise AssertionError("Expected result to be '%s' but was '%s'."
                                 % (expected_status, self_result))

    def status_should_not_contain(self, error_indicator):
	if error_indicator.lower() in self._result.lower():
	    raise AssertionError("Result contains '%s' indicating failure." 
				% (error_indicator))

    def prelim(self):
	self._result = "OK"

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
