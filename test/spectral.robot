*** Settings ***
Library           spectralRobotLibrary.py	${WORKDIR}
Library           OperatingSystem
Library           Collections
Library		  Process
Suite Setup       Start Spectral Unit Tests

*** Variables ***
${WORKDIR}	/tmp/spectral

*** Test Cases ***
Build Spectral
    [Documentation]	Build Spectral 
    [Tags]	Setup 
    Build Spectral
   

*** Keywords *** 
Build Spectral
    [Arguments]
    [Timeout]	10 seconds
    Build System	${WORKDIR}
    Status Should Not Contain	error
    Status Should Not Contain	failure



Start Spectral Unit Tests
    log	Starting Spectral Unit Testing

End Spectral Unit Tests
    log	Terminating Spectral Unit Testing



*** Settings ***
Suite Teardown    End Spectral Unit Tests
